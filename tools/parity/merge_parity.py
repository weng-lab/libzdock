#!/usr/bin/env python3
# Copyright (c) 2019-2026 Arjan van der Velde, Weng Lab
# SPDX-License-Identifier: BSD-2-Clause

"""Aggregate and compare parity-harness runs copied back from servers.

Reads the output directories written by ``tools/parity/run_parity.py``,
combines the shard counters, and deduplicates the per-line mismatch
reports into normalized signatures (record kind plus the divergence
position, with digits masked) so that millions of reports collapse into
the handful of underlying bugs. The merged structure is written as
``merged.json`` and the top signatures are printed.

Usage:

    python3 tools/parity/merge_parity.py RESULTS_DIR [MORE_DIRS ...]
        [--json PATH] [--top N]

    python3 tools/parity/merge_parity.py --compare A/merged.json
        B/merged.json [--top N]

The compare form is the fix-verification loop: after a parser fix,
re-run the corpus remotely, merge again, and diff the two merged files.
"""

from __future__ import annotations

import argparse
import datetime
import json
import re
import sys
from collections import Counter
from dataclasses import dataclass, field
from pathlib import Path
from typing import Optional

from common import COUNT_FIELDS, as_dict, as_int, as_list, as_str

# One report line: "<path>:<line>: <what>: <detail>"; the path and the
# what part may themselves contain colons, so known what-prefixes are
# matched explicitly after the line number.
REPORT_RE = re.compile(r"^(?P<path>.+):(?P<line>\d+): (?P<rest>.*)$")
WHAT_PREFIXES = (
    "new parser rejected: ",
    "UNEXPECTED baseline difference: ",
    "PDBRUN state mismatch: ",
    "field mismatch: ",
    "text mismatch: ",
    "type mismatch: ",
)
# Input-quality census classes: reported separately from parity
# mismatches and never part of any exit decision.
INPUT_WHAT_PREFIXES = ("coordinate unknown: ",)
INPUT_WHATS = ("coordinate unknown",)
QUOTE_RE = re.compile(r"^legacy '(.*?)' new '(.*)'$", re.DOTALL)
BASELINE_QUOTE_RE = re.compile(r"^baseline '(.*?)' modern '(.*)'$", re.DOTALL)


def norm(text: str) -> str:
    """Mask digits and collapse whitespace so variants share a signature."""
    masked = re.sub(r"\d+", "#", text)
    return re.sub(r"\s+", " ", masked).strip()


def truncate(text: str, limit: int) -> str:
    """Cut a signature component to a bounded length."""
    return text if len(text) <= limit else text[:limit] + "..."


def first_word(text: str) -> str:
    """Return the leading record name of a raw line, if any."""
    stripped = text.lstrip()
    word = stripped.split(" ", 1)[0] if stripped else ""
    return word[:12]


def divergence_signature(template: str, left: str, right: str) -> str:
    """Describe the first differing column between two rendered texts."""
    index = 0
    while index < len(left) and index < len(right) and left[index] == right[index]:
        index += 1
    window = 10
    lo = max(0, index - 4)
    return (
        f"{template}@col {index} "
        f"{truncate(norm(left[lo:index + window]), 24)} != "
        f"{truncate(norm(right[lo:index + window]), 24)}"
    )


def quoted_pair(detail: str, pattern: re.Pattern[str]) -> Optional[tuple[str, str]]:
    """Extract the (left, right) quoted pair from a report detail."""
    match = pattern.match(detail)
    if match is None:
        return None
    return match.group(1), match.group(2)


def rejected_signature(detail: str) -> str:
    """Signature for a modern rejection of a legacy-accepted line."""
    head, _, line = detail.partition(" | ")
    return f"{head} | {first_word(line)}"


def type_signature(detail: str) -> str:
    """Signature for a record-classification divergence."""
    head, _, line = detail.partition(" | ")
    return f"type {norm(head)} | {first_word(line)}"


def state_signature(detail: str) -> str:
    """Signature for a PDBRUN version-state divergence."""
    return f"state {norm(detail)}"


def divergence_of(
    template: str, detail: str, pattern: re.Pattern[str]
) -> Optional[str]:
    """Signature at the first differing column of a quoted text pair."""
    pair = quoted_pair(detail, pattern)
    if pair is None:
        return None
    return divergence_signature(template, pair[0], pair[1])


SPECIAL_SIGNATURES = {
    "new parser rejected": rejected_signature,
    "type mismatch": type_signature,
    "PDBRUN state mismatch": state_signature,
}


def signature_for(what: str, detail: str) -> str:
    """Map one report to its normalized dedup signature."""
    special = SPECIAL_SIGNATURES.get(what)
    if special is not None:
        return special(detail)
    if what == "text mismatch":
        signature = divergence_of("text", detail, QUOTE_RE)
        if signature is not None:
            return signature
    elif what == "UNEXPECTED baseline difference":
        signature = divergence_of("baseline", detail, BASELINE_QUOTE_RE)
        if signature is not None:
            return signature
    elif what.startswith("field "):
        return truncate(f"{what} {norm(detail)}", 120)
    return truncate(f"{what}: {norm(detail)}", 120)


def parse_report_line(raw: str) -> Optional[tuple[str, str, str, str]]:
    """Split one report line into (path, line, what, detail)."""
    match = REPORT_RE.match(raw)
    if match is None:
        return None
    rest = match.group("rest")
    for prefix in INPUT_WHAT_PREFIXES + WHAT_PREFIXES:
        if rest.startswith(prefix):
            return (
                match.group("path"),
                match.group("line"),
                prefix[:-2],
                rest[len(prefix) :],
            )
    return match.group("path"), match.group("line"), rest, ""


@dataclass
class SignatureEntry:
    """One normalized mismatch signature with its count and example."""

    signature: str
    count: int
    example: str


@dataclass
class MergedResult:
    """Aggregated counters and deduplicated signatures over run dirs."""

    sources: list[str]
    totals: dict[str, int] = field(default_factory=dict)
    report_lines_unparsed: int = 0
    files_with_mismatches: int = 0
    signatures: list[SignatureEntry] = field(default_factory=list)
    input_signatures: list[SignatureEntry] = field(default_factory=list)

    def to_json(self) -> dict[str, object]:
        """Return the JSON-ready structure written to merged.json."""
        return {
            "generated_at": datetime.datetime.now(datetime.timezone.utc).isoformat(
                timespec="seconds"
            ),
            "sources": self.sources,
            "totals": self.totals,
            "report_lines_unparsed": self.report_lines_unparsed,
            "files_with_mismatches": self.files_with_mismatches,
            "signature_count": len(self.signatures),
            "signatures": [
                {
                    "signature": entry.signature,
                    "count": entry.count,
                    "example": entry.example,
                }
                for entry in self.signatures
            ],
            "input_signatures": [
                {
                    "signature": entry.signature,
                    "count": entry.count,
                    "example": entry.example,
                }
                for entry in self.input_signatures
            ],
        }

    @classmethod
    def from_json(cls, data: dict[str, object]) -> MergedResult:
        """Rebuild a merged result from its merged.json structure."""

        def entries(key: str) -> list[SignatureEntry]:
            return [
                SignatureEntry(
                    as_str(entry.get("signature", "")),
                    as_int(entry.get("count", 0)),
                    as_str(entry.get("example", "")),
                )
                for entry in (
                    as_dict(item) for item in as_list(data.get(key, []))
                )
            ]

        return cls(
            sources=[as_str(entry) for entry in as_list(data.get("sources", []))],
            totals={
                as_str(key): as_int(value)
                for key, value in as_dict(data.get("totals", {})).items()
            },
            report_lines_unparsed=as_int(data.get("report_lines_unparsed", 0)),
            files_with_mismatches=as_int(data.get("files_with_mismatches", 0)),
            signatures=entries("signatures"),
            input_signatures=entries("input_signatures"),
        )


def load_run_totals(results_dir: Path) -> dict[str, int]:
    """Return the summed counters for one results directory."""
    summary_path = results_dir / "summary.json"
    if summary_path.is_file():
        summary = as_dict(json.loads(summary_path.read_text(encoding="utf-8")))
        return {
            as_str(key): as_int(value)
            for key, value in as_dict(summary.get("totals", {})).items()
        }
    totals: dict[str, int] = {}
    for status_path in sorted((results_dir / "shards").glob("*.status.json")):
        status = as_dict(json.loads(status_path.read_text(encoding="utf-8")))
        counters = as_dict(status.get("counters") or {})
        for field_name, value in counters.items():
            totals[as_str(field_name)] = totals.get(as_str(field_name), 0) + as_int(
                value
            )
    return totals


class SignatureAccumulator:
    """Collects normalized signatures, counts, and example locations.

    Parity mismatches and input-quality census reports are kept in
    separate tables: the census describes the input corpus, not parser
    divergence, and is never part of an exit decision.
    """

    def __init__(self) -> None:
        self.counts: Counter[str] = Counter()
        self.examples: dict[str, str] = {}
        self.mismatched_files: set[str] = set()
        self.input_counts: Counter[str] = Counter()
        self.input_examples: dict[str, str] = {}
        self.input_files: dict[str, set[str]] = {}
        self.unparsed_lines = 0

    def add_report(self, raw: str) -> None:
        """Consume one report line (or count it as unparsed noise)."""
        parsed = parse_report_line(raw.rstrip("\n"))
        if parsed is None:
            if raw.strip():
                self.unparsed_lines += 1
            return
        path, line, what, detail = parsed
        signature = signature_for(what, detail)
        location = f"{path}:{line}"
        if what in INPUT_WHATS:
            self.input_counts[signature] += 1
            self.input_examples.setdefault(signature, location)
            self.input_files.setdefault(signature, set()).add(path)
            return
        self.counts[signature] += 1
        self.examples.setdefault(signature, location)
        self.mismatched_files.add(path)

    def add_dir(self, results_dir: Path) -> None:
        """Consume every mismatch report of one results directory."""
        for report_path in sorted((results_dir / "shards").glob("*.mismatches.log")):
            with report_path.open(encoding="utf-8", errors="replace") as handle:
                for raw in handle:
                    self.add_report(raw)


def merge_dirs(dirs: list[Path]) -> MergedResult:
    """Merge one or more results directories into a merged result."""
    merged = MergedResult(sources=[str(directory) for directory in dirs])
    accumulator = SignatureAccumulator()
    for results_dir in dirs:
        for field_name, value in load_run_totals(results_dir).items():
            merged.totals[field_name] = merged.totals.get(field_name, 0) + value
        accumulator.add_dir(results_dir)
    merged.report_lines_unparsed = accumulator.unparsed_lines
    merged.files_with_mismatches = len(accumulator.mismatched_files)
    merged.signatures = [
        SignatureEntry(signature, count, accumulator.examples[signature])
        for signature, count in accumulator.counts.most_common()
    ]
    merged.input_signatures = [
        SignatureEntry(
            signature, count, accumulator.input_examples[signature]
        )
        for signature, count in accumulator.input_counts.most_common()
    ]
    return merged


def print_merge(merged: MergedResult, top: int) -> None:
    """Print the human-readable merge summary."""
    print(" | ".join(f"{name} {merged.totals.get(name, 0)}" for name in COUNT_FIELDS))
    print(
        f"unique signatures {len(merged.signatures)} | files with "
        f"mismatches {merged.files_with_mismatches} | unparsed report "
        f"lines {merged.report_lines_unparsed}"
    )
    if not merged.signatures:
        print("no mismatch reports found (quiet run or full parity)")
    for entry in merged.signatures[:top]:
        print(f"{entry.count:>9}  {truncate(entry.signature, 110)}")
        print(f"            e.g. {entry.example}")
    if merged.input_signatures:
        print(
            "input-quality census (separate from parser mismatches; never "
            "part of any exit decision):"
        )
        for entry in merged.input_signatures[:top]:
            print(f"{entry.count:>9}  {truncate(entry.signature, 110)}")
            print(f"            e.g. {entry.example}")


def print_total_deltas(left: MergedResult, right: MergedResult) -> None:
    """Print the totals that changed between two merged results."""
    print("totals (left -> right, delta):")
    keys = sorted(set(left.totals) | set(right.totals))
    changed = False
    for key in keys:
        before = left.totals.get(key, 0)
        after = right.totals.get(key, 0)
        if before != after:
            changed = True
            print(f"  {key}: {before} -> {after} ({after - before:+d})")
    if not changed:
        print("  unchanged")


def signature_counts(merged: MergedResult) -> Counter[str]:
    """Return signature -> count for one merged result."""
    return Counter({entry.signature: entry.count for entry in merged.signatures})


def compare(left: MergedResult, right: MergedResult, top: int) -> int:
    """Compare two merged results and print fix-verification deltas."""
    print_total_deltas(left, right)
    left_counts = signature_counts(left)
    right_counts = signature_counts(right)
    improved = {
        signature: right_counts[signature] - left_counts[signature]
        for signature in left_counts
        if right_counts[signature] != left_counts[signature]
    }
    introduced = {
        signature: delta
        for signature, delta in right_counts.items() - left_counts.items()
        if delta > 0
    }
    print(
        f"signatures: {len(left_counts)} -> {len(right_counts)} "
        f"({sum(1 for delta in improved.values() if delta < 0)} "
        f"improved, {len(introduced)} new/worse)"
    )
    for label, items in (("improved", improved), ("new/worse", introduced)):
        ordered = sorted(items.items(), key=lambda pair: pair[1], reverse=True)
        for signature, delta in ordered[:top]:
            print(f"  {label}: {delta:+d}  {truncate(signature, 110)}")
    left_input = Counter({entry.signature: entry.count for entry in left.input_signatures})
    right_input = Counter(
        {entry.signature: entry.count for entry in right.input_signatures}
    )
    input_delta = {
        signature: right_input[signature] - left_input[signature]
        for signature in set(left_input) | set(right_input)
        if right_input[signature] != left_input[signature]
    }
    if input_delta:
        print("input-quality census deltas (informational):")
        for signature in sorted(input_delta, key=lambda sig: input_delta[sig],
                                reverse=True)[:top]:
            print(f"  {input_delta[signature]:+d}  {truncate(signature, 110)}")
    return 0 if not introduced else 1


def parse_args(argv: Optional[list[str]] = None) -> argparse.Namespace:
    """Parse command-line arguments."""
    parser = argparse.ArgumentParser(
        description="Aggregate and compare parity-harness results."
    )
    parser.add_argument(
        "dirs", nargs="*", type=Path, help="results directories from run_parity.py"
    )
    parser.add_argument(
        "--compare",
        nargs=2,
        type=Path,
        metavar=("A", "B"),
        help="compare two merged.json files instead",
    )
    parser.add_argument(
        "--json",
        type=Path,
        default=None,
        help="where to write the merged structure "
        "(default: <first-dir>/merged.json)",
    )
    parser.add_argument(
        "--top",
        type=int,
        default=25,
        help="signatures to show in printed output "
        "(default: %(default)s; the JSON file keeps "
        "all)",
    )
    return parser.parse_args(argv)


def write_merged(merged: MergedResult, json_path: Path) -> None:
    """Write the merged result as JSON."""
    json_path.write_text(
        json.dumps(merged.to_json(), indent=2) + "\n", encoding="utf-8"
    )


def main(argv: Optional[list[str]] = None) -> int:
    """Merge or compare parity results and return the process exit code."""
    options = parse_args(argv)
    if options.compare is not None:
        left_path, right_path = options.compare
        left = MergedResult.from_json(json.loads(left_path.read_text(encoding="utf-8")))
        right = MergedResult.from_json(
            json.loads(right_path.read_text(encoding="utf-8"))
        )
        return compare(left, right, options.top)
    if not options.dirs:
        print(
            "error: give at least one results directory (or --compare A B)",
            file=sys.stderr,
        )
        return 2
    for directory in options.dirs:
        if (
            not (directory / "shards").is_dir()
            and not (directory / "summary.json").is_file()
        ):
            print(
                f"error: not a run_parity.py output directory: " f"{directory}",
                file=sys.stderr,
            )
            return 2
    merged = merge_dirs(options.dirs)
    json_path = (
        options.json if options.json is not None else options.dirs[0] / "merged.json"
    )
    write_merged(merged, json_path)
    print(
        f"merged {len(options.dirs)} director"
        f"{'y' if len(options.dirs) == 1 else 'ies'} -> {json_path}"
    )
    print_merge(merged, options.top)
    return 0


if __name__ == "__main__":
    sys.exit(main())
