#!/usr/bin/env python3
# Copyright (c) 2019-2026 Arjan van der Velde, Weng Lab
# SPDX-License-Identifier: BSD-2-Clause

"""Run the legacy-vs-modern parser parity harness over a large corpus.

Driver for ``bin/compare_parsers`` sized for corpora far larger than one
process should walk alone (hundreds of thousands of PDB files). It
enumerates the input tree deterministically, splits the file list into
shards, runs one harness process per shard across a worker pool, and
stores every shard's summary, mismatch report, and exit status under an
output directory together with a manifest. The output directory is
self-contained so it can be archived and copied back to the development
repository, where ``tools/parity/merge_parity.py`` aggregates and
compares runs.

The harness reads input files by the paths given in each shard list, so
the lists carry paths relative to the input root and the harness runs
with the input root as its working directory; mismatch reports then name
files identically no matter which server produced them.

Usage (on the remote server, from the bundle root):

    make compare-parsers
    python3 tools/parity/run_parity.py --input /path/to/input \\
        --output /path/to/parity-results --compress
"""

from __future__ import annotations

import argparse
import concurrent.futures
import datetime
import hashlib
import json
import os
import platform
import re
import socket
import subprocess
import sys
import time
from dataclasses import dataclass, field
from pathlib import Path
from typing import Optional

from common import COUNT_FIELDS, as_dict, as_float, as_int, as_list, as_str

SUMMARY_RE = re.compile(
    r"files (?P<files>\d+) \| lines (?P<lines>\d+) "
    r"\| type mismatches (?P<type_mismatches>\d+) "
    r"\| text mismatches (?P<text_mismatches>\d+) "
    r"\| field mismatches (?P<field_mismatches>\d+) "
    r"\| state mismatches (?P<state_mismatches>\d+) "
    r"\| rejections (?P<rejections>\d+)"
)
BASELINE_RE = re.compile(
    r"baseline: expected ZDATOM (?P<zdatom>\d+) "
    r"\| expected v5 token (?P<v5_token>\d+) "
    r"\| UNEXPECTED (?P<unexpected>\d+)"
)


def parse_args(argv: Optional[list[str]] = None) -> argparse.Namespace:
    """Parse command-line arguments."""
    repo_root = Path(__file__).resolve().parents[2]
    parser = argparse.ArgumentParser(
        description="Sharded legacy-vs-modern parser parity run."
    )
    parser.add_argument(
        "--input",
        required=True,
        type=Path,
        help="corpus root directory to enumerate",
    )
    parser.add_argument(
        "--binary",
        type=Path,
        default=repo_root / "bin" / "compare_parsers",
        help="path to the compare_parsers harness binary (default: %(default)s)",
    )
    parser.add_argument(
        "--output",
        type=Path,
        default=None,
        help="output directory (default: parity-results-<timestamp> in the "
        "current directory)",
    )
    parser.add_argument(
        "--jobs",
        type=int,
        default=None,
        help="parallel harness processes (default: one per CPU)",
    )
    parser.add_argument(
        "--shard-size",
        type=int,
        default=2000,
        help="input files per shard (default: %(default)s)",
    )
    parser.add_argument(
        "--max-files",
        type=int,
        default=0,
        help="enumerate at most N files, 0 for all (default: %(default)s)",
    )
    parser.add_argument(
        "--baseline",
        type=Path,
        default=None,
        help="optional baseline dump executable passed through to the harness",
    )
    parser.add_argument(
        "--all-extensions",
        action="store_true",
        help="visit every file, not only *.pdb",
    )
    parser.add_argument(
        "--quiet",
        action="store_true",
        help="suppress per-line mismatch reports",
    )
    parser.add_argument(
        "--resume",
        action="store_true",
        help="skip shards whose status file records a completed run",
    )
    parser.add_argument(
        "--compress",
        action="store_true",
        help="write <output>.tar.gz next to the output directory when finished",
    )
    return parser.parse_args(argv)


def sha256_of(path: Path) -> str:
    """Return the hex SHA-256 digest of a file."""
    digest = hashlib.sha256()
    with path.open("rb") as handle:
        for chunk in iter(lambda: handle.read(1 << 20), b""):
            digest.update(chunk)
    return digest.hexdigest()


def repo_version(repo_root: Path) -> Optional[str]:
    """Return git describe output for the repository, or None."""
    info_file = repo_root / "BUNDLE_INFO.txt"
    if info_file.is_file():
        text = info_file.read_text(encoding="utf-8").strip()
        return text or None
    if not (repo_root / ".git").exists():
        return None
    try:
        result = subprocess.run(
            ["git", "-C", str(repo_root), "describe", "--always", "--dirty"],
            capture_output=True,
            text=True,
            timeout=10,
            check=False,
        )
    except (OSError, subprocess.TimeoutExpired):
        return None
    if result.returncode != 0:
        return None
    return result.stdout.strip() or None


def enumerate_files(root: Path, all_extensions: bool, max_files: int) -> list[str]:
    """Return sorted paths under root, relative to root, POSIX style."""
    found: list[str] = []
    for dirpath, dirnames, filenames in os.walk(root):
        dirnames.sort()
        for name in sorted(filenames):
            if not all_extensions and not name.lower().endswith(".pdb"):
                continue
            relative = Path(dirpath, name).relative_to(root)
            found.append(relative.as_posix())
            if max_files and len(found) >= max_files:
                return found
    return found


def parse_summary(text: str) -> Optional[dict[str, int]]:
    """Parse the harness summary line into a counter dictionary."""
    match = SUMMARY_RE.search(text)
    if match is None:
        return None
    counters = {name: int(match.group(name)) for name in COUNT_FIELDS}
    baseline = BASELINE_RE.search(text)
    if baseline is not None:
        counters["baseline_expected_zdatom"] = int(baseline.group("zdatom"))
        counters["baseline_expected_v5_token"] = int(baseline.group("v5_token"))
        counters["baseline_unexpected"] = int(baseline.group("unexpected"))
    return counters


def shard_status_path(output: Path, index: int) -> Path:
    """Return the status file path for a shard index."""
    return output / "shards" / f"shard-{index:05d}.status.json"


@dataclass
class ShardStatus:
    """Recorded outcome of one harness shard invocation."""

    shard: int
    files: int
    returncode: int
    seconds: float
    counters: Optional[dict[str, int]]
    summary: str
    complete: bool

    def to_json(self) -> dict[str, object]:
        """Return the JSON-ready status written to disk."""
        return {
            "shard": self.shard,
            "files": self.files,
            "returncode": self.returncode,
            "seconds": self.seconds,
            "counters": self.counters,
            "summary": self.summary,
            "complete": self.complete,
        }

    @classmethod
    def from_json(cls, data: dict[str, object]) -> ShardStatus:
        """Rebuild a status from its on-disk JSON structure."""
        raw_counters = data.get("counters")
        counters = (
            {as_str(key): as_int(value) for key, value in as_dict(raw_counters).items()}
            if raw_counters is not None
            else None
        )
        return cls(
            shard=as_int(data.get("shard", 0)),
            files=as_int(data.get("files", 0)),
            returncode=as_int(data.get("returncode", 0)),
            seconds=as_float(data.get("seconds", 0.0)),
            counters=counters,
            summary=as_str(data.get("summary", "")),
            complete=data.get("complete") is True,
        )


def load_shard_status(output: Path, index: int) -> Optional[ShardStatus]:
    """Load one shard status from disk, or None when absent/invalid."""
    status_path = shard_status_path(output, index)
    if not status_path.is_file():
        return None
    try:
        data = json.loads(status_path.read_text(encoding="utf-8"))
    except (OSError, ValueError):
        return None
    if not isinstance(data, dict):
        return None
    return ShardStatus.from_json(data)


def completed_shard(output: Path, index: int) -> bool:
    """Return True when a shard already ran to completion (resume mode)."""
    status = load_shard_status(output, index)
    return status is not None and status.complete


def shard_command(
    binary: Path, list_path: Path, options: argparse.Namespace
) -> list[str]:
    """Build the harness command line for one shard."""
    command = [str(binary), "CORPUS", "--files-from", str(list_path)]
    if options.baseline is not None:
        command += ["--baseline", str(options.baseline)]
    if options.quiet:
        command.append("--quiet")
    return command


def count_list_files(list_path: Path) -> int:
    """Return the number of non-empty lines in a shard list."""
    return sum(
        1 for line in list_path.read_text(encoding="utf-8").splitlines() if line.strip()
    )


def write_status(status: ShardStatus, output: Path) -> ShardStatus:
    """Finalize a shard status from its recorded returncode and write it."""
    finished = ShardStatus(
        shard=status.shard,
        files=status.files,
        returncode=status.returncode,
        seconds=status.seconds,
        counters=status.counters,
        summary=status.summary,
        complete=status.returncode in (0, 1) and status.counters is not None,
    )
    shard_status_path(output, finished.shard).write_text(
        json.dumps(finished.to_json(), indent=2, sort_keys=True) + "\n",
        encoding="utf-8",
    )
    return finished


def write_manifest(
    output: Path,
    options: argparse.Namespace,
    input_root: Path,
    binary: Path,
    file_count: int,
    shard_count: int,
) -> None:
    """Write the run manifest describing tool and environment."""
    repo_root = Path(__file__).resolve().parents[2]
    manifest = {
        "generated_at": datetime.datetime.now(datetime.timezone.utc).isoformat(
            timespec="seconds"
        ),
        "hostname": socket.gethostname(),
        "platform": platform.platform(),
        "python": sys.version.split()[0],
        "cpu_count": os.cpu_count(),
        "tool_version": repo_version(repo_root),
        "binary": str(binary),
        "binary_sha256": sha256_of(binary),
        "baseline": str(options.baseline) if options.baseline is not None else None,
        "baseline_sha256": (
            sha256_of(options.baseline) if options.baseline is not None else None
        ),
        "input_root": str(input_root),
        "files": file_count,
        "shard_size": options.shard_size,
        "shards": shard_count,
        "jobs": options.jobs,
        "all_extensions": options.all_extensions,
        "quiet": options.quiet,
    }
    (output / "manifest.json").write_text(
        json.dumps(manifest, indent=2, sort_keys=True) + "\n", encoding="utf-8"
    )


def add_counters(total: dict[str, int], counters: dict[str, int]) -> None:
    """Accumulate one shard's counters into the totals."""
    for name, value in counters.items():
        total[name] = total.get(name, 0) + value


@dataclass
class RunSummary:
    """Aggregated counters and shard outcomes for the whole run."""

    totals: dict[str, int] = field(default_factory=dict)
    shards: int = 0
    shard_errors: list[int] = field(default_factory=list)
    incomplete: bool = False

    def to_json(self) -> dict[str, object]:
        """Return the JSON-ready summary written to disk."""
        return {
            "totals": self.totals,
            "shards": self.shards,
            "shard_errors": self.shard_errors,
            "incomplete": self.incomplete,
        }

    @classmethod
    def from_json(cls, data: dict[str, object]) -> RunSummary:
        """Rebuild a summary from its on-disk JSON structure."""
        return cls(
            totals={
                as_str(key): as_int(value)
                for key, value in as_dict(data.get("totals", {})).items()
            },
            shards=as_int(data.get("shards", 0)),
            shard_errors=[
                as_int(item) for item in as_list(data.get("shard_errors", []))
            ],
            incomplete=data.get("incomplete") is True,
        )


def aggregate(statuses: list[ShardStatus], shard_count: int) -> RunSummary:
    """Combine shard statuses into the run summary structure."""
    totals: dict[str, int] = {name: 0 for name in COUNT_FIELDS}
    for status in statuses:
        if status.counters is not None:
            add_counters(totals, status.counters)
    errors = [status.shard for status in statuses if not status.complete]
    return RunSummary(
        totals=totals,
        shards=shard_count,
        shard_errors=errors,
        incomplete=len(statuses) < shard_count or bool(errors),
    )


def print_status(index: int, total: int, status: ShardStatus) -> None:
    """Print one progress line for a finished shard."""
    counters = status.counters or {}
    interesting = " ".join(
        f"{name} {counters[name]}" for name in COUNT_FIELDS[2:] if counters.get(name)
    )
    detail = f" ({interesting})" if interesting else ""
    print(
        f"[{index + 1}/{total}] shard-{index:05d}: "
        f"{status.files} files, "
        f"{counters.get('lines', 0)} lines, "
        f"{status.seconds:.1f}s{detail}",
        flush=True,
    )


def compress_output(output: Path) -> Path:
    """Create a gzip tarball of the output directory and return its path."""
    archive = output.with_suffix(".tar.gz")
    subprocess.run(
        ["tar", "-czf", str(archive), "-C", str(output.parent), output.name],
        check=True,
        timeout=3600,
    )
    return archive


def validate_options(options: argparse.Namespace) -> Optional[str]:
    """Resolve and validate paths and defaults; return an error message."""
    options.input = options.input.resolve()
    options.binary = options.binary.resolve()
    if not options.input.is_dir():
        return f"input directory not found: {options.input}"
    if not options.binary.is_file():
        return (
            f"harness binary not found: {options.binary} "
            "(build it with 'make compare-parsers')"
        )
    if options.baseline is not None:
        options.baseline = options.baseline.resolve()
        if not options.baseline.is_file():
            return f"baseline executable not found: {options.baseline}"
    if options.jobs is None:
        options.jobs = os.cpu_count() or 1
    if options.shard_size < 1:
        return "--shard-size must be at least 1"
    return None


def write_shard_lists(output: Path, files: list[str], shard_size: int) -> list[int]:
    """Write one file list per shard and return the shard indexes."""
    shard_dir = output / "shards"
    shard_dir.mkdir(parents=True, exist_ok=True)
    starts = range(0, len(files), shard_size)
    for index, start in enumerate(starts):
        list_path = shard_dir / f"shard-{index:05d}.list"
        list_path.write_text(
            "\n".join(files[start : start + shard_size]) + "\n",
            encoding="utf-8",
        )
    return list(range(len(starts)))


def execute_shards(
    binary: Path,
    input_root: Path,
    output: Path,
    indexes: list[int],
    options: argparse.Namespace,
) -> bool:
    """Run the pending shards across the worker pool.

    Returns True when the run was interrupted; shard statuses are always
    left on disk for aggregation.
    """
    pending = [
        index
        for index in indexes
        if not (options.resume and completed_shard(output, index))
    ]
    skipped = len(indexes) - len(pending)
    if skipped:
        print(f"resume: skipping {skipped} completed shards", flush=True)
    pool = concurrent.futures.ThreadPoolExecutor(max_workers=options.jobs)
    try:
        futures = {
            pool.submit(
                run_and_record, binary, input_root, output, index, options
            ): index
            for index in pending
        }
        for future in concurrent.futures.as_completed(futures):
            status = future.result()
            print_status(status.shard, len(indexes), status)
        pool.shutdown(wait=True)
        return False
    except KeyboardInterrupt:
        pool.shutdown(wait=False, cancel_futures=True)
        print(
            "interrupted; use --resume on the same output directory to continue",
            file=sys.stderr,
        )
        return True


def run_and_record(
    binary: Path,
    input_root: Path,
    output: Path,
    index: int,
    options: argparse.Namespace,
) -> ShardStatus:
    """Run one shard, persist its status, and return it."""
    started = time.monotonic()
    shard_dir = output / "shards"
    list_path = shard_dir / f"shard-{index:05d}.list"
    summary_path = shard_dir / f"shard-{index:05d}.summary.txt"
    report_path = shard_dir / f"shard-{index:05d}.mismatches.log"
    with summary_path.open("wb") as out, report_path.open("wb") as err:
        process = subprocess.run(
            shard_command(binary, list_path, options),
            stdout=out,
            stderr=err,
            cwd=str(input_root),
            check=False,
        )
    summary_text = summary_path.read_text(encoding="utf-8", errors="replace")
    return write_status(
        ShardStatus(
            shard=index,
            files=count_list_files(list_path),
            returncode=process.returncode,
            seconds=round(time.monotonic() - started, 3),
            counters=parse_summary(summary_text),
            summary=summary_text.strip(),
            complete=False,
        ),
        output,
    )


def report_run(output: Path, summary: RunSummary, interrupted: bool) -> int:
    """Write the summary file and print the run report; return the code."""
    (output / "summary.json").write_text(
        json.dumps(summary.to_json(), indent=2, sort_keys=True) + "\n",
        encoding="utf-8",
    )
    if interrupted:
        return 130
    totals = summary.totals
    print(
        f"files {totals['files']} | lines {totals['lines']} | "
        f"type mismatches {totals['type_mismatches']} | "
        f"text mismatches {totals['text_mismatches']} | "
        f"field mismatches {totals['field_mismatches']} | "
        f"state mismatches {totals['state_mismatches']} | "
        f"rejections {totals['rejections']}"
    )
    if summary.shard_errors:
        print(
            f"error: {len(summary.shard_errors)} shard(s) failed: "
            f"{summary.shard_errors[:10]}",
            file=sys.stderr,
        )
    print(f"results in {output}; aggregate with tools/parity/merge_parity.py {output}")
    return 0 if not summary.shard_errors else 1


def collect_statuses(output: Path, shard_count: int) -> list[ShardStatus]:
    """Load every recorded shard status from disk, ordered by shard."""
    statuses = []
    for index in range(shard_count):
        status = load_shard_status(output, index)
        if status is not None:
            statuses.append(status)
    statuses.sort(key=lambda status: status.shard)
    return statuses


def load_shard_lists(output: Path) -> list[int]:
    """Return the shard indexes of the lists already on disk."""
    shard_dir = output / "shards"
    return sorted(
        int(path.stem.split("-")[1]) for path in shard_dir.glob("shard-*.list")
    )


def default_output() -> Path:
    """Return the timestamped default output directory name."""
    stamp = datetime.datetime.now().strftime("%Y%m%d-%H%M%S")
    return Path(f"parity-results-{stamp}")


def main(argv: Optional[list[str]] = None) -> int:
    """Run the sharded parity harness and return the process exit code."""
    options = parse_args(argv)
    error = validate_options(options)
    if error is not None:
        print(f"error: {error}", file=sys.stderr)
        return 2
    if options.output is None:
        options.output = default_output()
    output = options.output.resolve()

    if options.resume:
        # a resume reuses the shard lists of the original run verbatim;
        # re-deriving them could pair a different shard-size with stale
        # per-shard statuses
        indexes = load_shard_lists(output)
        if not indexes:
            print(
                f"error: no shard lists under {output / 'shards'} to resume",
                file=sys.stderr,
            )
            return 2
        file_count = sum(
            count_list_files(output / "shards" / f"shard-{index:05d}.list")
            for index in indexes
        )
    else:
        print(f"enumerating {options.input} ...", flush=True)
        files = enumerate_files(
            options.input, options.all_extensions, options.max_files
        )
        if not files:
            print("error: no input files found", file=sys.stderr)
            return 2
        indexes = write_shard_lists(output, files, options.shard_size)
        write_manifest(
            output, options, options.input, options.binary, len(files), len(indexes)
        )
        file_count = len(files)
    print(
        f"{file_count} files in {len(indexes)} shards, "
        f"{options.jobs} jobs; output {output}",
        flush=True,
    )

    interrupted = execute_shards(
        options.binary, options.input, output, indexes, options
    )
    statuses = collect_statuses(output, len(indexes))
    summary = aggregate(statuses, len(indexes))
    if interrupted:
        summary.incomplete = True
    code = report_run(output, summary, interrupted)
    if code == 130:
        return code

    if options.compress:
        archive = compress_output(output)
        print(f"archived {archive} ({archive.stat().st_size // (1 << 20)} MiB)")
    return code


if __name__ == "__main__":
    sys.exit(main())
