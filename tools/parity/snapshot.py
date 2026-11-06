#!/usr/bin/env python3
# Copyright (c) 2019-2026 Arjan van der Velde, Weng Lab
# SPDX-License-Identifier: BSD-2-Clause

"""Freeze, verify, and compare the legacy parser output for the corpus.

The legacy libpdb++ parser is being retired; these modes keep its
observable behavior runnable against the modern parser without building
any legacy code:

generate  run the bb0032c baseline dump over the input corpus and write
          its records as zstd-compressed JSON Lines in ~1 GiB chunks
          (never splitting an input file across chunks; every chunk is
          written to a temporary name and renamed when complete, so an
          interrupted run never leaves a partial chunk that looks
          finished), plus a manifest listing every chunk (SHA-256, size,
          record count) and every input file (path, hash, chunk) with
          the dump's commit, compiler, flags, and harness version -- all
          paths relative to the snapshot directory, so it can move.

verify    check every chunk against the manifest (existence, SHA-256,
          size, file and record counts, per-file paths and hashes) and
          report any chunk file the manifest does not know about.

compare   check the modern parser against the snapshot: stream each
          chunk, re-parse every frozen raw input line with the modern
          parser, and classify every record difference into the same
          expected-difference groups as the live bb0032c comparison
          (the pinned pre-fb720ad ZDATOM 80-column clamp, the
          pre-ea22558 unbounded v5 USER token reads, and the deliberate
          segid/element salvage gains); everything else -- including
          parsed-field differences -- is unexpected and fails the run.
          The compare is self-sufficient: the frozen raw lines are the
          modern parser's input, so the live corpus is never touched
          and the snapshot stays valid when files move or vanish.

The snapshot lives outside the repository (parity-private at first,
later /zdockdata on CephFS); no corpus data is ever committed.

Chunk layout (relative to the snapshot directory):

    manifest.json
    chunks/chunk-NNNNNN.jsonl.zst

Each chunk holds one JSON object per input file:

    {"path": "...", "sha256": "...", "records": [
        {"line": N, "type": T, "text": "...", "fields": {...}}, ...]}
"""

from __future__ import annotations

import argparse
import concurrent.futures
import datetime
import hashlib
import json
import os
import shutil
import subprocess
import sys
import tempfile
from pathlib import Path
from typing import Iterator, Optional

from common import as_dict, as_int, as_list, as_str
from run_parity import enumerate_files, sha256_of

SNAPSHOT_VERSION = 1
# zdock::pdb::RecordType::USER_GFX_LABEL by numeric value: the one
# record whose xyz fields the live comparison skips when the modern
# parse left them unset (the v5 matching-failure semantics).
USER_GFX_LABEL_TYPE = 43
XYZ_KEYS = ("xyz[0]", "xyz[1]", "xyz[2]")
# The v5 records whose unbounded %s reads the bounded reads changed.
V5_TOKEN_PREFIXES = ("USER  CNAME ", "USER  COLOR ", "USER  GFX COLOR ",
                     "USER  GFX FONT ")
# The record types whose renders carry a coordinate-record prefix
# (ATOM, HETATM, SIGATM by numeric value in zdock::pdb::RecordType).
ATOMISH_TYPES = (2, 14, 29)
ATOMISH_PREFIXES = ("ATOM", "HETA", "SIGA")


def repo_version() -> Optional[str]:
    """Return the harness source version (the repository commit)."""
    repo = Path(__file__).resolve().parents[2]
    try:
        result = subprocess.run(
            ["git", "-C", str(repo), "describe", "--always", "--dirty"],
            capture_output=True, text=True, timeout=10, check=False)
    except (OSError, subprocess.TimeoutExpired):
        return None
    if result.returncode != 0:
        return None
    return result.stdout.strip() or None


def read_meta(dump: Path) -> dict[str, str]:
    """Return the dump binary's recorded build metadata."""
    result = subprocess.run([str(dump), "--meta"], capture_output=True,
                            text=True, timeout=30, check=True)
    return {as_str(key): as_str(value)
            for key, value in json.loads(result.stdout).items()}


def parse_args(argv: Optional[list[str]] = None) -> argparse.Namespace:
    """Parse command-line arguments."""
    repo_root = Path(__file__).resolve().parents[2]
    parser = argparse.ArgumentParser(
        description="Freeze, verify, and compare legacy parser output.")
    sub = parser.add_subparsers(dest="command", required=True)

    gen = sub.add_parser("generate", help="snapshot the bb0032c dump")
    gen.add_argument("--input", required=True, type=Path,
                     help="corpus root directory")
    gen.add_argument("--dump", required=True, type=Path,
                     help="path to the bb0032c baseline dump binary")
    gen.add_argument("--output", required=True, type=Path,
                     help="snapshot directory to create")
    gen.add_argument("--zstdio", type=Path,
                     default=repo_root / "bin" / "zstdio",
                     help="zstd stream filter (default: %(default)s)")
    gen.add_argument("--jobs", type=int, default=4,
                     help="parallel dump processes (default: %(default)s)")
    gen.add_argument("--chunk-gb", type=float, default=1.0,
                     help="approximate compressed chunk size in GiB "
                          "(default: %(default)s)")
    gen.add_argument("--all-extensions", action="store_true",
                     help="visit every file, not only *.pdb")
    gen.add_argument("--max-files", type=int, default=0,
                     help="enumerate at most N files, 0 for all")

    ver = sub.add_parser("verify", help="check chunks against the manifest")
    ver.add_argument("--snapshot", required=True, type=Path,
                     help="snapshot directory")
    ver.add_argument("--zstdio", type=Path,
                     default=repo_root / "bin" / "zstdio",
                     help="zstd stream filter (default: %(default)s)")

    cmp = sub.add_parser("compare", help="compare the modern parser")
    cmp.add_argument("--snapshot", required=True, type=Path,
                     help="snapshot directory")
    cmp.add_argument("--binary", type=Path,
                     default=repo_root / "bin" / "compare_parsers",
                     help="modern harness binary (default: %(default)s)")
    cmp.add_argument("--zstdio", type=Path,
                     default=repo_root / "bin" / "zstdio",
                     help="zstd stream filter (default: %(default)s)")
    cmp.add_argument("--jobs", type=int, default=4,
                     help="parallel modern-dump processes "
                          "(default: %(default)s)")
    return parser.parse_args(argv)


# ---------------------------------------------------------------------------
# shared streaming helpers
# ---------------------------------------------------------------------------


def stream_chunk(snapshot: Path, zstdio: Path,
                 name: str) -> Iterator[dict[str, object]]:
    """Yield the per-file JSON objects of one chunk, decompressed."""
    with open(snapshot / name, "rb") as chunk_file:
        with subprocess.Popen([str(zstdio), "d"], stdin=chunk_file,
                              stdout=subprocess.PIPE,
                              stderr=subprocess.PIPE) as proc:
            assert proc.stdout is not None
            with proc.stdout:
                for raw in proc.stdout:
                    yield json.loads(raw)
            assert proc.stderr is not None
            stderr = proc.stderr.read().decode(errors="replace")
            proc.stderr.close()
            if proc.wait() != 0:
                raise RuntimeError(f"zstd decompressor failed for {name}: "
                                   f"{stderr}")


def modern_dump(binary: Path, path: str, raw_lines: list[str]
                ) -> dict[str, object]:
    """Run the modern dump over one file's frozen lines.

    The lines are written to a scratch file so the modern parser reads
    them exactly as it reads any input file (the pinned line-buffer
    semantics); the live corpus is never touched.
    """
    scratch = Path(tempfile.mkdtemp()) / "frozen.pdb"
    try:
        # the frozen lines carry only U+0000..U+00FF (the dump escapes
        # every non-ASCII byte), so latin-1 restores the original bytes;
        # an empty file stays empty (no phantom trailing line)
        body = "\n".join(raw_lines)
        if body:
            body += "\n"
        scratch.write_bytes(body.encode("latin-1"))
        result = subprocess.run([str(binary), "--dump-json", str(scratch)],
                                capture_output=True, timeout=600,
                                check=False)
    finally:
        shutil.rmtree(scratch.parent, ignore_errors=True)
    if result.returncode != 0:
        raise RuntimeError(
            f"modern dump failed for {path}: "
            f"{result.stderr.decode(errors='replace')[:200]}")
    return json.loads(result.stdout)


# ---------------------------------------------------------------------------
# generate
# ---------------------------------------------------------------------------


class ChunkWriter:
    """Streams JSON lines through zstd into ~1 GiB compressed chunks.

    Each chunk is written under a temporary name and renamed when
    complete, so an interrupted run never leaves a partial chunk that
    looks finished; the manifest, written last, is the completion
    record.
    """

    def __init__(self, chunks_dir: Path, zstdio: Path, chunk_bytes: int):
        self.chunks_dir = chunks_dir
        self.chunks_dir.mkdir(parents=True, exist_ok=True)
        self.zstdio = zstdio
        self.chunk_bytes = chunk_bytes
        self.index = 0
        self.records = 0
        self.files = 0
        self.process: Optional[subprocess.Popen] = None
        self.tmp_path: Optional[Path] = None
        self.completed: list[dict[str, object]] = []

    def _start(self) -> None:
        # the chunk file deliberately outlives this call: the zstd
        # subprocess writes into it until _finish closes both
        tmp = (self.chunks_dir /
               f".chunk-{self.index:06d}.jsonl.zst.tmp-{os.getpid()}")
        out = tmp.open("wb")
        self.process = subprocess.Popen([str(self.zstdio), "c"],
                                        stdin=subprocess.PIPE,
                                        stdout=out,
                                        stderr=subprocess.PIPE)
        self.tmp_path = tmp
        self.records = 0
        self.files = 0

    def _finish(self) -> None:
        process = self.process
        tmp_path = self.tmp_path
        assert process is not None and process.stdin is not None
        assert tmp_path is not None and process.stderr is not None
        process.stdin.close()
        stderr = process.stderr.read().decode(errors="replace")
        process.stderr.close()
        code = process.wait()
        if code != 0:
            raise RuntimeError(f"zstd compressor failed: {stderr}")
        final = self.chunks_dir / f"chunk-{self.index:06d}.jsonl.zst"
        os.replace(tmp_path, final)
        self.completed.append({
            "name": f"chunks/{final.name}",
            "sha256": sha256_of(final),
            "size": final.stat().st_size,
            "records": self.records,
            "files": self.files,
        })
        self.index += 1
        self.process = None
        self.tmp_path = None

    def add(self, line: bytes, records: int) -> None:
        """Append one file's JSON line, cutting the chunk when large."""
        if self.process is None:
            self._start()
        assert self.process is not None and self.process.stdin is not None
        self.process.stdin.write(line)
        self.records += records
        self.files += 1
        assert self.tmp_path is not None
        if self.tmp_path.stat().st_size >= self.chunk_bytes:
            self._finish()

    def close(self) -> None:
        """Complete the trailing chunk, if one is open."""
        if self.process is not None:
            self._finish()


def dump_bytes_line(dump: Path, path: str, data: bytes,
                    digest: str) -> tuple[bytes, int]:
    """Run the bb0032c dump over captured bytes; return the JSON line.

    The bytes are written to a scratch file so the frozen hash and the
    dumped records describe exactly the same captured content, even
    while the live corpus churns.
    """
    scratch = Path(tempfile.mkdtemp()) / "capture.pdb"
    try:
        scratch.write_bytes(data)
        result = subprocess.run([str(dump), "--json", str(scratch)],
                                capture_output=True, timeout=600,
                                check=False)
    finally:
        shutil.rmtree(scratch.parent, ignore_errors=True)
    if result.returncode != 0:
        raise RuntimeError(f"dump failed for {path}: "
                           f"{result.stderr.decode(errors='replace')[:200]}")
    parsed = json.loads(result.stdout)
    wrapped = json.dumps({
        "path": path,
        "sha256": digest,
        "records": parsed["records"],
    }, separators=(",", ":")).encode()
    return wrapped + b"\n", len(parsed["records"])


def prepare_generate(
        options: argparse.Namespace
) -> Optional[tuple[Path, Path, Path, Path, list[str]]]:
    """Validate the generate options.

    Returns (input root, dump, output, zstdio, files).
    """
    input_root = options.input.resolve()
    dump = options.dump.resolve()
    zstdio = options.zstdio.resolve()
    if not input_root.is_dir():
        print(f"error: input directory not found: {input_root}",
              file=sys.stderr)
        return None
    if not dump.is_file():
        print(f"error: dump binary not found: {dump}", file=sys.stderr)
        return None
    if not zstdio.is_file():
        print(f"error: zstdio binary not found: {zstdio} "
              f"(build it with 'make zstdio')", file=sys.stderr)
        return None
    output = options.output
    if output.exists():
        print(f"error: output directory already exists: {output}",
              file=sys.stderr)
        return None
    existing_parent = output.parent
    while not existing_parent.exists():
        existing_parent = existing_parent.parent
    free = shutil.disk_usage(existing_parent).free
    print(f"free space on the output filesystem: {free // (1 << 30)} GiB",
          flush=True)
    print(f"enumerating {input_root} ...", flush=True)
    files = enumerate_files(input_root, options.all_extensions,
                            options.max_files)
    if not files:
        print("error: no input files found", file=sys.stderr)
        return None
    return input_root, dump, output, zstdio, files


def run_generate(options: argparse.Namespace) -> int:
    """Snapshot the bb0032c dump output over the corpus."""
    prepared = prepare_generate(options)
    if prepared is None:
        return 2
    input_root, dump, output, zstdio, files = prepared
    print(f"{len(files)} files; hashing and dumping with {options.jobs} "
          f"jobs", flush=True)

    meta = read_meta(dump)
    print(f"dump meta: {json.dumps(meta)}", flush=True)

    missing: list[str] = []

    def work(path: str) -> Optional[tuple[str, str, bytes, int]]:
        """Capture one file once; skip files that vanish mid-run.

        The captured flag distinguishes a vanished file from a captured
        one; the caller reconciles the missing list after the sweep.
        """
        try:
            data = (input_root / path).read_bytes()
        except OSError:
            return None
        digest = hashlib.sha256(data).hexdigest()
        line, records = dump_bytes_line(dump, path, data, digest)
        return path, digest, line, records

    output.mkdir(parents=True, exist_ok=True)
    writer = ChunkWriter(output / "chunks", zstdio,
                         int(options.chunk_gb * (1 << 30)))
    manifest_files: list[dict[str, object]] = []
    total_records = 0
    pool = concurrent.futures.ThreadPoolExecutor(max_workers=options.jobs)
    interrupted = False
    try:
        for result in pool.map(work, files):
            if result is None:
                continue
            path, digest, line, records = result
            chunk_index = writer.index
            writer.add(line, records)
            manifest_files.append({
                "path": path,
                "sha256": digest,
                "chunk": chunk_index,
                "records": records,
            })
            total_records += records
            if len(manifest_files) % 5000 == 0:
                print(f"  {len(manifest_files)}/{len(files)} files, "
                      f"{total_records} records, "
                      f"{writer.index} chunks", flush=True)
        pool.shutdown(wait=True)
        writer.close()
        captured = {str(entry["path"]) for entry in manifest_files}
        missing = [path for path in files if path not in captured]
    except KeyboardInterrupt:
        pool.shutdown(wait=False, cancel_futures=True)
        interrupted = True
        print("interrupted; partial chunks remain as .tmp files and no "
              "manifest was written; remove the output directory to "
              "restart", file=sys.stderr)
        pool.shutdown(wait=True)

    if not interrupted:
        write_manifest(output, meta, dump, manifest_files, total_records,
                       writer.completed, missing)
        total_bytes = sum(as_int(chunk["size"])
                          for chunk in writer.completed)
        print(f"snapshot: {len(writer.completed)} chunks, "
              f"{total_bytes // (1 << 20)} MiB compressed, "
              f"{total_records} records, {len(manifest_files)} files, "
              f"{len(missing)} missing at generation")
        print(f"manifest: {output / 'manifest.json'}")
    return 130 if interrupted else 0


def write_manifest(output: Path, meta: dict[str, str], dump: Path,
                   manifest_files: list[dict[str, object]],
                   total_records: int, completed: list[dict[str, object]],
                   missing: list[str]) -> None:
    """Write the snapshot manifest atomically (the completion record)."""
    manifest = {
        "snapshot_version": SNAPSHOT_VERSION,
        "generated_at": datetime.datetime.now(
            datetime.timezone.utc).isoformat(timespec="seconds"),
        "commit": meta.get("commit", "unknown"),
        "compiler": meta.get("compiler", "unknown"),
        "flags": meta.get("flags", "unknown"),
        "harness_version": repo_version() or "unknown",
        "dump_sha256": sha256_of(dump),
        "files": len(manifest_files),
        "records": total_records,
        "chunks": completed,
        "file_index": manifest_files,
        "missing_at_generation": missing,
    }
    tmp_manifest = output / ".manifest.json.tmp"
    tmp_manifest.write_text(json.dumps(manifest), encoding="utf-8")
    os.replace(tmp_manifest, output / "manifest.json")


# ---------------------------------------------------------------------------
# verify
# ---------------------------------------------------------------------------


def run_verify(options: argparse.Namespace) -> int:
    """Check every chunk against the snapshot manifest."""
    snapshot = options.snapshot.resolve()
    manifest_path = snapshot / "manifest.json"
    if not manifest_path.is_file():
        print(f"error: no manifest in {snapshot}", file=sys.stderr)
        return 2
    manifest = json.loads(manifest_path.read_text(encoding="utf-8"))
    errors: list[str] = []

    expected_chunks = {as_str(chunk["name"]): chunk
                       for chunk in manifest["chunks"]}
    actual_chunks = {
        f"chunks/{path.name}"
        for path in (snapshot / "chunks").iterdir()
        if path.is_file() and not path.name.startswith(".")
    }
    for stray in sorted(actual_chunks - set(expected_chunks)):
        errors.append(f"chunk not in manifest: {stray}")
    for missing in sorted(set(expected_chunks) - actual_chunks):
        errors.append(f"manifest chunk missing: {missing}")

    file_index = {as_str(entry["path"]): as_dict(entry)
                  for entry in as_list(manifest["file_index"])}
    seen_files: dict[str, str] = {}
    total_records = 0

    for name in sorted(expected_chunks):
        chunk_path = snapshot / name
        if not chunk_path.is_file():
            continue  # already reported as missing
        expected = expected_chunks[name]
        digest = sha256_of(chunk_path)
        if digest != as_str(expected["sha256"]):
            errors.append(f"{name}: sha256 mismatch")
            continue
        size = chunk_path.stat().st_size
        if size != as_int(expected["size"]):
            errors.append(f"{name}: size {size} != "
                          f"{as_int(expected['size'])}")
            continue
        chunk_number = int(name.split("-")[1].split(".")[0])
        records = 0
        count = 0
        try:
            for entry in stream_chunk(snapshot, options.zstdio.resolve(),
                                      name):
                count += 1
                records += len(as_list(entry["records"]))
                path = as_str(entry["path"])
                seen_files[path] = name
                index_entry = file_index.get(path)
                if index_entry is None:
                    errors.append(f"{name}: file not in manifest: {path}")
                    continue
                if as_int(index_entry["chunk"]) != chunk_number:
                    errors.append(f"{name}: file assigned to another "
                                  f"chunk: {path}")
                if as_str(index_entry["sha256"]) != as_str(
                        entry["sha256"]):
                    errors.append(f"{name}: file hash mismatch: {path}")
                if as_int(index_entry["records"]) != len(
                        as_list(entry["records"])):
                    errors.append(f"{name}: file record count mismatch: "
                                  f"{path}")
        except (RuntimeError, ValueError) as error:
            errors.append(f"{name}: {error}")
            continue
        if records != as_int(expected["records"]):
            errors.append(f"{name}: records {records} != "
                          f"{as_int(expected['records'])}")
        if count != as_int(expected["files"]):
            errors.append(f"{name}: files {count} != "
                          f"{as_int(expected['files'])}")
        total_records += records

    for missing in sorted(set(file_index) - set(seen_files)):
        errors.append(f"manifest file never seen in a chunk: {missing}")

    if errors:
        for error in errors[:50]:
            print(f"verify: {error}", file=sys.stderr)
        print(f"verify FAILED: {len(errors)} problem(s); "
              f"{total_records} records reconciled before failure")
        return 1
    print(f"verify OK: {len(manifest['chunks'])} chunks, "
          f"{total_records} records, {len(file_index)} files, "
          f"commit {manifest['commit']}")
    return 0


# ---------------------------------------------------------------------------
# compare
# ---------------------------------------------------------------------------


class Counters:
    """Per-class record counters for one compare run."""

    def __init__(self) -> None:
        self.files = 0
        self.records = 0
        self.parity = 0
        self.expected_zdatom = 0
        self.expected_v5_token = 0
        self.expected_modern_gain = 0
        self.unexpected = 0
        self.field_mismatches = 0
        self.coordinate_unknowns = 0
        self.hash_mismatches = 0

    def as_dict(self) -> dict[str, int]:
        return dict(self.__dict__)

    def total_mismatches(self) -> int:
        """Every finding that fails a compare run."""
        return (self.unexpected + self.field_mismatches
                + self.hash_mismatches)


def classify(snap_type: int, snap_text: str, modern_type: int,
             modern_text: str) -> str:
    """Mirror the live comparison's expected-difference groups."""
    if modern_text == snap_text and modern_type == snap_type:
        return "parity"
    atomish = (modern_text.startswith(ATOMISH_PREFIXES)
               or modern_type in ATOMISH_TYPES)
    if atomish and snap_text and modern_text:
        # the pinned pre-fb720ad clamp cut ZDATOM renders at 80 columns
        modern80 = modern_text[:80]
        if (modern80 == snap_text
                or (len(modern_text) == 82 and len(snap_text) == 80
                    and modern_text[:80] == snap_text)):
            return "expected_zdatom"
    if (snap_text.startswith("UNKNOWN:")
        and (modern_type in ATOMISH_TYPES
             or modern_text.startswith(ATOMISH_PREFIXES))):
        # the deliberate segid/element salvage gain: legacy junk vs the
        # modern standard record the line announces
        return "expected_modern_gain"
    if (snap_text.startswith(V5_TOKEN_PREFIXES)
            or modern_text.startswith(V5_TOKEN_PREFIXES)):
        if len(snap_text) > 60 or len(modern_text) > 60:
            return "expected_v5_token"
    return "unexpected"


def fields_equal(snap_fields: dict[str, object],
                 modern_fields: dict[str, object],
                 modern_type: int) -> tuple[bool, list[str]]:
    """Compare parsed fields with the live comparison's semantics."""
    if (modern_type == USER_GFX_LABEL_TYPE
            and not all(key in modern_fields for key in XYZ_KEYS)):
        # the live comparison skips the xyz checks entirely when the
        # modern parse left them unset (the v5 matching failure)
        left = {key: value for key, value in snap_fields.items()
                if key not in XYZ_KEYS}
        right = {key: value for key, value in modern_fields.items()
                 if key not in XYZ_KEYS}
        if set(left) != set(right):
            return False, [f"field keys differ: snapshot-only "
                           f"{sorted(set(left) - set(right))} modern-only "
                           f"{sorted(set(right) - set(left))}"]
        return (left == right,
                [key for key in left if left[key] != right.get(key)])
    if set(snap_fields) != set(modern_fields):
        return False, [f"field keys differ: snapshot-only "
                       f"{sorted(set(snap_fields) - set(modern_fields))} "
                       f"modern-only "
                       f"{sorted(set(modern_fields) - set(snap_fields))}"]
    differing = [key for key in snap_fields
                 if snap_fields[key] != modern_fields[key]]
    return not differing, differing


def compare_file(path: str, snap_entry: dict[str, object],
                 modern: dict[str, object], counters: Counters,
                 examples: dict[str, list[str]]) -> None:
    """Compare one file's modern records against its snapshot records."""
    snap_records = [as_dict(record)
                    for record in as_list(snap_entry["records"])]
    modern_records = [as_dict(record)
                      for record in as_list(modern["records"])]
    if len(snap_records) != len(modern_records):
        counters.unexpected += 1
        bucket = examples.setdefault("record-count", [])
        if len(bucket) < 5:
            bucket.append(f"{path}: snapshot {len(snap_records)} modern "
                          f"{len(modern_records)}")
        return
    for snap_record, modern_record in zip(snap_records, modern_records):
        counters.records += 1
        snap_type = as_int(snap_record["type"])
        modern_type = as_int(modern_record["type"])
        snap_text = as_str(snap_record["text"])
        modern_text = as_str(modern_record["text"])
        klass = classify(snap_type, snap_text, modern_type, modern_text)
        setattr(counters, klass, getattr(counters, klass) + 1)
        if klass == "unexpected":
            bucket = examples.setdefault("unexpected", [])
            if len(bucket) < 5:
                bucket.append(
                    f"{path}:{as_int(snap_record['line'])}: snapshot "
                    f"type {snap_type} '{snap_text[:60]}' modern type "
                    f"{modern_type} '{modern_text[:60]}'")
            # the live comparison reports parsed-field differences
            # independently of the text class; mirror that when both
            # sides carry fields
            if snap_record["fields"] and modern_record["fields"]:
                equal, differing = fields_equal(snap_record["fields"],
                                                modern_record["fields"],
                                                modern_type)
                if not equal:
                    counters.field_mismatches += 1
                    bucket = examples.setdefault("field", [])
                    if len(bucket) < 5:
                        bucket.append(f"{path}:"
                                      f"{as_int(snap_record['line'])}: "
                                      f"{', '.join(differing[:4])}")
            continue
        if klass in ("parity", "expected_zdatom"):
            equal, differing = fields_equal(snap_record["fields"],
                                            modern_record["fields"],
                                            modern_type)
            if not equal:
                counters.field_mismatches += 1
                counters.unexpected += 1
                bucket = examples.setdefault("field", [])
                if len(bucket) < 5:
                    bucket.append(f"{path}:"
                                  f"{as_int(snap_record['line'])}: "
                                  f"{', '.join(differing[:4])}")


def prepare_compare(
        options: argparse.Namespace
) -> Optional[tuple[Path, Path, dict[str, object]]]:
    """Validate the compare options.

    Returns (snapshot, binary, manifest).
    """
    snapshot = options.snapshot.resolve()
    binary = options.binary.resolve()
    if not (snapshot / "manifest.json").is_file():
        print(f"error: no manifest in {snapshot}", file=sys.stderr)
        return None
    if not binary.is_file():
        print(f"error: modern binary not found: {binary}", file=sys.stderr)
        return None
    manifest = json.loads((snapshot / "manifest.json").read_text(
        encoding="utf-8"))
    return snapshot, binary, manifest


def run_compare(options: argparse.Namespace) -> int:
    """Compare the modern parser against the snapshot."""
    prepared = prepare_compare(options)
    if prepared is None:
        return 2
    snapshot, binary, manifest = prepared
    counters = Counters()
    examples: dict[str, list[str]] = {}
    errors: list[str] = []

    def work(snap_entry: dict[str, object]) -> dict[str, object]:
        """Run the modern dump over one file's frozen lines."""
        raw_lines = [
            as_str(record["raw"])
            for record in as_list(snap_entry["records"])
        ]
        return modern_dump(binary, as_str(snap_entry["path"]), raw_lines)

    pool = concurrent.futures.ThreadPoolExecutor(max_workers=options.jobs)
    pending: list[tuple[dict[str, object], concurrent.futures.Future]] = []
    window = max(options.jobs * 4, 8)
    done = False
    chunk_names = [as_str(chunk["name"])
                   for chunk in as_list(manifest["chunks"])]

    def stream() -> Iterator[dict[str, object]]:
        for name in chunk_names:
            yield from stream_chunk(snapshot, options.zstdio.resolve(),
                                    name)

    entries = stream()

    def consume(pending: list[tuple[dict[str, object],
                                     concurrent.futures.Future]]
                ) -> tuple[bool, Optional[str]]:
        """Consume one ordered result; return (interrupted, error)."""
        snap_entry, future = pending.pop(0)
        modern = future.result()
        counters.files += 1
        path = as_str(snap_entry["path"])
        counters.coordinate_unknowns += as_int(
            modern.get("rejectedCoordinates", 0))
        compare_file(path, snap_entry, modern, counters, examples)
        return False, None

    try:
        while True:
            while len(pending) < window and not done:
                item = next(entries, None)
                if item is None:
                    done = True
                    break
                pending.append((item, pool.submit(work, item)))
            if not pending:
                break
            interrupted, error = consume(pending)
            if error is not None:
                errors.append(error)
            if interrupted:
                break
    except KeyboardInterrupt:
        pool.shutdown(wait=False, cancel_futures=True)
        print("interrupted", file=sys.stderr)
        return 130
    pool.shutdown(wait=True)

    ok = (counters.total_mismatches() == 0 and not errors
          and counters.files == as_int(manifest["files"]))
    summary = {
        "snapshot": str(snapshot),
        "manifest_records": as_int(manifest["records"]),
        "manifest_files": as_int(manifest["files"]),
        "files_seen": counters.files,
        **counters.as_dict(),
        "errors": errors[:50],
        "examples": examples,
    }
    print(json.dumps(summary, indent=2))
    print("compare PASSED" if ok else "compare FAILED")
    return 0 if ok else 1


def main(argv: Optional[list[str]] = None) -> int:
    """Dispatch the snapshot subcommands."""
    options = parse_args(argv)
    if options.command == "generate":
        return run_generate(options)
    if options.command == "verify":
        return run_verify(options)
    return run_compare(options)


if __name__ == "__main__":
    sys.exit(main())
