<!--
Copyright (c) 2019-2026 Arjan van der Velde, Weng Lab
SPDX-License-Identifier: BSD-2-Clause
-->

# Remote parity-harness runs

How to run the legacy-vs-modern parser comparison (`tools/compare_parsers.cpp`,
the `compare-parsers` make target) over corpora far too large for one local
invocation — hundreds of thousands of PDB files — on another server, and how
to bring the results back for analysis.

## One-time setup

1. Build the bundle from this repository:

   ```console
   tools/parity/make_bundle.sh
   ```

   This writes `libzdock-parity-bundle-<commit>-<date>.tar.gz` (~25 MB): the
   full source tree including the pinned `contrib/` submodule checkouts and
   the `tools/parity/` scripts, with build output and caches excluded.

2. Copy it to the server that holds the corpus and unpack:

   ```console
   scp libzdock-parity-bundle-*.tar.gz server:~/
   ssh server
   tar xzf libzdock-parity-bundle-*.tar.gz && cd libzdock
   ```

## Remote run

The server needs only a C++ toolchain, GNU make, and Python 3 (stdlib only).

```console
make compare-parsers
python3 tools/parity/run_parity.py --input /path/to/input \
    --output /path/to/parity-results --compress
```

What the runner does: enumerates the input tree (`*.pdb` by default) in
sorted order, writes shard file lists (`--shard-size`, default 2000 files),
runs one harness process per shard across `--jobs` workers, and stores each
shard's summary, mismatch report, and exit status under the output
directory, plus a `manifest.json` (host, tool commit, binary SHA-256) and an
aggregated `summary.json`. An interrupted run resumes with `--resume` on the
same output directory (it reuses the original shard lists verbatim, so a
different `--shard-size` on resume cannot split already-recorded shards).
`--baseline <exe>` runs the optional bb0032c
comparison; `--quiet` drops per-line reports (counters only).

Mismatch reports are plain text with relative paths and can be large
(~3 MB per million lines on the local corpus); `--compress` archives the
output directory to `<output>.tar.gz` for the trip back.

## Copy back and analyze

```console
scp server:~/parity-results-*.tar.gz .
tar xzf parity-results-*.tar.gz
python3 tools/parity/merge_parity.py parity-results-<stamp>
```

The merge deduplicates every per-line report into normalized signatures
(record kind plus first-divergence column, digits masked), so millions of
reports collapse into the handful of underlying bugs; it writes
`merged.json` and prints the top signatures with one example location each.

After fixing parser bugs and re-running the corpus remotely, compare runs:

```console
python3 tools/parity/merge_parity.py --compare old/merged.json new/merged.json
```

This prints counter deltas and per-signature improvements/new introductions;
a nonzero exit means new or worsened signatures appeared.

## Output layout

| File | Contents |
| ---- | -------- |
| `manifest.json` | host, tool commit, binary SHA-256, options, file count |
| `summary.json` | aggregated counters, shard statuses, error list |
| `shards/shard-NNNNN.list` | input paths (relative) given to that shard |
| `shards/shard-NNNNN.summary.txt` | the harness summary line |
| `shards/shard-NNNNN.mismatches.log` | per-line mismatch reports (stderr) |
| `shards/shard-NNNNN.status.json` | exit code, seconds, parsed counters |
| `merged.json` | written locally by `merge_parity.py` |

Exit codes of the runner: 0 = all shards ran, 1 = shard errors, 2 = usage or
setup error, 130 = interrupted (resumable). Results stay out of the
repository, matching the rewrite plan's rule that mismatch reports are not
committed.
