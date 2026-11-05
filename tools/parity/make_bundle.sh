#!/bin/sh
# Copyright (c) 2019-2026 Arjan van der Velde, Weng Lab
# SPDX-License-Identifier: BSD-2-Clause

# Package the libzdock source tree for a remote parity-harness run.
#
# The bundle carries the full source (including the pinned contrib
# submodule checkouts, which git archive would omit) plus the
# tools/parity scripts, so a server only needs a C++ toolchain, GNU
# make, and Python 3 to build bin/compare_parsers and run
# tools/parity/run_parity.py over a local corpus. Build output, test
# results, caches, and any results directories are excluded.
#
# Usage: tools/parity/make_bundle.sh [-o BUNDLE.tar.gz]
#
# The bundle records its provenance in libzdock/BUNDLE_INFO.txt, which
# the remote runner copies into its manifest. Copy the bundle to the
# remote server, then:
#
#   tar xzf libzdock-parity-bundle-*.tar.gz
#   cd libzdock && make compare-parsers
#   python3 tools/parity/run_parity.py --input /corpus/input \
#       --output /corpus/parity-results --compress

set -eu

script_dir=$(CDPATH= cd -- "$(dirname -- "$0")" && pwd)
repo_root=$(CDPATH= cd -- "$script_dir/../.." && pwd)
output=""
while getopts o: option; do
    case $option in
        o) output=$OPTARG ;;
        *) echo "usage: $0 [-o BUNDLE.tar.gz]" >&2; exit 2 ;;
    esac
done

cd "$repo_root"
describe=$(git describe --always --dirty 2>/dev/null || echo unknown)
stamp=$(date +%Y%m%d-%H%M%S)
[ -n "$output" ] || output="libzdock-parity-bundle-${describe}-${stamp}.tar.gz"

stage=$(mktemp -d "${TMPDIR:-/tmp}/libzdock-bundle.XXXXXX")
trap 'rm -rf "$stage"' EXIT

mkdir "$stage/libzdock"
tar -cf - \
    --exclude=./.git \
    --exclude=./bin \
    --exclude=./lib \
    --exclude=./build \
    --exclude=./python/.venv \
    --exclude=./python/.uv-cache \
    --exclude=./python/.pylint.d \
    --exclude='./parity-results-*' \
    --exclude='./*.log' \
    --exclude=./compile_commands.json \
    --exclude=./.clangd \
    --exclude=./opencode.json \
    . | tar -xf - -C "$stage/libzdock"

{
    echo "repository: libzdock"
    echo "commit: $describe"
    echo "packed: $(date -u +%Y-%m-%dT%H:%M:%SZ)"
    echo "host: $(hostname 2>/dev/null || echo unknown)"
} > "$stage/libzdock/BUNDLE_INFO.txt"

tar -czf "$output" -C "$stage" libzdock
sha256sum "$output"
