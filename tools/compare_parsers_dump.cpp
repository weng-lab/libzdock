// Copyright (c) 2019-2026 Arjan van der Velde, Weng Lab
// SPDX-License-Identifier: BSD-2-Clause

// Per-line parser dump for the compare_parsers baseline mode. Reads a PDB
// file through operator>> exactly as production does and prints one line
// per record: "<type>\t<rendered text>". Works against any libpdb++
// revision, including the pre-rewrite bb0032c tree (build it there):
//
//   git worktree add /tmp/bb0032c bb0032c
//   make -C /tmp/bb0032c clean all
//   g++ -O2 -std=c++14 -I/tmp/bb0032c/contrib/eigen \
//       -I/tmp/bb0032c/include tools/compare_parsers_dump.cpp \
//       /tmp/bb0032c/lib/libzdock.a -o /tmp/baseline_dump

#include "pdb++.h"

#include <cstdio>
#include <cstring>
#include <string>

int main(int argc, char **argv) {
  if (argc < 2) {
    fprintf(stderr, "usage: %s <pdb-file>\n", argv[0]);
    return 2;
  }
  FILE *in = fopen(argv[1], "r");
  if (!in) {
    fprintf(stderr, "cannot open %s\n", argv[1]);
    return 2;
  }
  char buf[libpdb::PDB::BufLen * 4];
  while (fgets(buf, sizeof(buf), in)) {
    // strip the trailing newline (operator>> semantics)
    std::string line(buf);
    while (!line.empty() &&
           (line.back() == '\n' || line.back() == '\r'))
      line.pop_back();
    libpdb::PDB record(line.c_str());
    printf("%d\t%s\n", (int)record.type(), record.chars());
  }
  fclose(in);
  return 0;
}
