// Copyright (c) 2019-2026 Arjan van der Velde, Weng Lab
// SPDX-License-Identifier: BSD-2-Clause

#include "Split.hpp"
#include "Utils.hpp"
#include <fstream>
#include <iostream>
#include <unistd.h>

namespace zdock {

Split::Split(const std::string &zdockfn, const int chunksize,
             const std::string &prefix)
    : zdockfn_(zdockfn), chunksize_(chunksize), prefix_(prefix) {}

std::string Split::suffix(const uint64_t i) const {
  char ret[] = "aaaa";
  const char alpha[] = "abcdefghijklmnopqrstuvwxyz";
  constexpr uint64_t base = sizeof(alpha) - 1;
  constexpr uint64_t capacity = base * base * base * base;
  if (i >= capacity) {
    throw SplitException("Too many output chunks");
  }
  uint64_t x = i;
  for (size_t i = 0; i < sizeof(ret) - 1; ++i) {
    size_t v = 1;
    for (size_t j = i + 1; j < sizeof(ret) - 1; ++j) {
      v *= base;
    }
    if (x >= v) {
      ret[i] = alpha[x / v];
      x -= v * (x / v);
    }
  }
  return std::string(ret);
}

void Split::split() {
  const ZDOCK z(zdockfn_);
  const int chunksize = (-1 == chunksize_ ? z.npredictions() : chunksize_);
  if (chunksize <= 0) {
    throw SplitException("Chunk size must be greater than zero");
  }
  ZDOCK zz(z); // copy
  zz.predictions().clear();
  int i = 0, chunk = 0;
  const auto writeChunk = [&]() {
    const std::string ofn = prefix_ + suffix(chunk++);
    std::ofstream f(ofn);
    if (!f.is_open()) {
      throw SplitException("Error opening output file '" + ofn + "'.");
    }
    f << zz << std::endl;
    zz.predictions().clear();
  };
  for (const auto &v : z.predictions()) {
    zz.predictions().push_back(v);
    if (++i >= chunksize) {
      writeChunk();
      i = 0;
    }
  }
  if (!zz.predictions().empty()) {
    writeChunk();
  }
}

void usage(const std::string &cmd, const std::string &err = "") {
  // print error if any
  if ("" != err) {
    std::cerr << "Error: " << err << std::endl << std::endl;
  }
  // print usage
  std::cerr
      << "usage: " << cmd << " [options] <zdock output>\n\n"
      << "  -n <integer>    chunk size (defaults to input size)\n"
      << "  -p <string>     output filename prefix (defaults to \"zdsplit.\")\n"
      << std::endl;
}

} // namespace zdock

int main(int argc, char *argv[]) {
  std::string zdockfn;
  std::string prefix = "zdsplit.";
  int chunksize = -1;
  int c;
  try {
    while ((c = getopt(argc, argv, "hn:p:")) != -1) {
      switch (c) {
    case 'n':
      chunksize = zdock::Utils::parseInt(optarg);
      break;
    case 'p':
      prefix = std::string(optarg);
      break;
    case 'h': // usage
      zdock::usage(argv[0]);
      return 0;
    case '?':
      zdock::usage(argv[0]);
      return 1;
    default:
      return 1;
      }
    }
  } catch (const zdock::Exception &e) {
    zdock::usage(argv[0], e.what());
    return 1;
  }
  if (argc > optind) {
    zdockfn = argv[optind]; // zdock file
  } else {
    zdock::usage(argv[0], "No ZDOCK output file specified.");
    return 1;
  }
  try {
    const auto t1 = zdock::Utils::tic();
    zdock::Split p(zdockfn, chunksize, prefix);
    p.split();
    std::cerr << "duration: " << zdock::Utils::toc(t1) << " sec" << std::endl;
  } catch (const zdock::Exception &e) {
    // something went wrong
    zdock::usage(argv[0], e.what());
    return 1;
  }
}
