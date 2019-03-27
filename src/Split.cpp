/**
 * Copyright (c) 2019, Arjan van der Velde, Weng Lab
 * All rights reserved.
 *
 * Redistribution and use in source and binary forms, with or without
 * modification, are permitted provided that the following conditions are met:
 *
 * * Redistributions of source code must retain the above copyright notice, this
 *   list of conditions and the following disclaimer.
 *
 * * Redistributions in binary form must reproduce the above copyright notice,
 *   this list of conditions and the following disclaimer in the documentation
 *   and/or other materials provided with the distribution.
 *
 * THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND CONTRIBUTORS "AS IS"
 * AND ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT LIMITED TO, THE
 * IMPLIED WARRANTIES OF MERCHANTABILITY AND FITNESS FOR A PARTICULAR PURPOSE ARE
 * DISCLAIMED. IN NO EVENT SHALL THE COPYRIGHT HOLDER OR CONTRIBUTORS BE LIABLE
 * FOR ANY DIRECT, INDIRECT, INCIDENTAL, SPECIAL, EXEMPLARY, OR CONSEQUENTIAL
 * DAMAGES (INCLUDING, BUT NOT LIMITED TO, PROCUREMENT OF SUBSTITUTE GOODS OR
 * SERVICES; LOSS OF USE, DATA, OR PROFITS; OR BUSINESS INTERRUPTION) HOWEVER
 * CAUSED AND ON ANY THEORY OF LIABILITY, WHETHER IN CONTRACT, STRICT LIABILITY,
 * OR TORT (INCLUDING NEGLIGENCE OR OTHERWISE) ARISING IN ANY WAY OUT OF THE USE
 * OF THIS SOFTWARE, EVEN IF ADVISED OF THE POSSIBILITY OF SUCH DAMAGE.
 */

#include "Split.hpp"
#include "Utils.hpp"
#include <cassert>
#include <cmath>
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
  uint64_t x = i;
  for (size_t i = 0; i < sizeof(ret) - 1; ++i) {
    size_t v =
        static_cast<size_t>(std::pow(sizeof(alpha) - 1, sizeof(ret) - i - 2));
    if (x >= v) {
      ret[i] = alpha[x / v];
      x -= v * (x / v);
    }
  }
  assert(0 == x); // if not, x was too large to represent
  return std::string(ret);
}

void Split::split() {
  const ZDOCK z(zdockfn_);
  const int chunksize = (-1 == chunksize_ ? z.npredictions() : chunksize_);
  ZDOCK zz(z); // copy
  zz.predictions().clear();
  int i = 0, chunk = 0;
  for (const auto &v : z.predictions()) {
    zz.predictions().push_back(v);
    if (++i >= chunksize) {
      const std::string ofn = prefix_ + suffix(chunk++);
      std::ofstream f(ofn);
      if (f.is_open()) {
        f << zz << std::endl;
      } else {
        throw SplitException("Error opening output file '" + ofn + "'.");
      }
      zz.predictions().clear();
      i = 0;
    }
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
  while ((c = getopt(argc, argv, "hn:p:")) != -1) {
    switch (c) {
    case 'n':
      chunksize = std::stoi(optarg);
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

