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

#include "UnSplit.hpp"
#include "Utils.hpp"
#include <fstream>
#include <iostream>
#include <unistd.h>

namespace zdock {

UnSplit::UnSplit(const std::vector<std::string> &files) : files_(files) {}

void UnSplit::unsplit() {
  ZDOCK z(files_[0]);
  int i = 0;
  for (const auto &x : files_) {
    if (i++ > 0) {
      ZDOCK zz(x);
      for (const auto &p : zz.predictions()) {
        z.predictions().push_back(p);
      }
    }
  }
  std::cout << z << std::endl;
}

void usage(const std::string &cmd, const std::string &err = "") {
  // print error if any
  if ("" != err) {
    std::cerr << "Error: " << err << std::endl << std::endl;
  }
  // print usage
  std::cerr << "usage: " << cmd << " <zdock output> [file] [...]\n\n"
            << std::endl;
}

} // namespace zdock

int main(int argc, char *argv[]) {
  std::vector<std::string> files;
  int c;
  while ((c = getopt(argc, argv, "h")) != -1) {
    switch (c) {
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
  if (argc <= optind) {
    zdock::usage(argv[0], "No ZDOCK output file(s) specified.");
    return 1;
  }
  while (argc > optind) {
    files.push_back(argv[optind++]); // zdock file
  }
  try {
    const auto t1 = zdock::Utils::tic();
    zdock::UnSplit p(files);
    p.unsplit();
    std::cerr << "duration: " << zdock::Utils::toc(t1) << " sec" << std::endl;
  } catch (const zdock::Exception &e) {
    // something went wrong
    zdock::usage(argv[0], e.what());
    return 1;
  }
}

