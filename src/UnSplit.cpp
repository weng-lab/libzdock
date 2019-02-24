/**
 * Copyright 2019 Arjan van der Velde, vandervelde.ag [at] gmail
 * Licensed under the Apache License, Version 2.0 (the "License");
 * you may not use this file except in compliance with the License.
 * You may obtain a copy of the License at
 * http://www.apache.org/licenses/LICENSE-2.0
 *
 * Unless required by applicable law or agreed to in writing, software
 * distributed under the License is distributed on an "AS IS" BASIS,
 * WITHOUT WARRANTIES OR CONDITIONS OF ANY KIND, either express or implied.
 * See the License for the specific language governing permissions and
 * limitations under the License.
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

