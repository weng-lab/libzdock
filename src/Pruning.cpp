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

#include "Pruning.hpp"
#include <string>
#include <unistd.h>

namespace zdock {

Pruning::Pruning(const std::string &zdockfn, const std::string ligandfn,
                 const double cutoff)
    : zdockfn_(zdockfn), ligfn_(ligandfn), cutoff_(cutoff) {}

void Pruning::doPrune() {
  ZdockPruning p(zdockfn_, cutoff_, ligfn_);
  p.prune();
  std::cout << p.zdock() << std::endl;
}

void usage(const std::string &cmd, const std::string &err = "") {
  // print error if any
  if ("" != err) {
    std::cerr << "Error: " << err << std::endl << std::endl;
  }
  // print usage
  std::cerr << "usage: " << cmd << " [options] <zdock output>\n\n"
            << "  -c <double>     cutoff RMSD (defaults to 16.00)\n"
            << "  -l <filename>   ligand PDB filename; defaults to ligand in "
               "ZDOCK output\n"
            << std::endl;
}

} // namespace zdock

int main(int argc, char *argv[]) {
  std::string zdockfn, ligfn;
  double cutoff = 16.00;
  int c;
  while ((c = getopt(argc, argv, "h:c:l:")) != -1) {
    switch (c) {
    case 'c':
      cutoff = std::stod(optarg);
      break;
    case 'l':
      zdockfn = optarg;
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
    zdock::Pruning p(zdockfn, ligfn, cutoff);
    p.doPrune();
  } catch (const zdock::Exception &e) {
    // something went wrong
    throw zdock::ZdockPruningException(e.what());
  }
}

