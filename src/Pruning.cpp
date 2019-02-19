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
#include "PDB.hpp"
#include "Utils.hpp"
#include <cstdio>
#include <unistd.h>

namespace zdock {

char spinner() {
  const char s[4] = {'|', '/', '-', '\\'};
  static size_t i = 0;
  return s[(++i) % 4];
}

ZdockPruning::ZdockPruning(const std::string &zdockoutput, const double cutoff,
                           const std::string &ligandpdb)
    : zdock_(zdockoutput), cutoff_(cutoff), txl_(zdockoutput) {
  // ligand file name
  if (!zdock_.ismzdock()) {
    if ("" == ligandpdb) {
      ligfn_ = zdock_.ligand().filename;
      if ('/' != ligfn_[0]) { // relative
        ligfn_ = Utils::copath(zdockoutput, ligfn_);
      }
    } else {
      ligfn_ = ligandpdb;
    }
  }
}

void ZdockPruning::prune() {
  const auto v = zdock_.predictions(); // our copy
  const auto n = zdock_.npredictions();
  auto &preds = zdock_.predictions();              // our ref
  double min = std::numeric_limits<double>::max(); // big number
  std::vector<int> l(n, 0);
  int clusters = 0;

  // ZDOCK only for now
  if (zdock_.ismzdock()) {
    throw ZdockPruningException("M-ZDOCK not supported");
  }

  // read pdb file (CA only!)
  PDB ligpdb(ligfn_, [](const auto &r) {
    return Utils::trim_copy(r.atom.name) == "CA";
  });
  const double ligsize = ligpdb.matrix().cols();

  // pre-compute all poses
  std::vector<ZdockPruning::Matrix> poses;
  for (size_t i = 0; i < n; ++i) {
    poses.push_back(txl_.txLigand(ligpdb.matrix(), v[i]));
  }

  // find clusters
  zdock_.predictions().clear();
  char buf[100];
  const size_t interval = 100;
  for (size_t i = 0; i < n; ++i) {
    if (!(i % interval)) {
      std::snprintf(buf, sizeof(buf), "\r%c prediction: %ld, clusters: %d", spinner(), i, clusters);
      std::cerr << buf << std::flush;
    }
    if (!l.at(i)) {
      l[i] = clusters + 1;
      preds.push_back(v[i]);
      for (size_t j = i + 1; j < n; ++j) {
        if (!l.at(j)) {
          const double rmsd =
              std::sqrt((poses.at(i) - poses.at(j)).squaredNorm() / ligsize);
          min = std::min(min, rmsd);
          if (rmsd < cutoff_) {
            l[j] = clusters + 1;
          }
        }
      }
      clusters++;
    }
  }

  // copy out results
  clusters_ = l;
  ligsize_ = ligsize;
  nclusters_ = clusters;
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
    zdock::ZdockPruning p(zdockfn, cutoff, ligfn);
    p.prune();
    std::cout << p.zdock() << std::endl;
  } catch (const zdock::Exception &e) {
    // something went wrong
    throw zdock::ZdockPruningException(e.what());
  }
}

