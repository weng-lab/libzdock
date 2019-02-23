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

Pruning::Pruning(const std::string &zdockoutput, const double cutoff,
                 const std::string &structurefn)
    : zdock_(zdockoutput), cutoff_(cutoff), txl_(zdockoutput),
      txm_(zdockoutput) {

  // ligand file name
  if (!zdock_.ismzdock()) {
    // ZDOCK has "ligand"
    if ("" == structurefn) {
      strucfn_ = zdock_.ligand().filename;
      if ('/' != strucfn_[0]) { // relative
        strucfn_ = Utils::copath(zdockoutput, strucfn_);
      }
    } else {
      strucfn_ = structurefn;
    }
  } else {
    // M-ZDOCK only has "receptor"
    if ("" == structurefn) {
      strucfn_ = zdock_.receptor().filename;
      if ('/' != strucfn_[0]) { // relative
        strucfn_ = Utils::copath(zdockoutput, strucfn_);
      }
    } else {
      strucfn_ = structurefn;
    }
  }
}

void Pruning::prune() {
  const auto v = zdock_.predictions(); // our copy
  const auto n = zdock_.npredictions();
  auto &preds = zdock_.predictions();              // our ref
  double min = std::numeric_limits<double>::max(); // big number
  const bool ismzdock = zdock_.ismzdock();

  // ZDOCK only for now
  if (ismzdock) {
    std::cerr << "Pruning for M-ZDOCK" << std::endl;
  } else {
    std::cerr << "Pruning for ZDOCK" << std::endl;
  }

  // read pdb file (CA only!)
  PDB pdb(strucfn_,
          [](const auto &r) { return Utils::trim_copy(r.atom.name) == "CA"; });
  const double strucsize = pdb.matrix().cols();

  // pre-compute all poses
  std::vector<Pruning::Matrix> poses0, poses1;
  if (ismzdock) {
    for (size_t i = 0; i < n; ++i) {
      poses0.push_back(
          txm_.txMultimer(pdb.matrix(), v[i], 0)); // "left side" of "receptor"
      poses1.push_back(
          txm_.txMultimer(pdb.matrix(), v[i], 2)); // "right side" of "receptor"
    }
  } else {
    for (size_t i = 0; i < n; ++i) {
      poses0.push_back(txl_.txLigand(pdb.matrix(), v[i]));
    }
  }

  // find clusters
  zdock_.predictions().clear();
  std::vector<int> l(n, 0);
  int clusters = 0;
  char buf[100];
  const size_t interval = 100;
  for (size_t i = 0; i < n; ++i) {
    if (!(i % interval)) {
      std::snprintf(buf, sizeof(buf), "\r%c prediction: %ld, clusters: %d",
                    spinner(), i, clusters);
      std::cerr << buf << std::flush;
    }
    if (!l.at(i)) {
      l[i] = clusters + 1;
      preds.push_back(v[i]);
      for (size_t j = i + 1; j < n; ++j) {
        if (!l.at(j)) {
          double rmsd;
          if (ismzdock) {
            rmsd = std::min<double>(
                std::sqrt((poses0.at(i) - poses0.at(j)).squaredNorm() /
                          strucsize),
                std::sqrt((poses0.at(i) - poses1.at(j)).squaredNorm() /
                          strucsize));
          } else {
            rmsd = std::sqrt((poses0.at(i) - poses0.at(j)).squaredNorm() /
                             strucsize);
          }
          min = std::min(min, rmsd); // just for stats
          if (rmsd < cutoff_) {
            l[j] = clusters + 1;
          }
        }
      }
      clusters++;
    }
  }
  std::snprintf(buf, sizeof(buf), "\r%c prediction: %ld, clusters: %d", '-', n,
                clusters);
  std::cerr << buf << std::endl;

  // copy out results
  clusters_ = l;
  strucsize_ = strucsize;
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
      ligfn = optarg;
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
    zdock::Pruning p(zdockfn, cutoff, ligfn);
    p.prune();
    std::cout << p.zdock() << std::endl;
    std::cerr << "duration: " << zdock::Utils::toc(t1) << " sec" << std::endl;
  } catch (const zdock::Exception &e) {
    // something went wrong
    zdock::usage(argv[0], e.what());
    return 1;
  }
}

