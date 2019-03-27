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

#include "Pruning.hpp"
#include "PDB.hpp"
#include "Utils.hpp"
#include <cstdio>
#include <iomanip>
#include <unistd.h>

namespace zdock {

char spinner() {
  const char s[4] = {'|', '/', '-', '\\'};
  static size_t i = 0;
  return s[(++i) % 4];
}

Pruning::Pruning(const std::string &zdockoutput, const double cutoff,
                 const std::string &structurefn, const bool getclusters)
    : zdock_(zdockoutput), cutoff_(cutoff), txl_(zdockoutput),
      txm_(zdockoutput), getclusters_(getclusters) {

  // ligand file name
  if (!zdock_.ismzdock()) {
    // ZDOCK has "ligand"
    if ("" == structurefn) {
      strucfn_ = Utils::copath(zdockoutput, zdock_.ligand().filename);
    } else {
      strucfn_ = structurefn;
    }
  } else {
    // M-ZDOCK only has "receptor"
    if ("" == structurefn) {
      strucfn_ = Utils::copath(zdockoutput, zdock_.receptor().filename);
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

  // print some info on stderr
  if (ismzdock) {
    std::cerr << "Pruning for M-ZDOCK; cutoff: " << std::fixed
              << std::setprecision(2) << cutoff_ << std::endl;
  } else {
    std::cerr << "Pruning for ZDOCK; cutoff: " << std::fixed
              << std::setprecision(2) << cutoff_ << std::endl;
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
  int assigned = 0;
  char buf[100];
  const size_t interval = 100;
  for (size_t i = 0; i < n; ++i) {
    if (!(i % interval)) {
      std::snprintf(buf, sizeof(buf),
                    "\r%c prediction: %ld, clusters: %d (%.2f%%)", spinner(), i,
                    clusters, 100.0 * assigned / n);
      std::cerr << buf << std::flush;
    }
    if (!l.at(i)) {
      l[i] = clusters + 1;
      assigned++;
      if (getclusters_) {
        // create prediction object w/ cluster number as score
        auto tmppred = v[i];
        tmppred.score = static_cast<double>(l[i]);
        preds.push_back(tmppred);
      } else {
        preds.push_back(v[i]);
      }
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
            assigned++;
            if (getclusters_) {
              // create prediction object w/ cluster number as score
              auto tmppred = v[j];
              tmppred.score = static_cast<double>(l[j]);
              preds.push_back(tmppred);
            }
          }
        }
      }
      clusters++;
    }
  }
  std::snprintf(buf, sizeof(buf), "\r%c prediction: %ld, clusters: %d (%.2f%%)",
                '-', n, clusters, 100.0);
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
  std::cerr
      << "usage: " << cmd << " [options] <zdock output>\n\n"
      << "  -c <double>     cutoff RMSD (defaults to 16.00)\n"
      << "  -C              return all prediction, but with score replaced by\n"
      << "                  cluster number.\n"
      << "  -l <filename>   structure PDB filename; defaults to ligand in "
         "ZDOCK\n"
      << std::endl;
}

} // namespace zdock

int main(int argc, char *argv[]) {
  std::string zdockfn, ligfn;
  double cutoff = 16.00;
  bool getclusters = false;
  int c;
  while ((c = getopt(argc, argv, "hc:l:C")) != -1) {
    switch (c) {
    case 'c':
      cutoff = std::stod(optarg);
      break;
    case 'l':
      ligfn = optarg;
      break;
    case 'C':
      getclusters = true;
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
    zdock::Pruning p(zdockfn, cutoff, ligfn, getclusters);
    p.prune();
    std::cout << p.zdock() << std::endl;
    std::cerr << "duration: " << zdock::Utils::toc(t1) << " sec" << std::endl;
  } catch (const zdock::Exception &e) {
    // something went wrong
    zdock::usage(argv[0], e.what());
    return 1;
  }
}

