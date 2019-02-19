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

#include "ZdockPruning.hpp"
#include "Utils.hpp"
#include <cmath>
#include <iostream>
#include <limits>
//#include <map>

namespace e = Eigen;
namespace p = libpdb;

namespace zdock {

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

  // read pdb file
  PDB ligpdb(ligfn_);
  const double ligsize = ligpdb.matrix().cols();

  // pre-compute all poses
  std::vector<ZdockPruning::Matrix> poses;
  for (size_t i = 0; i < n; ++i) {
    poses.push_back(txl_.txLigand(ligpdb.matrix(), v[i]));
  }

  // find clusters
  zdock_.predictions().clear();
  for (size_t i = 0; i < n; ++i) {
    if (!(i % 100)) {
      std::cerr << i << '\t' << clusters << std::endl;
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

} // namespace zdock
