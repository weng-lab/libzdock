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

#include "FilterConstraints.hpp"
#include "Utils.hpp"
#include <cmath>
#include <iostream>
#include <limits>
#include <unistd.h>

namespace e = Eigen;
namespace p = libpdb;

namespace zdock {

FilterConstraints::FilterConstraints(const std::string &zdockoutput,
                                     const std::string &constraints,
                                     const std::string &receptorpdb,
                                     const std::string &ligandpdb)
    : zdock_(zdockoutput), txl_(zdockoutput), txm_(zdockoutput),
      confn_(constraints) {
  // receptor file name
  if ("" == receptorpdb) {
    recfn_ = Utils::copath(zdockoutput, zdock_.receptor().filename);
  } else {
    recfn_ = receptorpdb;
  }
  // ligand file name
  if (!zdock_.ismzdock()) {
    if ("" == ligandpdb) {
      ligfn_ = Utils::copath(zdockoutput, zdock_.ligand().filename);
    } else {
      ligfn_ = ligandpdb;
    }
  }
}

// m-zdock
void FilterConstraints::filterMZDOCKConstraints_() {
  // load constraints file
  Constraints ccc(confn_);

  // load full PDB files
  PDB structure(recfn_);

  // grab atoms for valid constraints
  PDB s0, s1, s2;
  std::vector<double> mindist;
  std::vector<double> maxdist;
  if (!ccc.constraints().size()) {
    throw ConstraintException("No constraints specified");
  }

  // make structures w/ selected atoms
  for (const auto &x : ccc.constraints()) {
    try {
      // these throw exceptions for bad constraints; we want to append
      // both_ l and r or _none_.
      const PDB::Record r1 = structure[x.recCoord];
      const PDB::Record r2 = structure[x.ligCoord]; // two copies
      const PDB::Record r3 = structure[x.ligCoord]; // two copies
      s0.append(r2);                                // one side
      s1.append(r1);                                // original
      s2.append(r3);                                // other side
      if (Constraint::MAX == x.constraintType) {
        maxdist.push_back(x.distance);
        mindist.push_back(-1.0); // negative distance
      } else {
        maxdist.push_back(std::numeric_limits<double>::max()); // big number
        mindist.push_back(x.distance);
      }
    } catch (const AtomNotFoundException e) {
      throw ConstraintException("Constraint Error: " + std::string(e.what()));
    }
  }

  // assess all poses
  const size_t n = maxdist.size();
  const auto v = zdock_.predictions(); // our copy
  auto &preds = zdock_.predictions();  // our ref
  preds.clear();
  e::Matrix<double, 2, e::Dynamic> m;
  m.resize(2, ccc.constraints().size());
  for (const Prediction &p : v) {
    // compare poses p0 and p2 to pose p1 ("middle" structure)
    FilterConstraints::Matrix p0 = txm_.txMultimer(s0.matrix(), p, 0);
    FilterConstraints::Matrix p1 = txm_.txMultimer(s1.matrix(), p, 1);
    FilterConstraints::Matrix p2 = txm_.txMultimer(s2.matrix(), p, 2);
    // row 0: dist p1 -> p0; row 1: dist p1 -> p2
    m << (p1 - p0.matrix()).colwise().squaredNorm().array().sqrt().matrix(),
        (p1 - p2.matrix()).colwise().squaredNorm().array().sqrt().matrix();
    bool accepted = true;
    for (size_t i = 0; i < n; ++i) {
      // actual filtering
      if (std::min<double>(m(0, i), m(1, i)) > maxdist[i] ||
          std::min<double>(m(0, i), m(1, i)) < mindist[i]) {
        accepted = false;
        break;
      }
    }
    if (accepted) {
      preds.push_back(p);
    }
  }
}

// zdock
void FilterConstraints::filterZDOCKConstraints_() {
  // load constraints file
  Constraints ccc(confn_);

  // load full PDB files
  PDB receptor(recfn_);
  PDB ligand(ligfn_);

  // grab atoms for valid constraints
  PDB ligatoms, recatoms;
  std::vector<double> mindist;
  std::vector<double> maxdist;
  if (!ccc.constraints().size()) {
    throw ConstraintException("No constraints specified");
  }

  // make structures w/ selected atoms
  for (const auto &x : ccc.constraints()) {
    try {
      // these throw exceptions for bad constraints; we want to append
      // both_ l and r or _none_.
      const PDB::Record r = receptor[x.recCoord];
      const PDB::Record l = ligand[x.ligCoord];
      recatoms.append(r);
      ligatoms.append(l);
      if (Constraint::MAX == x.constraintType) {
        maxdist.push_back(x.distance);
        mindist.push_back(-1.0); // negative distance
      } else {
        maxdist.push_back(std::numeric_limits<double>::max()); // big number
        mindist.push_back(x.distance);
      }
    } catch (const AtomNotFoundException e) {
      throw ConstraintException("Constraint Error: " + std::string(e.what()));
    }
  }

  // assess all poses
  const size_t n = maxdist.size();
  const auto v = zdock_.predictions(); // our copy
  auto &preds = zdock_.predictions();  // our ref
  preds.clear();
  for (const Prediction &p : v) {
    FilterConstraints::Matrix pose = txl_.txLigand(ligatoms.matrix(), p);
    e::Matrix<double, 1, e::Dynamic> m =
        (pose - recatoms.matrix()).colwise().squaredNorm().array().sqrt();
    bool accepted = true;
    for (size_t i = 0; i < n; ++i) {
      // actual filtering
      if (m(0, i) > maxdist[i] || m(0, i) < mindist[i]) {
        accepted = false;
        break;
      }
    }
    if (accepted) {
      preds.push_back(p);
    }
  }
}

void FilterConstraints::filter() {
  if (zdock_.ismzdock()) {
    std::cerr << "M-ZDOCK Constraints Filter" << std::endl;
    filterMZDOCKConstraints_();
  } else {
    std::cerr << "ZDOCK Constraints Filter" << std::endl;
    filterZDOCKConstraints_();
  }
}

void usage(const std::string &cmd, const std::string &err = "") {
  // print error if any
  if ("" != err) {
    std::cerr << "Error: " << err << std::endl << std::endl;
  }
  // print usage
  std::cerr << "usage: " << cmd
            << " [options] <zdock output> <constraints file>\n\n"
            << "  -r <filename>   receptor PDB filename; defaults to receptor "
               "in (M-)ZDOCK output\n"
            << "  -l <filename>   ligand PDB filename; defaults to ligand in "
               "ZDOCK output\n"
            << std::endl;
}

} // namespace zdock

int main(int argc, char *argv[]) {
  std::string zdockfn, confn, recfn, ligfn;
  int c;
  while ((c = getopt(argc, argv, "h:r:l:")) != -1) {
    switch (c) {
    case 'r':
      recfn = optarg;
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
  if (argc > optind + 1) {
    zdockfn = argv[optind];   // zdock file
    confn = argv[optind + 1]; // constraints file
  } else {
    if (argc > optind) {
      zdock::usage(argv[0], "No constraints file specified.");
    } else {
      zdock::usage(argv[0], "No ZDOCK output file specified.");
    }
    return 1;
  }
  try {
    const auto t1 = zdock::Utils::tic();
    zdock::FilterConstraints p(zdockfn, confn, recfn, ligfn);
    p.filter();
    std::cout << p.zdock() << std::endl;
    std::cerr << "duration: " << zdock::Utils::toc(t1) << " sec" << std::endl;
  } catch (const zdock::Exception &e) {
    // something went wrong
    zdock::usage(argv[0], e.what());
    return 1;
  }
}

