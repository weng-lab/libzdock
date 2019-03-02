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

#include "Centroids.hpp"
#include "TransformLigand.cpp"
#include "Utils.hpp"
#include "ZDOCK.hpp"
#include <cmath>
#include <unistd.h>

namespace p = libpdb;
namespace e = Eigen;

namespace zdock {
Centroids::Centroids(const std::string &zdockoutput, const std::string &ligand,
                     const size_t n)
    : zdockfn_(zdockoutput), ligandfn_(ligand), n_(n) {}

void Centroids::doCentroids() {
  std::string ligfn; // ligand file name
  ZDOCK z(zdockfn_); // zdock output parser

  // check n within range
  if (n_ < 1) {
    throw CentroidsException("Invalid prediction; valid range 1 - " +
                             std::to_string(z.npredictions()));
  }
  // limit to npreds
  n_ = std::min<size_t>(n_, z.npredictions());

  // figure out ligand file name
  if ("" == ligandfn_) {
    ligfn = Utils::copath(zdockfn_, z.ligand().filename);
  } else {
    ligfn = ligandfn_;
  }

  // load ligand
  PDB lig(ligfn), out;
  p::PDB x(p::PDB::HETATM);
  const TransformLigand txl(z);
  const e::Vector3d v = lig.centroid();
  for (size_t i = 0; i < n_; ++i) {
    const zdock::Prediction &pred = z.predictions()[i];
    e::Vector3d pose = txl.txLigand(v, pred);
    x.atom = templateAtom_;
    x.atom.serialNum = static_cast<int>(i) + 1;
    x.atom.residue.seqNum = static_cast<int>(i) + 1;
    x.atom.xyz[0] = pose(0);
    x.atom.xyz[1] = pose(1);
    x.atom.xyz[2] = pose(2);
    out.append(x);
  }
  for (const auto &x : out.atoms()) {
    std::cout << *x << '\n';
  }
}

void usage(const std::string &cmd, const std::string &err = "") {
  // print error if any
  if ("" != err) {
    std::cerr << "Error: " << err << std::endl << std::endl;
  }
  // print usage
  std::cerr << "usage: " << cmd << " [options] <zdock output>\n\n"
            << "  -n <integer>    number of centroids to generate (top-n) "
               "(defaults to 1; the top prediction)\n"
            << "  -l <filename>   ligand PDB filename; defaults to receptor in "
               "ZDOCK output\n"
            << std::endl;
}

} // namespace zdock

int main(int argc, char *argv[]) {
  std::string zdockfn, ligfn;
  int n = 1;
  int c;
  while ((c = getopt(argc, argv, "hn:l:")) != -1) {
    switch (c) {
    case 'n':
      n = std::stoi(optarg);
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
    zdock::Centroids ct(zdockfn, ligfn, n);
    ct.doCentroids();
    std::cerr << "duration: " << zdock::Utils::toc(t1) << " sec" << std::endl;
  } catch (const zdock::Exception &e) {
    // something went wrong
    zdock::usage(argv[0], e.what());
    return 1;
  }
}
