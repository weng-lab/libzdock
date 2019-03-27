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
                     const size_t n, const std::string &chain)
    : zdockfn_(zdockoutput), ligandfn_(ligand), chain_(chain), n_(n) {}

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
    x.atom.residue.chainId = chain_.c_str()[0];
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
            << "  -l <filename>   ligand PDB filename; defaults to receptor in ZDOCK output\n"
            << "  -c <char>       chain id to use for output (defaults to 'Z')\n"
            << std::endl;
}

} // namespace zdock

int main(int argc, char *argv[]) {
  std::string zdockfn, ligfn;
  std::string chain("Z");
  int n = 1;
  int c;
  while ((c = getopt(argc, argv, "hn:l:c:")) != -1) {
    switch (c) {
    case 'c':
      chain = std::string(optarg)[0];
      break;
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
    zdock::Centroids ct(zdockfn, ligfn, n, chain);
    ct.doCentroids();
    std::cerr << "duration: " << zdock::Utils::toc(t1) << " sec" << std::endl;
  } catch (const zdock::Exception &e) {
    // something went wrong
    zdock::usage(argv[0], e.what());
    return 1;
  }
}
