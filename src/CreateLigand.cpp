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

#include "CreateLigand.hpp"
#include "PDB.hpp"
#include "TransformLigand.hpp"
#include "Utils.hpp"
#include <string>
#include <unistd.h>

namespace p = libpdb;
namespace zdock {

CreateLigand::CreateLigand(const std::string &zdockoutput,
                           const std::string &ligand,
                           const std::string &receptor, const size_t n,
                           const bool cmplx, const bool allrecords)
    : zdockfn_(zdockoutput), ligandfn_(ligand), receptorfn_(receptor), n_(n),
      complex_(cmplx), allrecords_(allrecords) {}

void CreateLigand::doCreate() {
  std::string ligfn; // ligand file name
  std::string recfn; // receptor file name
  ZDOCK z(zdockfn_); // zdock output parser

  // check n within range
  if (n_ < 1 || n_ > z.npredictions()) {
    throw CreateLigandException("Invalid prediction; valid range 1 - " +
                                std::to_string(z.npredictions()));
  }
  const zdock::Prediction &pred = z.predictions()[n_ - 1];

  // figure out ligand file name
  if ("" == ligandfn_) {
    ligfn = Utils::copath(zdockfn_, z.ligand().filename);
  } else {
    ligfn = ligandfn_;
  }

  // load ligand
  PDB lig(ligfn);
  PDB rec;
  TransformLigand t(z);
  lig.setMatrix(t.txLigand(lig.matrix(), pred));

  if (complex_) {
    // figure out receptor file name
    if ("" == receptorfn_) {
      recfn = Utils::copath(zdockfn_, z.receptor().filename);
    } else {
      recfn = receptorfn_;
    }
    rec = PDB(recfn);
  }

  for (const auto &x : (allrecords_ ? lig.records() : lig.atoms())) {
    std::cout << *x << '\n'; // no flush
  }

  if (complex_) {
    for (const auto &x : (allrecords_ ? rec.records() : rec.atoms())) {
      std::cout << *x << '\n'; // no flush
    }
  }
}

void usage(const std::string &cmd, const std::string &err = "") {
  // print error if any
  if ("" != err) {
    std::cerr << "Error: " << err << std::endl << std::endl;
  }
  // print usage
  std::cerr
      << "usage: " << cmd << " [options] <zdock output>\n\n"
      << "  -n <integer>    index of prediction in ZDOCK file (defaults "
         "to 1; the top prediction)\n"
      << "  -c              create complex; by default only ligand is created\n"
      << "  -r <filename>   receptor PDB filename; defaults to receptor in "
         "ZDOCK output\n"
      << "  -l <filename>   ligand PDB filename; defaults to ligand in "
         "ZDOCK output\n"
      << "  -a              return all records (by default only ATOM and "
         "HETATM are returned)\n"

      << std::endl;
}

} // namespace zdock

int main(int argc, char *argv[]) {
  std::string zdockfn, ligfn, recfn;
  size_t n = 1;
  int c;
  bool cmplx = false;
  bool allrecords = false;
  while ((c = getopt(argc, argv, "achn:l:r:")) != -1) {
    switch (c) {
    case 'a':
      allrecords = true;
      break;
    case 'c':
      cmplx = true;
      break;
    case 'n': // prediction index
      n = std::stoi(optarg);
      break;
    case 'l': // alternative ligand
      ligfn = optarg;
      break;
    case 'r': // alternative ligand
      recfn = optarg;
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
    zdock::CreateLigand c(zdockfn, ligfn, recfn, n, cmplx, allrecords);
    c.doCreate();
  } catch (const zdock::Exception &e) {
    // something went wrong
    zdock::usage(argv[0], e.what());
    return 1;
  }
  return 0;
}

