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

#include "CreateMultimer.hpp"
#include "PDB.hpp"
#include "TransformMultimer.hpp"
#include "Utils.hpp"
#include <string>
#include <unistd.h>

namespace p = libpdb;
namespace zdock {

CreateMultimer::CreateMultimer(const std::string &zdockoutput,
                               const std::string &receptor, const size_t n,
                               const int mer)
    : zdockfn_(zdockoutput), receptorfn_(receptor), n_(n), mer_(mer) {}

void CreateMultimer::doCreate() {
  std::string recfn; // receptor file name
  ZDOCK z(zdockfn_); // zdock output parser

  // check n within range
  if (n_ < 1 || n_ > z.npredictions()) {
    throw CreateMultimerException("Invalid prediction; valid range 1 - " +
                                  std::to_string(z.npredictions()));
  }

  // check m within range
  if (mer_ >= z.symmetry() ||
      static_cast<int>(sizeof(TransformMultimer::CHAINS)) <= mer_) {
    throw CreateMultimerException("Invalid component; valid range 0 - " +
                                  std::to_string(z.symmetry() - 1));
  }

  // grab prediction
  const zdock::Prediction &pred = z.predictions()[n_ - 1];

  // figure out receptor file name
  if ("" == receptorfn_) {
    recfn = z.receptor().filename;
    if ('/' != recfn[0]) { // relative
      recfn = Utils::copath(zdockfn_, recfn);
    }
  } else {
    recfn = receptorfn_;
  }

  // load receptor
  PDB rec(recfn);
  TransformMultimer t(z);
  const auto m = rec.matrix();
  if (-1 == mer_) {
    int serial = 0;
    for (int i = 0; i < z.symmetry(); ++i) {
      rec.setMatrix(t.txMultimer(m, pred, i));
      for (const auto &x : rec.records()) {
        x->atom.serialNum = ++serial;
        x->atom.residue.chainId = TransformMultimer::CHAINS[i].c_str()[0];
        std::cout << *x << '\n'; // no flush
      }
    }
  } else {
    rec.setMatrix(t.txMultimer(m, pred, mer_));
    for (const auto &x : rec.records()) {
      x->atom.residue.chainId = TransformMultimer::CHAINS[mer_].c_str()[0];
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
      << "  -r <filename>   receptor PDB filename; defaults to receptor in "
         "ZDOCK output\n"
      << "  -m <mer>        component of multimer to output\n"
      << std::endl;
}

} // namespace zdock

int main(int argc, char *argv[]) {
  std::string zdockfn, recfn;
  size_t n = 1;
  int m = -1;
  int c;
  while ((c = getopt(argc, argv, "hn:r:m:")) != -1) {
    switch (c) {
    case 'm': // prediction index
      m = std::stoi(optarg);
      break;
    case 'n': // prediction index
      n = std::stoi(optarg);
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
    zdock::CreateMultimer c(zdockfn, recfn, n, m);
    c.doCreate();
  } catch (const zdock::Exception &e) {
    // something went wrong
    zdock::usage(argv[0], e.what());
    return 1;
  }
  return 0;
}

