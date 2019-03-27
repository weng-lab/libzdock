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
                               const int mer, const bool allrecords)
    : zdockfn_(zdockoutput), receptorfn_(receptor), n_(n), mer_(mer),
      allrecords_(allrecords) {}

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
    recfn = Utils::copath(zdockfn_, z.receptor().filename);
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
      for (const auto &x : (allrecords_ ? rec.records() : rec.atoms())) {
        if (p::PDB::ATOM == x->type() || p::PDB::HETATM == x->type()) {
          x->atom.serialNum = ++serial;
          x->atom.residue.chainId = TransformMultimer::CHAINS[i].c_str()[0];
        }
        std::cout << *x << '\n'; // no flush
      }
    }
  } else {
    rec.setMatrix(t.txMultimer(m, pred, mer_));
    for (const auto &x : (allrecords_ ? rec.records() : rec.atoms())) {
      if (p::PDB::ATOM == x->type() || p::PDB::HETATM == x->type()) {
        x->atom.residue.chainId = TransformMultimer::CHAINS[mer_].c_str()[0];
      }
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
  std::cerr << "usage: " << cmd << " [options] <zdock output>\n\n"
            << "  -n <integer>    index of prediction in M-ZDOCK file "
               "(defaults to 1; the top prediction)\n"
            << "  -r <filename>   receptor PDB filename; defaults to receptor "
               "in M-ZDOCK output\n"
            << "  -m <mer>        component of multimer to output (all if not "
               "specified)\n"
            << "  -a              return all records (by default only ATOM and "
               "HETATM are returned)\n"
            << std::endl;
}

} // namespace zdock

int main(int argc, char *argv[]) {
  std::string zdockfn, recfn;
  bool allrecords = false;
  size_t n = 1;
  int m = -1;
  int c;
  while ((c = getopt(argc, argv, "ahn:r:m:")) != -1) {
    switch (c) {
    case 'a':
      allrecords = true;
      break;
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
    zdock::CreateMultimer c(zdockfn, recfn, n, m, allrecords);
    c.doCreate();
  } catch (const zdock::Exception &e) {
    // something went wrong
    zdock::usage(argv[0], e.what());
    return 1;
  }
  return 0;
}

