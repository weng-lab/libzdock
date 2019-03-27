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

#include "Eigen/Dense"
#include "PDB.hpp"
#include "ZDOCK.hpp"
#include "TransformLigand.hpp"
#include "Test.hpp"

#include <string>

TEST_CASE("PDB rotation", "[prediction]") {
  const std::string _zdock = "2OOB/zdock.out.pruned";
  const std::string _ligand = "2OOB/ligand.pdb";
  const std::string _prefix = "2OOB/ligand.";

  // known values for tests
  const int nrows = 574;
  const int ncols = 3;
  const int npredictions = 292;
  const std::string ligandfn = "/tmp/tmpbOIZbH";
  const double ligtranslation[] = {12.926, -3.686, 27.546};
  const double ligrotation[] = {2.380392, 2.500101, -0.120466};

  // be sure these paths exist.
  REQUIRE_NOTHROW(test::getpath(_zdock));
  REQUIRE_NOTHROW(test::getpath(_ligand));

  const std::string zdock = test::getpath(_zdock);
  const std::string ligand = test::getpath(_ligand);
  
  SECTION("Load ZDOCK and PDB file") {
    zdock::ZDOCK z(zdock);
    zdock::PDB p(ligand);
    zdock::TransformLigand txl(z);

    SECTION("Count predictions") {
      REQUIRE(npredictions == z.npredictions());
      REQUIRE(ligandfn == z.ligand().filename);
      REQUIRE((ligtranslation[0] == z.ligand().translation[0] &&
               ligtranslation[1] == z.ligand().translation[1] &&
               ligtranslation[2] == z.ligand().translation[2]));
      REQUIRE((ligrotation[0] == z.ligand().rotation[0] &&
               ligrotation[1] == z.ligand().rotation[1] &&
               ligrotation[2] == z.ligand().rotation[2]));
    }

    SECTION("PDB dimensions") {
      REQUIRE(nrows == p.matrix().cols());
      REQUIRE(ncols == p.matrix().rows());
    }

    SECTION("PDB rotations") {
      const double epsilon = 1.5e-04;
      for (int i = 1; i < 11; ++i) {
        SECTION("Prediction " + std::to_string(i)) {
          REQUIRE_NOTHROW(test::getpath(_prefix + std::to_string(i) + ".pdb"));
          zdock::PDB l(test::getpath(_prefix + std::to_string(i) + ".pdb"));
          REQUIRE(nrows == l.matrix().cols());
          REQUIRE(ncols == l.matrix().rows());
          REQUIRE((txl.txLigand(p.matrix(), z.predictions()[i-1]) - l.matrix()).squaredNorm() < epsilon);
        }
      }
    }

  }

}
