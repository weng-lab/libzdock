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
#include "ZDOCK.hpp"
#include "Test.hpp"
#include <string>

TEST_CASE("ZDOCK parser functionality", "[ZDOCK]") {
  // data directory, relative to DATADIR
  const std::string _zdock = "ZDOCK";

  // == Regular M-ZDOCK ==
  SECTION("Load regular M-ZDOCK output file") {

    // known values
    const int npreds = 1500;
    const int symmetry = 24;
    const int boxsize = 81;
    const double spacing = 1.2;

    // input file
    REQUIRE_NOTHROW(test::getpath(_zdock + "/" + "mzdock.out"));
    zdock::ZDOCK z(test::getpath(_zdock + "/" + "mzdock.out"));

    SECTION("test ismzdock()") {
      REQUIRE(z.ismzdock());
      REQUIRE(!z.iszdock());
    }

    SECTION("receptor() _must_ throw exception") {
      REQUIRE_THROWS(z.receptor());
    }

    SECTION("ligand() _must_ throw exception") {
      REQUIRE_THROWS(z.ligand());
    }

    SECTION("structure() must _not_ throw exception") {
      REQUIRE_NOTHROW(z.structure());
    }

    SECTION("count predictions") {
      REQUIRE(npreds == z.npredictions());
    }

    SECTION("get symmetry") {
      REQUIRE(symmetry == z.symmetry());
    }

    SECTION("get box size") {
      REQUIRE(boxsize == z.boxsize());
    }

    SECTION("get grid spacing") {
      REQUIRE(spacing == z.spacing());
    }

    SECTION("should never be switched; this is mzdock") {
      REQUIRE(!z.isswitched());
    }

    SECTION("inspect structure definition") {
      const std::string fn = "/tmp/structure.pdb";
      const double rotation[] = {2.380392, 2.500101, 0.0};
      const double translation[] = {12.348, 1.190, 16.040};
      REQUIRE(fn == z.structure().filename);
      REQUIRE((z.structure().translation[0] == translation[0] &&
               z.structure().translation[1] == translation[1] &&
               z.structure().translation[2] == translation[2]));
      REQUIRE((z.structure().rotation[0] == rotation[0] &&
               z.structure().rotation[1] == rotation[1] &&
               z.structure().rotation[2] == rotation[2]));
    }

    SECTION("inspect prediction (10)") {
      const int translation[] = {3, 6, 0};
      const double rotation[] = {2.410160, 1.065020, 0.0};
      const double score = 18.06;
      const auto pred = z.predictions()[10];
      REQUIRE((translation[0] == pred.translation[0] &&
               translation[1] == pred.translation[1] &&
               translation[2] == pred.translation[2]));
      REQUIRE((rotation[0] == pred.rotation[0] &&
               rotation[1] == pred.rotation[1] &&
               rotation[2] == pred.rotation[2]));
      REQUIRE(score == pred.score);
    }

    SECTION("inspect prediction (200)") {
      const int translation[] = {4, 5, 0};
      const double rotation[] = {2.649080, 0.822633, 0.0};
      const double score = 14.14;
      const auto pred = z.predictions()[200];
      REQUIRE((translation[0] == pred.translation[0] &&
               translation[1] == pred.translation[1] &&
               translation[2] == pred.translation[2]));
      REQUIRE((rotation[0] == pred.rotation[0] &&
               rotation[1] == pred.rotation[1] &&
               rotation[2] == pred.rotation[2]));
      REQUIRE(score == pred.score);
    }
  }


  // == old style ZDOCK ==
  SECTION("Load old-style ZDOCK output file") {

    // known values
    const int npreds = 54000;
    const int boxsize = 128;
    const double spacing = 1.2;

    // input file
    REQUIRE_NOTHROW(test::getpath(_zdock + "/" + "2MTA.zd.out"));
    zdock::ZDOCK z(test::getpath(_zdock + "/" + "2MTA.zd.out"));

    SECTION("test ismzdock()") {
      REQUIRE(!z.ismzdock());
      REQUIRE(z.iszdock());
    }

    SECTION("receptor() must _not_ throw exception") {
      REQUIRE_NOTHROW(z.receptor());
    }

    SECTION("ligand() must _not_ throw exception") {
      REQUIRE_NOTHROW(z.ligand());
    }

    SECTION("structure() _must_ throw exception") {
      REQUIRE_THROWS(z.structure());
    }

    SECTION("count predictions") {
      REQUIRE(npreds == z.npredictions());
    }

    SECTION("symmetry() _must_ throw exception") {
      REQUIRE_THROWS(z.symmetry());
    }

    SECTION("get box size") {
      REQUIRE(boxsize == z.boxsize());
    }

    SECTION("get grid spacing") {
      REQUIRE(spacing == z.spacing());
    }

    SECTION("old-style ZDOCK output does not switch ligand/receptor") {
      REQUIRE(!z.isswitched());
    }

    SECTION("version should be 0 (for old-style)") {
      const int version = 0;
      REQUIRE(version == z.version());
    }

    SECTION("inspect receptor definition") {
      const std::string fn = "2MTA_r_u.pdb";
      const double rotation[] = {0.0, 0.0, 0.0};
      const double translation[] = {24.214, 2.719, 19.743};
      REQUIRE(fn == z.receptor().filename);
      REQUIRE((z.receptor().translation[0] == translation[0] &&
               z.receptor().translation[1] == translation[1] &&
               z.receptor().translation[2] == translation[2]));
      REQUIRE((z.receptor().rotation[0] == rotation[0] &&
               z.receptor().rotation[1] == rotation[1] &&
               z.receptor().rotation[2] == rotation[2]));
    }

    SECTION("inspect ligand definition") {
      const std::string fn = "2MTA_l_u.pdb";
      const double rotation[] = {-2.161817, 2.994959, 0.315678};
      const double translation[] = {5.128, 30.543, 21.422};
      REQUIRE(fn == z.ligand().filename);
      REQUIRE((z.ligand().translation[0] == translation[0] &&
               z.ligand().translation[1] == translation[1] &&
               z.ligand().translation[2] == translation[2]));
      REQUIRE((z.ligand().rotation[0] == rotation[0] &&
               z.ligand().rotation[1] == rotation[1] &&
               z.ligand().rotation[2] == rotation[2]));
    }

    SECTION("inspect prediction (10)") {
      const int translation[] = {12, 105, 122};
      const double rotation[] = {-1.675516, 1.515235, 1.068565};
      const double score = 1182.094;
      const auto pred = z.predictions()[10];
      REQUIRE((translation[0] == pred.translation[0] &&
               translation[1] == pred.translation[1] &&
               translation[2] == pred.translation[2]));
      REQUIRE((rotation[0] == pred.rotation[0] &&
               rotation[1] == pred.rotation[1] &&
               rotation[2] == pred.rotation[2]));
      REQUIRE(score == pred.score);
    }

    SECTION("inspect prediction (200)") {
      const int translation[] = {122, 125, 97};
      const double rotation[] = {-2.722714, 0.977824, 0.441928};
      const double score = 1052.015;
      const auto pred = z.predictions()[200];
      REQUIRE((translation[0] == pred.translation[0] &&
               translation[1] == pred.translation[1] &&
               translation[2] == pred.translation[2]));
      REQUIRE((rotation[0] == pred.rotation[0] &&
               rotation[1] == pred.rotation[1] &&
               rotation[2] == pred.rotation[2]));
      REQUIRE(score == pred.score);
    }
  }

  // == new style ZDOCK (no switching) ==
  SECTION("Load new-style ZDOCK output file, non-switched structures") {

    // known values
    const int npreds = 2000;
    const int boxsize = 144;
    const double spacing = 1.2;

    // input file
    REQUIRE_NOTHROW(test::getpath(_zdock + "/" + "6GWC.zd.out"));
    zdock::ZDOCK z(test::getpath(_zdock + "/" + "6GWC.zd.out"));

    SECTION("test ismzdock()") {
      REQUIRE(!z.ismzdock());
      REQUIRE(z.iszdock());
    }

    SECTION("receptor() must _not_ throw exception") {
      REQUIRE_NOTHROW(z.receptor());
    }

    SECTION("ligand() must _not_ throw exception") {
      REQUIRE_NOTHROW(z.ligand());
    }

    SECTION("structure() _must_ throw exception") {
      REQUIRE_THROWS(z.structure());
    }

    SECTION("count predictions") {
      REQUIRE(npreds == z.npredictions());
    }

    SECTION("symmetry() _must_ throw exception") {
      REQUIRE_THROWS(z.symmetry());
    }

    SECTION("get box size") {
      REQUIRE(boxsize == z.boxsize());
    }

    SECTION("get grid spacing") {
      REQUIRE(spacing == z.spacing());
    }

    SECTION("old-style ZDOCK output does not switch ligand/receptor") {
      REQUIRE(!z.isswitched());
    }

    SECTION("version should be 0 (for old-style)") {
      const int version = 1;
      REQUIRE(version == z.version());
    }

    SECTION("inspect receptor definition") {
      const std::string fn = "/tmp/tmpx09xwv";
      const double rotation[] = {-3.141593, 2.153910, -1.714959};
      const double translation[] = {22.977, -55.978, -37.418};
      REQUIRE(fn == z.receptor().filename);
      REQUIRE((z.receptor().translation[0] == translation[0] &&
               z.receptor().translation[1] == translation[1] &&
               z.receptor().translation[2] == translation[2]));
      REQUIRE((z.receptor().rotation[0] == rotation[0] &&
               z.receptor().rotation[1] == rotation[1] &&
               z.receptor().rotation[2] == rotation[2]));
    }

    SECTION("inspect ligand definition") {
      const std::string fn = "/tmp/tmpL7uaFH";
      const double rotation[] = {2.380392, 2.500101, -0.120466};
      const double translation[] = {9.821, 3.935, 2.207};
      REQUIRE(fn == z.ligand().filename);
      REQUIRE((z.ligand().translation[0] == translation[0] &&
               z.ligand().translation[1] == translation[1] &&
               z.ligand().translation[2] == translation[2]));
      REQUIRE((z.ligand().rotation[0] == rotation[0] &&
               z.ligand().rotation[1] == rotation[1] &&
               z.ligand().rotation[2] == rotation[2]));
    }

    SECTION("inspect prediction (80)") {
      const int translation[] = {132, 129, 7};
      const double rotation[] = {-0.523599, 2.343981, -1.040456};
      const double score = 1061.583;
      const auto pred = z.predictions()[80];
      REQUIRE((translation[0] == pred.translation[0] &&
               translation[1] == pred.translation[1] &&
               translation[2] == pred.translation[2]));
      REQUIRE((rotation[0] == pred.rotation[0] &&
               rotation[1] == pred.rotation[1] &&
               rotation[2] == pred.rotation[2]));
      REQUIRE(score == pred.score);
    }

    SECTION("inspect prediction (933)") {
      const int translation[] = {21, 13, 17};
      const double rotation[] = {0.000000, 1.569311, -3.003978};
      const double score = 855.715;
      const auto pred = z.predictions()[933];
      REQUIRE((translation[0] == pred.translation[0] &&
               translation[1] == pred.translation[1] &&
               translation[2] == pred.translation[2]));
      REQUIRE((rotation[0] == pred.rotation[0] &&
               rotation[1] == pred.rotation[1] &&
               rotation[2] == pred.rotation[2]));
      REQUIRE(score == pred.score);
    }
  }

  // == new style ZDOCK (_with_ switching) ==
  SECTION("Load new-style ZDOCK output file, switched structures") {

    // known values
    const int npreds = 2000;
    const int boxsize = 128;
    const double spacing = 1.2;

    // input file
    REQUIRE_NOTHROW(test::getpath(_zdock + "/" + "4EEW.zd.out"));
    zdock::ZDOCK z(test::getpath(_zdock + "/" + "4EEW.zd.out"));

    SECTION("test ismzdock()") {
      REQUIRE(!z.ismzdock());
      REQUIRE(z.iszdock());
    }

    SECTION("receptor() must _not_ throw exception") {
      REQUIRE_NOTHROW(z.receptor());
    }

    SECTION("ligand() must _not_ throw exception") {
      REQUIRE_NOTHROW(z.ligand());
    }

    SECTION("structure() _must_ throw exception") {
      REQUIRE_THROWS(z.structure());
    }

    SECTION("count predictions") {
      REQUIRE(npreds == z.npredictions());
    }

    SECTION("symmetry() _must_ throw exception") {
      REQUIRE_THROWS(z.symmetry());
    }

    SECTION("get box size") {
      REQUIRE(boxsize == z.boxsize());
    }

    SECTION("get grid spacing") {
      REQUIRE(spacing == z.spacing());
    }

    SECTION("old-style ZDOCK output does not switch ligand/receptor") {
      REQUIRE(z.isswitched());
    }

    SECTION("version should be 0 (for old-style)") {
      const int version = 1;
      REQUIRE(version == z.version());
    }

    SECTION("inspect receptor definition") {
      const std::string fn = "/tmp/tmprW_Py4";
      const double rotation[] = {-2.617994, 1.450119, -2.739050};
      const double translation[] = {-4.599, -21.244, -12.573};
      REQUIRE(fn == z.receptor().filename);
      REQUIRE((z.receptor().translation[0] == translation[0] &&
               z.receptor().translation[1] == translation[1] &&
               z.receptor().translation[2] == translation[2]));
      REQUIRE((z.receptor().rotation[0] == rotation[0] &&
               z.receptor().rotation[1] == rotation[1] &&
               z.receptor().rotation[2] == rotation[2]));
    }

    SECTION("inspect ligand definition") {
      const std::string fn = "/tmp/tmpWuoYDX";
      const double rotation[] = {2.380392, 2.500101, -0.120466};
      const double translation[] = {-37.470, -28.981, -1.051};
      REQUIRE(fn == z.ligand().filename);
      REQUIRE((z.ligand().translation[0] == translation[0] &&
               z.ligand().translation[1] == translation[1] &&
               z.ligand().translation[2] == translation[2]));
      REQUIRE((z.ligand().rotation[0] == rotation[0] &&
               z.ligand().rotation[1] == rotation[1] &&
               z.ligand().rotation[2] == rotation[2]));
    }

    SECTION("inspect prediction (1999)") {
      const int translation[] = {121, 99, 125};
      const double rotation[] = {-1.308997, 1.504028, 0.636812};;
      const double score = 565.418;
      const auto pred = z.predictions()[1999];
      REQUIRE((translation[0] == pred.translation[0] &&
               translation[1] == pred.translation[1] &&
               translation[2] == pred.translation[2]));
      REQUIRE((rotation[0] == pred.rotation[0] &&
               rotation[1] == pred.rotation[1] &&
               rotation[2] == pred.rotation[2]));
      REQUIRE(score == pred.score);
    }

    SECTION("inspect prediction (0)") {
      const int translation[] = {4, 34, 17};;
      const double rotation[] = {-3.141593, 0.899350, -2.069907};
      const double score = 1099.992;
      const auto pred = z.predictions()[0];
      REQUIRE((translation[0] == pred.translation[0] &&
               translation[1] == pred.translation[1] &&
               translation[2] == pred.translation[2]));
      REQUIRE((rotation[0] == pred.rotation[0] &&
               rotation[1] == pred.rotation[1] &&
               rotation[2] == pred.rotation[2]));
      REQUIRE(score == pred.score);
    }
  }
}
