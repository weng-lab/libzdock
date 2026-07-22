// Copyright (c) 2019-2026 Arjan van der Velde, Weng Lab
// SPDX-License-Identifier: BSD-2-Clause

#include "Constraints.hpp"
#include "Exception.hpp"
#include "Test.hpp"

#include <sstream>
#include <string>

TEST_CASE("Constraint records parse and format", "[constraints]") {
  // Invariant: valid records preserve both atom identities, distance, and type.
  std::istringstream input(
      "13 OE2 GLU A 5 101 OD1 ASP b 12 7.3 MIN\n"
      "14 CA ALA B 6 102 N GLY c 13 8.5\n");
  zdock::Constraint minimum;
  zdock::Constraint defaultMaximum;

  REQUIRE_NOTHROW(input >> minimum);
  REQUIRE_NOTHROW(input >> defaultMaximum);

  REQUIRE(minimum.recCoord.serialNum == 13);
  REQUIRE(minimum.recCoord.atomName == "OE2");
  REQUIRE(minimum.recCoord.resName == "GLU");
  REQUIRE(minimum.recCoord.chain == 'A');
  REQUIRE(minimum.recCoord.resNum == 5);
  REQUIRE(minimum.ligCoord.serialNum == 101);
  REQUIRE(minimum.ligCoord.atomName == "OD1");
  REQUIRE(minimum.ligCoord.resName == "ASP");
  REQUIRE(minimum.ligCoord.chain == 'b');
  REQUIRE(minimum.ligCoord.resNum == 12);
  REQUIRE(minimum.distance == Catch::Approx(7.3));
  REQUIRE(minimum.constraintType == zdock::Constraint::MIN);
  REQUIRE(defaultMaximum.constraintType == zdock::Constraint::MAX);

  std::ostringstream output;
  output << minimum;
  REQUIRE(output.str() == "13\tOE2\tGLU\tA\t5\t101\tOD1\tASP\tb\t12\t7.3\tMIN");
}

TEST_CASE("Constraint collections load valid files", "[constraints]") {
  // Invariant: a constraints file yields one ordered object per valid input line.
  const zdock::Constraints constraints(
      test::getpath("Constraints/valid.constraints"));

  REQUIRE(constraints.constraints().size() == 2);
  REQUIRE(constraints.constraints()[0].constraintType == zdock::Constraint::MIN);
  REQUIRE(constraints.constraints()[1].constraintType == zdock::Constraint::MAX);
  REQUIRE(constraints.constraints()[1].distance == Catch::Approx(8.5));
}

TEST_CASE("Malformed constraint files identify their line", "[constraints]") {
  // Invariant: malformed input fails atomically and reports its one-based line.
  REQUIRE_THROWS_WITH(
      zdock::Constraints(test::getpath("Constraints/invalid.constraints")),
      "Error reading constraints (line: 2)");
}

TEST_CASE("Malformed constraint records are rejected", "[constraints]") {
  // Invariant: the stream parser rejects records outside the documented grammar.
  std::istringstream input("not a constraint\n");
  zdock::Constraint constraint;

  REQUIRE_THROWS_WITH(input >> constraint, "Constraint format error");
}
