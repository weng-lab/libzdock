// Copyright (c) 2019-2026 Arjan van der Velde, Weng Lab
// SPDX-License-Identifier: BSD-2-Clause

#include "Exception.hpp"
#include "Test.hpp"
#include "ZDOCK.hpp"

#include <sstream>
#include <string>

TEST_CASE("ZDOCK rejects unavailable and malformed input", "[ZDOCK]") {
  // Invariant: unreadable files and structurally invalid headers never parse.
  REQUIRE_THROWS_AS(zdock::ZDOCK("/path/that/does/not/exist.zdock"),
                    zdock::ZDOCKInvalidFormat);
  REQUIRE_THROWS_AS(
      zdock::ZDOCK(test::getpath("ZDOCK/invalid_header.out")),
      zdock::ZDOCKInvalidFormat);
  REQUIRE_THROWS_AS(
      zdock::ZDOCK(test::getpath("ZDOCK/incomplete_translation.out")),
      zdock::ZDOCKInvalidFormat);
}

TEST_CASE("ZDOCK rejects mixed prediction formats", "[ZDOCK]") {
  // Invariant: one file cannot switch between ZDOCK and M-ZDOCK predictions.
  REQUIRE_THROWS_AS(
      zdock::ZDOCK(test::getpath("ZDOCK/mixed_after_mzdock.out")),
      zdock::ZDOCKInvalidFormat);
  REQUIRE_THROWS_AS(
      zdock::ZDOCK(test::getpath("ZDOCK/mixed_after_zdock.out")),
      zdock::ZDOCKInvalidFormat);
}

TEST_CASE("ZDOCK rejects malformed prediction rows", "[ZDOCK]") {
  // Invariant: once predictions begin, every remaining row must be a prediction.
  REQUIRE_THROWS_AS(
      zdock::ZDOCK(test::getpath("ZDOCK/invalid_prediction.out")),
      zdock::ZDOCKInvalidFormat);
}

TEST_CASE("M-ZDOCK requires at least threefold symmetry", "[ZDOCK]") {
  // Invariant: M-ZDOCK symmetry values below three are rejected during parsing.
  REQUIRE_THROWS_WITH(
      zdock::ZDOCK(test::getpath("ZDOCK/invalid_symmetry.out")),
      Catch::Matchers::ContainsSubstring(
          "M-ZDOCK symmetry cannot be less than 3"));
}

TEST_CASE("Const ZDOCK accessors enforce format boundaries", "[ZDOCK]") {
  // Invariant: const structure access follows the same format contract as mutable access.
  const zdock::ZDOCK regular(test::getpath("ZDOCK/6GWC.zd.out"));
  const zdock::ZDOCK multimer(test::getpath("ZDOCK/mzdock.out"));

  REQUIRE_NOTHROW(regular.receptor());
  REQUIRE_NOTHROW(regular.ligand());
  REQUIRE_THROWS_AS(regular.structure(), zdock::ZDOCKUnsupported);
  REQUIRE_THROWS_AS(multimer.receptor(), zdock::ZDOCKUnsupported);
  REQUIRE_THROWS_AS(multimer.ligand(), zdock::ZDOCKUnsupported);
  REQUIRE_NOTHROW(multimer.structure());
}

TEST_CASE("ZDOCK serialization preserves format-specific fields", "[ZDOCK]") {
  // Invariant: textual output retains each format's header and prediction shape.
  const zdock::ZDOCK modern(test::getpath("ZDOCK/6GWC.zd.out"));
  const zdock::ZDOCK fixed(test::getpath("ZDOCK/2MTA.zd.out"));
  const zdock::ZDOCK multimer(test::getpath("ZDOCK/mzdock.out"));
  std::ostringstream modernText;
  std::ostringstream fixedText;
  std::ostringstream multimerText;

  modernText << modern;
  fixedText << fixed;
  multimerText << multimer;

  REQUIRE(modernText.str().find("144\t1.2\t0\n") == 0);
  REQUIRE(fixedText.str().find("128\t1.2\n") == 0);
  REQUIRE(multimerText.str().find("81\t1.2\t24\n") == 0);
  REQUIRE(multimerText.str().find("4.916730\t0.938406\t78\t73\t21.22") !=
          std::string::npos);
}
