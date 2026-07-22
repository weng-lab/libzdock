// Copyright (c) 2019-2026 Arjan van der Velde, Weng Lab
// SPDX-License-Identifier: BSD-2-Clause

#include "Utils.hpp"
#include "Test.hpp"

namespace test {
const std::string getpath(const std::string &p) {
  if ("" == p) {
    return zdock::Utils::realpath(STR(DATADIR));
  }
  return zdock::Utils::realpath(std::string(STR(DATADIR)) + "/" + p);
}
} // namespace test

// be sure datadir exists
TEST_CASE("test data dir", "[pre]") {
  // Invariant: the configured test-data directory exists and is canonicalizable.
  REQUIRE_NOTHROW(test::getpath());
}

// keep this file empty
