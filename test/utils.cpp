// Copyright (c) 2019-2026 Arjan van der Velde, Weng Lab
// SPDX-License-Identifier: BSD-2-Clause

#include "Exception.hpp"
#include "Test.hpp"
#include "Utils.hpp"

TEST_CASE("Utils rejects empty companion paths", "[utils]") {
  // Invariant: companion-path resolution never indexes an empty filename.
  REQUIRE_THROWS_AS(zdock::Utils::copath(test::getpath(), ""),
                    zdock::PathException);
}
