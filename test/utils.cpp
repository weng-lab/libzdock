#include "Exception.hpp"
#include "Test.hpp"
#include "Utils.hpp"

TEST_CASE("Utils rejects empty companion paths", "[utils]") {
  // Invariant: companion-path resolution never indexes an empty filename.
  REQUIRE_THROWS_AS(zdock::Utils::copath(test::getpath(), ""),
                    zdock::PathException);
}
