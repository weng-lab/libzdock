#define CATCH_CONFIG_MAIN
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
  REQUIRE_NOTHROW(test::getpath());
}

// keep this file empty
