#pragma once

#define STR_EXPAND(tok) #tok
#define STR(tok) STR_EXPAND(tok)
#include "Utils.hpp"
#include "catch2/catch.hpp"

namespace test {
const std::string getpath(const std::string &p = "");
} // namespace test
