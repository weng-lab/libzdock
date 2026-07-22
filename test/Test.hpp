// Copyright (c) 2019-2026 Arjan van der Velde, Weng Lab
// SPDX-License-Identifier: BSD-2-Clause

#pragma once

#define STR_EXPAND(tok) #tok
#define STR(tok) STR_EXPAND(tok)
#include "Utils.hpp"
#include "catch_amalgamated.hpp"

namespace test {
const std::string getpath(const std::string &p = "");
} // namespace test
