#pragma once

#define STR_EXPAND(tok) #tok
#define STR(tok) STR_EXPAND(tok)
#include "catch2/catch.hpp"
