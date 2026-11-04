// Copyright (c) 2019-2026 Arjan van der Velde, Weng Lab
// SPDX-License-Identifier: BSD-2-Clause

#include "Test.hpp"
#include "pdb++.h"

#include <climits>
#include <cstring>
#include <limits>

namespace {
void require_bounded_record(const libpdb::PDB &record) {
  const char *text = record.chars();
  REQUIRE(text != nullptr);
  REQUIRE(std::strlen(text) <= libpdb::PDB::BufLen - 1);
}

void require_truncated_record(const libpdb::PDB &record) {
  const char *text = record.chars();
  REQUIRE(text != nullptr);
  REQUIRE(std::strlen(text) == libpdb::PDB::BufLen - 1);
}

template <size_t N> void fill_field(char (&field)[N]) {
  std::memset(field, 'x', N - 1);
  field[N - 1] = '\0';
}
} // namespace

TEST_CASE("PDBRUN v5 user records are bounded", "[pdb-format]") {
  // Invariant: every supported v5 USER record fits in PDB's fixed output buffer.
  libpdb::PDB::PdbrunOutputVersion(5);

  SECTION("background color") {
    libpdb::PDB record(libpdb::PDB::USER_BGCOLOR);
    record.userBgColor.rgb[0] = std::numeric_limits<double>::max();
    record.userBgColor.rgb[1] = std::numeric_limits<double>::max();
    record.userBgColor.rgb[2] = std::numeric_limits<double>::max();
    require_truncated_record(record);
  }

  SECTION("named color") {
    libpdb::PDB record(libpdb::PDB::USER_CNAME);
    fill_field(record.userCName.name);
    record.userCName.rgb[0] = std::numeric_limits<double>::max();
    record.userCName.rgb[1] = std::numeric_limits<double>::max();
    record.userCName.rgb[2] = std::numeric_limits<double>::max();
    require_truncated_record(record);
  }

  SECTION("color specification") {
    libpdb::PDB record(libpdb::PDB::USER_COLOR);
    fill_field(record.userColor.spec);
    record.userColor.rgb[0] = std::numeric_limits<double>::max();
    record.userColor.rgb[1] = std::numeric_limits<double>::max();
    record.userColor.rgb[2] = std::numeric_limits<double>::max();
    require_truncated_record(record);
  }

  SECTION("chain") {
    libpdb::PDB record(libpdb::PDB::USER_CHAIN);
    record.userChain.atom0 = INT_MIN;
    record.userChain.atom1 = INT_MAX;
    require_bounded_record(record);
  }

  SECTION("graphics color") {
    libpdb::PDB record(libpdb::PDB::USER_GFX_COLOR);
    fill_field(record.userGfxColor.spec);
    record.userGfxColor.rgb[0] = std::numeric_limits<double>::max();
    record.userGfxColor.rgb[1] = std::numeric_limits<double>::max();
    record.userGfxColor.rgb[2] = std::numeric_limits<double>::max();
    require_truncated_record(record);
  }

  SECTION("graphics font") {
    libpdb::PDB record(libpdb::PDB::USER_GFX_FONT);
    fill_field(record.userGfxFont.name);
    record.userGfxFont.size = INT_MIN;
    require_bounded_record(record);
  }

  SECTION("graphics label") {
    libpdb::PDB record(libpdb::PDB::USER_GFX_LABEL);
    record.userGfxLabel.xyz[0] = std::numeric_limits<double>::max();
    record.userGfxLabel.xyz[1] = std::numeric_limits<double>::max();
    record.userGfxLabel.xyz[2] = std::numeric_limits<double>::max();
    fill_field(record.userGfxLabel.text);
    require_truncated_record(record);
  }
}

TEST_CASE("Standard records are bounded at maximum width", "[pdb-format]") {
  // Invariant: full-width standard records fit the buffer exactly, and
  // overflowing numeric fields degrade to asterisks without exceeding it.
  libpdb::PDB::PdbrunOutputVersion(6);

  SECTION("full-width remark") {
    libpdb::PDB record(libpdb::PDB::REMARK);
    record.agrdes.num = 1;
    fill_field(record.agrdes.text);
    require_bounded_record(record);
    REQUIRE(std::strlen(record.chars()) == 79);
  }

  SECTION("anisou overflow") {
    libpdb::PDB record(libpdb::PDB::ANISOU);
    record.anisou.serialNum = 1;
    for (int &u : record.anisou.u) {
      u = 99999999;
    }
    require_bounded_record(record);
    REQUIRE(std::strstr(record.chars(), "*******") != nullptr);
  }

  SECTION("atom coordinate overflow") {
    libpdb::PDB record(libpdb::PDB::ATOM);
    record.atom.serialNum = 1;
    record.atom.xyz[0] = 99999.999;
    record.atom.xyz[1] = 99999.999;
    record.atom.xyz[2] = 99999.999;
    record.atom.occupancy = 99999.999;
    record.atom.tempFactor = 99999.999;
    require_bounded_record(record);
    REQUIRE(std::strstr(record.chars(), "********") != nullptr);
    REQUIRE(std::strstr(record.chars(), "******") != nullptr);
  }

  SECTION("master count overflow") {
    libpdb::PDB record(libpdb::PDB::MASTER);
    record.master.numCoordinate = 100000;
    require_bounded_record(record);
    REQUIRE(std::strstr(record.chars(), "*****") != nullptr);
  }

  SECTION("v6 background color overflow") {
    libpdb::PDB record(libpdb::PDB::USER_BGCOLOR);
    record.userBgColor.rgb[0] = 99999.999;
    record.userBgColor.rgb[1] = 99999.999;
    record.userBgColor.rgb[2] = 99999.999;
    require_bounded_record(record);
    REQUIRE(std::strstr(record.chars(), "*****") != nullptr);
  }
}
