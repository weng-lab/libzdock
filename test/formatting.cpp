#include "Test.hpp"
#include "pdb++.h"

#include <cstring>

namespace {
void require_bounded_record(const libpdb::PDB &record) {
  const char *text = record.chars();
  REQUIRE(text != nullptr);
  REQUIRE(std::strlen(text) <= libpdb::PDB::BufLen - 1);
}
} // namespace

TEST_CASE("PDBRUN v5 user records are bounded", "[pdb-format]") {
  libpdb::PDB::PdbrunOutputVersion(5);

  SECTION("background color") {
    libpdb::PDB record(libpdb::PDB::USER_BGCOLOR);
    record.userBgColor.rgb[0] = 0.1;
    record.userBgColor.rgb[1] = 0.2;
    record.userBgColor.rgb[2] = 0.3;
    require_bounded_record(record);
  }

  SECTION("named color") {
    libpdb::PDB record(libpdb::PDB::USER_CNAME);
    std::strcpy(record.userCName.name, "name");
    record.userCName.rgb[0] = 0.1;
    record.userCName.rgb[1] = 0.2;
    record.userCName.rgb[2] = 0.3;
    require_bounded_record(record);
  }

  SECTION("color specification") {
    libpdb::PDB record(libpdb::PDB::USER_COLOR);
    std::strcpy(record.userColor.spec, "spec");
    record.userColor.rgb[0] = 0.1;
    record.userColor.rgb[1] = 0.2;
    record.userColor.rgb[2] = 0.3;
    require_bounded_record(record);
  }

  SECTION("chain") {
    libpdb::PDB record(libpdb::PDB::USER_CHAIN);
    record.userChain.atom0 = 1;
    record.userChain.atom1 = 2;
    require_bounded_record(record);
  }

  SECTION("graphics color") {
    libpdb::PDB record(libpdb::PDB::USER_GFX_COLOR);
    std::strcpy(record.userGfxColor.spec, "spec");
    record.userGfxColor.rgb[0] = 0.1;
    record.userGfxColor.rgb[1] = 0.2;
    record.userGfxColor.rgb[2] = 0.3;
    require_bounded_record(record);
  }

  SECTION("graphics font") {
    libpdb::PDB record(libpdb::PDB::USER_GFX_FONT);
    std::strcpy(record.userGfxFont.name, "font");
    record.userGfxFont.size = 12;
    require_bounded_record(record);
  }

  SECTION("graphics label") {
    libpdb::PDB record(libpdb::PDB::USER_GFX_LABEL);
    record.userGfxLabel.xyz[0] = 1.0;
    record.userGfxLabel.xyz[1] = 2.0;
    record.userGfxLabel.xyz[2] = 3.0;
    std::strcpy(record.userGfxLabel.text, "label");
    require_bounded_record(record);
  }
}
