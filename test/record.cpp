// Copyright (c) 2019-2026 Arjan van der Velde, Weng Lab
// SPDX-License-Identifier: BSD-2-Clause

// Suite for the modern record class (Phase 2). At this stage it covers the
// RecordType enumeration only; parsing, writing, and parity tests with the
// legacy implementation arrive with the later Phase 2 commits.

#include "Test.hpp"
#include "Record.hpp"

#include <set>
#include <string>
#include <vector>

namespace {

using zdock::pdb::RecordType;

struct NamedType {
  RecordType type;
  const char *name;
};

// Every enumerator in declaration order with its expected spelling. The
// list mirrors the legacy libpdb::PDB::RecordType names one-for-one, so a
// renamed, reordered, or missing type fails here.
const std::vector<NamedType> &namedTypes() {
  static const std::vector<NamedType> types = {
      {RecordType::UNKNOWN, "UNKNOWN"},
      {RecordType::ANISOU, "ANISOU"},
      {RecordType::ATOM, "ATOM"},
      {RecordType::AUTHOR, "AUTHOR"},
      {RecordType::CAVEAT, "CAVEAT"},
      {RecordType::COMPND, "COMPND"},
      {RecordType::CONECT, "CONECT"},
      {RecordType::CRYST1, "CRYST1"},
      {RecordType::END, "END"},
      {RecordType::FORMUL, "FORMUL"},
      {RecordType::FTNOTE, "FTNOTE"},
      {RecordType::HEADER, "HEADER"},
      {RecordType::HELIX, "HELIX"},
      {RecordType::HET, "HET"},
      {RecordType::HETATM, "HETATM"},
      {RecordType::JRNL, "JRNL"},
      {RecordType::KEYWDS, "KEYWDS"},
      {RecordType::MASTER, "MASTER"},
      {RecordType::MDLTYP, "MDLTYP"},
      {RecordType::MTRIX, "MTRIX"},
      {RecordType::NUMMDL, "NUMMDL"},
      {RecordType::OBSLTE, "OBSLTE"},
      {RecordType::ORIGX, "ORIGX"},
      {RecordType::REMARK, "REMARK"},
      {RecordType::REVDAT, "REVDAT"},
      {RecordType::SCALE, "SCALE"},
      {RecordType::SEQRES, "SEQRES"},
      {RecordType::SHEET, "SHEET"},
      {RecordType::SIGATM, "SIGATM"},
      {RecordType::SIGUIJ, "SIGUIJ"},
      {RecordType::SITE, "SITE"},
      {RecordType::SOURCE, "SOURCE"},
      {RecordType::SPLIT, "SPLIT"},
      {RecordType::SPRSDE, "SPRSDE"},
      {RecordType::SSBOND, "SSBOND"},
      {RecordType::TER, "TER"},
      {RecordType::TITLE, "TITLE"},
      {RecordType::TURN, "TURN"},
      {RecordType::TVECT, "TVECT"},
      {RecordType::USER, "USER"},
      {RecordType::MODEL, "MODEL"},
      {RecordType::ENDMDL, "ENDMDL"},
      {RecordType::EXPDTA, "EXPDTA"},
      {RecordType::SYMDES, "SYMDES"},
      {RecordType::SYMOP, "SYMOP"},
      {RecordType::MTXDES, "MTXDES"},
      {RecordType::CMPDES, "CMPDES"},
      {RecordType::CMPONT, "CMPONT"},
      {RecordType::TRNSFM, "TRNSFM"},
      {RecordType::AGRDES, "AGRDES"},
      {RecordType::AGGRGT, "AGGRGT"},
      {RecordType::ZDATOM, "ZDATOM"},
      {RecordType::USER_PDBRUN, "USER_PDBRUN"},
      {RecordType::USER_EYEPOS, "USER_EYEPOS"},
      {RecordType::USER_ATPOS, "USER_ATPOS"},
      {RecordType::USER_WINDOW, "USER_WINDOW"},
      {RecordType::USER_FOCUS, "USER_FOCUS"},
      {RecordType::USER_VIEWPORT, "USER_VIEWPORT"},
      {RecordType::USER_BGCOLOR, "USER_BGCOLOR"},
      {RecordType::USER_ANGLE, "USER_ANGLE"},
      {RecordType::USER_DISTANCE, "USER_DISTANCE"},
      {RecordType::USER_FILE, "USER_FILE"},
      {RecordType::USER_MARKNAME, "USER_MARKNAME"},
      {RecordType::USER_MARK, "USER_MARK"},
      {RecordType::USER_CNAME, "USER_CNAME"},
      {RecordType::USER_COLOR, "USER_COLOR"},
      {RecordType::USER_RADIUS, "USER_RADIUS"},
      {RecordType::USER_OBJECT, "USER_OBJECT"},
      {RecordType::USER_ENDOBJ, "USER_ENDOBJ"},
      {RecordType::USER_CHAIN, "USER_CHAIN"},
      {RecordType::USER_GFX_BEGIN, "USER_GFX_BEGIN"},
      {RecordType::USER_GFX_END, "USER_GFX_END"},
      {RecordType::USER_GFX_COLOR, "USER_GFX_COLOR"},
      {RecordType::USER_GFX_NORMAL, "USER_GFX_NORMAL"},
      {RecordType::USER_GFX_VERTEX, "USER_GFX_VERTEX"},
      {RecordType::USER_GFX_FONT, "USER_GFX_FONT"},
      {RecordType::USER_GFX_TEXTPOS, "USER_GFX_TEXTPOS"},
      {RecordType::USER_GFX_LABEL, "USER_GFX_LABEL"},
      {RecordType::USER_GFX_MOVE, "USER_GFX_MOVE"},
      {RecordType::USER_GFX_DRAW, "USER_GFX_DRAW"},
      {RecordType::USER_GFX_MARKER, "USER_GFX_MARKER"},
      {RecordType::USER_GFX_POINT, "USER_GFX_POINT"},
  };
  return types;
}

} // namespace

TEST_CASE("Every record type renders its enumerator spelling",
          "[record-types]") {
  // Invariant: to_string maps each enumerator to its exact spelling, so
  // diagnostics and future error messages name record types unambiguously.
  for (const auto &named : namedTypes()) {
    CAPTURE(named.type);
    REQUIRE(zdock::pdb::to_string(named.type) == named.name);
  }
}

TEST_CASE("Record type names are distinct and values contiguous",
          "[record-types]") {
  // Invariant: the enumeration holds one distinct name per distinct value,
  // numbered contiguously from zero in declaration order, mirroring the
  // legacy enum the Phase 3 migration correlates against.
  const auto &types = namedTypes();
  std::set<std::string> names;
  for (std::size_t i = 0; i < types.size(); ++i) {
    names.insert(types[i].name);
    CAPTURE(types[i].name);
    REQUIRE(static_cast<int>(types[i].type) == static_cast<int>(i));
  }
  REQUIRE(names.size() == types.size());
}

TEST_CASE("Out-of-range record types render as invalid", "[record-types]") {
  // Invariant: values outside the enumeration render as "<invalid>" rather
  // than indexing a name table past its end.
  REQUIRE(zdock::pdb::to_string(static_cast<RecordType>(82)) == "<invalid>");
  REQUIRE(zdock::pdb::to_string(static_cast<RecordType>(200)) == "<invalid>");
}
