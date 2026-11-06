// Copyright (c) 2019-2026 Arjan van der Velde, Weng Lab
// SPDX-License-Identifier: BSD-2-Clause

#include "Record.hpp"

#include "Columns.hpp"

#include <cassert>
#include <cstring>
#include <utility>

namespace zdock::pdb {
namespace {

constexpr std::size_t kJunkWidth = 72; // the pinned unknown-junk field

// Record-name dispatch over the upper-cased first four characters
// (space-padded), mirroring the legacy first-4-char classification. The
// USER subtypes dispatch on the PDBRUN version and land with the USER
// records; for now every USER line classifies as USER.
RecordType classifyLine(std::string_view line) {
  struct KeyType {
    const char *key;
    RecordType type;
  };
  static constexpr KeyType kKeys[] = {
      {"ANIS", RecordType::ANISOU},   {"ATOM", RecordType::ATOM},
      {"AUTH", RecordType::AUTHOR},   {"CAVE", RecordType::CAVEAT},
      {"COMP", RecordType::COMPND},   {"CONE", RecordType::CONECT},
      {"CRYS", RecordType::CRYST1},   {"END ", RecordType::END},
      {"ENDM", RecordType::ENDMDL},   {"FORM", RecordType::FORMUL},
      {"FTNO", RecordType::FTNOTE},   {"HEAD", RecordType::HEADER},
      {"HELI", RecordType::HELIX},    {"HET ", RecordType::HET},
      {"HETA", RecordType::HETATM},   {"JRNL", RecordType::JRNL},
      {"KEYW", RecordType::KEYWDS},   {"MAST", RecordType::MASTER},
      {"MDLT", RecordType::MDLTYP},   {"MTRI", RecordType::MTRIX},
      {"NUMM", RecordType::NUMMDL},   {"OBSL", RecordType::OBSLTE},
      {"ORIG", RecordType::ORIGX},    {"REMA", RecordType::REMARK},
      {"REVD", RecordType::REVDAT},   {"SCAL", RecordType::SCALE},
      {"SEQR", RecordType::SEQRES},   {"SHEE", RecordType::SHEET},
      {"SIGA", RecordType::SIGATM},   {"SIGU", RecordType::SIGUIJ},
      {"SITE", RecordType::SITE},     {"SOUR", RecordType::SOURCE},
      {"SPLI", RecordType::SPLIT},    {"SPRS", RecordType::SPRSDE},
      {"SSBO", RecordType::SSBOND},   {"TER ", RecordType::TER},
      {"TITL", RecordType::TITLE},    {"TURN", RecordType::TURN},
      {"TVEC", RecordType::TVECT},    {"USER", RecordType::USER},
      {"MODE", RecordType::MODEL},    {"EXPD", RecordType::EXPDTA},
      {"SYMD", RecordType::SYMDES},   {"SYMO", RecordType::SYMOP},
      {"MTXD", RecordType::MTXDES},   {"CMPD", RecordType::CMPDES},
      {"CMPO", RecordType::CMPONT},   {"TRNS", RecordType::TRNSFM},
      {"AGRD", RecordType::AGRDES},   {"AGGR", RecordType::AGGRGT},
  };
  char key[4];
  std::size_t i = 0;
  for (; i < 4; ++i) {
    if (i >= line.size() || line[i] == '\0' || line[i] == '\n')
      break;
    const char c = line[i];
    key[i] = (c >= 'a' && c <= 'z') ? static_cast<char>(c - 'a' + 'A') : c;
  }
  for (; i < 4; ++i)
    key[i] = ' ';
  for (const auto &entry : kKeys)
    if (std::memcmp(key, entry.key, 4) == 0)
      return entry.type;
  return RecordType::UNKNOWN;
}

// The USER subtype dispatch for one PDBRUN version, mirroring the legacy
// pdbrun5Type/pdbrun6Type matchers over the text after "USER  ": the
// first subtype letter selects the group and the remaining characters are
// matched case-insensitively against fragments that carry their own
// separating blanks.
RecordType classifyUserSubtype(std::string_view text, int pdbrunVersion) {
  const auto matches = [&text](std::size_t offset, const char *literal) {
    const std::size_t len = std::strlen(literal);
    if (offset + len > text.size())
      return false;
    for (std::size_t i = 0; i < len; ++i) {
      const char c = text[offset + i];
      const char l = literal[i];
      if (c != l && (c < 'a' || c > 'z' || c - 'a' + 'A') != l)
        return false;
    }
    return true;
  };
  const auto boundary = [&text](std::size_t offset) {
    return offset >= text.size() || text[offset] == '\0' ||
           text[offset] == '\n' || text[offset] == ' ';
  };
  if (pdbrunVersion != 5 && pdbrunVersion != 6) {
    // the legacy default branch recognizes only PDBRUN
    if (matches(1, "DBRUN "))
      return RecordType::USER_PDBRUN;
    return RecordType::USER;
  }
  switch (text.empty() ? '\0' : text[0]) {
  case 'A':
  case 'a':
    if (matches(1, "NGLE "))
      return RecordType::USER_ANGLE;
    if (matches(1, "TPOS "))
      return RecordType::USER_ATPOS;
    break;
  case 'B':
  case 'b':
    if (matches(1, "GCOLOR "))
      return RecordType::USER_BGCOLOR;
    break;
  case 'C':
  case 'c':
    if (matches(1, "HAIN "))
      return RecordType::USER_CHAIN;
    if (matches(1, "NAME "))
      return RecordType::USER_CNAME;
    if (matches(1, "OLOR "))
      return RecordType::USER_COLOR;
    break;
  case 'D':
  case 'd':
    if (matches(1, "ISTANCE "))
      return RecordType::USER_DISTANCE;
    break;
  case 'E':
  case 'e':
    if (pdbrunVersion == 6 && matches(1, "NDOBJ") && boundary(6))
      return RecordType::USER_ENDOBJ;
    if (pdbrunVersion == 5 && matches(1, "NDOBJ "))
      return RecordType::USER_ENDOBJ;
    if (matches(1, "YEPOS "))
      return RecordType::USER_EYEPOS;
    break;
  case 'F':
  case 'f':
    if (matches(1, "ILE "))
      return RecordType::USER_FILE;
    if (matches(1, "OCUS "))
      return RecordType::USER_FOCUS;
    break;
  case 'G':
  case 'g':
    if (pdbrunVersion == 6) {
      if (text.size() < 4 || (text[1] != 'F' && text[1] != 'f') ||
          text[2] != 'X' || text[3] != ' ')
        break;
      switch (text[4]) {
      case 'B':
      case 'b':
        if (matches(5, "EGIN "))
          return RecordType::USER_GFX_BEGIN;
        break;
      case 'C':
      case 'c':
        if (matches(5, "OLOR "))
          return RecordType::USER_GFX_COLOR;
        break;
      case 'E':
      case 'e':
        if (matches(5, "ND") && boundary(7))
          return RecordType::USER_GFX_END;
        break;
      case 'F':
      case 'f':
        if (matches(5, "ONT "))
          return RecordType::USER_GFX_FONT;
        break;
      case 'L':
      case 'l':
        if (matches(5, "ABEL "))
          return RecordType::USER_GFX_LABEL;
        break;
      case 'N':
      case 'n':
        if (matches(5, "ORMAL "))
          return RecordType::USER_GFX_NORMAL;
        break;
      case 'T':
      case 't':
        if (matches(5, "EXTPOS "))
          return RecordType::USER_GFX_TEXTPOS;
        break;
      case 'V':
      case 'v':
        if (matches(5, "ERTEX "))
          return RecordType::USER_GFX_VERTEX;
        break;
      }
    } else {
      if (!matches(1, "FX "))
        break;
      if (matches(4, "COLOR "))
        return RecordType::USER_GFX_COLOR;
      if (matches(4, "DRAW "))
        return RecordType::USER_GFX_DRAW;
      if (matches(4, "FONT "))
        return RecordType::USER_GFX_FONT;
      if (matches(4, "LABEL "))
        return RecordType::USER_GFX_LABEL;
      if (matches(4, "MARKER "))
        return RecordType::USER_GFX_MARKER;
      if (matches(4, "MOVE "))
        return RecordType::USER_GFX_MOVE;
      if (matches(4, "POINT "))
        return RecordType::USER_GFX_POINT;
    }
    break;
  case 'M':
  case 'm':
    if (matches(1, "ARKNAME ")) // the longer name must win
      return RecordType::USER_MARKNAME;
    if (matches(1, "ARK "))
      return RecordType::USER_MARK;
    break;
  case 'O':
  case 'o':
    if (pdbrunVersion == 6) {
      if (matches(1, "BJECT") && boundary(6))
        return RecordType::USER_OBJECT;
    } else {
      if (matches(1, "BJECT "))
        return RecordType::USER_OBJECT;
    }
    break;
  case 'P':
  case 'p':
    if (matches(1, "DBRUN "))
      return RecordType::USER_PDBRUN;
    break;
  case 'R':
  case 'r':
    if (matches(1, "ADIUS "))
      return RecordType::USER_RADIUS;
    break;
  case 'V':
  case 'v':
    if (matches(1, "IEWPORT "))
      return RecordType::USER_VIEWPORT;
    break;
  case 'W':
  case 'w':
    if (matches(1, "INDOW "))
      return RecordType::USER_WINDOW;
    break;
  }
  return RecordType::USER;
}

// Layout tables for the simple families. Read and write walks share the
// columns; the continuation columns render %D (zero as blanks) and the
// text columns left-justify at their precision, exactly as the legacy
// write formats did.
constexpr ColumnSpec kAuthorCols[] = {
    {8, 2, ColumnKind::Int, 0, false, ' ', true},
    {10, 69, ColumnKind::Str, 0, true, ' ', false, 69},
};
constexpr ColumnSpec kTitleCols[] = {
    {8, 2, ColumnKind::Int, 0, false, ' ', true},
    {10, 70, ColumnKind::Str, 0, true, ' ', false, 70},
};
constexpr ColumnSpec kCompndCols[] = {
    {7, 3, ColumnKind::Int, 0, false, ' ', true},
    {10, 70, ColumnKind::Str, 0, true, ' ', false, 70},
};
constexpr ColumnSpec kKeywdsCols[] = {
    {8, 2, ColumnKind::Int, 0, false, ' ', true},
    {10, 69, ColumnKind::Str, 0, true, ' ', false, 69},
};
constexpr ColumnSpec kSourceCols[] = {
    {7, 3, ColumnKind::Int, 0, false, ' ', true},
    {10, 69, ColumnKind::Str, 0, true, ' ', false, 69},
};
constexpr ColumnSpec kMdltypCols[] = {
    {8, 2, ColumnKind::Int, 0, false, ' ', true},
    {10, 70, ColumnKind::Str, 0, true, ' ', false, 70},
};

// Standard ATOM/HETATM/SIGATM layout (read and write columns coincide;
// the write styles ride along). The Literal runs are the format's literal
// spaces, which the legacy walk matched one by one -- a non-space there is
// a hard mismatch driving the fallback chain.
constexpr ColumnSpec kAtomCols[] = {
    {0, 6, ColumnKind::Skip},
    {6, 5, ColumnKind::Int},
    {11, 1, ColumnKind::Literal, 0, false, ' '},
    {12, 4, ColumnKind::Str, 0, true, ' ', false, -1},
    {16, 1, ColumnKind::Char},
    {17, 4, ColumnKind::Str, 0, true, ' ', false, -1},
    {21, 1, ColumnKind::Char},
    {22, 4, ColumnKind::Int},
    {26, 1, ColumnKind::Char},
    {27, 3, ColumnKind::Literal, 0, false, ' '},
    {30, 8, ColumnKind::Real, 3},
    {38, 8, ColumnKind::Real, 3},
    {46, 8, ColumnKind::Real, 3},
    {54, 6, ColumnKind::Real, 2},
    {60, 6, ColumnKind::Real, 2},
    {66, 10, ColumnKind::Literal, 0, false, ' '},
    {76, 2, ColumnKind::Str},
    {78, 2, ColumnKind::Str},
};

// The distinct 82-column ZDATOM layout (mark_sur extension): type,
// surface, radius, segment id, and charge instead of occupancy,
// tempFactor, element, and charge. The declared record length is 82.
constexpr ColumnSpec kZdatomCols[] = {
    {0, 6, ColumnKind::Skip},
    {6, 5, ColumnKind::Int},
    {11, 1, ColumnKind::Literal, 0, false, ' '},
    {12, 4, ColumnKind::Str, 0, true, ' ', false, -1},
    {16, 1, ColumnKind::Char},
    {17, 4, ColumnKind::Str, 0, true, ' ', false, -1},
    {21, 1, ColumnKind::Char},
    {22, 4, ColumnKind::Int},
    {26, 1, ColumnKind::Char},
    {27, 3, ColumnKind::Literal, 0, false, ' '},
    {30, 8, ColumnKind::Real, 3},
    {38, 8, ColumnKind::Real, 3},
    {46, 8, ColumnKind::Real, 3},
    {54, 1, ColumnKind::Literal, 0, false, ' '},
    {55, 2, ColumnKind::Int},
    {57, 4, ColumnKind::Literal, 0, false, ' '},
    {61, 2, ColumnKind::Int},
    {63, 1, ColumnKind::Literal, 0, false, ' '},
    {64, 4, ColumnKind::Real, 2},
    {68, 4, ColumnKind::Literal, 0, false, ' '},
    {72, 4, ColumnKind::Str},
    {76, 1, ColumnKind::Literal, 0, false, ' '},
    {77, 5, ColumnKind::Real, 2},
};

// ---------------------------------------------------------------------------
// Residue-heavy and structural families (Phase 2, commit 6). Read and
// write columns coincide except where noted (REVDAT, OBSLTE, SPRSDE,
// SPLIT, MODEL have dual tables).
// ---------------------------------------------------------------------------

constexpr ColumnSpec kHelixCols[] = {
    {0, 7, ColumnKind::Skip},
    {7, 3, ColumnKind::Int, 0, false, ' ', true},
    {10, 1, ColumnKind::Literal, 0, false, ' '},
    {11, 3, ColumnKind::Str},
    {14, 1, ColumnKind::Literal, 0, false, ' '},
    {15, 4, ColumnKind::Str, 0, true},
    {19, 1, ColumnKind::Char},
    {20, 1, ColumnKind::Literal, 0, false, ' '},
    {21, 4, ColumnKind::Int},
    {25, 1, ColumnKind::Char},
    {26, 1, ColumnKind::Literal, 0, false, ' '},
    {27, 4, ColumnKind::Str, 0, true},
    {31, 1, ColumnKind::Char},
    {32, 1, ColumnKind::Literal, 0, false, ' '},
    {33, 4, ColumnKind::Int},
    {37, 1, ColumnKind::Char},
    {38, 2, ColumnKind::Int, 0, false, ' ', true},
    {40, 30, ColumnKind::Str, 0, true},
    {70, 1, ColumnKind::Literal, 0, false, ' '},
    {71, 5, ColumnKind::Int, 0, false, ' ', true},
};

constexpr ColumnSpec kSheetCols[] = {
    {0, 6, ColumnKind::Skip},
    {6, 4, ColumnKind::Int, 0, false, ' ', true},
    {10, 1, ColumnKind::Literal, 0, false, ' '},
    {11, 3, ColumnKind::Str},
    {14, 2, ColumnKind::Int},
    {16, 1, ColumnKind::Literal, 0, false, ' '},
    {17, 4, ColumnKind::Str, 0, true},
    {21, 1, ColumnKind::Char},
    {22, 4, ColumnKind::Int},
    {26, 1, ColumnKind::Char},
    {27, 1, ColumnKind::Literal, 0, false, ' '},
    {28, 4, ColumnKind::Str, 0, true},
    {32, 1, ColumnKind::Char},
    {33, 4, ColumnKind::Int},
    {37, 1, ColumnKind::Char},
    {38, 2, ColumnKind::Int},
    {40, 1, ColumnKind::Literal, 0, false, ' '},
    {41, 4, ColumnKind::Str, 0, true},
    {45, 4, ColumnKind::Str, 0, true},
    {49, 1, ColumnKind::Char},
    {50, 4, ColumnKind::Int, 0, false, ' ', true},
    {54, 1, ColumnKind::Char},
    {55, 1, ColumnKind::Literal, 0, false, ' '},
    {56, 4, ColumnKind::Str, 0, true},
    {60, 4, ColumnKind::Str, 0, true},
    {64, 1, ColumnKind::Char},
    {65, 4, ColumnKind::Int, 0, false, ' ', true},
    {69, 1, ColumnKind::Char},
};

constexpr ColumnSpec kSiteCols[] = {
    {0, 7, ColumnKind::Skip},
    {7, 3, ColumnKind::Int},
    {10, 1, ColumnKind::Literal, 0, false, ' '},
    {11, 3, ColumnKind::Str},
    {14, 1, ColumnKind::Literal, 0, false, ' '},
    {15, 2, ColumnKind::Int},
    {17, 1, ColumnKind::Literal, 0, false, ' '},
    {18, 4, ColumnKind::Str, 0, true},
    {22, 1, ColumnKind::Char},
    {23, 4, ColumnKind::Int, 0, false, ' ', true},
    {27, 1, ColumnKind::Char},
    {28, 1, ColumnKind::Literal, 0, false, ' '},
    {29, 4, ColumnKind::Str, 0, true},
    {33, 1, ColumnKind::Char},
    {34, 4, ColumnKind::Int, 0, false, ' ', true},
    {38, 1, ColumnKind::Char},
    {39, 1, ColumnKind::Literal, 0, false, ' '},
    {40, 4, ColumnKind::Str, 0, true},
    {44, 1, ColumnKind::Char},
    {45, 4, ColumnKind::Int, 0, false, ' ', true},
    {49, 1, ColumnKind::Char},
    {50, 1, ColumnKind::Literal, 0, false, ' '},
    {51, 4, ColumnKind::Str, 0, true},
    {55, 1, ColumnKind::Char},
    {56, 4, ColumnKind::Int, 0, false, ' ', true},
    {60, 1, ColumnKind::Char},
};

constexpr ColumnSpec kSsbondCols[] = {
    {0, 7, ColumnKind::Skip},
    {7, 3, ColumnKind::Int, 0, false, ' ', true},
    {10, 1, ColumnKind::Literal, 0, false, ' '},
    {11, 4, ColumnKind::Str, 0, true},
    {15, 1, ColumnKind::Char},
    {16, 1, ColumnKind::Literal, 0, false, ' '},
    {17, 4, ColumnKind::Int},
    {21, 1, ColumnKind::Char},
    {22, 3, ColumnKind::Literal, 0, false, ' '},
    {25, 4, ColumnKind::Str, 0, true},
    {29, 1, ColumnKind::Char},
    {30, 1, ColumnKind::Literal, 0, false, ' '},
    {31, 4, ColumnKind::Int, 0, false, ' ', true},
    {35, 1, ColumnKind::Char},
    {36, 23, ColumnKind::Skip},
    {59, 6, ColumnKind::Int, 0, false, ' ', true},
    {65, 1, ColumnKind::Literal, 0, false, ' '},
    {66, 6, ColumnKind::Int, 0, false, ' ', true},
    {72, 1, ColumnKind::Literal, 0, false, ' '},
    {73, 5, ColumnKind::Real, 2},
};

constexpr ColumnSpec kTerCols[] = {
    {0, 6, ColumnKind::Skip},
    {6, 5, ColumnKind::Int},
    {11, 6, ColumnKind::Skip},
    {17, 4, ColumnKind::Str, 0, true},
    {21, 1, ColumnKind::Char},
    {22, 4, ColumnKind::Int},
    {26, 1, ColumnKind::Char},
};

constexpr ColumnSpec kTurnCols[] = {
    {0, 7, ColumnKind::Skip},
    {7, 3, ColumnKind::Int, 0, false, ' ', true},
    {10, 1, ColumnKind::Literal, 0, false, ' '},
    {11, 3, ColumnKind::Str},
    {14, 1, ColumnKind::Literal, 0, false, ' '},
    {15, 4, ColumnKind::Str, 0, true},
    {19, 1, ColumnKind::Char},
    {20, 4, ColumnKind::Int},
    {24, 1, ColumnKind::Char},
    {25, 1, ColumnKind::Literal, 0, false, ' '},
    {26, 4, ColumnKind::Str, 0, true},
    {30, 1, ColumnKind::Char},
    {31, 4, ColumnKind::Int},
    {35, 1, ColumnKind::Char},
    {36, 4, ColumnKind::Skip},
    {40, 30, ColumnKind::Str, 0, true},
};

constexpr ColumnSpec kCmpontCols[] = {
    {0, 7, ColumnKind::Skip},
    {7, 3, ColumnKind::Int},
    {10, 1, ColumnKind::Literal, 0, false, ' '},
    {11, 4, ColumnKind::Str},
    {15, 1, ColumnKind::Char},
    {16, 1, ColumnKind::Literal, 0, false, ' '},
    {17, 4, ColumnKind::Int},
    {21, 1, ColumnKind::Char},
    {22, 1, ColumnKind::Literal, 0, false, ' '},
    {23, 4, ColumnKind::Str},
    {27, 1, ColumnKind::Char},
    {28, 1, ColumnKind::Literal, 0, false, ' '},
    {29, 4, ColumnKind::Int},
    {33, 1, ColumnKind::Char},
};

constexpr ColumnSpec kConectCols[] = {
    {0, 6, ColumnKind::Skip},
    {6, 5, ColumnKind::Int},
    {11, 5, ColumnKind::Int, 0, false, ' ', true},
    {16, 5, ColumnKind::Int, 0, false, ' ', true},
    {21, 5, ColumnKind::Int, 0, false, ' ', true},
    {26, 5, ColumnKind::Int, 0, false, ' ', true},
    {31, 5, ColumnKind::Int, 0, false, ' ', true},
    {36, 5, ColumnKind::Int, 0, false, ' ', true},
    {41, 5, ColumnKind::Int, 0, false, ' ', true},
    {46, 5, ColumnKind::Int, 0, false, ' ', true},
    {51, 5, ColumnKind::Int, 0, false, ' ', true},
    {56, 5, ColumnKind::Int, 0, false, ' ', true},
};

constexpr ColumnSpec kCryst1Cols[] = {
    {0, 6, ColumnKind::Skip},
    {6, 9, ColumnKind::Real, 3},
    {15, 9, ColumnKind::Real, 3},
    {24, 9, ColumnKind::Real, 3},
    {33, 7, ColumnKind::Real, 2},
    {40, 7, ColumnKind::Real, 2},
    {47, 7, ColumnKind::Real, 2},
    {54, 1, ColumnKind::Literal, 0, false, ' '},
    {55, 11, ColumnKind::Str, 0, true},
    {66, 4, ColumnKind::Int},
};

constexpr ColumnSpec kMtrixCols[] = {
    {0, 5, ColumnKind::Skip},
    {5, 1, ColumnKind::Int},
    {6, 1, ColumnKind::Literal, 0, false, ' '},
    {7, 3, ColumnKind::Int},
    {10, 10, ColumnKind::Real, 6},
    {20, 10, ColumnKind::Real, 6},
    {30, 10, ColumnKind::Real, 6},
    {40, 5, ColumnKind::Skip},
    {45, 10, ColumnKind::Real, 5},
    {55, 3, ColumnKind::Literal, 0, false, ' '},
    {58, 2, ColumnKind::Int, 0, false, ' ', true},
};

constexpr ColumnSpec kOrigxScaleCols[] = {
    {0, 5, ColumnKind::Skip},
    {5, 1, ColumnKind::Int},
    {6, 4, ColumnKind::Skip},
    {10, 10, ColumnKind::Real, 6},
    {20, 10, ColumnKind::Real, 6},
    {30, 10, ColumnKind::Real, 6},
    {40, 5, ColumnKind::Skip},
    {45, 10, ColumnKind::Real, 5},
};

constexpr ColumnSpec kSymopCols[] = {
    {0, 5, ColumnKind::Skip},
    {5, 1, ColumnKind::Int},
    {6, 1, ColumnKind::Literal, 0, false, ' '},
    {7, 3, ColumnKind::Int},
    {10, 10, ColumnKind::Real, 6},
    {20, 10, ColumnKind::Real, 6},
    {30, 10, ColumnKind::Real, 6},
    {40, 5, ColumnKind::Skip},
    {45, 10, ColumnKind::Real, 5},
};

constexpr ColumnSpec kRevdatReadCols[] = {
    {0, 7, ColumnKind::Skip},
    {7, 3, ColumnKind::Int, 0, false, ' ', true},
    {10, 2, ColumnKind::Int, 0, false, ' ', true},
    {12, 1, ColumnKind::Literal, 0, false, ' '},
    {13, 9, ColumnKind::Str, 0, true},
    {22, 1, ColumnKind::Literal, 0, false, ' '},
    {23, 4, ColumnKind::Str},
    {27, 4, ColumnKind::Literal, 0, false, ' '},
    {31, 1, ColumnKind::Int},
    {32, 7, ColumnKind::Literal, 0, false, ' '},
    {39, 6, ColumnKind::Str, 0, true},
    {45, 1, ColumnKind::Literal, 0, false, ' '},
    {46, 6, ColumnKind::Str, 0, true},
    {52, 1, ColumnKind::Literal, 0, false, ' '},
    {53, 6, ColumnKind::Str, 0, true},
    {59, 1, ColumnKind::Literal, 0, false, ' '},
    {60, 6, ColumnKind::Str, 0, true},
};

constexpr ColumnSpec kRevdatWriteCols[] = {
    {0, 7, ColumnKind::Skip},
    {7, 3, ColumnKind::Int, 0, false, ' ', true},
    {10, 2, ColumnKind::Int, 0, false, ' ', true},
    {12, 1, ColumnKind::Literal, 0, false, ' '},
    {13, 9, ColumnKind::Str, 0, true},
    {22, 1, ColumnKind::Literal, 0, false, ' '},
    {23, 4, ColumnKind::Str},
    {27, 4, ColumnKind::Literal, 0, false, ' '},
    {31, 1, ColumnKind::Int},
    {32, 7, ColumnKind::Literal, 0, false, ' '},
    {39, 7, ColumnKind::Str, 0, true},
    {46, 7, ColumnKind::Str, 0, true},
    {53, 7, ColumnKind::Str, 0, true},
    {60, 6, ColumnKind::Str, 0, true},
};

constexpr ColumnSpec kSeqresCols[] = {
    {0, 6, ColumnKind::Skip},
    {6, 4, ColumnKind::Int},
    {10, 1, ColumnKind::Literal, 0, false, ' '},
    {11, 1, ColumnKind::Char},
    {12, 1, ColumnKind::Literal, 0, false, ' '},
    {13, 4, ColumnKind::Int},
    {17, 2, ColumnKind::Literal, 0, false, ' '},
    {19, 4, ColumnKind::Str, 0, true},
    {23, 4, ColumnKind::Str, 0, true},
    {27, 4, ColumnKind::Str, 0, true},
    {31, 4, ColumnKind::Str, 0, true},
    {35, 4, ColumnKind::Str, 0, true},
    {39, 4, ColumnKind::Str, 0, true},
    {43, 4, ColumnKind::Str, 0, true},
    {47, 4, ColumnKind::Str, 0, true},
    {51, 4, ColumnKind::Str, 0, true},
    {55, 4, ColumnKind::Str, 0, true},
    {59, 4, ColumnKind::Str, 0, true},
    {63, 4, ColumnKind::Str, 0, true},
    {67, 4, ColumnKind::Str, 0, true},
    {71, 4, ColumnKind::Str, 0, true},
};

constexpr ColumnSpec kSplitReadCols[] = {
    {0, 8, ColumnKind::Skip},
    {8, 2, ColumnKind::Int, 0, false, ' ', true},
    {10, 1, ColumnKind::Literal, 0, false, ' '},
    {11, 4, ColumnKind::Str},
    {15, 1, ColumnKind::Literal, 0, false, ' '},
    {16, 4, ColumnKind::Str},
    {20, 1, ColumnKind::Literal, 0, false, ' '},
    {21, 4, ColumnKind::Str},
    {25, 1, ColumnKind::Literal, 0, false, ' '},
    {26, 4, ColumnKind::Str},
    {30, 1, ColumnKind::Literal, 0, false, ' '},
    {31, 4, ColumnKind::Str},
    {35, 1, ColumnKind::Literal, 0, false, ' '},
    {36, 4, ColumnKind::Str},
    {40, 1, ColumnKind::Literal, 0, false, ' '},
    {41, 4, ColumnKind::Str},
    {45, 1, ColumnKind::Literal, 0, false, ' '},
    {46, 4, ColumnKind::Str},
    {50, 1, ColumnKind::Literal, 0, false, ' '},
    {51, 4, ColumnKind::Str},
    {55, 1, ColumnKind::Literal, 0, false, ' '},
    {56, 4, ColumnKind::Str},
    {60, 1, ColumnKind::Literal, 0, false, ' '},
    {61, 4, ColumnKind::Str},
    {65, 1, ColumnKind::Literal, 0, false, ' '},
    {66, 4, ColumnKind::Str},
    {70, 1, ColumnKind::Literal, 0, false, ' '},
    {71, 4, ColumnKind::Str},
    {75, 1, ColumnKind::Literal, 0, false, ' '},
    {76, 4, ColumnKind::Str},
};

constexpr ColumnSpec kSplitWriteCols[] = {
    {0, 8, ColumnKind::Skip},
    {8, 2, ColumnKind::Int, 0, false, ' ', true},
    {10, 1, ColumnKind::Literal, 0, false, ' '},
    {11, 5, ColumnKind::Str, 0, true},
    {16, 5, ColumnKind::Str, 0, true},
    {21, 5, ColumnKind::Str, 0, true},
    {26, 5, ColumnKind::Str, 0, true},
    {31, 5, ColumnKind::Str, 0, true},
    {36, 5, ColumnKind::Str, 0, true},
    {41, 5, ColumnKind::Str, 0, true},
    {46, 5, ColumnKind::Str, 0, true},
    {51, 5, ColumnKind::Str, 0, true},
    {56, 5, ColumnKind::Str, 0, true},
    {61, 5, ColumnKind::Str, 0, true},
    {66, 5, ColumnKind::Str, 0, true},
    {71, 5, ColumnKind::Str, 0, true},
    {76, 4, ColumnKind::Str, 0, true},
};

constexpr ColumnSpec kObslteReadCols[] = {
    {0, 8, ColumnKind::Skip},
    {8, 2, ColumnKind::Int, 0, false, ' ', true},
    {10, 1, ColumnKind::Literal, 0, false, ' '},
    {11, 9, ColumnKind::Str, 0, true},
    {20, 1, ColumnKind::Literal, 0, false, ' '},
    {21, 4, ColumnKind::Str},
    {25, 6, ColumnKind::Skip},
    {31, 4, ColumnKind::Str},
    {35, 1, ColumnKind::Literal, 0, false, ' '},
    {36, 4, ColumnKind::Str},
    {40, 1, ColumnKind::Literal, 0, false, ' '},
    {41, 4, ColumnKind::Str},
    {45, 1, ColumnKind::Literal, 0, false, ' '},
    {46, 4, ColumnKind::Str},
    {50, 1, ColumnKind::Literal, 0, false, ' '},
    {51, 4, ColumnKind::Str},
    {55, 1, ColumnKind::Literal, 0, false, ' '},
    {56, 4, ColumnKind::Str},
    {60, 1, ColumnKind::Literal, 0, false, ' '},
    {61, 4, ColumnKind::Str},
    {65, 1, ColumnKind::Literal, 0, false, ' '},
    {66, 4, ColumnKind::Str},
    {70, 1, ColumnKind::Literal, 0, false, ' '},
    {71, 4, ColumnKind::Str},
};

constexpr ColumnSpec kObslteWriteCols[] = {
    {0, 8, ColumnKind::Skip},
    {8, 2, ColumnKind::Int, 0, false, ' ', true},
    {10, 1, ColumnKind::Literal, 0, false, ' '},
    {11, 9, ColumnKind::Str, 0, true},
    {20, 1, ColumnKind::Literal, 0, false, ' '},
    {21, 10, ColumnKind::Str, 0, true},
    {31, 5, ColumnKind::Str, 0, true},
    {36, 5, ColumnKind::Str, 0, true},
    {41, 5, ColumnKind::Str, 0, true},
    {46, 5, ColumnKind::Str, 0, true},
    {51, 5, ColumnKind::Str, 0, true},
    {56, 5, ColumnKind::Str, 0, true},
    {61, 5, ColumnKind::Str, 0, true},
    {66, 5, ColumnKind::Str, 0, true},
    {71, 4, ColumnKind::Str, 0, true},
};

constexpr ColumnSpec kMasterCols[] = {
    {0, 10, ColumnKind::Skip},
    {10, 5, ColumnKind::Int},
    {15, 5, ColumnKind::Int},
    {20, 5, ColumnKind::Int},
    {25, 5, ColumnKind::Int},
    {30, 5, ColumnKind::Int},
    {35, 5, ColumnKind::Int},
    {40, 5, ColumnKind::Int},
    {45, 5, ColumnKind::Int},
    {50, 5, ColumnKind::Int},
    {55, 5, ColumnKind::Int},
    {60, 5, ColumnKind::Int},
    {65, 5, ColumnKind::Int},
};

constexpr ColumnSpec kAggrgtCols[] = {
    {0, 7, ColumnKind::Skip},
    {7, 3, ColumnKind::Int},
    {10, 1, ColumnKind::Literal, 0, false, ' '},
    {11, 3, ColumnKind::Int},
    {14, 1, ColumnKind::Literal, 0, false, ' '},
    {15, 3, ColumnKind::Int},
    {18, 1, ColumnKind::Literal, 0, false, ' '},
    {19, 3, ColumnKind::Int},
    {22, 1, ColumnKind::Literal, 0, false, ' '},
    {23, 3, ColumnKind::Int},
    {26, 1, ColumnKind::Literal, 0, false, ' '},
    {27, 3, ColumnKind::Int},
    {30, 1, ColumnKind::Literal, 0, false, ' '},
    {31, 3, ColumnKind::Int},
    {34, 1, ColumnKind::Literal, 0, false, ' '},
    {35, 3, ColumnKind::Int},
    {38, 1, ColumnKind::Literal, 0, false, ' '},
    {39, 3, ColumnKind::Int},
    {42, 1, ColumnKind::Literal, 0, false, ' '},
    {43, 3, ColumnKind::Int},
    {46, 1, ColumnKind::Literal, 0, false, ' '},
    {47, 3, ColumnKind::Int},
    {50, 1, ColumnKind::Literal, 0, false, ' '},
    {51, 3, ColumnKind::Int},
    {54, 1, ColumnKind::Literal, 0, false, ' '},
    {55, 3, ColumnKind::Int},
    {58, 1, ColumnKind::Literal, 0, false, ' '},
    {59, 3, ColumnKind::Int},
    {62, 1, ColumnKind::Literal, 0, false, ' '},
    {63, 3, ColumnKind::Int},
    {66, 1, ColumnKind::Literal, 0, false, ' '},
    {67, 3, ColumnKind::Int},
    {70, 1, ColumnKind::Literal, 0, false, ' '},
    {71, 3, ColumnKind::Int},
};

constexpr ColumnSpec kTrnsfmCols[] = {
    {0, 7, ColumnKind::Skip},
    {7, 3, ColumnKind::Int},
    {10, 1, ColumnKind::Literal, 0, false, ' '},
    {11, 3, ColumnKind::Int},
    {14, 1, ColumnKind::Literal, 0, false, ' '},
    {15, 3, ColumnKind::Int},
};

constexpr ColumnSpec kTvectCols[] = {
    {0, 7, ColumnKind::Skip},
    {7, 3, ColumnKind::Int, 0, false, ' ', true},
    {10, 10, ColumnKind::Real, 5},
    {20, 10, ColumnKind::Real, 5},
    {30, 10, ColumnKind::Real, 5},
    {40, 30, ColumnKind::Str, 0, true},
};

constexpr ColumnSpec kHetCols[] = {
    {0, 7, ColumnKind::Skip},
    {7, 4, ColumnKind::Str, 0, true},
    {11, 1, ColumnKind::Literal, 0, false, ' '},
    {12, 1, ColumnKind::Char},
    {13, 4, ColumnKind::Int},
    {17, 1, ColumnKind::Char},
    {18, 2, ColumnKind::Literal, 0, false, ' '},
    {20, 5, ColumnKind::Int},
    {25, 5, ColumnKind::Skip},
    {30, 40, ColumnKind::Str, 0, true},
};

constexpr ColumnSpec kFormulCols[] = {
    {0, 8, ColumnKind::Skip},
    {8, 2, ColumnKind::Int, 0, false, ' ', true},
    {10, 2, ColumnKind::Literal, 0, false, ' '},
    {12, 4, ColumnKind::Str, 0, true},
    {16, 2, ColumnKind::Int, 0, false, ' ', true},
    {18, 1, ColumnKind::Char},
    {19, 51, ColumnKind::Str, 0, true},
};

constexpr ColumnSpec kHeaderCols[] = {
    {0, 10, ColumnKind::Skip},
    {10, 40, ColumnKind::Str, 0, true},
    {50, 9, ColumnKind::Str},
    {59, 3, ColumnKind::Literal, 0, false, ' '},
    {62, 4, ColumnKind::Str},
};

// HEADER writes a different geometry than it reads: the legacy writer
// renders the date left-justified in 11 columns followed by a literal
// space (so a short or spilled date sits at the left of its field),
// while the reader slices a 9-column date. Well-formed 9-character
// dates render identically either way; an overlong classification that
// spills into the date columns does not.
constexpr ColumnSpec kHeaderWriteCols[] = {
    {0, 10, ColumnKind::Skip},
    {10, 40, ColumnKind::Str, 0, true},
    {50, 11, ColumnKind::Str, 0, true},
    {61, 1, ColumnKind::Literal, 0, false, ' '},
    {62, 4, ColumnKind::Str, 0, true},
};

constexpr ColumnSpec kCaveatCols[] = {
    {0, 8, ColumnKind::Skip},
    {8, 2, ColumnKind::Int, 0, false, ' ', true},
    {10, 1, ColumnKind::Literal, 0, false, ' '},
    {11, 4, ColumnKind::Str},
    {15, 4, ColumnKind::Literal, 0, false, ' '},
    {19, 60, ColumnKind::Str, 0, true},
};

constexpr ColumnSpec kJrnlCols[] = {
    {0, 12, ColumnKind::Skip},
    {12, 67, ColumnKind::Str, 0, true},
};

// The description family reads identically but writes in two styles:
// FTNOTE and REMARK render %3D (zero as blanks) and a left-justified
// text, while AGRDES, SYMDES, MTXDES, and CMPDES render a plain %3d and
// a right-justified 68-column text.
constexpr ColumnSpec kAgrdesCols[] = {
    {0, 7, ColumnKind::Skip},
    {7, 3, ColumnKind::Int, 0, false, ' ', true},
    {10, 1, ColumnKind::Literal, 0, false, ' '},
    {11, 68, ColumnKind::Str, 0, true},
};

constexpr ColumnSpec kAgrdesPlainWriteCols[] = {
    {0, 7, ColumnKind::Skip},
    {7, 3, ColumnKind::Int},
    {10, 1, ColumnKind::Literal, 0, false, ' '},
    {11, 68, ColumnKind::Str},
};

// ANISOU and SIGUIJ share the atom-layout identification columns and the
// six %7D thermal tensors (zero rendered as blanks on write).
constexpr ColumnSpec kAnisouCols[] = {
    {0, 6, ColumnKind::Skip},
    {6, 5, ColumnKind::Int},
    {11, 1, ColumnKind::Literal, 0, false, ' '},
    {12, 4, ColumnKind::Str, 0, true},
    {16, 1, ColumnKind::Char},
    {17, 4, ColumnKind::Str, 0, true},
    {21, 1, ColumnKind::Char},
    {22, 4, ColumnKind::Int},
    {26, 1, ColumnKind::Char},
    {27, 1, ColumnKind::Literal, 0, false, ' '},
    {28, 7, ColumnKind::Int},
    {35, 7, ColumnKind::Int},
    {42, 7, ColumnKind::Int},
    {49, 7, ColumnKind::Int},
    {56, 7, ColumnKind::Int},
    {63, 7, ColumnKind::Int},
};

// SIGUIJ writes the same layout as ANISOU except that the u columns use
// the zero-as-blank emitter (%7D), as the legacy write table does.
constexpr ColumnSpec kSiguijWriteCols[] = {
    {0, 6, ColumnKind::Skip},
    {6, 5, ColumnKind::Int},
    {11, 1, ColumnKind::Literal, 0, false, ' '},
    {12, 4, ColumnKind::Str, 0, true},
    {16, 1, ColumnKind::Char},
    {17, 4, ColumnKind::Str, 0, true},
    {21, 1, ColumnKind::Char},
    {22, 4, ColumnKind::Int},
    {26, 1, ColumnKind::Char},
    {27, 1, ColumnKind::Literal, 0, false, ' '},
    {28, 7, ColumnKind::Int, 0, false, ' ', true},
    {35, 7, ColumnKind::Int, 0, false, ' ', true},
    {42, 7, ColumnKind::Int, 0, false, ' ', true},
    {49, 7, ColumnKind::Int, 0, false, ' ', true},
    {56, 7, ColumnKind::Int, 0, false, ' ', true},
    {63, 7, ColumnKind::Int, 0, false, ' ', true},
};

// MODEL reads its number from columns 7-14 (non-standard layout used by
// the repository fixtures) but writes %5d at columns 10-14; the pinned
// quirk round-trips because both parse to the same value.
constexpr ColumnSpec kModelReadCols[] = {
    {0, 6, ColumnKind::Skip},
    {6, 8, ColumnKind::Int},
};

constexpr ColumnSpec kModelWriteCols[] = {
    {0, 9, ColumnKind::Skip},
    {9, 5, ColumnKind::Int},
};

constexpr ColumnSpec kNummdlCols[] = {
    {0, 10, ColumnKind::Skip},
    {10, 4, ColumnKind::Int},
};

constexpr ColumnSpec kNummdlWriteCols[] = {
    {0, 10, ColumnKind::Skip},
    {10, 4, ColumnKind::Int, 0, true},
};

constexpr ColumnSpec kExpdtaCols[] = {
    {0, 8, ColumnKind::Skip},
    {8, 2, ColumnKind::Int, 0, false, ' ', true},
    {10, 69, ColumnKind::Str, 0, true},
};

// Continuation-plus-text families: one typed int and one string, bound to
// the payload members directly (no offsets, no varargs). The record-name
// region before the first column is skipped, its width taken from the
// first column's offset.
template <typename Payload>
bool parseContText(ColumnReader &reader, Record &record,
                   const ColumnSpec *cols, int Payload::*continuation,
                   std::string Payload::*text) {
  reader.skipChars(cols[0].offset);
  std::int64_t value = 0;
  if (!reader.readInt(cols[0].width, &value))
    return false;
  Payload &payload = record.as<Payload>();
  payload.*continuation = static_cast<int>(value);
  payload.*text = reader.readString(cols[1].width);
  return true;
}

template <typename Payload>
void writeContText(const Record &record, ColumnWriter &writer,
                   const ColumnSpec *cols, int Payload::*continuation,
                   std::string Payload::*text) {
  writer.writeInt(record.as<Payload>().*continuation, cols[0]);
  writer.writeString(record.as<Payload>().*text, cols[1]);
}

// Asserts the write-side prefix exactly fills the region before the first
// column of its family.
void assertPrefix(const char *prefix, const ColumnSpec *cols) {
  assert(std::strlen(prefix) ==
         (cols[0].kind == ColumnKind::Skip ? cols[0].width : cols[0].offset));
}

// The atom payload binder: reads the standard layout columns in table
// order into the typed members. Occupancy and tempFactor are set only
// when their columns parsed or defaulted.
bool parseAtomColumns(ColumnReader &reader, Record &record,
                      const ColumnSpec *cols) {
  reader.skipChars(cols[0].width);
  std::int64_t serial = 0;
  if (!reader.readInt(cols[1].width, &serial))
    return false;
  AtomRecord &atom = record.as<AtomRecord>();
  atom.serialNum = static_cast<int>(serial);
  if (!reader.matchLiterals(cols[2].fill, cols[2].width))
    return false;
  atom.name = reader.readString(cols[3].width);
  atom.altLoc = reader.readChars(cols[4].width)[0];
  atom.residue.name = reader.readString(cols[5].width);
  atom.residue.chainId = reader.readChars(cols[6].width)[0];
  std::int64_t seq = 0;
  if (!reader.readInt(cols[7].width, &seq))
    return false;
  atom.residue.seqNum = static_cast<int>(seq);
  atom.residue.insertCode = reader.readChars(cols[8].width)[0];
  if (!reader.matchLiterals(cols[9].fill, cols[9].width))
    return false;
  double value = 0.0;
  if (!reader.readReal(cols[10].width, &value))
    return false;
  atom.xyz[0] = value;
  if (!reader.readReal(cols[11].width, &value))
    return false;
  atom.xyz[1] = value;
  if (!reader.readReal(cols[12].width, &value))
    return false;
  atom.xyz[2] = value;
  if (!reader.readReal(cols[13].width, &value))
    return false;
  atom.occupancy = value;
  if (!reader.readReal(cols[14].width, &value))
    return false;
  atom.tempFactor = value;
  if (!reader.matchLiterals(cols[15].fill, cols[15].width))
    return false;
  atom.element = reader.readString(cols[16].width);
  atom.charge = reader.readString(cols[17].width);
  return true;
}

// The ZDATOM payload binder: fills the mark_sur extension on a fresh
// record, so a failed standard attempt leaves nothing behind.
bool parseZdatomColumns(ColumnReader &reader, Record &record,
                        const ColumnSpec *cols) {
  reader.skipChars(cols[0].width);
  std::int64_t serial = 0;
  if (!reader.readInt(cols[1].width, &serial))
    return false;
  AtomRecord &atom = record.as<AtomRecord>();
  atom.serialNum = static_cast<int>(serial);
  if (!reader.matchLiterals(cols[2].fill, cols[2].width))
    return false;
  atom.name = reader.readString(cols[3].width);
  atom.altLoc = reader.readChars(cols[4].width)[0];
  atom.residue.name = reader.readString(cols[5].width);
  atom.residue.chainId = reader.readChars(cols[6].width)[0];
  std::int64_t seq = 0;
  if (!reader.readInt(cols[7].width, &seq))
    return false;
  atom.residue.seqNum = static_cast<int>(seq);
  atom.residue.insertCode = reader.readChars(cols[8].width)[0];
  if (!reader.matchLiterals(cols[9].fill, cols[9].width))
    return false;
  double value = 0.0;
  if (!reader.readReal(cols[10].width, &value))
    return false;
  atom.xyz[0] = value;
  if (!reader.readReal(cols[11].width, &value))
    return false;
  atom.xyz[1] = value;
  if (!reader.readReal(cols[12].width, &value))
    return false;
  atom.xyz[2] = value;
  if (!reader.matchLiterals(cols[13].fill, cols[13].width))
    return false;
  Zdatom &zdatom = atom.zdatom.emplace();
  std::int64_t small = 0;
  if (!reader.readInt(cols[14].width, &small))
    return false;
  zdatom.type = static_cast<int>(small);
  if (!reader.matchLiterals(cols[15].fill, cols[15].width))
    return false;
  if (!reader.readInt(cols[16].width, &small))
    return false;
  zdatom.surface = static_cast<int>(small);
  if (!reader.matchLiterals(cols[17].fill, cols[17].width))
    return false;
  if (!reader.readReal(cols[18].width, &value))
    return false;
  zdatom.rad = value;
  if (!reader.matchLiterals(cols[19].fill, cols[19].width))
    return false;
  zdatom.segid = reader.readString(cols[20].width);
  if (!reader.matchLiterals(cols[21].fill, cols[21].width))
    return false;
  if (!reader.readReal(cols[22].width, &value))
    return false;
  zdatom.chg = value;
  return true;
}

void writeAtomColumns(const AtomRecord &atom, ColumnWriter &writer,
                      const ColumnSpec *cols) {
  writer.writeInt(atom.serialNum, cols[1]);
  writer.writeString(atom.name, cols[3]);
  writer.writeChar(atom.altLoc, cols[4]);
  writer.writeString(atom.residue.name, cols[5]);
  writer.writeChar(atom.residue.chainId, cols[6]);
  writer.writeInt(atom.residue.seqNum, cols[7]);
  writer.writeChar(atom.residue.insertCode, cols[8]);
  writer.writeReal(atom.xyz[0], cols[10]);
  writer.writeReal(atom.xyz[1], cols[11]);
  writer.writeReal(atom.xyz[2], cols[12]);
  writer.writeReal(atom.occupancy.value_or(0.0), cols[13]);
  writer.writeReal(atom.tempFactor.value_or(0.0), cols[14]);
  writer.writeString(atom.element, cols[16]);
  writer.writeString(atom.charge, cols[17]);
}

void writeZdatomColumns(const AtomRecord &atom, ColumnWriter &writer,
                        const ColumnSpec *cols) {
  const Zdatom &zdatom = *atom.zdatom;
  writer.writeInt(atom.serialNum, cols[1]);
  writer.writeString(atom.name, cols[3]);
  writer.writeChar(atom.altLoc, cols[4]);
  writer.writeString(atom.residue.name, cols[5]);
  writer.writeChar(atom.residue.chainId, cols[6]);
  writer.writeInt(atom.residue.seqNum, cols[7]);
  writer.writeChar(atom.residue.insertCode, cols[8]);
  writer.writeReal(atom.xyz[0], cols[10]);
  writer.writeReal(atom.xyz[1], cols[11]);
  writer.writeReal(atom.xyz[2], cols[12]);
  writer.writeInt(zdatom.type, cols[14]);
  writer.writeInt(zdatom.surface, cols[16]);
  writer.writeReal(zdatom.rad, cols[18]);
  writer.writeString(zdatom.segid, cols[20]);
  writer.writeReal(zdatom.chg, cols[22]);
}

// Small typed helpers for the residue-heavy binders.
bool readResidue(ColumnReader &reader, const ColumnSpec *cols,
                 std::size_t nameIndex, Residue *residue) {
  residue->name = reader.readString(cols[nameIndex].width);
  residue->chainId = reader.readChars(cols[nameIndex + 1].width)[0];
  std::int64_t seq = 0;
  if (!reader.readInt(cols[nameIndex + 2].width, &seq))
    return false;
  residue->seqNum = static_cast<int>(seq);
  residue->insertCode = reader.readChars(cols[nameIndex + 3].width)[0];
  return true;
}

void writeResidue(const Residue &residue, ColumnWriter &writer,
                  const ColumnSpec *cols, std::size_t nameIndex) {
  writer.writeString(residue.name, cols[nameIndex]);
  writer.writeChar(residue.chainId, cols[nameIndex + 1]);
  writer.writeInt(residue.seqNum, cols[nameIndex + 2]);
  writer.writeChar(residue.insertCode, cols[nameIndex + 3]);
}

#define RETURN_ON_LITERAL(reader, cols, index)                                \
  do {                                                                        \
    if (!(reader).matchLiterals((cols)[index].fill, (cols)[index].width))     \
      return false;                                                                         \
  } while (false)

#define RETURN_ON_INT(reader, cols, index, out)                               \
  do {                                                                        \
    std::int64_t macroValue = 0;                                              \
    if (!(reader).readInt((cols)[index].width, &macroValue))                  \
      return false;                                                           \
    *(out) = static_cast<int>(macroValue);                                    \
  } while (false)

#define RETURN_ON_REAL(reader, cols, index, out)                              \
  do {                                                                        \
    double macroValue = 0.0;                                                  \
    if (!(reader).readReal((cols)[index].width, &macroValue))                 \
      return false;                                                           \
    *(out) = macroValue;                                                      \
  } while (false)

// The HELIX residue groups carry the separator literal between the chain
// id and the sequence number (%4s%c %4d%c), so the binder addresses the
// five columns explicitly.
bool parseHelixColumns(ColumnReader &reader, Record &record,
                       const ColumnSpec *cols) {
  reader.skipChars(cols[0].width);
  HelixRecord &h = record.as<HelixRecord>();
  RETURN_ON_INT(reader, cols, 1, &h.serialNum);
  RETURN_ON_LITERAL(reader, cols, 2);
  h.id = reader.readString(cols[3].width);
  RETURN_ON_LITERAL(reader, cols, 4);
  for (std::size_t i = 0; i < 2; ++i) {
    if (i > 0) // the separator between the two residue groups
      RETURN_ON_LITERAL(reader, cols, 10);
    const std::size_t base = 5 + i * 6; // five columns plus the separator
    h.residues[i].name = reader.readString(cols[base].width);
    h.residues[i].chainId = reader.readChars(cols[base + 1].width)[0];
    RETURN_ON_LITERAL(reader, cols, base + 2);
    std::int64_t seq = 0;
    if (!reader.readInt(cols[base + 3].width, &seq))
      return false;
    h.residues[i].seqNum = static_cast<int>(seq);
    h.residues[i].insertCode = reader.readChars(cols[base + 4].width)[0];
  }
  RETURN_ON_INT(reader, cols, 16, &h.type);
  h.comment = reader.readString(cols[17].width);
  RETURN_ON_LITERAL(reader, cols, 18);
  RETURN_ON_INT(reader, cols, 19, &h.length);
  return true;
}

void writeHelixColumns(const HelixRecord &h, ColumnWriter &writer,
                       const ColumnSpec *cols) {
  writer.writeInt(h.serialNum, cols[1]);
  writer.writeString(h.id, cols[3]);
  for (std::size_t i = 0; i < 2; ++i) {
    const std::size_t base = 5 + i * 6; // five columns plus the separator
    writer.writeString(h.residues[i].name, cols[base]);
    writer.writeChar(h.residues[i].chainId, cols[base + 1]);
    writer.writeInt(h.residues[i].seqNum, cols[base + 3]);
    writer.writeChar(h.residues[i].insertCode, cols[base + 4]);
  }
  writer.writeInt(h.type, cols[16]);
  writer.writeString(h.comment, cols[17]);
  writer.writeInt(h.length, cols[19]);
}

bool parseSheetColumns(ColumnReader &reader, Record &record,
                       const ColumnSpec *cols) {
  reader.skipChars(cols[0].width);
  SheetRecord &s = record.as<SheetRecord>();
  RETURN_ON_INT(reader, cols, 1, &s.strandNum);
  RETURN_ON_LITERAL(reader, cols, 2);
  s.id = reader.readString(cols[3].width);
  RETURN_ON_INT(reader, cols, 4, &s.count);
  RETURN_ON_LITERAL(reader, cols, 5);
  if (!readResidue(reader, cols, 6, &s.residues[0]))
    return false;
  RETURN_ON_LITERAL(reader, cols, 10);
  if (!readResidue(reader, cols, 11, &s.residues[1]))
    return false;
  RETURN_ON_INT(reader, cols, 15, &s.sense);
  RETURN_ON_LITERAL(reader, cols, 16);
  s.atoms[0].name = reader.readString(cols[17].width);
  if (!readResidue(reader, cols, 18, &s.atoms[0].residue))
    return false;
  RETURN_ON_LITERAL(reader, cols, 22);
  s.atoms[1].name = reader.readString(cols[23].width);
  if (!readResidue(reader, cols, 24, &s.atoms[1].residue))
    return false;
  return true;
}

void writeSheetColumns(const SheetRecord &s, ColumnWriter &writer,
                       const ColumnSpec *cols) {
  writer.writeInt(s.strandNum, cols[1]);
  writer.writeString(s.id, cols[3]);
  writer.writeInt(s.count, cols[4]);
  writeResidue(s.residues[0], writer, cols, 6);
  writeResidue(s.residues[1], writer, cols, 11);
  writer.writeInt(s.sense, cols[15]);
  writer.writeString(s.atoms[0].name, cols[17]);
  writeResidue(s.atoms[0].residue, writer, cols, 18);
  writer.writeString(s.atoms[1].name, cols[23]);
  writeResidue(s.atoms[1].residue, writer, cols, 24);
}

bool parseSiteColumns(ColumnReader &reader, Record &record,
                      const ColumnSpec *cols) {
  reader.skipChars(cols[0].width);
  SiteRecord &s = record.as<SiteRecord>();
  RETURN_ON_INT(reader, cols, 1, &s.seqNum);
  RETURN_ON_LITERAL(reader, cols, 2);
  s.id = reader.readString(cols[3].width);
  RETURN_ON_LITERAL(reader, cols, 4);
  RETURN_ON_INT(reader, cols, 5, &s.count);
  for (std::size_t i = 0; i < 4; ++i) {
    const std::size_t base = 6 + i * 5;
    RETURN_ON_LITERAL(reader, cols, base);
    if (!readResidue(reader, cols, base + 1, &s.residues[i]))
      return false;
  }
  return true;
}

void writeSiteColumns(const SiteRecord &s, ColumnWriter &writer,
                      const ColumnSpec *cols) {
  writer.writeInt(s.seqNum, cols[1]);
  writer.writeString(s.id, cols[3]);
  writer.writeInt(s.count, cols[5]);
  for (std::size_t i = 0; i < 4; ++i)
    writeResidue(s.residues[i], writer, cols, 7 + i * 5);
}

// SSBOND residue groups carry the separator literal between the chain id
// and the sequence number, as HELIX does.
bool parseSsbondColumns(ColumnReader &reader, Record &record,
                        const ColumnSpec *cols) {
  reader.skipChars(cols[0].width);
  SsbondRecord &s = record.as<SsbondRecord>();
  RETURN_ON_INT(reader, cols, 1, &s.seqNum);
  RETURN_ON_LITERAL(reader, cols, 2);
  for (std::size_t i = 0; i < 2; ++i) {
    if (i > 0) // the separator between the two residue groups
      RETURN_ON_LITERAL(reader, cols, 8);
    const std::size_t base = 3 + i * 6; // five columns plus the separator
    s.residues[i].name = reader.readString(cols[base].width);
    s.residues[i].chainId = reader.readChars(cols[base + 1].width)[0];
    RETURN_ON_LITERAL(reader, cols, base + 2);
    std::int64_t seq = 0;
    if (!reader.readInt(cols[base + 3].width, &seq))
      return false;
    s.residues[i].seqNum = static_cast<int>(seq);
    s.residues[i].insertCode = reader.readChars(cols[base + 4].width)[0];
  }
  reader.skipChars(cols[14].width);
  RETURN_ON_INT(reader, cols, 15, &s.sym1);
  RETURN_ON_LITERAL(reader, cols, 16);
  RETURN_ON_INT(reader, cols, 17, &s.sym2);
  RETURN_ON_LITERAL(reader, cols, 18);
  RETURN_ON_REAL(reader, cols, 19, &s.length);
  return true;
}

void writeSsbondColumns(const SsbondRecord &s, ColumnWriter &writer,
                        const ColumnSpec *cols) {
  writer.writeInt(s.seqNum, cols[1]);
  for (std::size_t i = 0; i < 2; ++i) {
    const std::size_t base = 3 + i * 6; // five columns plus the separator
    writer.writeString(s.residues[i].name, cols[base]);
    writer.writeChar(s.residues[i].chainId, cols[base + 1]);
    writer.writeInt(s.residues[i].seqNum, cols[base + 3]);
    writer.writeChar(s.residues[i].insertCode, cols[base + 4]);
  }
  writer.writeInt(s.sym1, cols[15]);
  writer.writeInt(s.sym2, cols[17]);
  writer.writeReal(s.length, cols[19]);
}

bool parseTerColumns(ColumnReader &reader, Record &record,
                     const ColumnSpec *cols) {
  reader.skipChars(cols[0].width);
  TerRecord &t = record.as<TerRecord>();
  RETURN_ON_INT(reader, cols, 1, &t.serialNum);
  reader.skipChars(cols[2].width);
  if (!readResidue(reader, cols, 3, &t.residue))
    return false;
  return true;
}

void writeTerColumns(const TerRecord &t, ColumnWriter &writer,
                     const ColumnSpec *cols) {
  writer.writeInt(t.serialNum, cols[1]);
  writeResidue(t.residue, writer, cols, 3);
}

bool parseTurnColumns(ColumnReader &reader, Record &record,
                      const ColumnSpec *cols) {
  reader.skipChars(cols[0].width);
  TurnRecord &t = record.as<TurnRecord>();
  RETURN_ON_INT(reader, cols, 1, &t.seqNum);
  RETURN_ON_LITERAL(reader, cols, 2);
  t.id = reader.readString(cols[3].width);
  RETURN_ON_LITERAL(reader, cols, 4);
  if (!readResidue(reader, cols, 5, &t.residues[0]))
    return false;
  RETURN_ON_LITERAL(reader, cols, 9);
  if (!readResidue(reader, cols, 10, &t.residues[1]))
    return false;
  reader.skipChars(cols[14].width);
  t.comment = reader.readString(cols[15].width);
  return true;
}

void writeTurnColumns(const TurnRecord &t, ColumnWriter &writer,
                      const ColumnSpec *cols) {
  writer.writeInt(t.seqNum, cols[1]);
  writer.writeString(t.id, cols[3]);
  writeResidue(t.residues[0], writer, cols, 5);
  writeResidue(t.residues[1], writer, cols, 10);
  writer.writeString(t.comment, cols[15]);
}

// CMPONT residue groups carry the separator literal between the chain id
// and the sequence number, as HELIX does.
bool parseCmpontColumns(ColumnReader &reader, Record &record,
                        const ColumnSpec *cols) {
  reader.skipChars(cols[0].width);
  CmpontRecord &c = record.as<CmpontRecord>();
  RETURN_ON_INT(reader, cols, 1, &c.seqNum);
  RETURN_ON_LITERAL(reader, cols, 2);
  for (std::size_t i = 0; i < 2; ++i) {
    if (i > 0) // the separator between the two residue groups
      RETURN_ON_LITERAL(reader, cols, 8);
    const std::size_t base = 3 + i * 6; // five columns plus the separator
    c.residues[i].name = reader.readString(cols[base].width);
    c.residues[i].chainId = reader.readChars(cols[base + 1].width)[0];
    RETURN_ON_LITERAL(reader, cols, base + 2);
    std::int64_t seq = 0;
    if (!reader.readInt(cols[base + 3].width, &seq))
      return false;
    c.residues[i].seqNum = static_cast<int>(seq);
    c.residues[i].insertCode = reader.readChars(cols[base + 4].width)[0];
  }
  return true;
}

void writeCmpontColumns(const CmpontRecord &c, ColumnWriter &writer,
                        const ColumnSpec *cols) {
  writer.writeInt(c.seqNum, cols[1]);
  for (std::size_t i = 0; i < 2; ++i) {
    const std::size_t base = 3 + i * 6; // five columns plus the separator
    writer.writeString(c.residues[i].name, cols[base]);
    writer.writeChar(c.residues[i].chainId, cols[base + 1]);
    writer.writeInt(c.residues[i].seqNum, cols[base + 3]);
    writer.writeChar(c.residues[i].insertCode, cols[base + 4]);
  }
}

bool parseConectColumns(ColumnReader &reader, Record &record,
                        const ColumnSpec *cols) {
  reader.skipChars(cols[0].width);
  ConectRecord &c = record.as<ConectRecord>();
  RETURN_ON_INT(reader, cols, 1, &c.serialNum);
  // columns 2-5: covalent[0..3]; columns 6-7: bonds[0].hydrogen[0..1];
  // column 8: bonds[0].salt; columns 9-10: bonds[1].hydrogen[0..1];
  // column 11: bonds[1].salt
  for (std::size_t i = 2; i <= 5; ++i)
    RETURN_ON_INT(reader, cols, i, &c.covalent[i - 2]);
  RETURN_ON_INT(reader, cols, 6, &c.bonds[0].hydrogen[0]);
  RETURN_ON_INT(reader, cols, 7, &c.bonds[0].hydrogen[1]);
  RETURN_ON_INT(reader, cols, 8, &c.bonds[0].salt);
  RETURN_ON_INT(reader, cols, 9, &c.bonds[1].hydrogen[0]);
  RETURN_ON_INT(reader, cols, 10, &c.bonds[1].hydrogen[1]);
  RETURN_ON_INT(reader, cols, 11, &c.bonds[1].salt);
  return true;
}

void writeConectColumns(const ConectRecord &c, ColumnWriter &writer,
                        const ColumnSpec *cols) {
  writer.writeInt(c.serialNum, cols[1]);
  for (std::size_t i = 0; i < 4; ++i)
    writer.writeInt(c.covalent[i], cols[2 + i]);
  writer.writeInt(c.bonds[0].hydrogen[0], cols[6]);
  writer.writeInt(c.bonds[0].hydrogen[1], cols[7]);
  writer.writeInt(c.bonds[0].salt, cols[8]);
  writer.writeInt(c.bonds[1].hydrogen[0], cols[9]);
  writer.writeInt(c.bonds[1].hydrogen[1], cols[10]);
  writer.writeInt(c.bonds[1].salt, cols[11]);
}

bool parseCryst1Columns(ColumnReader &reader, Record &record,
                        const ColumnSpec *cols) {
  reader.skipChars(cols[0].width);
  Cryst1Record &c = record.as<Cryst1Record>();
  RETURN_ON_REAL(reader, cols, 1, &c.a);
  RETURN_ON_REAL(reader, cols, 2, &c.b);
  RETURN_ON_REAL(reader, cols, 3, &c.c);
  RETURN_ON_REAL(reader, cols, 4, &c.alpha);
  RETURN_ON_REAL(reader, cols, 5, &c.beta);
  RETURN_ON_REAL(reader, cols, 6, &c.gamma);
  RETURN_ON_LITERAL(reader, cols, 7);
  c.spaceGrp = reader.readString(cols[8].width);
  RETURN_ON_INT(reader, cols, 9, &c.z);
  return true;
}

void writeCryst1Columns(const Cryst1Record &c, ColumnWriter &writer,
                        const ColumnSpec *cols) {
  writer.writeReal(c.a, cols[1]);
  writer.writeReal(c.b, cols[2]);
  writer.writeReal(c.c, cols[3]);
  writer.writeReal(c.alpha, cols[4]);
  writer.writeReal(c.beta, cols[5]);
  writer.writeReal(c.gamma, cols[6]);
  writer.writeString(c.spaceGrp, cols[8]);
  writer.writeInt(c.z, cols[9]);
}

bool parseMtrixColumns(ColumnReader &reader, Record &record,
                       const ColumnSpec *cols) {
  reader.skipChars(cols[0].width);
  MtrixRecord &m = record.as<MtrixRecord>();
  RETURN_ON_INT(reader, cols, 1, &m.rowNum);
  RETURN_ON_LITERAL(reader, cols, 2);
  RETURN_ON_INT(reader, cols, 3, &m.serialNum);
  RETURN_ON_REAL(reader, cols, 4, &m.m1);
  RETURN_ON_REAL(reader, cols, 5, &m.m2);
  RETURN_ON_REAL(reader, cols, 6, &m.m3);
  reader.skipChars(cols[7].width);
  RETURN_ON_REAL(reader, cols, 8, &m.v);
  RETURN_ON_LITERAL(reader, cols, 9);
  RETURN_ON_INT(reader, cols, 10, &m.given);
  return true;
}

void writeMtrixColumns(const MtrixRecord &m, ColumnWriter &writer,
                       const ColumnSpec *cols) {
  writer.writeInt(m.rowNum, cols[1]);
  writer.writeInt(m.serialNum, cols[3]);
  writer.writeReal(m.m1, cols[4]);
  writer.writeReal(m.m2, cols[5]);
  writer.writeReal(m.m3, cols[6]);
  writer.writeReal(m.v, cols[8]);
  writer.writeInt(m.given, cols[10]);
}

bool parseOrigxScaleColumns(ColumnReader &reader, Record &record,
                            const ColumnSpec *cols, bool scale) {
  reader.skipChars(cols[0].width);
  int rowNum = 0;
  RETURN_ON_INT(reader, cols, 1, &rowNum);
  reader.skipChars(cols[2].width);
  double v1 = 0.0, v2 = 0.0, v3 = 0.0, v4 = 0.0;
  RETURN_ON_REAL(reader, cols, 3, &v1);
  RETURN_ON_REAL(reader, cols, 4, &v2);
  RETURN_ON_REAL(reader, cols, 5, &v3);
  reader.skipChars(cols[6].width);
  RETURN_ON_REAL(reader, cols, 7, &v4);
  if (scale) {
    ScaleRecord &s = record.as<ScaleRecord>();
    s.rowNum = rowNum;
    s.s1 = v1;
    s.s2 = v2;
    s.s3 = v3;
    s.u = v4;
  } else {
    OrigxRecord &o = record.as<OrigxRecord>();
    o.rowNum = rowNum;
    o.o1 = v1;
    o.o2 = v2;
    o.o3 = v3;
    o.t = v4;
  }
  return true;
}

void writeOrigxScaleColumns(const Record &record, ColumnWriter &writer,
                            const ColumnSpec *cols, bool scale) {
  if (scale) {
    const ScaleRecord &s = record.as<ScaleRecord>();
    writer.writeInt(s.rowNum, cols[1]);
    writer.writeReal(s.s1, cols[3]);
    writer.writeReal(s.s2, cols[4]);
    writer.writeReal(s.s3, cols[5]);
    writer.writeReal(s.u, cols[7]);
  } else {
    const OrigxRecord &o = record.as<OrigxRecord>();
    writer.writeInt(o.rowNum, cols[1]);
    writer.writeReal(o.o1, cols[3]);
    writer.writeReal(o.o2, cols[4]);
    writer.writeReal(o.o3, cols[5]);
    writer.writeReal(o.t, cols[7]);
  }
}

bool parseSymopColumns(ColumnReader &reader, Record &record,
                       const ColumnSpec *cols) {
  reader.skipChars(cols[0].width);
  SymopRecord &s = record.as<SymopRecord>();
  RETURN_ON_INT(reader, cols, 1, &s.rowNum);
  RETURN_ON_LITERAL(reader, cols, 2);
  RETURN_ON_INT(reader, cols, 3, &s.serialNum);
  RETURN_ON_REAL(reader, cols, 4, &s.s1);
  RETURN_ON_REAL(reader, cols, 5, &s.s2);
  RETURN_ON_REAL(reader, cols, 6, &s.s3);
  reader.skipChars(cols[7].width);
  RETURN_ON_REAL(reader, cols, 8, &s.t);
  return true;
}

void writeSymopColumns(const SymopRecord &s, ColumnWriter &writer,
                       const ColumnSpec *cols) {
  writer.writeInt(s.rowNum, cols[1]);
  writer.writeInt(s.serialNum, cols[3]);
  writer.writeReal(s.s1, cols[4]);
  writer.writeReal(s.s2, cols[5]);
  writer.writeReal(s.s3, cols[6]);
  writer.writeReal(s.t, cols[8]);
}

bool parseRevdatColumns(ColumnReader &reader, Record &record,
                        const ColumnSpec *cols) {
  reader.skipChars(cols[0].width);
  RevdatRecord &r = record.as<RevdatRecord>();
  RETURN_ON_INT(reader, cols, 1, &r.modification);
  RETURN_ON_INT(reader, cols, 2, &r.continuation);
  RETURN_ON_LITERAL(reader, cols, 3);
  r.timestamp = reader.readString(cols[4].width);
  RETURN_ON_LITERAL(reader, cols, 5);
  r.id = reader.readString(cols[6].width);
  RETURN_ON_LITERAL(reader, cols, 7);
  RETURN_ON_INT(reader, cols, 8, &r.modType);
  RETURN_ON_LITERAL(reader, cols, 9);
  r.record[0] = reader.readString(cols[10].width);
  RETURN_ON_LITERAL(reader, cols, 11);
  r.record[1] = reader.readString(cols[12].width);
  RETURN_ON_LITERAL(reader, cols, 13);
  r.record[2] = reader.readString(cols[14].width);
  RETURN_ON_LITERAL(reader, cols, 15);
  r.record[3] = reader.readString(cols[16].width);
  return true;
}

void writeRevdatColumns(const RevdatRecord &r, ColumnWriter &writer,
                        const ColumnSpec *cols) {
  writer.writeInt(r.modification, cols[1]);
  writer.writeInt(r.continuation, cols[2]);
  writer.writeString(r.timestamp, cols[4]);
  writer.writeString(r.id, cols[6]);
  writer.writeInt(r.modType, cols[8]);
  writer.writeString(r.record[0], cols[10]);
  writer.writeString(r.record[1], cols[11]);
  writer.writeString(r.record[2], cols[12]);
  writer.writeString(r.record[3], cols[13]);
}

bool parseSeqresColumns(ColumnReader &reader, Record &record,
                        const ColumnSpec *cols) {
  reader.skipChars(cols[0].width);
  SeqresRecord &s = record.as<SeqresRecord>();
  RETURN_ON_INT(reader, cols, 1, &s.serialNum);
  RETURN_ON_LITERAL(reader, cols, 2);
  s.chainId = reader.readChars(cols[3].width)[0];
  RETURN_ON_LITERAL(reader, cols, 4);
  RETURN_ON_INT(reader, cols, 5, &s.count);
  RETURN_ON_LITERAL(reader, cols, 6);
  for (std::size_t i = 0; i < 13; ++i)
    s.names[i] = reader.readString(cols[7 + i].width);
  return true;
}

void writeSeqresColumns(const SeqresRecord &s, ColumnWriter &writer,
                        const ColumnSpec *cols) {
  writer.writeInt(s.serialNum, cols[1]);
  writer.writeChar(s.chainId, cols[3]);
  writer.writeInt(s.count, cols[5]);
  for (std::size_t i = 0; i < 13; ++i)
    writer.writeString(s.names[i], cols[7 + i]);
}

bool parseSplitColumns(ColumnReader &reader, Record &record,
                       const ColumnSpec *cols) {
  reader.skipChars(cols[0].width);
  SplitRecord &s = record.as<SplitRecord>();
  RETURN_ON_INT(reader, cols, 1, &s.continuation);
  RETURN_ON_LITERAL(reader, cols, 2);
  for (std::size_t i = 0; i < 14; ++i) {
    s.idMap[i] = reader.readString(cols[3 + i * 2].width);
    if (i + 1 < 14) // the separator literal between the id columns
      RETURN_ON_LITERAL(reader, cols, 4 + i * 2);
  }
  return true;
}

void writeSplitColumns(const SplitRecord &s, ColumnWriter &writer,
                       const ColumnSpec *wcols) {
  writer.writeInt(s.continuation, wcols[1]);
  for (std::size_t i = 0; i < 13; ++i)
    writer.writeString(s.idMap[i], wcols[3 + i]);
  writer.writeString(s.idMap[13], wcols[16]);
}

bool parseObslteSprsdeColumns(ColumnReader &reader, Record &record,
                              const ColumnSpec *cols, bool sprsde) {
  reader.skipChars(cols[0].width);
  int continuation = 0;
  RETURN_ON_INT(reader, cols, 1, &continuation);
  RETURN_ON_LITERAL(reader, cols, 2);
  const std::string timestamp = reader.readString(cols[3].width);
  RETURN_ON_LITERAL(reader, cols, 4);
  const std::string id = reader.readString(cols[5].width);
  reader.skipChars(cols[6].width);
  std::string ids[9];
  for (std::size_t i = 0; i < 9; ++i) {
    ids[i] = reader.readString(cols[7 + i * 2].width);
    if (i + 1 < 9)
      RETURN_ON_LITERAL(reader, cols, 8 + i * 2);
  }
  if (sprsde) {
    SprsdeRecord &s = record.as<SprsdeRecord>();
    s.continuation = continuation;
    s.timestamp = timestamp;
    s.id = id;
    for (std::size_t i = 0; i < 9; ++i)
      s.supersede[i] = ids[i];
  } else {
    ObslteRecord &o = record.as<ObslteRecord>();
    o.continuation = continuation;
    o.timestamp = timestamp;
    o.oldId = id;
    for (std::size_t i = 0; i < 9; ++i)
      o.idMap[i] = ids[i];
  }
  return true;
}

void writeObslteSprsdeColumns(const Record &record, ColumnWriter &writer,
                              const ColumnSpec *wcols, bool sprsde) {
  if (sprsde) {
    const SprsdeRecord &s = record.as<SprsdeRecord>();
    writer.writeInt(s.continuation, wcols[1]);
    writer.writeString(s.timestamp, wcols[3]);
    writer.writeString(s.id, wcols[5]);
    for (std::size_t i = 0; i < 8; ++i)
      writer.writeString(s.supersede[i], wcols[6 + i]);
    writer.writeString(s.supersede[8], wcols[14]);
  } else {
    const ObslteRecord &o = record.as<ObslteRecord>();
    writer.writeInt(o.continuation, wcols[1]);
    writer.writeString(o.timestamp, wcols[3]);
    writer.writeString(o.oldId, wcols[5]);
    for (std::size_t i = 0; i < 8; ++i)
      writer.writeString(o.idMap[i], wcols[6 + i]);
    writer.writeString(o.idMap[8], wcols[14]);
  }
}

bool parseMasterColumns(ColumnReader &reader, Record &record,
                        const ColumnSpec *cols) {
  reader.skipChars(cols[0].width);
  MasterRecord &m = record.as<MasterRecord>();
  RETURN_ON_INT(reader, cols, 1, &m.numRemark);
  RETURN_ON_INT(reader, cols, 2, &m.numFtnote);
  RETURN_ON_INT(reader, cols, 3, &m.numHet);
  RETURN_ON_INT(reader, cols, 4, &m.numHelix);
  RETURN_ON_INT(reader, cols, 5, &m.numSheet);
  RETURN_ON_INT(reader, cols, 6, &m.numTurn);
  RETURN_ON_INT(reader, cols, 7, &m.numSite);
  RETURN_ON_INT(reader, cols, 8, &m.numTransform);
  RETURN_ON_INT(reader, cols, 9, &m.numCoordinate);
  RETURN_ON_INT(reader, cols, 10, &m.numTer);
  RETURN_ON_INT(reader, cols, 11, &m.numConect);
  RETURN_ON_INT(reader, cols, 12, &m.numSeqres);
  return true;
}

void writeMasterColumns(const MasterRecord &m, ColumnWriter &writer,
                        const ColumnSpec *cols) {
  writer.writeInt(m.numRemark, cols[1]);
  writer.writeInt(m.numFtnote, cols[2]);
  writer.writeInt(m.numHet, cols[3]);
  writer.writeInt(m.numHelix, cols[4]);
  writer.writeInt(m.numSheet, cols[5]);
  writer.writeInt(m.numTurn, cols[6]);
  writer.writeInt(m.numSite, cols[7]);
  writer.writeInt(m.numTransform, cols[8]);
  writer.writeInt(m.numCoordinate, cols[9]);
  writer.writeInt(m.numTer, cols[10]);
  writer.writeInt(m.numConect, cols[11]);
  writer.writeInt(m.numSeqres, cols[12]);
}

bool parseAggrgtColumns(ColumnReader &reader, Record &record,
                        const ColumnSpec *cols) {
  reader.skipChars(cols[0].width);
  AggrgtRecord &a = record.as<AggrgtRecord>();
  RETURN_ON_INT(reader, cols, 1, &a.serialNum);
  RETURN_ON_LITERAL(reader, cols, 2);
  RETURN_ON_INT(reader, cols, 3, &a.numComponents);
  for (std::size_t i = 0; i < 14; ++i) {
    // the separator literal before every component column
    RETURN_ON_LITERAL(reader, cols, 4 + i * 2);
    RETURN_ON_INT(reader, cols, 5 + i * 2, &a.cmpontSerialNums[i]);
  }
  return true;
}

void writeAggrgtColumns(const AggrgtRecord &a, ColumnWriter &writer,
                        const ColumnSpec *cols) {
  writer.writeInt(a.serialNum, cols[1]);
  writer.writeInt(a.numComponents, cols[3]);
  for (std::size_t i = 0; i < 14; ++i)
    writer.writeInt(a.cmpontSerialNums[i], cols[5 + i * 2]);
}

bool parseTrnsfmColumns(ColumnReader &reader, Record &record,
                        const ColumnSpec *cols) {
  reader.skipChars(cols[0].width);
  TrnsfmRecord &t = record.as<TrnsfmRecord>();
  RETURN_ON_INT(reader, cols, 1, &t.resultSerialNum);
  RETURN_ON_LITERAL(reader, cols, 2);
  RETURN_ON_INT(reader, cols, 3, &t.applySerialNum);
  RETURN_ON_LITERAL(reader, cols, 4);
  RETURN_ON_INT(reader, cols, 5, &t.sourceSerialNum);
  return true;
}

void writeTrnsfmColumns(const TrnsfmRecord &t, ColumnWriter &writer,
                        const ColumnSpec *cols) {
  writer.writeInt(t.resultSerialNum, cols[1]);
  writer.writeInt(t.applySerialNum, cols[3]);
  writer.writeInt(t.sourceSerialNum, cols[5]);
}

bool parseTvectColumns(ColumnReader &reader, Record &record,
                       const ColumnSpec *cols) {
  reader.skipChars(cols[0].width);
  TvectRecord &t = record.as<TvectRecord>();
  RETURN_ON_INT(reader, cols, 1, &t.serialNum);
  RETURN_ON_REAL(reader, cols, 2, &t.t1);
  RETURN_ON_REAL(reader, cols, 3, &t.t2);
  RETURN_ON_REAL(reader, cols, 4, &t.t3);
  t.comment = reader.readString(cols[5].width);
  return true;
}

void writeTvectColumns(const TvectRecord &t, ColumnWriter &writer,
                       const ColumnSpec *cols) {
  writer.writeInt(t.serialNum, cols[1]);
  writer.writeReal(t.t1, cols[2]);
  writer.writeReal(t.t2, cols[3]);
  writer.writeReal(t.t3, cols[4]);
  writer.writeString(t.comment, cols[5]);
}

bool parseHetColumns(ColumnReader &reader, Record &record,
                     const ColumnSpec *cols) {
  reader.skipChars(cols[0].width);
  HetRecord &h = record.as<HetRecord>();
  h.hetGrp.name = reader.readString(cols[1].width);
  RETURN_ON_LITERAL(reader, cols, 2);
  h.hetGrp.chainId = reader.readChars(cols[3].width)[0];
  RETURN_ON_INT(reader, cols, 4, &h.hetGrp.seqNum);
  h.hetGrp.insertCode = reader.readChars(cols[5].width)[0];
  RETURN_ON_LITERAL(reader, cols, 6);
  RETURN_ON_INT(reader, cols, 7, &h.numAtoms);
  reader.skipChars(cols[8].width);
  h.text = reader.readString(cols[9].width);
  return true;
}

void writeHetColumns(const HetRecord &h, ColumnWriter &writer,
                     const ColumnSpec *cols) {
  writer.writeString(h.hetGrp.name, cols[1]);
  writer.writeChar(h.hetGrp.chainId, cols[3]);
  writer.writeInt(h.hetGrp.seqNum, cols[4]);
  writer.writeChar(h.hetGrp.insertCode, cols[5]);
  writer.writeInt(h.numAtoms, cols[7]);
  writer.writeString(h.text, cols[9]);
}

bool parseFormulColumns(ColumnReader &reader, Record &record,
                        const ColumnSpec *cols) {
  reader.skipChars(cols[0].width);
  FormulRecord &f = record.as<FormulRecord>();
  RETURN_ON_INT(reader, cols, 1, &f.component);
  RETURN_ON_LITERAL(reader, cols, 2);
  f.hetId = reader.readString(cols[3].width);
  RETURN_ON_INT(reader, cols, 4, &f.continuation);
  f.exclude = reader.readChars(cols[5].width)[0];
  f.formula = reader.readString(cols[6].width);
  return true;
}

void writeFormulColumns(const FormulRecord &f, ColumnWriter &writer,
                        const ColumnSpec *cols) {
  writer.writeInt(f.component, cols[1]);
  writer.writeString(f.hetId, cols[3]);
  writer.writeInt(f.continuation, cols[4]);
  writer.writeChar(f.exclude, cols[5]);
  writer.writeString(f.formula, cols[6]);
}

bool parseHeaderColumns(ColumnReader &reader, Record &record,
                        const ColumnSpec *cols) {
  reader.skipChars(cols[0].width);
  HeaderRecord &h = record.as<HeaderRecord>();
  h.classification = reader.readString(cols[1].width);
  h.timestamp = reader.readString(cols[2].width);
  RETURN_ON_LITERAL(reader, cols, 3);
  h.id = reader.readString(cols[4].width);
  return true;
}

void writeHeaderColumns(const HeaderRecord &h, ColumnWriter &writer,
                        const ColumnSpec *cols) {
  writer.writeString(h.classification, cols[1]);
  writer.writeString(h.timestamp, cols[2]);
  writer.writeString(h.id, cols[4]);
}

bool parseCaveatColumns(ColumnReader &reader, Record &record,
                        const ColumnSpec *cols) {
  reader.skipChars(cols[0].width);
  CaveatRecord &c = record.as<CaveatRecord>();
  RETURN_ON_INT(reader, cols, 1, &c.continuation);
  RETURN_ON_LITERAL(reader, cols, 2);
  c.id = reader.readString(cols[3].width);
  RETURN_ON_LITERAL(reader, cols, 4);
  c.text = reader.readString(cols[5].width);
  return true;
}

void writeCaveatColumns(const CaveatRecord &c, ColumnWriter &writer,
                        const ColumnSpec *cols) {
  writer.writeInt(c.continuation, cols[1]);
  writer.writeString(c.id, cols[3]);
  writer.writeString(c.text, cols[5]);
}

bool parseJrnlColumns(ColumnReader &reader, Record &record,
                      const ColumnSpec *cols) {
  reader.skipChars(cols[0].width);
  record.as<JrnlRecord>().text = reader.readString(cols[1].width);
  return true;
}

void writeJrnlColumns(const JrnlRecord &j, ColumnWriter &writer,
                      const ColumnSpec *cols) {
  writer.writeString(j.text, cols[1]);
}

bool parseAgrdesColumns(ColumnReader &reader, Record &record,
                        const ColumnSpec *cols) {
  reader.skipChars(cols[0].width);
  AgrdesRecord &a = record.as<AgrdesRecord>();
  RETURN_ON_INT(reader, cols, 1, &a.num);
  RETURN_ON_LITERAL(reader, cols, 2);
  a.text = reader.readString(cols[3].width);
  return true;
}

void writeAgrdesColumns(const AgrdesRecord &a, ColumnWriter &writer,
                        const ColumnSpec *cols) {
  writer.writeInt(a.num, cols[1]);
  writer.writeString(a.text, cols[3]);
}

bool parseAnisouColumns(ColumnReader &reader, Record &record,
                        const ColumnSpec *cols) {
  reader.skipChars(cols[0].width);
  AnisouRecord &a = record.as<AnisouRecord>();
  std::int64_t serial = 0;
  if (!reader.readInt(cols[1].width, &serial))
    return false;
  a.serialNum = static_cast<int>(serial);
  RETURN_ON_LITERAL(reader, cols, 2);
  a.name = reader.readString(cols[3].width);
  a.altLoc = reader.readChars(cols[4].width)[0];
  a.residue.name = reader.readString(cols[5].width);
  a.residue.chainId = reader.readChars(cols[6].width)[0];
  std::int64_t seq = 0;
  if (!reader.readInt(cols[7].width, &seq)) {
      return false;
  }
  a.residue.seqNum = static_cast<int>(seq);
  a.residue.insertCode = reader.readChars(cols[8].width)[0];
  RETURN_ON_LITERAL(reader, cols, 9);
  for (std::size_t i = 0; i < 6; ++i)
    RETURN_ON_INT(reader, cols, 10 + i, &a.u[i]);
  return true;
}

void writeAnisouColumns(const AnisouRecord &a, ColumnWriter &writer,
                        const ColumnSpec *cols) {
  writer.writeInt(a.serialNum, cols[1]);
  writer.writeString(a.name, cols[3]);
  writer.writeChar(a.altLoc, cols[4]);
  writer.writeString(a.residue.name, cols[5]);
  writer.writeChar(a.residue.chainId, cols[6]);
  writer.writeInt(a.residue.seqNum, cols[7]);
  writer.writeChar(a.residue.insertCode, cols[8]);
  for (std::size_t i = 0; i < 6; ++i)
    writer.writeInt(a.u[i], cols[10 + i]);
}

bool parseModelColumns(ColumnReader &reader, Record &record,
                       const ColumnSpec *cols) {
  reader.skipChars(cols[0].width);
  std::int64_t num = 0;
  if (!reader.readInt(cols[1].width, &num))
    return false;
  record.as<ModelRecord>().num = static_cast<int>(num);
  return true;
}

void writeModelColumns(const ModelRecord &m, ColumnWriter &writer,
                       const ColumnSpec *cols) {
  writer.writeInt(m.num, cols[1]);
}

bool parseNummdlColumns(ColumnReader &reader, Record &record,
                        const ColumnSpec *cols) {
  reader.skipChars(cols[0].width);
  std::int64_t number = 0;
  if (!reader.readInt(cols[1].width, &number))
    return false;
  record.as<NummdlRecord>().number = static_cast<int>(number);
  return true;
}

void writeNummdlColumns(const NummdlRecord &n, ColumnWriter &writer,
                        const ColumnSpec *cols) {
  writer.writeInt(n.number, cols[1]);
}

bool parseExpdtaColumns(ColumnReader &reader, Record &record,
                        const ColumnSpec *cols) {
  reader.skipChars(cols[0].width);
  ExpdtaRecord &e = record.as<ExpdtaRecord>();
  RETURN_ON_INT(reader, cols, 1, &e.continuation);
  e.technique = reader.readString(cols[2].width);
  return true;
}

void writeExpdtaColumns(const ExpdtaRecord &e, ColumnWriter &writer,
                        const ColumnSpec *cols) {
  writer.writeInt(e.continuation, cols[1]);
  writer.writeString(e.technique, cols[2]);
}


// ---------------------------------------------------------------------------
// USER records (Phase 2, commit 7). The PDBRUN layouts are version-dependent
// (v5 and v6 tables); the v5 libc scanf paths use whitespace-delimited
// tokens instead of fixed columns, with dedicated bounded parsers.
// ---------------------------------------------------------------------------

constexpr ColumnSpec kUserPdbrunCols[] = {
    {0, 12, ColumnKind::Skip},
    {12, 1, ColumnKind::Literal, 0, false, ' '},
    {13, 2, ColumnKind::Int},
};
constexpr ColumnSpec kUserEyePosCols[] = {
    {0, 12, ColumnKind::Skip},
    {12, 1, ColumnKind::Literal, 0, false, ' '},
    {13, 9, ColumnKind::Real, 3},
    {22, 1, ColumnKind::Literal, 0, false, ' '},
    {23, 9, ColumnKind::Real, 3},
    {32, 1, ColumnKind::Literal, 0, false, ' '},
    {33, 9, ColumnKind::Real, 3},
};
constexpr ColumnSpec kUserWindowCols[] = {
    {0, 12, ColumnKind::Skip},
    {12, 1, ColumnKind::Literal, 0, false, ' '},
    {13, 9, ColumnKind::Real, 3},
    {22, 1, ColumnKind::Literal, 0, false, ' '},
    {23, 9, ColumnKind::Real, 3},
    {32, 1, ColumnKind::Literal, 0, false, ' '},
    {33, 9, ColumnKind::Real, 3},
    {42, 1, ColumnKind::Literal, 0, false, ' '},
    {43, 9, ColumnKind::Real, 3},
    {52, 1, ColumnKind::Literal, 0, false, ' '},
    {53, 9, ColumnKind::Real, 3},
    {62, 1, ColumnKind::Literal, 0, false, ' '},
    {63, 9, ColumnKind::Real, 3},
};
constexpr ColumnSpec kUserFocusCols[] = {
    {0, 11, ColumnKind::Skip},
    {11, 2, ColumnKind::Literal, 0, false, ' '},
    {13, 9, ColumnKind::Real, 3},
};
constexpr ColumnSpec kUserViewportCols[] = {
    {0, 14, ColumnKind::Skip},
    {14, 1, ColumnKind::Literal, 0, false, ' '},
    {15, 9, ColumnKind::Real, 3},
    {24, 1, ColumnKind::Literal, 0, false, ' '},
    {25, 9, ColumnKind::Real, 3},
    {34, 1, ColumnKind::Literal, 0, false, ' '},
    {35, 9, ColumnKind::Real, 3},
    {44, 1, ColumnKind::Literal, 0, false, ' '},
    {45, 9, ColumnKind::Real, 3},
};
constexpr ColumnSpec kUserBgColorCols[] = {
    {0, 13, ColumnKind::Skip},
    {13, 1, ColumnKind::Literal, 0, false, ' '},
    {14, 5, ColumnKind::Real, 3},
    {19, 1, ColumnKind::Literal, 0, false, ' '},
    {20, 5, ColumnKind::Real, 3},
    {25, 1, ColumnKind::Literal, 0, false, ' '},
    {26, 5, ColumnKind::Real, 3},
};
constexpr ColumnSpec kUserRadiusCols[] = {
    {0, 12, ColumnKind::Skip},
    {12, 1, ColumnKind::Literal, 0, false, ' '},
    {13, 7, ColumnKind::Real, 3},
};
constexpr ColumnSpec kUserObjectCols[] = {
    {0, 12, ColumnKind::Skip},
    {12, 1, ColumnKind::Literal, 0, false, ' '},
    {13, 3, ColumnKind::Int},
};
constexpr ColumnSpec kUserAngle6Cols[] = {
    {0, 11, ColumnKind::Skip},
    {11, 1, ColumnKind::Literal, 0, false, ' '},
    {12, 6, ColumnKind::Int},
    {18, 1, ColumnKind::Literal, 0, false, ' '},
    {19, 6, ColumnKind::Int},
    {24, 1, ColumnKind::Literal, 0, false, ' '},
    {25, 6, ColumnKind::Int},
    {30, 1, ColumnKind::Literal, 0, false, ' '},
    {31, 6, ColumnKind::Int},
    {36, 1, ColumnKind::Literal, 0, false, ' '},
    {37, 9, ColumnKind::Real, 3},
};
constexpr ColumnSpec kUserAngle5Cols[] = {
    {0, 11, ColumnKind::Skip},
    {11, 1, ColumnKind::Literal, 0, false, ' '},
    {12, 2, ColumnKind::Int},
    {14, 1, ColumnKind::Literal, 0, false, ' '},
    {15, 6, ColumnKind::Int},
    {21, 1, ColumnKind::Literal, 0, false, ' '},
    {22, 6, ColumnKind::Int},
    {27, 1, ColumnKind::Literal, 0, false, ' '},
    {28, 6, ColumnKind::Int},
    {33, 1, ColumnKind::Literal, 0, false, ' '},
    {34, 6, ColumnKind::Int},
    {39, 1, ColumnKind::Literal, 0, false, ' '},
    {40, 8, ColumnKind::Real, 3},
};
constexpr ColumnSpec kUserDistance6Cols[] = {
    {0, 14, ColumnKind::Skip},
    {14, 1, ColumnKind::Literal, 0, false, ' '},
    {15, 6, ColumnKind::Int},
    {21, 1, ColumnKind::Literal, 0, false, ' '},
    {22, 6, ColumnKind::Int},
    {27, 1, ColumnKind::Literal, 0, false, ' '},
    {28, 9, ColumnKind::Real, 3},
};
constexpr ColumnSpec kUserDistance5Cols[] = {
    {0, 14, ColumnKind::Skip},
    {14, 1, ColumnKind::Literal, 0, false, ' '},
    {15, 2, ColumnKind::Int},
    {17, 1, ColumnKind::Literal, 0, false, ' '},
    {18, 6, ColumnKind::Int},
    {24, 1, ColumnKind::Literal, 0, false, ' '},
    {25, 6, ColumnKind::Int},
    {31, 1, ColumnKind::Literal, 0, false, ' '},
    {32, 8, ColumnKind::Real, 3},
};
constexpr ColumnSpec kUserFile6Cols[] = {
    {0, 10, ColumnKind::Skip},
    {10, 1, ColumnKind::Literal, 0, false, ' '},
    {11, 4, ColumnKind::Int},
    {15, 1, ColumnKind::Literal, 0, false, ' '},
    {16, 56, ColumnKind::Str, 0, true},
};
constexpr ColumnSpec kUserFile5Cols[] = {
    {0, 10, ColumnKind::Skip},
    {10, 1, ColumnKind::Literal, 0, false, ' '},
    {11, 61, ColumnKind::Str, 0, true},
};
constexpr ColumnSpec kUserMarkname6Cols[] = {
    {0, 14, ColumnKind::Skip},
    {14, 1, ColumnKind::Literal, 0, false, ' '},
    {15, 57, ColumnKind::Str, 0, true},
};
constexpr ColumnSpec kUserMark6Cols[] = {
    {0, 10, ColumnKind::Skip},
    {10, 1, ColumnKind::Literal, 0, false, ' '},
    {11, 57, ColumnKind::Str, 0, true},
};
constexpr ColumnSpec kUserCName6Cols[] = {
    {0, 11, ColumnKind::Skip},
    {11, 1, ColumnKind::Literal, 0, false, ' '},
    {12, 5, ColumnKind::Real, 3},
    {17, 1, ColumnKind::Literal, 0, false, ' '},
    {18, 5, ColumnKind::Real, 3},
    {23, 1, ColumnKind::Literal, 0, false, ' '},
    {24, 5, ColumnKind::Real, 3},
    {29, 1, ColumnKind::Literal, 0, false, ' '},
    {30, 38, ColumnKind::Str, 0, true},
};
constexpr ColumnSpec kUserChain6Cols[] = {
    {0, 11, ColumnKind::Skip},
    {11, 1, ColumnKind::Literal, 0, false, ' '},
    {12, 6, ColumnKind::Int},
    {18, 1, ColumnKind::Literal, 0, false, ' '},
    {19, 6, ColumnKind::Int},
};
constexpr ColumnSpec kUserGfxBegin6Cols[] = {
    {0, 15, ColumnKind::Skip},
    {15, 1, ColumnKind::Literal, 0, false, ' '},
    {16, 32, ColumnKind::Str, 0, true},
};
constexpr ColumnSpec kUserGfxColor6Cols[] = {
    {0, 15, ColumnKind::Skip},
    {15, 1, ColumnKind::Literal, 0, false, ' '},
    {16, 5, ColumnKind::Real, 3},
    {21, 1, ColumnKind::Literal, 0, false, ' '},
    {22, 5, ColumnKind::Real, 3},
    {27, 1, ColumnKind::Literal, 0, false, ' '},
    {28, 5, ColumnKind::Real, 3},
    {33, 1, ColumnKind::Literal, 0, false, ' '},
    {34, 38, ColumnKind::Str, 0, true},
};
constexpr ColumnSpec kUserGfxPoint6Cols[] = {
    {0, 16, ColumnKind::Skip},
    {16, 1, ColumnKind::Literal, 0, false, ' '},
    {17, 9, ColumnKind::Real, 3},
    {26, 1, ColumnKind::Literal, 0, false, ' '},
    {27, 9, ColumnKind::Real, 3},
    {36, 1, ColumnKind::Literal, 0, false, ' '},
    {37, 9, ColumnKind::Real, 3},
};
constexpr ColumnSpec kUserGfxTextPos6Cols[] = {
    {0, 17, ColumnKind::Skip},
    {17, 1, ColumnKind::Literal, 0, false, ' '},
    {18, 9, ColumnKind::Real, 3},
    {27, 1, ColumnKind::Literal, 0, false, ' '},
    {28, 9, ColumnKind::Real, 3},
    {37, 1, ColumnKind::Literal, 0, false, ' '},
    {38, 9, ColumnKind::Real, 3},
};
constexpr ColumnSpec kUserGfxFont6Cols[] = {
    {0, 14, ColumnKind::Skip},
    {14, 1, ColumnKind::Literal, 0, false, ' '},
    {15, 3, ColumnKind::Int},
    {18, 1, ColumnKind::Literal, 0, false, ' '},
    {19, 53, ColumnKind::Str, 0, true},
};
constexpr ColumnSpec kUserGfxLabel6Cols[] = {
    {0, 15, ColumnKind::Skip},
    {15, 1, ColumnKind::Literal, 0, false, ' '},
    {16, 56, ColumnKind::Str, 0, true},
};
constexpr ColumnSpec kUserGfxMove5Cols[] = {
    {0, 14, ColumnKind::Skip},
    {14, 1, ColumnKind::Literal, 0, false, ' '},
    {15, 9, ColumnKind::Real, 3},
    {24, 1, ColumnKind::Literal, 0, false, ' '},
    {25, 9, ColumnKind::Real, 3},
    {34, 1, ColumnKind::Literal, 0, false, ' '},
    {35, 9, ColumnKind::Real, 3},
};
constexpr ColumnSpec kUserGfxMarker5Cols[] = {
    {0, 16, ColumnKind::Skip},
    {16, 1, ColumnKind::Literal, 0, false, ' '},
    {17, 9, ColumnKind::Real, 3},
    {26, 1, ColumnKind::Literal, 0, false, ' '},
    {27, 9, ColumnKind::Real, 3},
    {36, 1, ColumnKind::Literal, 0, false, ' '},
    {37, 9, ColumnKind::Real, 3},
};
constexpr ColumnSpec kUserGfxPoint5Cols[] = {
    {0, 15, ColumnKind::Skip},
    {15, 1, ColumnKind::Literal, 0, false, ' '},
    {16, 9, ColumnKind::Real, 3},
    {25, 1, ColumnKind::Literal, 0, false, ' '},
    {26, 9, ColumnKind::Real, 3},
    {35, 1, ColumnKind::Literal, 0, false, ' '},
    {36, 9, ColumnKind::Real, 3},
};
constexpr ColumnSpec kUserCols[] = {
    {0, 4, ColumnKind::Skip},
    {4, 2, ColumnKind::Str},
    {6, 66, ColumnKind::Str, 0, true},
};


// Whitespace-delimited token and number scanners for the v5 libc paths.
// A token is bounded by its destination size (the %38s/%53s bounds); the
// numbers parse with the shared fixed-format grammar on each token, so a
// missing or non-numeric token keeps the member at its default and the
// record still parses (libc sscanf matching-failure semantics: assigned
// fields keep their values, the rest stay untouched, the call succeeds).
bool scanWhitespaceToken(std::string_view &rest, std::string_view &token,
                         int bound) {
  std::size_t i = 0;
  while (i < rest.size() && rest[i] == ' ')
    ++i;
  const std::size_t start = i;
  while (i < rest.size() && rest[i] != ' ' && rest[i] != '\0' &&
         rest[i] != '\n' && i - start < static_cast<std::size_t>(bound))
    ++i;
  if (i == start)
    return false;
  token = rest.substr(start, i - start);
  rest = rest.substr(i);
  return true;
}

bool scanWhitespaceReal(std::string_view &rest, double *value) {
  std::string_view token;
  if (!scanWhitespaceToken(rest, token, 32))
    return false;
  while (!token.empty() &&
         std::isspace(static_cast<unsigned char>(token.back())))
    token.remove_suffix(1);
  return parseFixedReal(token, value);
}

bool scanWhitespaceInt(std::string_view &rest, std::int64_t *value) {
  std::string_view token;
  if (!scanWhitespaceToken(rest, token, 20))
    return false;
  while (!token.empty() &&
         std::isspace(static_cast<unsigned char>(token.back())))
    token.remove_suffix(1);
  return parseFixedInt(token, value);
}

// The v5 libc USER records parse with whitespace-delimited tokens after
// their fixed subtype prefix; a matching failure leaves the remaining
// members at their defaults while the record still parses.
bool parseUserBgColor5(std::string_view line, Record &record) {
  std::string_view rest = line.substr(13); // after "USER  BGCOLOR"
  UserBgColorRecord &b = record.as<UserBgColorRecord>();
  double rgb[3] = {0.0, 0.0, 0.0};
  for (int i = 0; i < 3; ++i)
    if (!scanWhitespaceReal(rest, &rgb[i]))
      break; // matching failure: keep the defaults
  b.rgb[0] = rgb[0];
  b.rgb[1] = rgb[1];
  b.rgb[2] = rgb[2];
  return true;
}

bool parseUserCName5(std::string_view line, Record &record, bool gfx) {
  const std::size_t prefix =
      gfx ? 15 : 11; // "USER  GFX COLOR" vs "USER  CNAME"/"USER  COLOR"
  std::string_view rest = line.substr(std::min(prefix, line.size()));
  std::string_view token;
  std::string name;
  if (scanWhitespaceToken(rest, token, 38))
    name = std::string(token);
  double rgb[3] = {0.0, 0.0, 0.0};
  for (int i = 0; i < 3; ++i) {
    if (!scanWhitespaceReal(rest, &rgb[i]))
      break; // matching failure: keep the defaults
  }
  if (record.holds<UserCNameRecord>()) {
    UserCNameRecord &c = record.as<UserCNameRecord>();
    c.name = name;
    c.rgb = std::array<double, 3>{rgb[0], rgb[1], rgb[2]};
  } else {
    UserColorRecord &c = record.as<UserColorRecord>();
    c.spec = name;
    c.rgb = std::array<double, 3>{rgb[0], rgb[1], rgb[2]};
  }
  return true;
}

bool parseUserChain5(std::string_view line, Record &record) {
  std::string_view rest = line.substr(std::min<std::size_t>(11, line.size()));
  UserChainRecord &c = record.as<UserChainRecord>();
  std::int64_t atom0 = 0;
  std::int64_t atom1 = 0;
  if (scanWhitespaceInt(rest, &atom0))
    c.atom0 = static_cast<int>(atom0);
  if (scanWhitespaceInt(rest, &atom1))
    c.atom1 = static_cast<int>(atom1);
  return true;
}

bool parseUserGfxFont5(std::string_view line, Record &record) {
  std::string_view rest = line.substr(std::min<std::size_t>(15, line.size()));
  std::string_view token;
  std::string name;
  if (scanWhitespaceToken(rest, token, 53))
    name = std::string(token);
  UserGfxFontRecord &f = record.as<UserGfxFontRecord>();
  f.name = name;
  std::int64_t size = 0;
  if (scanWhitespaceInt(rest, &size))
    f.size = static_cast<int>(size);
  return true;
}

bool parseUserGfxLabel5(std::string_view line, Record &record) {
  std::string_view rest = line.substr(std::min<std::size_t>(15, line.size()));
  UserGfxLabelRecord &l = record.as<UserGfxLabelRecord>();
  double xyz[3] = {0.0, 0.0, 0.0};
  for (int i = 0; i < 3; ++i) {
    if (!scanWhitespaceReal(rest, &xyz[i]))
      return true; // matching failure: everything keeps its default
  }
  l.xyz = std::array<double, 3>{xyz[0], xyz[1], xyz[2]};
  // the scanset takes the rest of the line verbatim, bounded to 26
  std::size_t i = 0;
  while (i < rest.size() && rest[i] == ' ')
    ++i;
  std::size_t end = i;
  while (end < rest.size() && rest[end] != '\n' && end - i < 26)
    ++end;
  l.text = std::string(rest.substr(i, end - i));
  return true;
}


// Fixed-width and minimal field renderers for the USER appends; the
// minimal forms reproduce the libc %d / %.Nf conversions the v5 write
// formats use (no width, no padding).
std::string renderFixedReal(double value, int width, int places) {
  char buf[48];
  std::memset(buf, ' ', sizeof(buf));
  ColumnWriter writer(buf, sizeof(buf));
  writer.writeReal(value, ColumnSpec{0, static_cast<std::uint8_t>(width),
                                     ColumnKind::Real,
                                     static_cast<std::uint8_t>(places)});
  return std::string(buf, static_cast<std::size_t>(width));
}

std::string renderFixedInt(std::int64_t value, int width, bool zeroAsBlank,
                           bool leftJustify) {
  char buf[48];
  std::memset(buf, ' ', sizeof(buf));
  ColumnWriter writer(buf, sizeof(buf));
  writer.writeInt(value, ColumnSpec{0, static_cast<std::uint8_t>(width),
                                    ColumnKind::Int, 0, leftJustify, ' ',
                                    zeroAsBlank});
  return std::string(buf, static_cast<std::size_t>(width));
}

std::string renderMinimalReal(double value, int places) {
  std::string rendered = renderFixedReal(value, 24, places);
  const std::size_t first = rendered.find_first_not_of(' ');
  return first == std::string::npos ? "" : rendered.substr(first);
}

std::string renderMinimalInt(std::int64_t value) {
  std::string rendered = renderFixedInt(value, 24, false, false);
  const std::size_t first = rendered.find_first_not_of(' ');
  return first == std::string::npos ? "" : rendered.substr(first);
}

std::string renderFixedStr(std::string_view text, int width, bool leftJustify) {
  char buf[128];
  std::memset(buf, ' ', sizeof(buf));
  ColumnWriter writer(buf, sizeof(buf));
  writer.writeString(text, ColumnSpec{0, static_cast<std::uint8_t>(width),
                                      ColumnKind::Str, 0, leftJustify});
  return std::string(buf, static_cast<std::size_t>(width));
}

// The graphics primitive name for a captured USER GFX BEGIN text, mirroring
// the legacy getGfxType matchers.
GfxType gfxTypeOf(std::string_view text) {
  const auto equals = [&text](std::size_t offset, const char *literal) {
    const std::size_t len = std::strlen(literal);
    if (offset + len > text.size())
      return false;
    for (std::size_t i = 0; i < len; ++i) {
      const char c = text[offset + i];
      const char l = literal[i];
      if (c != l && (c < 'a' || c > 'z' || c - 'a' + 'A') != l)
        return false;
    }
    return true;
  };
  switch (text.empty() ? '\0' : text[0]) {
  case 'L':
  case 'l':
    if (equals(1, "INE-LOOP"))
      return GfxType::LINE_LOOP;
    if (equals(1, "INE-STRIP"))
      return GfxType::LINE_STRIP;
    if (equals(1, "INES"))
      return GfxType::LINES;
    break;
  case 'M':
  case 'm':
    if (equals(1, "ARKERS"))
      return GfxType::MARKERS;
    break;
  case 'P':
  case 'p':
    if (equals(1, "OINTS"))
      return GfxType::POINTS;
    if (equals(1, "OLYGON"))
      return GfxType::POLYGON;
    break;
  case 'Q':
  case 'q':
    if (equals(1, "UAD-STRIP"))
      return GfxType::QUAD_STRIP;
    if (equals(1, "UADS"))
      return GfxType::QUADS;
    break;
  case 'T':
  case 't':
    if (equals(1, "RIANGLE-FAN"))
      return GfxType::TRIANGLE_FAN;
    if (equals(1, "RIANGLE-STRIP"))
      return GfxType::TRIANGLE_STRIP;
    if (equals(1, "RIANGLES"))
      return GfxType::TRIANGLES;
    break;
  }
  return GfxType::UNKNOWN;
}

std::string gfxName(GfxType type) {
  switch (type) {
  case GfxType::POINTS:
    return "POINTS";
  case GfxType::MARKERS:
    return "MARKERS";
  case GfxType::LINES:
    return "LINES";
  case GfxType::LINE_STRIP:
    return "LINE-STRIP";
  case GfxType::LINE_LOOP:
    return "LINE-LOOP";
  case GfxType::TRIANGLES:
    return "TRIANGLES";
  case GfxType::TRIANGLE_STRIP:
    return "TRIANGLE-STRIP";
  case GfxType::TRIANGLE_FAN:
    return "TRIANGLE-FAN";
  case GfxType::QUADS:
    return "QUADS";
  case GfxType::QUAD_STRIP:
    return "QUAD-STRIP";
  case GfxType::POLYGON:
    return "POLYGON";
  case GfxType::UNKNOWN:
    return "UNKNOWN";
  }
  return "UNKNOWN";
}

// Fixed-column USER binders; every field reads through the shared parser
// and the version selects the table.
bool parseUserFixedRealFields(ColumnReader &reader, double *values,
                              const ColumnSpec *cols, std::size_t first,
                              std::size_t count, std::size_t stride) {
  for (std::size_t i = 0; i < count; ++i) {
    if (i > 0)
      RETURN_ON_LITERAL(reader, cols, first + i * stride - 1);
    RETURN_ON_REAL(reader, cols, first + i * stride, &values[i]);
  }
  return true;
}

bool parseUserPdbrun(ColumnReader &reader, Record &record,
                     const ColumnSpec *cols) {
  reader.skipChars(cols[0].width);
  RETURN_ON_LITERAL(reader, cols, 1);
  std::int64_t version = 0;
  if (!reader.readInt(cols[2].width, &version))
    return false;
  record.as<UserPdbrunRecord>().version = static_cast<int>(version);
  return true;
}

bool parseUserEyePos(ColumnReader &reader, Record &record,
                     const ColumnSpec *cols) {
  reader.skipChars(cols[0].width);
  RETURN_ON_LITERAL(reader, cols, 1);
  UserEyePosRecord &e = record.as<UserEyePosRecord>();
  double xyz[3] = {0.0, 0.0, 0.0};
  if (!parseUserFixedRealFields(reader, xyz, cols, 2, 3, 2))
    return false;
  e.xyz = std::array<double, 3>{xyz[0], xyz[1], xyz[2]};
  return true;
}

bool parseUserWindow(ColumnReader &reader, Record &record,
                     const ColumnSpec *cols) {
  reader.skipChars(cols[0].width);
  RETURN_ON_LITERAL(reader, cols, 1);
  UserWindowRecord &w = record.as<UserWindowRecord>();
  double values[6] = {0.0, 0.0, 0.0, 0.0, 0.0, 0.0};
  if (!parseUserFixedRealFields(reader, values, cols, 2, 6, 2))
    return false;
  w.left = values[0];
  w.right = values[1];
  w.bottom = values[2];
  w.top = values[3];
  w.hither = values[4];
  w.yon = values[5];
  return true;
}

bool parseUserFocus(ColumnReader &reader, Record &record,
                    const ColumnSpec *cols) {
  reader.skipChars(cols[0].width);
  RETURN_ON_LITERAL(reader, cols, 1);
  double focus = 0.0;
  RETURN_ON_REAL(reader, cols, 2, &focus);
  record.as<UserFocusRecord>().focus = focus;
  return true;
}

bool parseUserViewport(ColumnReader &reader, Record &record,
                       const ColumnSpec *cols) {
  reader.skipChars(cols[0].width);
  RETURN_ON_LITERAL(reader, cols, 1);
  UserViewportRecord &v = record.as<UserViewportRecord>();
  double values[4] = {0.0, 0.0, 0.0, 0.0};
  if (!parseUserFixedRealFields(reader, values, cols, 2, 4, 2))
    return false;
  v.xmin = values[0];
  v.xmax = values[1];
  v.ymin = values[2];
  v.ymax = values[3];
  return true;
}

bool parseUserBgColor6(ColumnReader &reader, Record &record,
                       const ColumnSpec *cols) {
  reader.skipChars(cols[0].width);
  RETURN_ON_LITERAL(reader, cols, 1);
  UserBgColorRecord &b = record.as<UserBgColorRecord>();
  double rgb[3] = {0.0, 0.0, 0.0};
  if (!parseUserFixedRealFields(reader, rgb, cols, 2, 3, 2))
    return false;
  b.rgb = std::array<double, 3>{rgb[0], rgb[1], rgb[2]};
  return true;
}

bool parseUserRadius(ColumnReader &reader, Record &record,
                     const ColumnSpec *cols) {
  reader.skipChars(cols[0].width);
  RETURN_ON_LITERAL(reader, cols, 1);
  double radius = 0.0;
  RETURN_ON_REAL(reader, cols, 2, &radius);
  record.as<UserRadiusRecord>().radius = radius;
  return true;
}

bool parseUserAngle(ColumnReader &reader, Record &record, const ColumnSpec *cols,
                    bool v5) {
  reader.skipChars(cols[0].width);
  RETURN_ON_LITERAL(reader, cols, 1);
  UserAngleRecord &a = record.as<UserAngleRecord>();
  std::int64_t value = 0;
  std::size_t base = 2;
  if (v5) {
    RETURN_ON_INT(reader, cols, 2, &value);
    a.which = static_cast<int>(value);
    RETURN_ON_LITERAL(reader, cols, 3);
    base = 4;
  } else {
    base = 2;
  }
  RETURN_ON_INT(reader, cols, base, &a.atom0);
  RETURN_ON_LITERAL(reader, cols, base + 1);
  RETURN_ON_INT(reader, cols, base + 2, &a.atom1);
  RETURN_ON_LITERAL(reader, cols, base + 3);
  RETURN_ON_INT(reader, cols, base + 4, &a.atom2);
  RETURN_ON_LITERAL(reader, cols, base + 5);
  RETURN_ON_INT(reader, cols, base + 6, &a.atom3);
  RETURN_ON_LITERAL(reader, cols, base + 7);
  double angle = 0.0;
  RETURN_ON_REAL(reader, cols, base + 8, &angle);
  a.angle = angle;
  return true;
}

bool parseUserDistance(ColumnReader &reader, Record &record,
                       const ColumnSpec *cols, bool v5) {
  reader.skipChars(cols[0].width);
  RETURN_ON_LITERAL(reader, cols, 1);
  UserDistanceRecord &d = record.as<UserDistanceRecord>();
  std::int64_t value = 0;
  std::size_t base = 2;
  if (v5) {
    RETURN_ON_INT(reader, cols, 2, &value);
    d.which = static_cast<int>(value);
    RETURN_ON_LITERAL(reader, cols, 3);
    base = 4;
  } else {
    base = 2;
  }
  RETURN_ON_INT(reader, cols, base, &d.atom0);
  RETURN_ON_LITERAL(reader, cols, base + 1);
  RETURN_ON_INT(reader, cols, base + 2, &d.atom1);
  RETURN_ON_LITERAL(reader, cols, base + 3);
  double distance = 0.0;
  RETURN_ON_REAL(reader, cols, base + 4, &distance);
  d.distance = distance;
  return true;
}

bool parseUserFile(ColumnReader &reader, Record &record, const ColumnSpec *cols,
                   bool v5) {
  reader.skipChars(cols[0].width);
  RETURN_ON_LITERAL(reader, cols, 1);
  UserFileRecord &f = record.as<UserFileRecord>();
  if (v5) {
    f.filename = reader.readString(cols[2].width);
    return true;
  }
  std::int64_t model = 0;
  RETURN_ON_INT(reader, cols, 2, &model);
  f.model = static_cast<int>(model);
  RETURN_ON_LITERAL(reader, cols, 3);
  f.filename = reader.readString(cols[4].width);
  return true;
}

bool parseUserMarkname(ColumnReader &reader, Record &record,
                       const ColumnSpec *cols) {
  reader.skipChars(cols[0].width);
  RETURN_ON_LITERAL(reader, cols, 1);
  record.as<UserMarknameRecord>().markname =
      reader.readString(cols[2].width);
  return true;
}

bool parseUserRgbName(ColumnReader &reader, Record &record,
                      const ColumnSpec *cols) {
  reader.skipChars(cols[0].width);
  RETURN_ON_LITERAL(reader, cols, 1);
  double rgb[3] = {0.0, 0.0, 0.0};
  if (!parseUserFixedRealFields(reader, rgb, cols, 2, 3, 2))
    return false;
  RETURN_ON_LITERAL(reader, cols, 7);
  const std::string name = reader.readString(cols[8].width);
  if (record.holds<UserCNameRecord>()) {
    UserCNameRecord &c = record.as<UserCNameRecord>();
    c.rgb = std::array<double, 3>{rgb[0], rgb[1], rgb[2]};
    c.name = name;
  } else {
    UserColorRecord &c = record.as<UserColorRecord>();
    c.rgb = std::array<double, 3>{rgb[0], rgb[1], rgb[2]};
    c.spec = name;
  }
  return true;
}

bool parseUserChain6(ColumnReader &reader, Record &record,
                     const ColumnSpec *cols) {
  reader.skipChars(cols[0].width);
  RETURN_ON_LITERAL(reader, cols, 1);
  UserChainRecord &c = record.as<UserChainRecord>();
  std::int64_t atom0 = 0;
  std::int64_t atom1 = 0;
  RETURN_ON_INT(reader, cols, 2, &atom0);
  RETURN_ON_LITERAL(reader, cols, 3);
  RETURN_ON_INT(reader, cols, 4, &atom1);
  c.atom0 = static_cast<int>(atom0);
  c.atom1 = static_cast<int>(atom1);
  return true;
}

bool parseUserObject(ColumnReader &reader, Record &record,
                     const ColumnSpec *cols) {
  reader.skipChars(cols[0].width);
  RETURN_ON_LITERAL(reader, cols, 1);
  std::int64_t model = 0;
  RETURN_ON_INT(reader, cols, 2, &model);
  record.as<UserObjectRecord>().model = static_cast<int>(model);
  return true;
}

bool parseUserGfxBegin6(ColumnReader &reader, Record &record,
                        const ColumnSpec *cols) {
  reader.skipChars(cols[0].width);
  RETURN_ON_LITERAL(reader, cols, 1);
  UserGfxBeginRecord &b = record.as<UserGfxBeginRecord>();
  b.unknown = reader.readString(cols[2].width);
  b.primitive = gfxTypeOf(b.unknown);
  return true;
}

bool parseUserGfxColor6(ColumnReader &reader, Record &record,
                        const ColumnSpec *cols) {
  reader.skipChars(cols[0].width);
  RETURN_ON_LITERAL(reader, cols, 1);
  UserColorRecord &c = record.as<UserColorRecord>();
  double rgb[3] = {0.0, 0.0, 0.0};
  if (!parseUserFixedRealFields(reader, rgb, cols, 2, 3, 2))
    return false;
  RETURN_ON_LITERAL(reader, cols, 7);
  c.rgb = std::array<double, 3>{rgb[0], rgb[1], rgb[2]};
  c.spec = reader.readString(cols[8].width);
  return true;
}

bool parseUserGfxPoint6(ColumnReader &reader, Record &record,
                        const ColumnSpec *cols) {
  reader.skipChars(cols[0].width);
  RETURN_ON_LITERAL(reader, cols, 1);
  UserGfxPointRecord &p = record.as<UserGfxPointRecord>();
  double xyz[3] = {0.0, 0.0, 0.0};
  if (!parseUserFixedRealFields(reader, xyz, cols, 2, 3, 2))
    return false;
  p.xyz = std::array<double, 3>{xyz[0], xyz[1], xyz[2]};
  return true;
}

bool parseUserGfxFont6(ColumnReader &reader, Record &record,
                       const ColumnSpec *cols) {
  reader.skipChars(cols[0].width);
  RETURN_ON_LITERAL(reader, cols, 1);
  UserGfxFontRecord &f = record.as<UserGfxFontRecord>();
  std::int64_t size = 0;
  RETURN_ON_INT(reader, cols, 2, &size);
  f.size = static_cast<int>(size);
  RETURN_ON_LITERAL(reader, cols, 3);
  f.name = reader.readString(cols[4].width);
  return true;
}

bool parseUserGfxLabel6(ColumnReader &reader, Record &record,
                        const ColumnSpec *cols) {
  reader.skipChars(cols[0].width);
  RETURN_ON_LITERAL(reader, cols, 1);
  record.as<UserGfxLabelRecord>().text = reader.readString(cols[2].width);
  return true;
}

bool parseUserGeneric(ColumnReader &reader, Record &record,
                      const ColumnSpec *cols) {
  reader.skipChars(cols[0].width);
  UserRecord &u = record.as<UserRecord>();
  u.subtype = reader.readString(cols[1].width);
  u.text = reader.readString(cols[2].width);
  return true;
}
} // namespace

std::optional<Record> Record::parse(std::string_view line,
                                    ParseContext &context) {
  (void)&context; // USER PDBRUN state handling lands with the USER records
  RecordType type = classifyLine(line);
  if (type == RecordType::UNKNOWN)
    return unknownFrom(line);
  if (type == RecordType::USER)
    type = classifyUserSubtype(line.substr(std::min<std::size_t>(6, line.size())),
                               context.pdbrunVersion);

  switch (type) {
  case RecordType::END:
  case RecordType::ENDMDL:
    return Record(type);
  case RecordType::AUTHOR: {
    ColumnReader reader(line);
    Record record(type);
    if (!parseContText(reader, record, kAuthorCols, &AuthorRecord::continuation,
                       &AuthorRecord::authors))
      return std::nullopt;
    return record;
  }
  case RecordType::TITLE: {
    ColumnReader reader(line);
    Record record(type);
    if (!parseContText(reader, record, kTitleCols, &TitleRecord::continuation,
                       &TitleRecord::text))
      return std::nullopt;
    return record;
  }
  case RecordType::COMPND: {
    ColumnReader reader(line);
    Record record(type);
    if (!parseContText(reader, record, kCompndCols, &CompndRecord::continuation,
                       &CompndRecord::compound))
      return std::nullopt;
    return record;
  }
  case RecordType::KEYWDS: {
    ColumnReader reader(line);
    Record record(type);
    if (!parseContText(reader, record, kKeywdsCols, &KeywdsRecord::continuation,
                       &KeywdsRecord::keywds))
      return std::nullopt;
    return record;
  }
  case RecordType::SOURCE: {
    ColumnReader reader(line);
    Record record(type);
    if (!parseContText(reader, record, kSourceCols, &SourceRecord::continuation,
                       &SourceRecord::name))
      return std::nullopt;
    return record;
  }
  case RecordType::MDLTYP: {
    ColumnReader reader(line);
    Record record(type);
    if (!parseContText(reader, record, kMdltypCols, &MdltypRecord::continuation,
                       &MdltypRecord::comment))
      return std::nullopt;
    return record;
  }
  case RecordType::ATOM:
  case RecordType::HETATM:
  case RecordType::SIGATM: {
    { // the standard 80-column layout
      ColumnReader reader(line);
      Record record(type);
      if (parseAtomColumns(reader, record, kAtomCols))
        return record;
    }
    { // the ZDATOM fallback: a fresh record, so the failed standard
      // attempt leaves no partially written fields behind
      ColumnReader reader(line);
      Record record(type);
      if (parseZdatomColumns(reader, record, kZdatomCols))
        return record;
    }
    return std::nullopt; // the caller falls back to the junk record
  }
  case RecordType::HELIX: {
    ColumnReader reader(line);
    Record record(type);
    if (!parseHelixColumns(reader, record, kHelixCols))
      return std::nullopt;
    return record;
  }
  case RecordType::SHEET: {
    ColumnReader reader(line);
    Record record(type);
    if (!parseSheetColumns(reader, record, kSheetCols))
      return std::nullopt;
    return record;
  }
  case RecordType::SITE: {
    ColumnReader reader(line);
    Record record(type);
    if (!parseSiteColumns(reader, record, kSiteCols))
      return std::nullopt;
    return record;
  }
  case RecordType::SSBOND: {
    ColumnReader reader(line);
    Record record(type);
    if (!parseSsbondColumns(reader, record, kSsbondCols))
      return std::nullopt;
    return record;
  }
  case RecordType::TER: {
    ColumnReader reader(line);
    Record record(type);
    if (!parseTerColumns(reader, record, kTerCols))
      return std::nullopt;
    return record;
  }
  case RecordType::TURN: {
    ColumnReader reader(line);
    Record record(type);
    if (!parseTurnColumns(reader, record, kTurnCols))
      return std::nullopt;
    return record;
  }
  case RecordType::CMPONT: {
    ColumnReader reader(line);
    Record record(type);
    if (!parseCmpontColumns(reader, record, kCmpontCols))
      return std::nullopt;
    return record;
  }
  case RecordType::CONECT: {
    ColumnReader reader(line);
    Record record(type);
    if (!parseConectColumns(reader, record, kConectCols))
      return std::nullopt;
    return record;
  }
  case RecordType::CRYST1: {
    ColumnReader reader(line);
    Record record(type);
    if (!parseCryst1Columns(reader, record, kCryst1Cols))
      return std::nullopt;
    return record;
  }
  case RecordType::MTRIX: {
    ColumnReader reader(line);
    Record record(type);
    if (!parseMtrixColumns(reader, record, kMtrixCols))
      return std::nullopt;
    return record;
  }
  case RecordType::ORIGX:
  case RecordType::SCALE: {
    ColumnReader reader(line);
    Record record(type);
    if (!parseOrigxScaleColumns(reader, record, kOrigxScaleCols,
                                type == RecordType::SCALE))
      return std::nullopt;
    return record;
  }
  case RecordType::SYMOP: {
    ColumnReader reader(line);
    Record record(type);
    if (!parseSymopColumns(reader, record, kSymopCols))
      return std::nullopt;
    return record;
  }
  case RecordType::REVDAT: {
    ColumnReader reader(line);
    Record record(type);
    if (!parseRevdatColumns(reader, record, kRevdatReadCols))
      return std::nullopt;
    return record;
  }
  case RecordType::SEQRES: {
    ColumnReader reader(line);
    Record record(type);
    if (!parseSeqresColumns(reader, record, kSeqresCols))
      return std::nullopt;
    return record;
  }
  case RecordType::SPLIT: {
    ColumnReader reader(line);
    Record record(type);
    if (!parseSplitColumns(reader, record, kSplitReadCols))
      return std::nullopt;
    return record;
  }
  case RecordType::OBSLTE:
  case RecordType::SPRSDE: {
    ColumnReader reader(line);
    Record record(type);
    if (!parseObslteSprsdeColumns(reader, record, kObslteReadCols,
                                  type == RecordType::SPRSDE))
      return std::nullopt;
    return record;
  }
  case RecordType::MASTER: {
    ColumnReader reader(line);
    Record record(type);
    if (!parseMasterColumns(reader, record, kMasterCols))
      return std::nullopt;
    return record;
  }
  case RecordType::AGGRGT: {
    ColumnReader reader(line);
    Record record(type);
    if (!parseAggrgtColumns(reader, record, kAggrgtCols))
      return std::nullopt;
    return record;
  }
  case RecordType::TRNSFM: {
    ColumnReader reader(line);
    Record record(type);
    if (!parseTrnsfmColumns(reader, record, kTrnsfmCols))
      return std::nullopt;
    return record;
  }
  case RecordType::TVECT: {
    ColumnReader reader(line);
    Record record(type);
    if (!parseTvectColumns(reader, record, kTvectCols))
      return std::nullopt;
    return record;
  }
  case RecordType::HET: {
    ColumnReader reader(line);
    Record record(type);
    if (!parseHetColumns(reader, record, kHetCols))
      return std::nullopt;
    return record;
  }
  case RecordType::FORMUL: {
    ColumnReader reader(line);
    Record record(type);
    if (!parseFormulColumns(reader, record, kFormulCols))
      return std::nullopt;
    return record;
  }
  case RecordType::HEADER: {
    ColumnReader reader(line);
    Record record(type);
    if (!parseHeaderColumns(reader, record, kHeaderCols))
      return std::nullopt;
    return record;
  }
  case RecordType::CAVEAT: {
    ColumnReader reader(line);
    Record record(type);
    if (!parseCaveatColumns(reader, record, kCaveatCols))
      return std::nullopt;
    return record;
  }
  case RecordType::JRNL: {
    ColumnReader reader(line);
    Record record(type);
    if (!parseJrnlColumns(reader, record, kJrnlCols))
      return std::nullopt;
    return record;
  }
  case RecordType::MODEL: {
    ColumnReader reader(line);
    Record record(type);
    if (!parseModelColumns(reader, record, kModelReadCols))
      return std::nullopt;
    return record;
  }
  case RecordType::ANISOU:
  case RecordType::SIGUIJ: {
    ColumnReader reader(line);
    Record record(type);
    const bool ok = parseAnisouColumns(reader, record, kAnisouCols);
    if (!ok)
      return std::nullopt;
    return record;
  }
  case RecordType::NUMMDL: {
    ColumnReader reader(line);
    Record record(type);
    if (!parseNummdlColumns(reader, record, kNummdlCols))
      return std::nullopt;
    return record;
  }
  case RecordType::EXPDTA: {
    ColumnReader reader(line);
    Record record(type);
    if (!parseExpdtaColumns(reader, record, kExpdtaCols))
      return std::nullopt;
    return record;
  }
  case RecordType::FTNOTE:
  case RecordType::REMARK:
  case RecordType::SYMDES:
  case RecordType::MTXDES:
  case RecordType::CMPDES:
  case RecordType::AGRDES: {
    ColumnReader reader(line);
    Record record(type);
    if (!parseAgrdesColumns(reader, record, kAgrdesCols))
      return std::nullopt;
    return record;
  }
  case RecordType::USER:
  case RecordType::USER_PDBRUN:
  case RecordType::USER_EYEPOS:
  case RecordType::USER_ATPOS:
  case RecordType::USER_WINDOW:
  case RecordType::USER_FOCUS:
  case RecordType::USER_VIEWPORT:
  case RecordType::USER_BGCOLOR:
  case RecordType::USER_ANGLE:
  case RecordType::USER_DISTANCE:
  case RecordType::USER_FILE:
  case RecordType::USER_MARKNAME:
  case RecordType::USER_MARK:
  case RecordType::USER_CNAME:
  case RecordType::USER_COLOR:
  case RecordType::USER_RADIUS:
  case RecordType::USER_OBJECT:
  case RecordType::USER_ENDOBJ:
  case RecordType::USER_CHAIN:
  case RecordType::USER_GFX_BEGIN:
  case RecordType::USER_GFX_END:
  case RecordType::USER_GFX_COLOR:
  case RecordType::USER_GFX_NORMAL:
  case RecordType::USER_GFX_VERTEX:
  case RecordType::USER_GFX_FONT:
  case RecordType::USER_GFX_TEXTPOS:
  case RecordType::USER_GFX_LABEL:
  case RecordType::USER_GFX_MOVE:
  case RecordType::USER_GFX_DRAW:
  case RecordType::USER_GFX_MARKER:
  case RecordType::USER_GFX_POINT: {
    const bool v5 = context.pdbrunVersion < 6;
    {
      ColumnReader reader(line);
      Record record(type);
      bool ok = false;
      switch (type) {
      case RecordType::USER:
        ok = parseUserGeneric(reader, record, kUserCols);
        break;
      case RecordType::USER_PDBRUN:
        ok = parseUserPdbrun(reader, record, kUserPdbrunCols);
        if (ok)
          context.pdbrunVersion = record.as<UserPdbrunRecord>().version;
        break;
      case RecordType::USER_EYEPOS:
      case RecordType::USER_ATPOS:
        ok = parseUserEyePos(reader, record, kUserEyePosCols);
        break;
      case RecordType::USER_WINDOW:
        ok = parseUserWindow(reader, record, kUserWindowCols);
        break;
      case RecordType::USER_FOCUS:
        ok = parseUserFocus(reader, record, kUserFocusCols);
        break;
      case RecordType::USER_VIEWPORT:
        ok = parseUserViewport(reader, record, kUserViewportCols);
        break;
      case RecordType::USER_BGCOLOR:
        if (v5)
          ok = parseUserBgColor5(line, record);
        else
          ok = parseUserBgColor6(reader, record, kUserBgColorCols);
        break;
      case RecordType::USER_ANGLE:
        ok = parseUserAngle(reader, record,
                            v5 ? kUserAngle5Cols : kUserAngle6Cols, v5);
        break;
      case RecordType::USER_DISTANCE:
        ok = parseUserDistance(reader, record,
                               v5 ? kUserDistance5Cols : kUserDistance6Cols,
                               v5);
        break;
      case RecordType::USER_FILE:
        ok = parseUserFile(reader, record,
                           v5 ? kUserFile5Cols : kUserFile6Cols, v5);
        break;
      case RecordType::USER_MARKNAME:
        if (v5) {
          ok = false; // v5 has no MARKNAME reader; it never classifies
        } else {
          ok = parseUserMarkname(reader, record, kUserMarkname6Cols);
        }
        break;
      case RecordType::USER_MARK:
        if (v5) {
          ok = false;
        } else {
          ok = parseUserMarkname(reader, record, kUserMark6Cols);
        }
        break;
      case RecordType::USER_CNAME:
        if (v5)
          ok = parseUserCName5(line, record, false);
        else
          ok = parseUserRgbName(reader, record, kUserCName6Cols);
        break;
      case RecordType::USER_COLOR:
        if (v5)
          ok = parseUserCName5(line, record, false);
        else
          ok = parseUserRgbName(reader, record, kUserCName6Cols);
        break;
      case RecordType::USER_RADIUS:
        ok = parseUserRadius(reader, record, kUserRadiusCols);
        break;
      case RecordType::USER_OBJECT:
      case RecordType::USER_ENDOBJ:
        if (v5)
          ok = parseUserObject(reader, record, kUserObjectCols);
        else
          ok = true; // v6 renders the subtype with no fields
        break;
      case RecordType::USER_CHAIN:
        if (v5)
          ok = parseUserChain5(line, record);
        else
          ok = parseUserChain6(reader, record, kUserChain6Cols);
        break;
      case RecordType::USER_GFX_BEGIN:
        if (v5)
          ok = false; // v5 has no GFX BEGIN reader; it never classifies
        else
          ok = parseUserGfxBegin6(reader, record, kUserGfxBegin6Cols);
        break;
      case RecordType::USER_GFX_END:
        ok = true; // v6: the subtype renders with no fields
        break;
      case RecordType::USER_GFX_COLOR:
        if (v5)
          ok = parseUserCName5(line, record, true);
        else
          ok = parseUserGfxColor6(reader, record, kUserGfxColor6Cols);
        break;
      case RecordType::USER_GFX_NORMAL:
      case RecordType::USER_GFX_VERTEX:
        if (v5)
          ok = false; // v5 has no readers; these never classify
        else
          ok = parseUserGfxPoint6(reader, record, kUserGfxPoint6Cols);
        break;
      case RecordType::USER_GFX_TEXTPOS:
        if (v5)
          ok = false;
        else
          ok = parseUserGfxPoint6(reader, record, kUserGfxTextPos6Cols);
        break;
      case RecordType::USER_GFX_FONT:
        if (v5)
          ok = parseUserGfxFont5(line, record);
        else
          ok = parseUserGfxFont6(reader, record, kUserGfxFont6Cols);
        break;
      case RecordType::USER_GFX_LABEL:
        if (v5)
          ok = parseUserGfxLabel5(line, record);
        else
          ok = parseUserGfxLabel6(reader, record, kUserGfxLabel6Cols);
        break;
      case RecordType::USER_GFX_MOVE:
      case RecordType::USER_GFX_DRAW:
      case RecordType::USER_GFX_MARKER:
      case RecordType::USER_GFX_POINT:
        if (v5) {
          const ColumnSpec *cols = type == RecordType::USER_GFX_MOVE
                                       ? kUserGfxMove5Cols
                                       : type == RecordType::USER_GFX_DRAW
                                             ? kUserGfxMove5Cols
                                             : type == RecordType::USER_GFX_MARKER
                                                   ? kUserGfxMarker5Cols
                                                   : kUserGfxPoint5Cols;
          ok = parseUserGfxPoint6(reader, record, cols);
        } else {
          ok = true; // v6 has no writers and no readers: fields stay unset
        }
        break;
      default:
        ok = false;
        break;
      }
      if (ok) {
        return record;
      }
      // the pinned fallback: a failed USER record parses with the generic
      // USER columns, keeping the subtype text
      {
        ColumnReader reader(line);
        Record generic(RecordType::USER);
        if (parseUserGeneric(reader, generic, kUserCols))
          return generic;
      }
      return std::nullopt;
    }
  }
  default:
    // Families whose layout tables land in later Phase 2 commits are not
    // parsed yet.
    return std::nullopt;
  }
}

Record Record::unknownFrom(std::string_view line) {
  Record record(RecordType::UNKNOWN);
  std::string junk(kJunkWidth, ' ');
  ColumnWriter writer(junk);
  writer.writeString(std::string(line), ColumnSpec{0, kJunkWidth, ColumnKind::Str});
  record.as<UnknownRecord>().junk = std::move(junk);
  return record;
}


// The USER record renderer: version-dependent prefixes, fixed-width and
// minimal field renders, the count-zero suppressions, the GFX BEGIN
// primitive substitution, and the v5 label quoting. Trailing blanks are
// trimmed by the caller's convention.
std::string userChars(const Record &record, const WriteContext &context) {
  const bool v5 = context.pdbrunVersion < 6;
  std::string out;
  const auto sep = [&out]() { out += ' '; };
  switch (record.type()) {
  case RecordType::USER: {
    const UserRecord &u = record.as<UserRecord>();
    out += "USER";
    out += renderFixedStr(u.subtype, 2, true);
    out += renderFixedStr(u.text, 66, true);
    break;
  }
  case RecordType::USER_PDBRUN: {
    const UserPdbrunRecord &p = record.as<UserPdbrunRecord>();
    out += "USER  PDBRUN ";
    out += renderFixedInt(p.version, 2, false, false);
    break;
  }
  case RecordType::USER_EYEPOS: {
    const UserEyePosRecord &e = record.as<UserEyePosRecord>();
    out += "USER  EYEPOS ";
    for (int i = 0; i < 3; ++i) {
      if (i > 0)
        sep();
      out += renderFixedReal(e.xyz[i], 9, 3);
    }
    break;
  }
  case RecordType::USER_ATPOS: {
    const UserEyePosRecord &e = record.as<UserEyePosRecord>();
    out += "USER  ATPOS  ";
    for (int i = 0; i < 3; ++i) {
      if (i > 0)
        sep();
      out += renderFixedReal(e.xyz[i], 9, 3);
    }
    break;
  }
  case RecordType::USER_WINDOW: {
    const UserWindowRecord &w = record.as<UserWindowRecord>();
    out += "USER  WINDOW ";
    const double values[6] = {w.left, w.right, w.bottom, w.top, w.hither,
                              w.yon};
    for (int i = 0; i < 6; ++i) {
      if (i > 0)
        sep();
      out += renderFixedReal(values[i], 9, 3);
    }
    break;
  }
  case RecordType::USER_FOCUS: {
    const UserFocusRecord &f = record.as<UserFocusRecord>();
    out += v5 ? "USER  FOCUS " : "USER  FOCUS  ";
    out += renderFixedReal(f.focus, 9, 3);
    break;
  }
  case RecordType::USER_VIEWPORT: {
    const UserViewportRecord &v = record.as<UserViewportRecord>();
    out += "USER  VIEWPORT ";
    const double values[4] = {v.xmin, v.xmax, v.ymin, v.ymax};
    for (int i = 0; i < 4; ++i) {
      if (i > 0)
        sep();
      out += renderFixedReal(values[i], 9, 3);
    }
    break;
  }
  case RecordType::USER_BGCOLOR: {
    const UserBgColorRecord &b = record.as<UserBgColorRecord>();
    out += "USER  BGCOLOR ";
    for (int i = 0; i < 3; ++i) {
      if (i > 0)
        sep();
      out += v5 ? renderMinimalReal(b.rgb[i], 3)
                : renderFixedReal(b.rgb[i], 5, 3);
    }
    break;
  }
  case RecordType::USER_ANGLE: {
    const UserAngleRecord &a = record.as<UserAngleRecord>();
    out += "USER  ANGLE ";
    if (v5) {
      out += renderFixedInt(a.which.value_or(0), 2, false, false);
      sep();
    }
    out += renderFixedInt(a.atom0, 6, false, false);
    sep();
    out += renderFixedInt(a.atom1, 6, false, false);
    sep();
    out += renderFixedInt(a.atom2, 6, false, false);
    sep();
    out += renderFixedInt(a.atom3, 6, false, false);
    sep();
    out += renderFixedReal(a.angle, v5 ? 8 : 9, 3);
    break;
  }
  case RecordType::USER_DISTANCE: {
    const UserDistanceRecord &d = record.as<UserDistanceRecord>();
    out += "USER  DISTANCE ";
    if (v5) {
      out += renderFixedInt(d.which.value_or(0), 2, false, false);
      sep();
    }
    out += renderFixedInt(d.atom0, 6, false, false);
    sep();
    out += renderFixedInt(d.atom1, 6, false, false);
    sep();
    out += renderFixedReal(d.distance, v5 ? 8 : 9, 3);
    break;
  }
  case RecordType::USER_FILE: {
    const UserFileRecord &f = record.as<UserFileRecord>();
    out += "USER  FILE ";
    if (v5) {
      out += renderFixedStr(f.filename, 61, true);
    } else {
      out += renderFixedInt(f.model.value_or(0), 4, false, false);
      sep();
      out += renderFixedStr(f.filename, 56, true);
    }
    break;
  }
  case RecordType::USER_MARKNAME:
    if (v5)
      return ""; // the v5 writer cannot emit MARKNAME
    out += "USER  MARKNAME ";
    out += renderFixedStr(record.as<UserMarknameRecord>().markname, 57, true);
    break;
  case RecordType::USER_MARK:
    if (v5)
      return ""; // the v5 writer cannot emit MARK
    out += "USER  MARK ";
    out += renderFixedStr(record.as<UserMarknameRecord>().markname, 57, true);
    break;
  case RecordType::USER_CNAME: {
    const UserCNameRecord &c = record.as<UserCNameRecord>();
    out += "USER  CNAME ";
    if (v5) {
      // the v5 libc %-.38s is a precision: truncate, never pad
      out += c.name.substr(0, 38);
      sep();
      out += renderMinimalReal(c.rgb[0], 3);
      sep();
      out += renderMinimalReal(c.rgb[1], 3);
      sep();
      out += renderMinimalReal(c.rgb[2], 3);
    } else {
      out += renderFixedReal(c.rgb[0], 5, 3);
      sep();
      out += renderFixedReal(c.rgb[1], 5, 3);
      sep();
      out += renderFixedReal(c.rgb[2], 5, 3);
      sep();
      out += renderFixedStr(c.name, 42, true);
    }
    break;
  }
  case RecordType::USER_COLOR: {
    const UserColorRecord &c = record.as<UserColorRecord>();
    out += "USER  COLOR ";
    if (v5) {
      out += c.spec.substr(0, 38);
      sep();
      out += renderMinimalReal(c.rgb[0], 3);
      sep();
      out += renderMinimalReal(c.rgb[1], 3);
      sep();
      out += renderMinimalReal(c.rgb[2], 3);
    } else {
      out += renderFixedReal(c.rgb[0], 5, 3);
      sep();
      out += renderFixedReal(c.rgb[1], 5, 3);
      sep();
      out += renderFixedReal(c.rgb[2], 5, 3);
      sep();
      out += renderFixedStr(c.spec, 42, true);
    }
    break;
  }
  case RecordType::USER_RADIUS: {
    const UserRadiusRecord &r = record.as<UserRadiusRecord>();
    out += "USER  RADIUS ";
    out += renderFixedReal(r.radius, 7, 3);
    break;
  }
  case RecordType::USER_OBJECT:
    if (v5) {
      out += "USER  OBJECT ";
      out += renderFixedInt(record.as<UserObjectRecord>().model, 3, false,
                            false);
    } else {
      out += "USER  OBJECT";
    }
    break;
  case RecordType::USER_ENDOBJ:
    if (v5) {
      out += "USER  ENDOBJ ";
      out += renderFixedInt(record.as<UserObjectRecord>().model, 3, false,
                            false);
    } else {
      out += "USER  ENDOBJ";
    }
    break;
  case RecordType::USER_CHAIN: {
    const UserChainRecord &c = record.as<UserChainRecord>();
    out += "USER  CHAIN ";
    if (v5) {
      out += renderMinimalInt(c.atom0);
      sep();
      out += renderMinimalInt(c.atom1);
    } else {
      out += renderFixedInt(c.atom0, 6, false, false);
      sep();
      out += renderFixedInt(c.atom1, 6, false, false);
    }
    break;
  }
  case RecordType::USER_GFX_BEGIN: {
    const UserGfxBeginRecord &b = record.as<UserGfxBeginRecord>();
    out += "USER  GFX BEGIN ";
    out += renderFixedStr(b.primitive == GfxType::UNKNOWN
                              ? b.unknown
                              : gfxName(b.primitive),
                          32, true);
    break;
  }
  case RecordType::USER_GFX_END:
    out += "USER  GFX END";
    break;
  case RecordType::USER_GFX_COLOR: {
    const UserColorRecord &c = record.as<UserColorRecord>();
    out += "USER  GFX COLOR ";
    if (v5) {
      out += c.spec.substr(0, 38); // %-.38s precision: truncate, no pad
      sep();
      out += renderMinimalReal(c.rgb[0], 3);
      sep();
      out += renderMinimalReal(c.rgb[1], 3);
      sep();
      out += renderMinimalReal(c.rgb[2], 3);
    } else {
      out += renderFixedReal(c.rgb[0], 5, 3);
      sep();
      out += renderFixedReal(c.rgb[1], 5, 3);
      sep();
      out += renderFixedReal(c.rgb[2], 5, 3);
      sep();
      out += c.spec; // the legacy v6 write uses a plain %s
    }
    break;
  }
  case RecordType::USER_GFX_NORMAL:
  case RecordType::USER_GFX_VERTEX:
  case RecordType::USER_GFX_TEXTPOS: {
    const UserGfxPointRecord &p = record.as<UserGfxPointRecord>();
    out += record.type() == RecordType::USER_GFX_NORMAL
               ? "USER  GFX NORMAL "
               : record.type() == RecordType::USER_GFX_VERTEX
                     ? "USER  GFX VERTEX "
                     : "USER  GFX TEXTPOS ";
    for (int i = 0; i < 3; ++i) {
      if (i > 0)
        sep();
      out += renderFixedReal(p.xyz[i], 9, 3);
    }
    break;
  }
  case RecordType::USER_GFX_FONT: {
    const UserGfxFontRecord &f = record.as<UserGfxFontRecord>();
    out += "USER  GFX FONT ";
    if (v5) {
      out += f.name.substr(0, 53); // %-.53s precision: truncate, no pad
      sep();
      out += renderMinimalInt(f.size);
    } else {
      out += renderFixedInt(f.size, 3, false, false);
      sep();
      out += renderFixedStr(f.name, 53, true);
    }
    break;
  }
  case RecordType::USER_GFX_LABEL: {
    const UserGfxLabelRecord &l = record.as<UserGfxLabelRecord>();
    out += "USER  GFX LABEL ";
    if (v5) {
      const std::array<double, 3> xyz = l.xyz.value_or(
          std::array<double, 3>{0.0, 0.0, 0.0});
      out += renderFixedReal(xyz[0], 9, 3);
      sep();
      out += renderFixedReal(xyz[1], 9, 3);
      sep();
      out += renderFixedReal(xyz[2], 9, 3);
      sep();
      out += '"';
      out += l.text.substr(0, 24); // %-.24s precision: truncate, no pad
      out += '"';
    } else {
      out += renderFixedStr(l.text, 56, true);
    }
    break;
  }
  case RecordType::USER_GFX_MOVE:
  case RecordType::USER_GFX_DRAW:
  case RecordType::USER_GFX_MARKER:
  case RecordType::USER_GFX_POINT:
    if (!v5)
      return ""; // the v6 writer cannot emit these
    {
      const UserGfxPointRecord &p = record.as<UserGfxPointRecord>();
      out += record.type() == RecordType::USER_GFX_MOVE
                 ? "USER  GFX MOVE "
                 : record.type() == RecordType::USER_GFX_DRAW
                       ? "USER  GFX DRAW "
                       : record.type() == RecordType::USER_GFX_MARKER
                             ? "USER  GFX MARKER "
                             : "USER  GFX POINT ";
      for (int i = 0; i < 3; ++i) {
        if (i > 0)
          sep();
        out += renderFixedReal(p.xyz[i], 9, 3);
      }
    }
    break;
  default:
    return "";
  }
  return out;
}

std::string Record::chars() const { return chars(WriteContext{}); }

std::string Record::chars(const WriteContext &context) const {
  char buf[83];
  const std::size_t length = write_into(buf, context);
  return std::string(buf, length);
}

std::size_t Record::write_into(char (&out)[83],
                               const WriteContext &context) const {
  if (type_ == RecordType::USER || type_ >= RecordType::USER_PDBRUN) {
    // the USER records render sequentially (version-dependent prefixes,
    // minimal and fixed-width fields); the layouts stay within the
    // 80-column record length
    const std::string text = userChars(*this, context);
    const std::size_t length = std::min(text.size(), sizeof(out) - 1);
    std::memcpy(out, text.data(), length);
    std::size_t trimmed = length;
    while (trimmed > 1 &&
           std::isspace(static_cast<unsigned char>(out[trimmed - 1])))
      --trimmed;
    out[trimmed] = '\0';
    return trimmed;
  }
  // The declared record length: 80 columns for every standard record, 82
  // for a ZDATOM-marked atom; the buffer accommodates the widest layout
  // plus a terminator.
  const std::size_t length =
      isAtom() && asAtom().zdatom.has_value() ? 82 : 80;
  static_assert(82 + 1 <= sizeof(out), "record buffer too small");
  std::memset(out, ' ', length);
  ColumnWriter writer(out, length);
  switch (type_) {
  case RecordType::UNKNOWN: {
    const UnknownRecord &payload = as<UnknownRecord>();
    writer.writeLiteral(0, "UNKNOWN:  ??");
    writer.writeString(payload.junk,
                       ColumnSpec{12, 6, ColumnKind::Str, 0, true, ' ', false, 6});
    writer.writeLiteral(18, "??");
    break;
  }
  case RecordType::END:
    writer.writeLiteral(0, "END");
    break;
  case RecordType::ENDMDL:
    writer.writeLiteral(0, "ENDMDL");
    break;
  case RecordType::ATOM:
  case RecordType::HETATM:
  case RecordType::SIGATM:
    writer.writeLiteral(0, type_ == RecordType::ATOM     ? "ATOM  "
                           : type_ == RecordType::HETATM ? "HETATM"
                                                         : "SIGATM");
    if (asAtom().zdatom.has_value()) {
      writeZdatomColumns(asAtom(), writer, kZdatomCols);
    } else {
      writeAtomColumns(asAtom(), writer, kAtomCols);
    }
    break;
  case RecordType::AUTHOR:
    assertPrefix("AUTHOR  ", kAuthorCols);
    writer.writeLiteral(0, "AUTHOR  ");
    writeContText(*this, writer, kAuthorCols, &AuthorRecord::continuation,
                  &AuthorRecord::authors);
    break;
  case RecordType::TITLE:
    assertPrefix("TITLE   ", kTitleCols);
    writer.writeLiteral(0, "TITLE   ");
    writeContText(*this, writer, kTitleCols, &TitleRecord::continuation,
                  &TitleRecord::text);
    break;
  case RecordType::COMPND:
    assertPrefix("COMPND ", kCompndCols);
    writer.writeLiteral(0, "COMPND ");
    writeContText(*this, writer, kCompndCols, &CompndRecord::continuation,
                  &CompndRecord::compound);
    break;
  case RecordType::KEYWDS:
    assertPrefix("KEYWDS  ", kKeywdsCols);
    writer.writeLiteral(0, "KEYWDS  ");
    writeContText(*this, writer, kKeywdsCols, &KeywdsRecord::continuation,
                  &KeywdsRecord::keywds);
    break;
  case RecordType::SOURCE:
    assertPrefix("SOURCE ", kSourceCols);
    writer.writeLiteral(0, "SOURCE ");
    writeContText(*this, writer, kSourceCols, &SourceRecord::continuation,
                  &SourceRecord::name);
    break;
  case RecordType::MDLTYP:
    assertPrefix("MDLTYP  ", kMdltypCols);
    writer.writeLiteral(0, "MDLTYP  ");
    writeContText(*this, writer, kMdltypCols, &MdltypRecord::continuation,
                  &MdltypRecord::comment);
    break;
  case RecordType::HELIX:
    assertPrefix("HELIX  ", kHelixCols);
    writer.writeLiteral(0, "HELIX  ");
    writeHelixColumns(as<HelixRecord>(), writer, kHelixCols);
    break;
  case RecordType::SHEET:
    assertPrefix("SHEET ", kSheetCols);
    writer.writeLiteral(0, "SHEET ");
    writeSheetColumns(as<SheetRecord>(), writer, kSheetCols);
    break;
  case RecordType::SITE:
    assertPrefix("SITE   ", kSiteCols);
    writer.writeLiteral(0, "SITE   ");
    writeSiteColumns(as<SiteRecord>(), writer, kSiteCols);
    break;
  case RecordType::SSBOND:
    assertPrefix("SSBOND ", kSsbondCols);
    writer.writeLiteral(0, "SSBOND ");
    writeSsbondColumns(as<SsbondRecord>(), writer, kSsbondCols);
    break;
  case RecordType::TER:
    assertPrefix("TER   ", kTerCols);
    writer.writeLiteral(0, "TER   ");
    writeTerColumns(as<TerRecord>(), writer, kTerCols);
    break;
  case RecordType::TURN:
    assertPrefix("TURN   ", kTurnCols);
    writer.writeLiteral(0, "TURN   ");
    writeTurnColumns(as<TurnRecord>(), writer, kTurnCols);
    break;
  case RecordType::CMPONT:
    assertPrefix("CMPONT ", kCmpontCols);
    writer.writeLiteral(0, "CMPONT ");
    writeCmpontColumns(as<CmpontRecord>(), writer, kCmpontCols);
    break;
  case RecordType::CONECT:
    assertPrefix("CONECT", kConectCols);
    writer.writeLiteral(0, "CONECT");
    writeConectColumns(as<ConectRecord>(), writer, kConectCols);
    break;
  case RecordType::CRYST1:
    assertPrefix("CRYST1", kCryst1Cols);
    writer.writeLiteral(0, "CRYST1");
    writeCryst1Columns(as<Cryst1Record>(), writer, kCryst1Cols);
    break;
  case RecordType::MTRIX:
    assertPrefix("MTRIX", kMtrixCols);
    writer.writeLiteral(0, "MTRIX");
    writeMtrixColumns(as<MtrixRecord>(), writer, kMtrixCols);
    break;
  case RecordType::ORIGX:
  case RecordType::SCALE:
    assertPrefix("ORIGX", kOrigxScaleCols);
    writer.writeLiteral(0, type_ == RecordType::ORIGX ? "ORIGX" : "SCALE");
    writeOrigxScaleColumns(*this, writer, kOrigxScaleCols,
                           type_ == RecordType::SCALE);
    break;
  case RecordType::SYMOP:
    assertPrefix("SYMOP", kSymopCols);
    writer.writeLiteral(0, "SYMOP");
    writeSymopColumns(as<SymopRecord>(), writer, kSymopCols);
    break;
  case RecordType::REVDAT:
    assertPrefix("REVDAT ", kRevdatWriteCols);
    writer.writeLiteral(0, "REVDAT ");
    writeRevdatColumns(as<RevdatRecord>(), writer, kRevdatWriteCols);
    break;
  case RecordType::SEQRES:
    assertPrefix("SEQRES", kSeqresCols);
    writer.writeLiteral(0, "SEQRES");
    writeSeqresColumns(as<SeqresRecord>(), writer, kSeqresCols);
    break;
  case RecordType::SPLIT:
    assertPrefix("SPLIT   ", kSplitWriteCols);
    writer.writeLiteral(0, "SPLIT   ");
    writeSplitColumns(as<SplitRecord>(), writer, kSplitWriteCols);
    break;
  case RecordType::OBSLTE:
  case RecordType::SPRSDE:
    assertPrefix("OBSLTE  ", kObslteWriteCols);
    writer.writeLiteral(0, type_ == RecordType::OBSLTE ? "OBSLTE  "
                                                       : "SPRSDE  ");
    writeObslteSprsdeColumns(*this, writer, kObslteWriteCols,
                             type_ == RecordType::SPRSDE);
    break;
  case RecordType::MASTER:
    assertPrefix("MASTER    ", kMasterCols);
    writer.writeLiteral(0, "MASTER    ");
    writeMasterColumns(as<MasterRecord>(), writer, kMasterCols);
    break;
  case RecordType::AGGRGT:
    assertPrefix("AGGRGT ", kAggrgtCols);
    writer.writeLiteral(0, "AGGRGT ");
    writeAggrgtColumns(as<AggrgtRecord>(), writer, kAggrgtCols);
    break;
  case RecordType::TRNSFM:
    assertPrefix("TRNSFM ", kTrnsfmCols);
    writer.writeLiteral(0, "TRNSFM ");
    writeTrnsfmColumns(as<TrnsfmRecord>(), writer, kTrnsfmCols);
    break;
  case RecordType::TVECT:
    assertPrefix("TVECT  ", kTvectCols);
    writer.writeLiteral(0, "TVECT  ");
    writeTvectColumns(as<TvectRecord>(), writer, kTvectCols);
    break;
  case RecordType::HET:
    assertPrefix("HET    ", kHetCols);
    writer.writeLiteral(0, "HET    ");
    writeHetColumns(as<HetRecord>(), writer, kHetCols);
    break;
  case RecordType::FORMUL:
    assertPrefix("FORMUL  ", kFormulCols);
    writer.writeLiteral(0, "FORMUL  ");
    writeFormulColumns(as<FormulRecord>(), writer, kFormulCols);
    break;
  case RecordType::HEADER:
    assertPrefix("HEADER    ", kHeaderWriteCols);
    writer.writeLiteral(0, "HEADER    ");
    writeHeaderColumns(as<HeaderRecord>(), writer, kHeaderWriteCols);
    break;
  case RecordType::CAVEAT:
    assertPrefix("CAVEAT  ", kCaveatCols);
    writer.writeLiteral(0, "CAVEAT  ");
    writeCaveatColumns(as<CaveatRecord>(), writer, kCaveatCols);
    break;
  case RecordType::JRNL:
    assertPrefix("JRNL        ", kJrnlCols);
    writer.writeLiteral(0, "JRNL        ");
    writeJrnlColumns(as<JrnlRecord>(), writer, kJrnlCols);
    break;
  case RecordType::MODEL:
    assertPrefix("MODEL    ", kModelWriteCols);
    writer.writeLiteral(0, "MODEL    ");
    writeModelColumns(as<ModelRecord>(), writer, kModelWriteCols);
    break;
  case RecordType::ANISOU:
  case RecordType::SIGUIJ:
    // the u columns write differently: legacy ANISOU renders %7d (a
    // parsed zero is a right-justified 0) while legacy SIGUIJ renders
    // %7D (a parsed zero is a blank field)
    if (type_ == RecordType::ANISOU) {
      assertPrefix("ANISOU", kAnisouCols);
      writer.writeLiteral(0, "ANISOU");
      writeAnisouColumns(as<AnisouRecord>(), writer, kAnisouCols);
    } else {
      assertPrefix("SIGUIJ", kSiguijWriteCols);
      writer.writeLiteral(0, "SIGUIJ");
      writeAnisouColumns(as<AnisouRecord>(), writer, kSiguijWriteCols);
    }
    break;
  case RecordType::NUMMDL:
    assertPrefix("NUMMDL    ", kNummdlWriteCols);
    writer.writeLiteral(0, "NUMMDL    ");
    writeNummdlColumns(as<NummdlRecord>(), writer, kNummdlWriteCols);
    break;
  case RecordType::EXPDTA:
    assertPrefix("EXPDTA  ", kExpdtaCols);
    writer.writeLiteral(0, "EXPDTA  ");
    writeExpdtaColumns(as<ExpdtaRecord>(), writer, kExpdtaCols);
    break;
  case RecordType::FTNOTE:
  case RecordType::REMARK:
  case RecordType::SYMDES:
  case RecordType::MTXDES:
  case RecordType::CMPDES:
  case RecordType::AGRDES: {
    static constexpr const char *kNames[] = {"FTNOTE ", "REMARK ", "SYMDES ",
                                             "MTXDES ", "CMPDES ", "AGRDES "};
    const bool blankZero =
        type_ == RecordType::FTNOTE || type_ == RecordType::REMARK;
    const ColumnSpec *wcols =
        blankZero ? kAgrdesCols : kAgrdesPlainWriteCols;
    const std::size_t base = static_cast<std::size_t>(
        (type_ == RecordType::FTNOTE    ? 0
         : type_ == RecordType::REMARK  ? 1
         : type_ == RecordType::SYMDES  ? 2
         : type_ == RecordType::MTXDES  ? 3
         : type_ == RecordType::CMPDES  ? 4
                                        : 5));
    assertPrefix(kNames[base], wcols);
    writer.writeLiteral(0, kNames[base]);
    writeAgrdesColumns(as<AgrdesRecord>(), writer, wcols);
    break;
  }
  default:
    // Families whose write tables land in later Phase 2 commits format as
    // empty.
    out[0] = '\0';
    return 0;
  }
  std::size_t trimmed = length;
  while (trimmed > 1 &&
         std::isspace(static_cast<unsigned char>(out[trimmed - 1])))
    --trimmed;
  out[trimmed] = '\0';
  return trimmed;
}

} // namespace zdock::pdb
