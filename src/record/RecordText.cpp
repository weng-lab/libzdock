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

constexpr ColumnSpec kSprsdeReadCols[] = {
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

constexpr ColumnSpec kSprsdeWriteCols[] = {
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
    std::int64_t value = 0;                                                   \
    if (!(reader).readInt((cols)[index].width, &value))                       \
      return false;                                                                         \
    *(out) = static_cast<int>(value);                                         \
  } while (false)

#define RETURN_ON_REAL(reader, cols, index, out)                              \
  do {                                                                        \
    double value = 0.0;                                                       \
    if (!(reader).readReal((cols)[index].width, &value))                      \
      return false;                                                                         \
    *(out) = value;                                                           \
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

} // namespace

std::optional<Record> Record::parse(std::string_view line,
                                    ParseContext &context) {
  (void)&context; // USER PDBRUN state handling lands with the USER records
  const RecordType type = classifyLine(line);
  if (type == RecordType::UNKNOWN)
    return unknownFrom(line);

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

std::string Record::chars() const { return chars(WriteContext{}); }

std::string Record::chars(const WriteContext &context) const {
  char buf[83];
  const std::size_t length = write_into(buf, context);
  return std::string(buf, length);
}

std::size_t Record::write_into(char (&out)[83],
                               const WriteContext & /*context*/) const {
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
    writer.writeLiteral(0, "ATOM  ");
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
    assertPrefix("HEADER    ", kHeaderCols);
    writer.writeLiteral(0, "HEADER    ");
    writeHeaderColumns(as<HeaderRecord>(), writer, kHeaderCols);
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
