// Copyright (c) 2019-2026 Arjan van der Velde, Weng Lab
// SPDX-License-Identifier: BSD-2-Clause

// Suite for the modern record class (Phase 2). At this stage it covers the
// RecordType enumeration only; parsing, writing, and parity tests with the
// legacy implementation arrive with the later Phase 2 commits.

#include "Test.hpp"
#include "Columns.hpp"
#include "FieldReader.hpp"
#include "Record.hpp"
#include "RecordClassifier.hpp"
#include "pdb++.h"

#include <algorithm>
#include <cmath>
#include <cstring>
#include <filesystem>
#include <fstream>
#include <random>
#include <set>
#include <string>
#include <utility>
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

// ---------------------------------------------------------------------------
// Bounded column reader and writer (Phase 2 section 4.2, commit 2)
// ---------------------------------------------------------------------------

namespace {

using libpdb::PDB;
using zdock::pdb::ColumnKind;
using zdock::pdb::ColumnReader;
using zdock::pdb::ColumnSpec;
using zdock::pdb::ColumnWriter;

// The legacy field reader and the new column reader shaped identically so
// one read-format walk can drive either in lockstep.
struct LegacyWalk {
  libpdb::FieldReader reader;
  char text[128] = {};

  bool literal(char c) { return reader.matchLiteral(c); }
  bool integer(int width, std::int64_t *value) {
    int v = 0;
    const bool ok = reader.readInt(width, &v);
    *value = v;
    return ok;
  }
  bool real(int width, double *value) { return reader.readReal(width, value); }
  void str(int width) { reader.readString(width, text); }
  void chars(int width) { reader.readChars(width, text); }
  void skip(int width) { reader.skipChars(width); }
};

struct ModernWalk {
  zdock::pdb::ColumnReader reader;

  bool literal(char c) { return reader.matchLiteral(c); }
  bool integer(int width, std::int64_t *value) {
    return reader.readInt(width, value);
  }
  bool real(int width, double *value) { return reader.readReal(width, value); }
  void str(int width) { (void)reader.readString(width); }
  void chars(int width) { (void)reader.readChars(width); }
  void skip(int width) { reader.skipChars(width); }
};

// One numeric field sample: which kind, the value the walk stored, and the
// exact field text both readers consumed.
struct FieldSample {
  bool isInt;
  std::int64_t integerValue = 0;
  double realValue = 0.0;
  std::string text;
};

// Strips the trailing blanks both readers strip before converting.
std::string strippedField(const std::string &text) {
  std::string stripped = text;
  while (!stripped.empty() &&
         std::isspace(static_cast<unsigned char>(stripped.back())))
    stripped.pop_back();
  return stripped;
}

// Walks a legacy read-format string over the walker, mirroring the
// production conversion order and width-bounded slicing; numeric samples
// are appended in conversion order. Returns false on the first hard
// mismatch (literal or numeric), as the production parse does.
template <typename W>
bool walkFormat(W &walk, const char *fmt, const std::string &line,
                std::vector<FieldSample> *samples) {
  std::size_t eol = 0;
  while (eol < line.size() && line[eol] != '\0' && line[eol] != '\n')
    ++eol;
  std::size_t pos = 0;
  for (; *fmt != '\0'; fmt++) {
    if (*fmt != '%') {
      const bool ok = walk.literal(*fmt);
      if (ok && pos < eol && line[pos] == *fmt)
        ++pos;
      if (!ok)
        return false;
      continue;
    }
    int width = 0;
    for (fmt++; std::isdigit(static_cast<unsigned char>(*fmt)); fmt++)
      width = width * 10 + *fmt - '0';
    if (width == 0)
      width = 1;
    const std::size_t start = pos;
    const std::size_t take =
        std::min<std::size_t>(width, pos < eol ? eol - pos : 0);
    switch (*fmt) {
    case 'd': {
      std::int64_t value = 0;
      if (!walk.integer(width, &value))
        return false;
      pos += take;
      samples->push_back({true, value, 0.0, line.substr(start, take)});
      break;
    }
    case 'f': {
      double value = 0.0;
      if (!walk.real(width, &value))
        return false;
      pos += take;
      samples->push_back({false, 0, value, line.substr(start, take)});
      break;
    }
    case 's':
      walk.str(width);
      pos += take;
      break;
    case 'c':
      walk.chars(width);
      pos += take;
      break;
    case ' ':
      walk.skip(width);
      pos += take;
      break;
    default:
      return false;
    }
  }
  return true;
}

// The legacy read-format tables (upstream .i data, verbatim).
char const *const kReadFormats[PDB::NUM_TYPES] = {
#include "read_format.i"
};

char const *const kPdbrun5Read[PDB::NUM_USER] = {
#include "pdbrun5_read.i"
};

char const *const kPdbrun6Read[PDB::NUM_USER] = {
#include "pdbrun6_read.i"
};

const char *readFormat(PDB::RecordType type, int pdbrunVersion) {
  if (type < PDB::USER_PDBRUN)
    return kReadFormats[type];
  const int index = type - PDB::USER_PDBRUN;
  return pdbrunVersion < 6 ? kPdbrun5Read[index] : kPdbrun6Read[index];
}

// Characterization fixture lines with the PDBRUN input version they were
// pinned under; the same lines the golden suite exercises.
const std::vector<std::pair<std::string, int>> &differentialLines() {
  static const std::vector<std::pair<std::string, int>> lines = {
      // standard records (input version 6)
      {"ATOM    352  N   MET b   1      20.321  -9.325  38.937  1.00 15.91"
       "           N",
       6},
      {"HETATM 9999  O   HOH A 999       1.500  -2.250   3.750  0.50 20.00"
       "           O",
       6},
      {"SIGATM    1  N   ALA A   1       0.020   0.020   0.020  0.00  1.00",
       6},
      {"ANISOU    1  N   ALA A   1     2120   1957   1835     -6     10"
       "     -9",
       6},
      {"SIGUIJ    1  N   ALA A   1     2120   1957   1835     -6     10"
       "     -9",
       6},
      {"AUTHOR   1VELDE, A.V.D.", 6},
      {"CAVEAT   1 1LYZ    SOME WARNING TEXT", 6},
      {"COMPND   1LYSOZYME", 6},
      {"JRNL        TO BE PUBLISHED", 6},
      {"KEYWDS   1LYSOZYME", 6},
      {"SOURCE   1HUMAN", 6},
      {"EXPDTA   1X-RAY DIFFRACTION", 6},
      {"MDLTYP   1TETRAMER", 6},
      {"TITLE    1LYSOZYME", 6},
      {"CRYST1   78.900   78.900   85.700  90.00  90.00  90.00 P 21 21 21"
       "    8",
       6},
      {"CONECT    1    2    3    4    5    6    7    8", 6},
      {"FORMUL   2  HOH  1*(H2 O)", 6},
      {"HEADER    TRANSFERASE24-JAN-95   1LYZ", 6},
      {"HELIX    1  H1 ALA A   10  LEU A   15  1                        "
       "             6",
       6},
      {"HET    NAG  A   1      14     N-ACETYL-D-GLUCOSAMINE", 6},
      {"MASTER        2    0    1    1    1    1    1    0  126    1   10"
       "    1",
       6},
      {"MODEL        1", 6},
      {"NUMMDL    2", 6},
      {"OBSLTE   1 24-JAN-95 1LYZ      2LYZ", 6},
      {"ORIGX1      1.000000  0.000000  0.000000        0.00000", 6},
      {"REVDAT   1   29-MAR-95 1LYZ    1       REVDAT", 6},
      {"SCALE1      0.017400  0.000000  0.000000        0.00000", 6},
      {"SEQRES   1 A    2  GLY  ALA", 6},
      {"SHEET    1  S1 5 ALA A  10  LEU A  15  0  N  ALA A  10   O  LEU"
       " A  15",
       6},
      {"SITE     1 CAT  2 GLY A  10  ALA A  15  VAL A  20  LEU A  25", 6},
      {"SPLIT    1 1ABZ 2CDZ", 6},
      {"SPRSDE   1 24-JAN-95 1LYZ      2LYZ", 6},
      {"SSBOND   1 CYS A   10    CYS A   25                       1555"
       "   1555  2.05",
       6},
      {"TER       1      ALA A  10", 6},
      {"TURN     1  T1 ALA A  10  LEU A  15     TURN", 6},
      {"TVECT    1   0.10000   0.20000   0.30000VECTOR", 6},
      {"CMPONT   1  ALAA    1   ALAA    2", 6},
      {"TRNSFM   1   1   1", 6},
      {"AGRDES   1                                             AN AGGREGATE"
       " DESCRIPTION",
       6},
      {"SYMDES   1                                               A SYMMETRY"
       " DESCRIPTION",
       6},
      {"MTXDES   1                                                 A MATRIX"
       " DESCRIPTION",
       6},
      {"CMPDES   1                                                A COMPONENT"
       " DESCRIPTION",
       6},
      {"AGGRGT   1   2   3   4   0   0   0   0   0   0   0   0   0   0   0"
       "   0",
       6},
      {"SYMOP1   1  1.000000  0.000000  0.000000        0.00000", 6},
      {"END", 6},
      {"ENDMDL", 6},
      {"USER  XX SOMETHING", 6},
      {"AUTHOR   1AAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAA"
       "AAAAAAAAAAAAA",
       6},
      {"REMARK   1 "
       "RRRRRRRRRRRRRRRRRRRRRRRRRRRRRRRRRRRRRRRRRRRRRRRRRRRRRRRRRRRRRRRRRR"
       "RR",
       6},
      {"TITLE    1"
       "TTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTT"
       "TT",
       6},
      // PDBRUN v6 USER records
      {"USER  PDBRUN  6", 6},
      {"USER  EYEPOS     1.000     2.000     3.000", 6},
      {"USER  ATPOS      4.000     5.000     6.000", 6},
      {"USER  WINDOW     1.000     2.000     3.000     4.000     5.000"
       "     6.000",
       6},
      {"USER  FOCUS      1.500", 6},
      {"USER  VIEWPORT     1.000     2.000     3.000     4.000", 6},
      {"USER  BGCOLOR 0.100 0.200 0.300", 6},
      {"USER  ANGLE      1      2      3      4   109.500", 6},
      {"USER  DISTANCE      1      2     3.500", 6},
      {"USER  FILE    1 /tmp/foo.pdb", 6},
      {"USER  MARKNAME MARK1", 6},
      {"USER  MARK MARK1", 6},
      {"USER  CNAME 1.000 0.000 0.000 RED", 6},
      {"USER  COLOR 0.000 1.000 0.000 GREEN", 6},
      {"USER  RADIUS   1.250", 6},
      {"USER  OBJECT", 6},
      {"USER  ENDOBJ", 6},
      {"USER  CHAIN      1     10", 6},
      {"USER  GFX BEGIN LINES", 6},
      {"USER  GFX END", 6},
      {"USER  GFX COLOR 1.000 0.000 0.000 RED", 6},
      {"USER  GFX NORMAL     0.000     0.000     1.000", 6},
      {"USER  GFX VERTEX     1.000     2.000     3.000", 6},
      {"USER  GFX FONT  10 Helvetica", 6},
      {"USER  GFX TEXTPOS     1.000     2.000     3.000", 6},
      {"USER  GFX LABEL LABEL TEXT", 6},
      // PDBRUN v5 USER records
      {"USER  PDBRUN  5", 5},
      {"USER  EYEPOS     1.000     2.000     3.000", 5},
      {"USER  ATPOS      4.000     5.000     6.000", 5},
      {"USER  WINDOW     1.000     2.000     3.000     4.000     5.000"
       "     6.000",
       5},
      {"USER  FOCUS     1.500", 5},
      {"USER  VIEWPORT     1.000     2.000     3.000     4.000", 5},
      {"USER  BGCOLOR 0.100 0.200 0.300", 5},
      {"USER  ANGLE  1      1      2      3      4  109.500", 5},
      {"USER  DISTANCE  1      1      2    3.500", 5},
      {"USER  FILE /tmp/foo.pdb", 5},
      {"USER  OBJECT   1", 5},
      {"USER  ENDOBJ   1", 5},
      {"USER  CHAIN 1 10", 5},
      {"USER  RADIUS   1.250", 5},
      {"USER  CNAME RED 0.100 0.200 0.300", 5},
      {"USER  COLOR RED 1.000 0.000 0.000", 5},
      {"USER  GFX COLOR RED 1.000 0.000 0.000", 5},
      {"USER  GFX FONT Helvetica 10", 5},
      {"USER  GFX MOVE     1.000     2.000     3.000", 5},
      {"USER  GFX DRAW     1.000     2.000     3.000", 5},
      {"USER  GFX MARKER     1.000     2.000     3.000", 5},
      {"USER  GFX POINT     1.000     2.000     3.000", 5},
      {"USER  MARKNAME MARK1", 5},
      {"USER  MARK MARK1", 5},
      {"USER  GFX BEGIN LINES", 5},
      {"USER  GFX END", 5},
      {"USER  GFX NORMAL     0.000     0.000     1.000", 5},
      {"USER  GFX VERTEX     1.000     2.000     3.000", 5},
      {"USER  GFX TEXTPOS     1.000     2.000     3.000", 5},
      // characterization extras: CRLF, charge, continuation, fallback
      {"ATOM    352  N   MET b   1      20.321  -9.325  38.937  1.00 15.91"
       "           N\r\n",
       6},
      {"ATOM      1  CA  ALA A   1       1.000   2.000   3.000  1.00 "
       "10.00\r\n",
       6},
      {"TER       1      ALA A  10\r\n", 6},
      {"SOME UNKNOWN JUNK TEXT\r\n", 6},
      {"ATOM      1  CA  ALA A   1       1.000   2.000   3.000  1     1 "
       "1.23     ALA  1.50",
       6},
      {"ATOM      1  CA  ALA A   1      12.345  -6.789  10.000  1     1 "
       "1.23     ALA  1.50",
       6},
      {"ATOM      1  CA  ALA A   1       1.000   2.000   3.000  1     1 "
       "1.23     ALA  0.25",
       6},
      {"ATOM      1  CA  ALA A   1       1.000   2.000   3.000  1     1 "
       "1.23     ALA -0.45",
       6},
      {"ATOM      1  CA  ALA A   1       1.000   2.000   3.000  1     1 "
       "1.23     ALA  -.15",
       6},
      {"ATOM      1  CA  ALA A   1       1.000   2.000   3.000  1     1 "
       "1.23     ALA -12.3",
       6},
      {"ATOM      1  N   ASN A   1      15.639 -11.284  -9.203  1     0 "
       "1.63          -.15",
       6},
      {"TITLE    2 SOME TEXT", 6},
      {"AUTHOR   2 SMITH, J.", 6},
      {"AUTHOR    0 SMITH, J.", 6},
      {"COMPND  10 TEN", 6},
      {"ATOM", 6},
      {"ATOM   1X", 6},
      {"atom      1  N   ALA A   1       1.000   2.000   3.000  1.00 "
       "10.00",
       6},
      {std::string(300, 'X'), 6},
      {std::string(120, 'X'), 6},
  };
  return lines;
}

// The three-way field agreement for one line: the legacy reader path, the
// new shared parser, and strtod/strtol on the same field text must agree
// bit-for-bit, and the match-count bookkeeping must agree as well.
void checkDifferentialLine(const std::string &line, int pdbrunVersion) {
  const PDB::RecordType type =
      libpdb::RecordClassifier::getType(line.c_str(), pdbrunVersion);
  const char *fmt = readFormat(type, pdbrunVersion);
  if (fmt == nullptr)
    return; // UNKNOWN, END, ENDMDL carry no fields

  LegacyWalk legacy{libpdb::FieldReader(line.c_str())};
  ModernWalk modern{ColumnReader(line)};
  std::vector<FieldSample> legacyFields;
  std::vector<FieldSample> modernFields;
  const bool legacyOk = walkFormat(legacy, fmt, line, &legacyFields);
  const bool modernOk = walkFormat(modern, fmt, line, &modernFields);

  CAPTURE(line);
  CAPTURE(pdbrunVersion);
  REQUIRE(legacyOk == modernOk);
  REQUIRE(legacy.reader.nmatch() == modern.reader.nmatch());
  REQUIRE(legacyFields.size() == modernFields.size());

  for (std::size_t i = 0; i < legacyFields.size(); ++i) {
    const FieldSample &l = legacyFields[i];
    const FieldSample &m = modernFields[i];
    CAPTURE(i);
    CAPTURE(l.text);
    // the field text as both parsers see it: trailing blanks stripped
    const std::string stripped = strippedField(l.text);
    if (l.isInt) {
      // secondary leg: strtol on the same field text
      char *end = nullptr;
      const long value =
          stripped.empty() ? 0L : std::strtol(stripped.c_str(), &end, 10);
      CAPTURE(l.integerValue);
      CAPTURE(m.integerValue);
      REQUIRE((stripped.empty() || *end == '\0'));
      REQUIRE(l.integerValue == m.integerValue);
      REQUIRE(static_cast<std::int64_t>(value) == m.integerValue);
    } else {
      // secondary leg: strtod on the same field text
      char *end = nullptr;
      const double value =
          stripped.empty() ? 0.0 : std::strtod(stripped.c_str(), &end);
      CAPTURE(l.realValue);
      CAPTURE(m.realValue);
      REQUIRE((stripped.empty() || *end == '\0'));
      REQUIRE(std::memcmp(&l.realValue, &m.realValue, sizeof(double)) == 0);
      REQUIRE(std::memcmp(&value, &m.realValue, sizeof(double)) == 0);
    }
  }
}

// Generates one field text in the supported fixed-format form: random
// signs, digit runs, point positions, blank runs, and all-blank fields.
std::string randomField(std::mt19937 &rng, int width, bool allowPoint) {
  std::uniform_int_distribution<int> digit(0, 9);
  std::uniform_int_distribution<int> flip(0, 99);
  if (flip(rng) < 8) // all-blank field
    return std::string(static_cast<std::size_t>(width), ' ');
  std::uniform_int_distribution<int> lenDist(1, width);
  int len = lenDist(rng);
  bool negative = flip(rng) < 25 && len >= 2;
  if (negative)
    --len;
  std::string core;
  if (negative)
    core += '-';
  int intDigits = len;
  int fracDigits = 0;
  if (allowPoint && len >= 2 && flip(rng) < 50) {
    std::uniform_int_distribution<int> split(0, len - 1);
    intDigits = split(rng);
    fracDigits = len - 1 - intDigits;
  }
  for (int i = 0; i < intDigits; i++)
    core += static_cast<char>('0' + digit(rng));
  if (allowPoint &&
      (fracDigits > 0 ||
       (intDigits > 0 && flip(rng) < 15 && core.size() < static_cast<std::size_t>(width))))
    core += '.';
  for (int i = 0; i < fracDigits; i++)
    core += static_cast<char>('0' + digit(rng));
  const std::size_t blanks = static_cast<std::size_t>(width) - core.size();
  std::uniform_int_distribution<std::size_t> lead(0, blanks);
  const std::size_t leading = lead(rng);
  return std::string(leading, ' ') + core +
         std::string(blanks - leading, ' ');
}

} // namespace

TEST_CASE("Column reader yields defaults at end-of-line", "[record-columns]") {
  // Invariant: fields starting at or past end-of-line yield the pinned
  // defaults (0, 0.0, empty string, spaces) and do not count as started.
  ColumnReader empty("");
  std::int64_t i = -1;
  double d = -1.0;
  REQUIRE(empty.readInt(5, &i));
  REQUIRE(i == 0);
  REQUIRE(empty.readReal(6, &d));
  REQUIRE(d == 0.0);
  REQUIRE(empty.nmatch() == 0);
  REQUIRE(empty.readString(4).empty());
  REQUIRE(empty.readChars(3) == "   ");

  ColumnReader shortLine("AB");
  REQUIRE(shortLine.readString(2) == "AB"); // consume the whole line
  i = -1;
  REQUIRE(shortLine.readInt(3, &i)); // the next field starts past EOL
  REQUIRE(i == 0);
  REQUIRE(shortLine.nmatch() == 1);
  REQUIRE(shortLine.readChars(4) == "    ");
}

TEST_CASE("Column reader counts started fields including skipped ones",
          "[record-columns]") {
  // Invariant: every field started on the line counts toward nmatch,
  // including ignored fields, and all-blank numeric fields parse as zero
  // while counting as started.
  ColumnReader r("ABCDEFGH");
  r.skipChars(3);
  REQUIRE(r.nmatch() == 1);
  REQUIRE(r.readString(5) == "DEFGH");
  REQUIRE(r.nmatch() == 2);

  ColumnReader blanks("      ");
  std::int64_t i = -1;
  double d = -1.0;
  REQUIRE(blanks.readInt(3, &i));
  REQUIRE(i == 0);
  REQUIRE(blanks.readReal(3, &d));
  REQUIRE(d == 0.0);
  REQUIRE(blanks.nmatch() == 2);
}

TEST_CASE("Column reader hard-mismatches invalid numeric fields",
          "[record-columns]") {
  // Invariant: a numeric field that does not fully match the fixed-format
  // grammar is a hard mismatch, and the accumulated prefix value is
  // stored exactly as the legacy strtol/strtod path stored it.
  ColumnReader partial("  1X");
  std::int64_t i = -1;
  REQUIRE_FALSE(partial.readInt(5, &i));
  REQUIRE(i == 1);

  ColumnReader point("12.5");
  i = -1;
  REQUIRE_FALSE(point.readInt(4, &i));
  REQUIRE(i == 12);

  ColumnReader loneSign("-");
  double d = -1.0;
  REQUIRE_FALSE(loneSign.readReal(1, &d));

  ColumnReader lonePoint(".");
  REQUIRE_FALSE(lonePoint.readReal(1, &d));

  ColumnReader innerBlank("1 5");
  REQUIRE_FALSE(innerBlank.readReal(3, &d));

  ColumnReader midCr("1\r5");
  REQUIRE_FALSE(midCr.readReal(3, &d));

  ColumnReader hexInt("0x1A");
  i = -1;
  REQUIRE_FALSE(hexInt.readInt(4, &i));
  REQUIRE(i == 0);
}

TEST_CASE("Column reader accepts the full supported numeric form",
          "[record-columns]") {
  // Invariant: signed, zero-less, point-led, and trailing-point forms
  // convert to the exact decimal value; a trailing carriage return is
  // stripped as the legacy reader stripped it.
  struct Case {
    const char *field;
    double value;
  };
  const std::vector<Case> cases = {
      {"  1.50", 1.5},   {"-.15", -0.15},   {"-0.00", -0.0},
      {"1.", 1.0},       {".5", 0.5},       {"+.5", 0.5},
      {"00012", 12.0},   {"9999999999", 9999999999.0},
      {"1.50\r", 1.5},   {"  0  ", 0.0},
  };
  for (const auto &c : cases) {
    CAPTURE(c.field);
    ColumnReader r(c.field);
    double d = -1.0;
    REQUIRE(r.readReal(static_cast<int>(std::strlen(c.field)), &d));
    REQUIRE(std::memcmp(&d, &c.value, sizeof(double)) == 0);
  }

  ColumnReader eighteen("999999999999999999"); // at the Int digit cap
  std::int64_t i = 0;
  REQUIRE(eighteen.readInt(18, &i));
  REQUIRE(i == 999999999999999999);

  ColumnReader nineteen("9999999999999999999"); // past the Int digit cap
  REQUIRE_FALSE(nineteen.readInt(19, &i));

  ColumnReader wideReal("99999999999"); // eleven digits: past the Real cap
  double d = 0.0;
  REQUIRE_FALSE(wideReal.readReal(11, &d));
}

TEST_CASE("Column reader pins the deliberate strtod divergences",
          "[record-columns]") {
  // Invariant: exponent, tab, inf/nan, and hex-float forms are rejected
  // even though legacy strtod accepted them; fixed-format PDB never
  // contains them, and the rejection is deliberate and pinned.
  struct Divergence {
    const char *field;
    bool real;
  };
  const std::vector<Divergence> divergences = {
      {"1e5", true},   {"\t1.5", true}, {"inf", true},
      {"nan", true},   {"0x1A", true},  {"1e5", false},
      {"\t12", false}, {"0x1A", false},
  };
  for (const auto &dv : divergences) {
    CAPTURE(dv.field);
    ColumnReader r(dv.field);
    std::int64_t i = 0;
    double d = 0.0;
    if (dv.real)
      REQUIRE_FALSE(r.readReal(static_cast<int>(std::strlen(dv.field)), &d));
    else
      REQUIRE_FALSE(r.readInt(static_cast<int>(std::strlen(dv.field)), &i));
  }
}

TEST_CASE("Column reader captures chars verbatim and strings trimmed",
          "[record-columns]") {
  // Invariant: %c-style fields capture up to the width verbatim with
  // space padding, and %s-style fields strip trailing blanks only.
  ColumnReader chars("ABCD");
  REQUIRE(chars.readChars(3) == "ABC");

  ColumnReader whole("ABCD");
  REQUIRE(whole.readString(4) == "ABCD");

  ColumnReader padded(" AB  ");
  REQUIRE(padded.readString(5) == " AB");

  ColumnReader cr("VELDE, A.V.D.\r");
  REQUIRE(cr.readString(14) == "VELDE, A.V.D.");

  ColumnReader literals("AXXX");
  REQUIRE(literals.matchLiteral('A'));
  REQUIRE(literals.readString(3) == "XXX");

  ColumnReader mismatch("BXXX");
  REQUIRE_FALSE(mismatch.matchLiteral('A'));

  ColumnReader atEol("");
  REQUIRE(atEol.matchLiteral('A')); // tolerated past end-of-line
}

TEST_CASE("Column writer reproduces the pinned integer emitters",
          "[record-columns]") {
  // Invariant: justification, blank/zero fill, %D zero-as-blank, sign
  // placement, and '*' overflow match the legacy sprintf semantics.
  std::string line(16, ' ');
  ColumnWriter w(line);
  zdock::pdb::ColumnSpec spec{};

  spec = {0, 5, ColumnKind::Int};
  w.writeInt(42, spec);
  REQUIRE(line.substr(0, 5) == "   42");

  spec = {0, 3, ColumnKind::Int, 0, false, ' ', true}; // %3D
  w.writeInt(0, spec);
  REQUIRE(line.substr(0, 3) == "   ");

  spec = {0, 5, ColumnKind::Int, 0, false, '0', false}; // %05d
  w.writeInt(42, spec);
  REQUIRE(line.substr(0, 5) == "00042");
  w.writeInt(-42, spec);
  REQUIRE(line.substr(0, 5) == "-0042");

  spec = {0, 5, ColumnKind::Int, 0, true}; // %-5d
  w.writeInt(42, spec);
  REQUIRE(line.substr(0, 5) == "42   ");

  spec = {0, 3, ColumnKind::Int}; // %3d overflow
  w.writeInt(12345, spec);
  REQUIRE(line.substr(0, 3) == "***");
  w.writeInt(-12345, spec);
  REQUIRE(line.substr(0, 3) == "***");

  spec = {0, 5, ColumnKind::Int, 0, false, ' ', true}; // %5D
  w.writeInt(5, spec);
  REQUIRE(line.substr(0, 5) == "    5");
  w.writeInt(-5, spec);
  REQUIRE(line.substr(0, 5) == "   -5");
}

TEST_CASE("Column writer reproduces the pinned real emitters",
          "[record-columns]") {
  // Invariant: fixed-point rendering with the legacy rounding, sign
  // placement, and '*' overflow conventions.
  std::string line(16, ' ');
  ColumnWriter w(line);
  zdock::pdb::ColumnSpec spec{};

  spec = {0, 8, ColumnKind::Real, 3}; // %8.3f
  w.writeReal(1.2344, spec);
  REQUIRE(line.substr(0, 8) == "   1.234");
  w.writeReal(-1.5, spec);
  REQUIRE(line.substr(0, 8) == "  -1.500");

  spec = {0, 4, ColumnKind::Real, 2}; // %4.2f overflow
  w.writeReal(123.45, spec);
  REQUIRE(line.substr(0, 4) == "****");

  spec = {0, 6, ColumnKind::Real, 2}; // %6.2f
  w.writeReal(-0.45, spec);
  REQUIRE(line.substr(0, 6) == " -0.45");

  spec = {0, 5, ColumnKind::Real, 2}; // %5.2f full-width negative
  w.writeReal(-1.5, spec);
  REQUIRE(line.substr(0, 5) == "-1.50");
}

TEST_CASE("Column writer reproduces the pinned string and char emitters",
          "[record-columns]") {
  // Invariant: precision truncation, justification, NUL as blank, and
  // '*' overflow match the legacy sprintf semantics.
  std::string line(16, ' ');
  ColumnWriter w(line);
  zdock::pdb::ColumnSpec spec{};

  spec = {0, 6, ColumnKind::Str, 0, true, ' ', false, 6}; // %-6.6s
  w.writeString("AB", spec);
  REQUIRE(line.substr(0, 6) == "AB    ");

  spec = {0, 2, ColumnKind::Str}; // %2s overflow
  w.writeString("ABC", spec);
  REQUIRE(line.substr(0, 2) == "**");

  spec = {0, 4, ColumnKind::Str}; // %4s right-justified
  w.writeString("AB", spec);
  REQUIRE(line.substr(0, 4) == "  AB");

  spec = {0, 3, ColumnKind::Char}; // %3c
  w.writeChar('X', spec);
  REQUIRE(line.substr(0, 3) == "  X");

  spec = {0, 2, ColumnKind::Char, 0, true}; // %-2c
  w.writeChar('\0', spec);
  REQUIRE(line.substr(0, 2) == "  ");
}

TEST_CASE("Three-way field agreement across the characterization corpus",
          "[record-differential]") {
  // Invariant: for every numeric field of every characterization fixture
  // line (under its pinned PDBRUN version), the legacy reader path, the
  // new shared fixed-format parser, and strtod/strtol on the same field
  // text agree bit-for-bit, and the match-count bookkeeping agrees.
  for (const auto &fixture : differentialLines()) {
    checkDifferentialLine(fixture.first, fixture.second);
  }
}

TEST_CASE("Three-way field agreement across the data files",
          "[record-differential]") {
  // Invariant: for every numeric field of every line under test/data,
  // parsed with the real read-format tables, the legacy reader path, the
  // new shared parser, and strtod/strtol agree bit-for-bit.
  std::vector<std::string> files;
  for (const auto &entry :
       std::filesystem::recursive_directory_iterator(test::getpath()))
    if (entry.is_regular_file())
      files.push_back(entry.path().string());
  std::sort(files.begin(), files.end());
  for (const auto &file : files) {
    std::ifstream in(file);
    REQUIRE(in);
    std::string line;
    while (std::getline(in, line))
      checkDifferentialLine(line, 6);
  }
}

TEST_CASE("Three-way field agreement on randomized supported fields",
          "[record-differential]") {
  // Invariant: for randomized supported-form field texts (random signs,
  // digit runs, point positions, blank runs, all-blank), the legacy
  // reader path, the new shared parser, and strtod/strtol agree
  // bit-for-bit, so a parser change can never silently alter stored data.
  std::mt19937 rng(20261104);
  for (int i = 0; i < 200000; ++i) {
    const int width = 1 + static_cast<int>(rng() % 10);
    const std::string text = randomField(rng, width, true);
    CAPTURE(i);
    CAPTURE(text);
    libpdb::FieldReader legacy(text.c_str());
    ColumnReader modern(text);
    double lv = -1.0;
    double mv = -1.0;
    REQUIRE(legacy.readReal(width, &lv));
    REQUIRE(modern.readReal(width, &mv));
    REQUIRE(std::memcmp(&lv, &mv, sizeof(double)) == 0);
    const std::string stripped = strippedField(text);
    char *end = nullptr;
    const double reference =
        stripped.empty() ? 0.0 : std::strtod(stripped.c_str(), &end);
    REQUIRE((stripped.empty() || *end == '\0'));
    REQUIRE(std::memcmp(&reference, &mv, sizeof(double)) == 0);
  }
  for (int i = 0; i < 100000; ++i) {
    const int width = 1 + static_cast<int>(rng() % 9);
    const std::string text = randomField(rng, width, false);
    CAPTURE(i);
    CAPTURE(text);
    libpdb::FieldReader legacy(text.c_str());
    ColumnReader modern(text);
    int lv = -1;
    std::int64_t mv = -1;
    REQUIRE(legacy.readInt(width, &lv));
    REQUIRE(modern.readInt(width, &mv));
    REQUIRE(static_cast<std::int64_t>(lv) == mv);
    const std::string stripped = strippedField(text);
    char *end = nullptr;
    const long reference =
        stripped.empty() ? 0L : std::strtol(stripped.c_str(), &end, 10);
    REQUIRE((stripped.empty() || *end == '\0'));
    REQUIRE(static_cast<std::int64_t>(reference) == mv);
  }
}
