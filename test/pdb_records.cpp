// Copyright (c) 2019-2026 Arjan van der Velde, Weng Lab
// SPDX-License-Identifier: BSD-2-Clause

// Characterization suite for the legacy libpdb++ implementation.
//
// Each test pins externally observable behavior of the 1992-format parser and
// formatter (including documented quirks) so that the modern rewrite in this
// repository can be validated against it record by record. Behavior pinned
// here is intentional; deliberate deviations require updating this suite in
// the same change.

#include "Test.hpp"
#include "pdb++.h"

#include <cstring>
#include <string>
#include <vector>

namespace {

using libpdb::PDB;

std::string spaces(int n) { return std::string(n, ' '); }

struct RecordFixture {
  PDB::RecordType type;
  std::string line;
};

// Canonical lines: exactly what parsing the line and re-emitting it with
// chars() produces (trailing blanks trimmed by the formatter).
const std::vector<RecordFixture> &standardFixtures() {
  static const std::vector<RecordFixture> fixtures = {
      {PDB::ATOM,
       "ATOM    352  N   MET b   1      20.321  -9.325  38.937  1.00 15.91"
       "           N"},
      {PDB::HETATM,
       "HETATM 9999  O   HOH A 999       1.500  -2.250   3.750  0.50 20.00"
       "           O"},
      {PDB::SIGATM,
       "SIGATM    1  N   ALA A   1       0.020   0.020   0.020  0.00  1.00"},
      {PDB::ANISOU,
       "ANISOU    1  N   ALA A   1     2120   1957   1835     -6     10"
       "     -9"},
      {PDB::SIGUIJ,
       "SIGUIJ    1  N   ALA A   1     2120   1957   1835     -6     10"
       "     -9"},
      {PDB::AUTHOR, "AUTHOR   1VELDE, A.V.D."},
      {PDB::CAVEAT, "CAVEAT   1 1LYZ    SOME WARNING TEXT"},
      {PDB::COMPND, "COMPND   1LYSOZYME"},
      {PDB::JRNL, "JRNL        TO BE PUBLISHED"},
      {PDB::KEYWDS, "KEYWDS   1LYSOZYME"},
      {PDB::SOURCE, "SOURCE   1HUMAN"},
      {PDB::EXPDTA, "EXPDTA   1X-RAY DIFFRACTION"},
      {PDB::MDLTYP, "MDLTYP   1TETRAMER"},
      {PDB::TITLE, "TITLE    1LYSOZYME"},
      {PDB::CRYST1,
       "CRYST1   78.900   78.900   85.700  90.00  90.00  90.00 P 21 21 21"
       "    8"},
      {PDB::CONECT, "CONECT    1    2    3    4    5    6    7    8"},
      {PDB::FORMUL, "FORMUL   2  HOH  1*(H2 O)"},
      {PDB::HEADER, "HEADER    TRANSFERASE24-JAN-95   1LYZ"},
      {PDB::HELIX,
       "HELIX    1  H1 ALA A   10  LEU A   15  1" + spaces(31) + "    6"},
      {PDB::HET, "HET    NAG  A   1      14     N-ACETYL-D-GLUCOSAMINE"},
      {PDB::MASTER,
       "MASTER        2    0    1    1    1    1    1    0  126    1   10"
       "    1"},
      {PDB::MODEL, "MODEL        1"},
      {PDB::NUMMDL, "NUMMDL    2"},
      {PDB::OBSLTE, "OBSLTE   1 24-JAN-95 1LYZ      2LYZ"},
      {PDB::ORIGX,
       "ORIGX1      1.000000  0.000000  0.000000        0.00000"},
      {PDB::REVDAT, "REVDAT   1   29-MAR-95 1LYZ    1       REVDAT"},
      {PDB::SCALE,
       "SCALE1      0.017400  0.000000  0.000000        0.00000"},
      {PDB::SEQRES, "SEQRES   1 A    2  GLY  ALA"},
      {PDB::SHEET,
       "SHEET    1  S1 5 ALA A  10  LEU A  15  0  N  ALA A  10   O  LEU"
       " A  15"},
      {PDB::SITE,
       "SITE     1 CAT  2 GLY A  10  ALA A  15  VAL A  20  LEU A  25"},
      {PDB::SPLIT, "SPLIT    1 1ABZ 2CDZ"},
      {PDB::SPRSDE, "SPRSDE   1 24-JAN-95 1LYZ      2LYZ"},
      {PDB::SSBOND,
       "SSBOND   1 CYS A   10    CYS A   25 " + spaces(23) +
           "  1555   1555  2.05"},
      {PDB::TER, "TER       1      ALA A  10"},
      {PDB::TURN, "TURN     1  T1 ALA A  10  LEU A  15     TURN"},
      {PDB::TVECT, "TVECT    1   0.10000   0.20000   0.30000VECTOR"},
      {PDB::CMPONT, "CMPONT   1  ALAA    1   ALAA    2"},
      {PDB::TRNSFM, "TRNSFM   1   1   1"},
      {PDB::AGRDES,
       "AGRDES   1 " + spaces(44) + "AN AGGREGATE DESCRIPTION"},
      {PDB::SYMDES,
       "SYMDES   1 " + spaces(46) + "A SYMMETRY DESCRIPTION"},
      {PDB::MTXDES,
       "MTXDES   1 " + spaces(48) + "A MATRIX DESCRIPTION"},
      {PDB::CMPDES,
       "CMPDES   1 " + spaces(45) + "A COMPONENT DESCRIPTION"},
      {PDB::AGGRGT,
       "AGGRGT   1   2   3   4   0   0   0   0   0   0   0   0   0   0   0"
       "   0"},
      {PDB::SYMOP,
       "SYMOP1   1  1.000000  0.000000  0.000000        0.00000"},
      {PDB::END, "END"},
      {PDB::ENDMDL, "ENDMDL"},
      {PDB::USER, "USER  XX SOMETHING"},
      {PDB::AUTHOR, "AUTHOR   1" + std::string(69, 'A')},
      {PDB::REMARK, "REMARK   1 " + std::string(68, 'R')},
      {PDB::TITLE, "TITLE    1" + std::string(70, 'T')},
  };
  return fixtures;
}

void roundTrip(const RecordFixture &fixture) {
  CAPTURE(fixture.line);
  PDB record(fixture.line.c_str());
  REQUIRE(record.type() == fixture.type);
  REQUIRE(std::string(record.chars()) == fixture.line);
}

// PDBRUN v6 USER records; the input version must be 6 (the default) for
// these to classify as their specific record types.
const std::vector<RecordFixture> &pdbrun6Fixtures() {
  static const std::vector<RecordFixture> fixtures = {
      {PDB::USER_PDBRUN, "USER  PDBRUN  6"},
      {PDB::USER_EYEPOS, "USER  EYEPOS     1.000     2.000     3.000"},
      {PDB::USER_ATPOS, "USER  ATPOS      4.000     5.000     6.000"},
      {PDB::USER_WINDOW,
       "USER  WINDOW     1.000     2.000     3.000     4.000     5.000"
       "     6.000"},
      {PDB::USER_FOCUS, "USER  FOCUS      1.500"},
      {PDB::USER_VIEWPORT,
       "USER  VIEWPORT     1.000     2.000     3.000     4.000"},
      {PDB::USER_BGCOLOR, "USER  BGCOLOR 0.100 0.200 0.300"},
      {PDB::USER_ANGLE, "USER  ANGLE      1      2      3      4   109.500"},
      {PDB::USER_DISTANCE, "USER  DISTANCE      1      2     3.500"},
      {PDB::USER_FILE, "USER  FILE    1 /tmp/foo.pdb"},
      {PDB::USER_MARKNAME, "USER  MARKNAME MARK1"},
      {PDB::USER_MARK, "USER  MARK MARK1"},
      {PDB::USER_CNAME, "USER  CNAME 1.000 0.000 0.000 RED"},
      {PDB::USER_COLOR, "USER  COLOR 0.000 1.000 0.000 GREEN"},
      {PDB::USER_RADIUS, "USER  RADIUS   1.250"},
      {PDB::USER_OBJECT, "USER  OBJECT"},
      {PDB::USER_ENDOBJ, "USER  ENDOBJ"},
      {PDB::USER_CHAIN, "USER  CHAIN      1     10"},
      {PDB::USER_GFX_BEGIN, "USER  GFX BEGIN LINES"},
      {PDB::USER_GFX_END, "USER  GFX END"},
      {PDB::USER_GFX_COLOR, "USER  GFX COLOR 1.000 0.000 0.000 RED"},
      {PDB::USER_GFX_NORMAL,
       "USER  GFX NORMAL     0.000     0.000     1.000"},
      {PDB::USER_GFX_VERTEX,
       "USER  GFX VERTEX     1.000     2.000     3.000"},
      {PDB::USER_GFX_FONT, "USER  GFX FONT  10 Helvetica"},
      {PDB::USER_GFX_TEXTPOS,
       "USER  GFX TEXTPOS     1.000     2.000     3.000"},
      {PDB::USER_GFX_LABEL, "USER  GFX LABEL LABEL TEXT"},
  };
  return fixtures;
}

// PDBRUN v5 USER records; the input version must be 5 while parsing these.
// Records without a v5 reader (MARKNAME, MARK, and most GFX records)
// degrade to plain USER records carrying the remainder of the line.
const std::vector<RecordFixture> &pdbrun5Fixtures() {
  static const std::vector<RecordFixture> fixtures = {
      {PDB::USER_PDBRUN, "USER  PDBRUN  5"},
      {PDB::USER_EYEPOS, "USER  EYEPOS     1.000     2.000     3.000"},
      {PDB::USER_ATPOS, "USER  ATPOS      4.000     5.000     6.000"},
      {PDB::USER_WINDOW,
       "USER  WINDOW     1.000     2.000     3.000     4.000     5.000"
       "     6.000"},
      {PDB::USER_FOCUS, "USER  FOCUS     1.500"},
      {PDB::USER_VIEWPORT,
       "USER  VIEWPORT     1.000     2.000     3.000     4.000"},
      {PDB::USER_BGCOLOR, "USER  BGCOLOR 0.100 0.200 0.300"},
      {PDB::USER_ANGLE,
       "USER  ANGLE  1      1      2      3      4  109.500"},
      {PDB::USER_DISTANCE, "USER  DISTANCE  1      1      2    3.500"},
      {PDB::USER_FILE, "USER  FILE /tmp/foo.pdb"},
      {PDB::USER_OBJECT, "USER  OBJECT   1"},
      {PDB::USER_ENDOBJ, "USER  ENDOBJ   1"},
      {PDB::USER_CHAIN, "USER  CHAIN 1 10"},
      {PDB::USER_RADIUS, "USER  RADIUS   1.250"},
      {PDB::USER_CNAME, "USER  CNAME RED 0.100 0.200 0.300"},
      {PDB::USER_COLOR, "USER  COLOR RED 1.000 0.000 0.000"},
      {PDB::USER_GFX_COLOR, "USER  GFX COLOR RED 1.000 0.000 0.000"},
      {PDB::USER_GFX_FONT, "USER  GFX FONT Helvetica 10"},
      {PDB::USER_GFX_MOVE, "USER  GFX MOVE     1.000     2.000     3.000"},
      {PDB::USER_GFX_DRAW, "USER  GFX DRAW     1.000     2.000     3.000"},
      {PDB::USER_GFX_MARKER,
       "USER  GFX MARKER     1.000     2.000     3.000"},
      {PDB::USER_GFX_POINT, "USER  GFX POINT     1.000     2.000     3.000"},
      // v5 has no dedicated reader/writer: these parse and format as USER.
      {PDB::USER, "USER  MARKNAME MARK1"},
      {PDB::USER, "USER  MARK MARK1"},
      {PDB::USER, "USER  GFX BEGIN LINES"},
      {PDB::USER, "USER  GFX END"},
      {PDB::USER, "USER  GFX NORMAL     0.000     0.000     1.000"},
      {PDB::USER, "USER  GFX VERTEX     1.000     2.000     3.000"},
      {PDB::USER, "USER  GFX TEXTPOS     1.000     2.000     3.000"},
  };
  return fixtures;
}

} // namespace

TEST_CASE("Standard records round-trip through parse and format",
          "[pdb-records]") {
  // Invariant: every standard record line is reproduced byte-for-byte by a
  // parse/format cycle, with trailing blanks trimmed by the formatter.
  PDB::PdbrunInputVersion(6);
  PDB::PdbrunOutputVersion(6);
  for (const auto &fixture : standardFixtures()) {
    roundTrip(fixture);
  }
}

TEST_CASE("PDBRUN v6 user records round-trip", "[pdb-records]") {
  // Invariant: v6 USER records keep their specific types and text across a
  // parse/format cycle under the default PDBRUN input version.
  PDB::PdbrunInputVersion(6);
  PDB::PdbrunOutputVersion(6);
  for (const auto &fixture : pdbrun6Fixtures()) {
    roundTrip(fixture);
  }
  PDB::PdbrunInputVersion(6);
  PDB::PdbrunOutputVersion(6);
}

TEST_CASE("PDBRUN v5 user records round-trip", "[pdb-records]") {
  // Invariant: v5 USER records round-trip under v5, and records without a
  // v5 reader degrade to plain USER records carrying the line remainder.
  PDB::PdbrunInputVersion(5);
  PDB::PdbrunOutputVersion(5);
  for (const auto &fixture : pdbrun5Fixtures()) {
    roundTrip(fixture);
  }
  PDB::PdbrunInputVersion(6);
  PDB::PdbrunOutputVersion(6);
}

TEST_CASE("v5 formatting suppresses records the v5 writer cannot emit",
          "[pdb-records]") {
  // Invariant: USER MARKNAME formats as an empty string under v5 output.
  PDB::PdbrunOutputVersion(5);
  PDB record(PDB::USER_MARKNAME);
  std::strncpy(record.userMarkname.markname, "MARK1", 6);
  REQUIRE(std::string(record.chars()) == "");
  PDB::PdbrunOutputVersion(6);
}

TEST_CASE("v5 GFX LABEL quoting degenerates across cycles",
          "[pdb-records]") {
  // Invariant (pinned legacy quirk): the v5 reader captures the quotes into
  // the label text and the v5 writer wraps the text in quotes again, so each
  // parse/format cycle adds one pair of quotes; no stable form exists.
  PDB::PdbrunInputVersion(5);
  PDB::PdbrunOutputVersion(5);
  PDB label("USER  GFX LABEL     1.000     2.000     3.000 \"LABEL TEXT\"");
  REQUIRE(label.type() == PDB::USER_GFX_LABEL);
  REQUIRE(label.userGfxLabel.xyz[0] == 1.0);
  REQUIRE(std::string(label.userGfxLabel.text) == "\"LABEL TEXT\"");
  REQUIRE(std::string(label.chars()) ==
          "USER  GFX LABEL     1.000     2.000     3.000 \"\"LABEL TEXT\"\"");
  PDB::PdbrunInputVersion(6);
  PDB::PdbrunOutputVersion(6);
}

TEST_CASE("v5 libc %s reads are bounded to their destinations",
          "[pdb-records]") {
  // Invariant: the four v5 libc sscanf paths whose %s conversions had no
  // width (USER CNAME/COLOR/GFX COLOR name/spec[39] and USER GFX FONT
  // name[54]) now read at most destination-size-minus-one characters, so a
  // long token cannot run past the record object (an AddressSanitizer
  // stack-buffer-overflow before the bound); the numeric conversions that
  // follow no longer match, and their members keep their defaults.
  const std::string token(200, 'A');
  struct LongTokenCase {
    std::string line;
    PDB::RecordType type;
    int bound; // characters captured into the destination
  };
  const std::vector<LongTokenCase> cases = {
      {"USER  CNAME " + token, PDB::USER_CNAME, 38},
      {"USER  COLOR " + token, PDB::USER_COLOR, 38},
      {"USER  GFX COLOR " + token, PDB::USER_GFX_COLOR, 38},
      {"USER  GFX FONT " + token, PDB::USER_GFX_FONT, 53},
  };
  PDB::PdbrunInputVersion(5);
  PDB::PdbrunOutputVersion(5);
  for (const auto &c : cases) {
    CAPTURE(c.line);
    PDB record(c.line.c_str());
    REQUIRE(record.type() == c.type);
    const char *text = nullptr;
    if (record.type() == PDB::USER_CNAME) {
      text = record.userCName.name;
      REQUIRE(record.userCName.rgb[0] == 0.0);
      REQUIRE(record.userCName.rgb[1] == 0.0);
      REQUIRE(record.userCName.rgb[2] == 0.0);
    } else if (record.type() == PDB::USER_COLOR) {
      text = record.userColor.spec;
      REQUIRE(record.userColor.rgb[0] == 0.0);
      REQUIRE(record.userColor.rgb[1] == 0.0);
      REQUIRE(record.userColor.rgb[2] == 0.0);
    } else if (record.type() == PDB::USER_GFX_COLOR) {
      text = record.userGfxColor.spec;
      REQUIRE(record.userGfxColor.rgb[0] == 0.0);
      REQUIRE(record.userGfxColor.rgb[1] == 0.0);
      REQUIRE(record.userGfxColor.rgb[2] == 0.0);
    } else {
      text = record.userGfxFont.name;
      REQUIRE(record.userGfxFont.size == 0);
    }
    REQUIRE(std::string(text) == std::string(c.bound, 'A'));
  }

  // Short tokens parse exactly as before the bound: the whole token is
  // captured and the numeric conversions still match.
  PDB cname("USER  CNAME RED 0.100 0.200 0.300");
  REQUIRE(cname.type() == PDB::USER_CNAME);
  REQUIRE(std::string(cname.userCName.name) == "RED");
  REQUIRE(cname.userCName.rgb[0] == 0.1);
  PDB font("USER  GFX FONT Helvetica 10");
  REQUIRE(font.type() == PDB::USER_GFX_FONT);
  REQUIRE(std::string(font.userGfxFont.name) == "Helvetica");
  REQUIRE(font.userGfxFont.size == 10);
  PDB::PdbrunInputVersion(6);
  PDB::PdbrunOutputVersion(6);
}

TEST_CASE("Atom records parse their fixed columns", "[pdb-records]") {
  // Invariant: ATOM fields land in the struct members at their wwPDB
  // columns; string fields are right-trimmed on read (" N  " -> " N").
  PDB record("ATOM    352  N   MET b   1      20.321  -9.325  38.937  1.00"
             " 15.91           N");
  REQUIRE(record.type() == PDB::ATOM);
  REQUIRE(record.atom.serialNum == 352);
  REQUIRE(std::string(record.atom.name) == " N");
  REQUIRE(record.atom.altLoc == ' ');
  REQUIRE(std::string(record.atom.residue.name) == "MET");
  REQUIRE(record.atom.residue.chainId == 'b');
  REQUIRE(record.atom.residue.seqNum == 1);
  REQUIRE(record.atom.residue.insertCode == ' ');
  REQUIRE(record.atom.xyz[0] == 20.321);
  REQUIRE(record.atom.xyz[1] == -9.325);
  REQUIRE(record.atom.xyz[2] == 38.937);
  REQUIRE(record.atom.occupancy == 1.0);
  REQUIRE(record.atom.tempFactor == 15.91);
  REQUIRE(std::string(record.atom.element) == " N");
  REQUIRE(record.atom.charge[0] == '\0');
  REQUIRE_FALSE(record.atom.iszdatom);
}

TEST_CASE("Record classification is case-insensitive", "[pdb-records]") {
  // Invariant: the first four characters are upper-cased before dispatch.
  PDB record("atom      1  N   ALA A   1       1.000   2.000   3.000  1.00"
             " 10.00");
  REQUIRE(record.type() == PDB::ATOM);
  REQUIRE(record.atom.serialNum == 1);
}

TEST_CASE("Malformed atom lines fall back to unknown records",
          "[pdb-records]") {
  // Invariant: a hard field mismatch falls back ATOM -> ZDATOM -> UNKNOWN,
  // and the unknown content is right-aligned junk in a 72-char field.
  PDB record("ATOM   1X");
  REQUIRE(record.type() == PDB::UNKNOWN);
  REQUIRE(std::string(record.chars()) == "UNKNOWN:  ??      ??");
}

TEST_CASE("Truncated records parse with field defaults", "[pdb-records]") {
  // Invariant: a record that ends early yields zero defaults instead of a
  // mismatch, so short ATOM lines stay ATOM records.
  PDB record("ATOM");
  REQUIRE(record.type() == PDB::ATOM);
  REQUIRE(record.atom.serialNum == 0);
  REQUIRE(record.atom.xyz[0] == 0.0);
  REQUIRE(record.atom.occupancy == 0.0);
  REQUIRE_FALSE(record.atom.iszdatom);
}

TEST_CASE("ZDOCK atom records fall back to the marksur layout",
          "[pdb-records]") {
  // Invariant: an ATOM line that fails the standard layout parses with the
  // ZDATOM extension columns, keeps the ATOM record type, and round-trips
  // with the segment id right-justified by the writer.
  const std::string line =
      "ATOM      1  CA  ALA A   1       1.000   2.000   3.000  1     1 1.23"
      "     ALA  1.50";
  PDB record(line.c_str());
  CAPTURE(line);
  REQUIRE(record.type() == PDB::ATOM);
  REQUIRE(record.atom.iszdatom);
  REQUIRE(record.atom.serialNum == 1);
  REQUIRE(record.atom.xyz[0] == 1.0);
  REQUIRE(record.atom.type == 1);
  REQUIRE(record.atom.surface == 1);
  REQUIRE(record.atom.rad == 1.23);
  REQUIRE(std::string(record.atom.segid) == " ALA");
  REQUIRE(record.atom.chg == 1.5);
  // Quirk: the failed standard-layout attempt already wrote occupancy from
  // its first %6f field, and the ZDATOM layout never resets it.
  REQUIRE(record.atom.occupancy == 1.0);
  // The ZDATOM layout keeps all 82 columns on output; the charge renders
  // in full at columns 78-82.
  REQUIRE(std::string(record.chars()) ==
          "ATOM      1  CA  ALA A   1       1.000   2.000   3.000  1     1"
          " 1.23     ALA  1.50");
}

TEST_CASE("MODEL records tolerate the non-standard column layout",
          "[pdb-records]") {
  // Invariant: the model number is read from columns 7-14, matching the
  // non-standard layout used by repository fixtures.
  PDB record("MODEL        1");
  REQUIRE(record.type() == PDB::MODEL);
  REQUIRE(record.model.num == 1);
}

TEST_CASE("Full-width ZDATOM record stays within the record buffer",
          "[pdb-records]") {
  // Invariant: a ZDATOM record whose layout renders 82 columns writes no
  // byte past the record buffer, and its text keeps all 82 columns
  // including the full charge field at columns 78-82.
  const std::string line = "ATOM      1  CA  ALA A   1      12.345  -6.789"
                           "  10.000  1     1 1.23     ALA  1.50";
  REQUIRE(line.size() == 82);
  PDB record(line.c_str());
  REQUIRE(record.type() == PDB::ATOM);
  REQUIRE(record.atom.iszdatom);
  REQUIRE(record.atom.chg == 1.5);
  const std::string text = record.chars();
  REQUIRE(text.size() == 82);
  REQUIRE(text == "ATOM      1  CA  ALA A   1      12.345  -6.789  10.000 "
                  " 1     1 1.23     ALA  1.50");
}

TEST_CASE("ZDATOM charges survive a write and re-read", "[pdb-records]") {
  // Invariant: full-width ZDATOM records render the charge across all
  // five columns 78-82 of the 82-column layout and the rendered text
  // parses back to the same charge; a charge too wide for the field
  // degrades to the formatter's asterisk convention and the line no
  // longer parses (the zero-less "-.15" input form normalizes to "-0.15").
  const std::string base = "ATOM      1  CA  ALA A   1       1.000   2.000"
                           "   3.000  1     1 1.23     ALA ";
  REQUIRE(base.size() == 77);
  struct ChargeCase {
    std::string line;
    std::string rendered;
    double charge;
    bool parses;
  };
  const std::vector<ChargeCase> cases = {
      {base + " 1.50", base + " 1.50", 1.5, true},
      {base + " 0.25", base + " 0.25", 0.25, true},
      {base + "-0.45", base + "-0.45", -0.45, true},
      {base + " -.15", base + "-0.15", -0.15, true},
      {base + "-12.3", base + "*****", 0.0, false},
  };
  for (const auto &c : cases) {
    CAPTURE(c.line);
    REQUIRE(c.line.size() == 82);
    PDB record(c.line.c_str());
    REQUIRE(record.type() == PDB::ATOM);
    REQUIRE(record.atom.iszdatom);
    REQUIRE(std::string(record.chars()) == c.rendered);
    PDB reread(record.chars());
    REQUIRE(reread.type() == (c.parses ? PDB::ATOM : PDB::UNKNOWN));
    if (c.parses) {
      REQUIRE(reread.atom.iszdatom);
      REQUIRE(reread.atom.chg == c.charge);
    }
  }
}

TEST_CASE("PDBRUN versions are global mutable parse state", "[pdb-records]") {
  // Invariant (pinned legacy quirk): parsing a USER PDBRUN record mutates
  // process-global input version state observed by every later parse.
  PDB::PdbrunInputVersion(6);
  REQUIRE(PDB::PdbrunInputVersion() == 6);
  PDB record("USER  PDBRUN  5");
  REQUIRE(record.type() == PDB::USER_PDBRUN);
  REQUIRE(record.userPdbrun.version == 5);
  REQUIRE(PDB::PdbrunInputVersion() == 5);
  PDB::PdbrunInputVersion(6);
  REQUIRE(PDB::PdbrunInputVersion() == 6);
}

TEST_CASE("Formatting a PDBRUN record mutates global output version",
          "[pdb-records]") {
  // Invariant (pinned legacy quirk): chars() on USER_PDBRUN assigns the
  // process-global output version as a side effect.
  PDB::PdbrunOutputVersion(6);
  PDB record(PDB::USER_PDBRUN);
  record.userPdbrun.version = 5;
  (void)record.chars();
  REQUIRE(PDB::PdbrunOutputVersion() == 5);
  PDB::PdbrunOutputVersion(6);
}

TEST_CASE("chars() returns a shared static buffer", "[pdb-records]") {
  // Invariant (pinned legacy quirk): successive chars() calls reuse one
  // buffer, so previously returned pointers change content.
  PDB atom("ATOM    352  N   MET b   1      20.321  -9.325  38.937  1.00"
           " 15.91           N");
  PDB hetatm("HETATM 9999  O   HOH A 999       1.500  -2.250   3.750  0.50"
             " 20.00           O");
  const char *stale = atom.chars();
  REQUIRE(std::string(stale).substr(0, 6) == "ATOM  ");
  (void)hetatm.chars();
  REQUIRE(std::string(stale).substr(0, 6) == "HETATM");
}

TEST_CASE("Fresh records carry documented defaults", "[pdb-records]") {
  // Invariant: type() byte-clears the record; ATOM occupancy defaults to
  // 1.0 and everything else to zero.
  PDB atom(PDB::ATOM);
  REQUIRE(atom.type() == PDB::ATOM);
  REQUIRE(atom.atom.occupancy == 1.0);
  REQUIRE(atom.atom.serialNum == 0);
  REQUIRE(atom.atom.xyz[0] == 0.0);

  PDB unknown;
  REQUIRE(unknown.type() == PDB::UNKNOWN);
  REQUIRE(unknown.unknown.junk[0] == '\0');
}

TEST_CASE("Records compare by type and bytes", "[pdb-records]") {
  // Invariant: equality requires the same record type and identical bytes.
  PDB a(PDB::ATOM);
  a.atom.serialNum = 7;
  PDB b(PDB::ATOM);
  b.atom.serialNum = 7;
  REQUIRE(a == b);
  b.atom.serialNum = 8;
  REQUIRE(a != b);
}

TEST_CASE("isalpha() matches parsed CA atoms", "[pdb-records]") {
  // Invariant: parsed atom names are right-trimmed, so a CA atom stores
  // " CA" and the 5-byte comparison against the " CA" literal succeeds;
  // blank-padded hand-built names do not match.
  PDB parsed("ATOM      1  CA  ALA A   1       1.000   2.000   3.000  1.00"
             " 10.00");
  REQUIRE(parsed.type() == PDB::ATOM);
  REQUIRE(std::string(parsed.atom.name) == " CA");
  REQUIRE(parsed.isalpha());

  PDB padded(PDB::ATOM);
  std::memcpy(padded.atom.name, " CA ", 4);
  REQUIRE_FALSE(padded.isalpha());
}

TEST_CASE("Overlong input lines fail the record stream", "[pdb-records]") {
  // Invariant: operator>> reads at most 4*BufLen characters per line; a
  // longer line fails the stream and leaves the record untouched.
  std::istringstream in(std::string(400, 'X') + "\n");
  PDB record;
  in >> record;
  REQUIRE(in.fail());
  REQUIRE(record.type() == PDB::UNKNOWN);

  PDB parsed((std::string(120, 'X')).c_str());
  REQUIRE(parsed.type() == PDB::UNKNOWN);
  REQUIRE(std::string(parsed.chars()) == "UNKNOWN:  ??******??");
}

TEST_CASE("Far-overlong unknown lines stay bounded", "[pdb-records]") {
  // Invariant: an unknown record line far beyond the field width (the
  // closed stack-overflow todo's 300-character case) writes nothing past
  // the record buffer — the junk field renders 72 overflow characters
  // and the text is the bounded "UNKNOWN:  ??******??" form.
  PDB parsed((std::string(300, 'X')).c_str());
  REQUIRE(parsed.type() == PDB::UNKNOWN);
  const std::string text = parsed.chars();
  REQUIRE(text.size() == 20);
  REQUIRE(text == "UNKNOWN:  ??******??");

  // The same line through operator>> (an upload path) is bounded too.
  std::istringstream in((std::string(300, 'X')) + "\n");
  PDB streamed;
  in >> streamed;
  REQUIRE(in);
  REQUIRE(streamed.type() == PDB::UNKNOWN);
  REQUIRE(std::string(streamed.chars()) == "UNKNOWN:  ??******??");
}

TEST_CASE("Fixed-width sscanf semantics", "[pdb-records]") {
  // Invariant: the custom sscanf reads exact-width fields, yields zero
  // defaults past end-of-line, and reports -1 only on hard mismatches.
  int i = -1, j = -1;
  double d = -1.0;
  char s[16];
  char c[16];

  CHECK(PDB::sscanf("  12  ", "%5d", &i) == 1);
  CHECK(i == 12);

  CHECK(PDB::sscanf("", "%5d", &i) == 0);
  CHECK(i == 0);

  CHECK(PDB::sscanf("   1X", "%5d", &i) == -1);

  CHECK(PDB::sscanf("  1.50", "%6f", &d) == 1);
  CHECK(d == 1.5);

  CHECK(PDB::sscanf("", "%6f", &d) == 0);
  CHECK(d == 0.0);

  CHECK(PDB::sscanf(" AB  ", "%5s", s) == 1);
  CHECK(std::string(s) == " AB");

  CHECK(PDB::sscanf("", "%3c", c) == 0);
  CHECK(std::string(c, 3) == "   ");

  CHECK(PDB::sscanf("ABCDEFGH", "%3 %5s", s) == 2);
  CHECK(std::string(s) == "DEFGH");

  CHECK(PDB::sscanf("12345678", "%4d%4d", &i, &j) == 2);
  CHECK(i == 1234);
  CHECK(j == 5678);

  CHECK(PDB::sscanf("AXXX", "A%3s", s) == 1);
  CHECK(std::string(s) == "XXX");

  CHECK(PDB::sscanf("BXXX", "A%3s", s) == -1);

  // Literals past end-of-line are tolerated; missing fields count as no
  // match rather than a mismatch.
  CHECK(PDB::sscanf("AB", "AB%3s", s) == 0);
  CHECK(std::string(s) == "");
}

TEST_CASE("CRLF line endings", "[pdb-records]") {
  // Invariant (characterization): a trailing \r is ordinary content to the
  // fixed-width reader. String and numeric fields strip it (isspace), so
  // records whose last field is a string round-trip to the LF-only form;
  // a \r landing on a format literal is a hard mismatch, so short records
  // fall back to UNKNOWN; a \r landing in a %c field is captured verbatim
  // into the record and only disappears from the rendered text through the
  // trailing-blank trim.
  PDB author("AUTHOR   1VELDE, A.V.D.\r\n");
  REQUIRE(author.type() == PDB::AUTHOR);
  REQUIRE(std::string(author.author.authors) == "VELDE, A.V.D.");
  REQUIRE(std::string(author.chars()) == "AUTHOR   1VELDE, A.V.D.");

  // Full-width atom: the \r falls in the charge field and is stripped, so
  // the record parses and renders exactly like its LF-only form.
  PDB atom("ATOM    352  N   MET b   1      20.321  -9.325  38.937  1.00"
           " 15.91           N\r\n");
  REQUIRE(atom.type() == PDB::ATOM);
  REQUIRE(atom.atom.xyz[0] == 20.321);
  REQUIRE(atom.atom.occupancy == 1.0);
  REQUIRE(atom.atom.tempFactor == 15.91);
  REQUIRE(std::string(atom.chars()) ==
          "ATOM    352  N   MET b   1      20.321  -9.325  38.937  1.00"
          " 15.91           N");

  // Short atom: the \r sits where the format's trailing literal spaces
  // expect content, so the standard and ZDATOM layouts hard-fail and the
  // line degrades to an UNKNOWN junk record.
  PDB shortAtom("ATOM      1  CA  ALA A   1       1.000   2.000   3.000"
                "  1.00 10.00\r\n");
  REQUIRE(shortAtom.type() == PDB::UNKNOWN);
  REQUIRE(std::string(shortAtom.chars()) == "UNKNOWN:  ??    AT??");

  // The \r is captured verbatim by the insert-code field (0x0D) and is
  // hidden on output only by the trailing-blank trim.
  PDB ter("TER       1      ALA A  10\r\n");
  REQUIRE(ter.type() == PDB::TER);
  REQUIRE(ter.ter.residue.insertCode == '\r');
  REQUIRE(std::string(ter.chars()) == "TER       1      ALA A  10");

  // Unknown junk: the \r\n count toward the 72-character junk field.
  PDB junk("SOME UNKNOWN JUNK TEXT\r\n");
  REQUIRE(junk.type() == PDB::UNKNOWN);
  REQUIRE(std::string(junk.chars()) == "UNKNOWN:  ??      ??");

  // operator>> strips one trailing \r (Windows uploads read identically).
  std::istringstream in("AUTHOR   1VELDE, A.V.D.\r\n");
  PDB streamed;
  in >> streamed;
  REQUIRE(in);
  REQUIRE(streamed.type() == PDB::AUTHOR);
  REQUIRE(std::string(streamed.chars()) == "AUTHOR   1VELDE, A.V.D.");

  // The strip changes what a stream-fed short record sees: the \r no
  // longer lands on a format literal, so a short ATOM stays an ATOM and a
  // TER no longer captures a carriage return in its insert-code field.
  std::istringstream shortIn("ATOM      1  CA  ALA A   1       1.000   2.000"
                             "   3.000  1.00 10.00\r\n");
  PDB streamedAtom;
  shortIn >> streamedAtom;
  REQUIRE(shortIn);
  REQUIRE(streamedAtom.type() == PDB::ATOM);
  REQUIRE(streamedAtom.atom.serialNum == 1);
  REQUIRE(streamedAtom.atom.occupancy == 1.0);
  REQUIRE(streamedAtom.atom.tempFactor == 10.0);

  std::istringstream terIn("TER       1      ALA A  10\r\n");
  PDB streamedTer;
  terIn >> streamedTer;
  REQUIRE(terIn);
  REQUIRE(streamedTer.type() == PDB::TER);
  REQUIRE(streamedTer.ter.residue.insertCode == ' ');
}

TEST_CASE("Real mark_sur charge written as -.15", "[pdb-records]") {
  // Invariant (characterization): the zero-less charge form found in real
  // mark_sur data parses to the same value as "-0.15", and the writer
  // normalizes it to "-0.15" across the full 82-column ZDATOM layout, so
  // the value survives a write and re-read even though the literal form
  // does not survive byte-for-byte.
  const std::string line = "ATOM      1  N   ASN A   1      15.639 -11.284"
                           "  -9.203  1     0 1.63          -.15";
  PDB record(line.c_str());
  REQUIRE(record.type() == PDB::ATOM);
  REQUIRE(record.atom.iszdatom);
  REQUIRE(record.atom.chg == -0.15);
  REQUIRE(record.atom.type == 1);
  REQUIRE(record.atom.surface == 0);
  REQUIRE(record.atom.rad == 1.63);
  const std::string text = record.chars();
  REQUIRE(text.size() == 82);
  REQUIRE(text == "ATOM      1  N   ASN A   1      15.639 -11.284  -9.203"
                  "  1     0 1.63         -0.15");
  PDB reread(text.c_str());
  REQUIRE(reread.type() == PDB::ATOM);
  REQUIRE(reread.atom.iszdatom);
  REQUIRE(reread.atom.chg == -0.15);
}

TEST_CASE("Continued String records", "[pdb-records]") {
  // Invariant (characterization): continuation numbers render right-
  // justified with zero as blanks (%D), each String record's continuation
  // column has its own width (2 for AUTHOR/TITLE, 3 for COMPND), leading
  // blanks of the text survive because only trailing blanks are stripped,
  // and continued lines round-trip byte-for-byte.
  PDB title("TITLE    2 SOME TEXT");
  REQUIRE(title.type() == PDB::TITLE);
  REQUIRE(title.title.continuation == 2);
  REQUIRE(std::string(title.title.text) == " SOME TEXT");
  REQUIRE(std::string(title.chars()) == "TITLE    2 SOME TEXT");

  PDB author("AUTHOR   2 SMITH, J.");
  REQUIRE(author.type() == PDB::AUTHOR);
  REQUIRE(author.author.continuation == 2);
  REQUIRE(std::string(author.chars()) == "AUTHOR   2 SMITH, J.");

  // COMPND's continuation column is three wide.
  PDB compnd("COMPND  10 TEN");
  REQUIRE(compnd.type() == PDB::COMPND);
  REQUIRE(compnd.compnd.continuation == 10);
  REQUIRE(std::string(compnd.chars()) == "COMPND  10 TEN");

  // Continuation zero renders as blanks and re-parses as zero.
  PDB blankCont("AUTHOR    0 SMITH, J.");
  REQUIRE(blankCont.type() == PDB::AUTHOR);
  REQUIRE(blankCont.author.continuation == 0);
  REQUIRE(std::string(blankCont.chars()) == "AUTHOR    0 SMITH, J.");
}

TEST_CASE("Fixed-width sprintf semantics", "[pdb-records]") {
  // Invariant: the custom sprintf right/left-justifies within the requested
  // width, renders %D zero as blanks, and overflows to asterisks.
  char buf[PDB::BufLen];

  CHECK(PDB::sprintf(buf, "%5d", 42) == 5);
  CHECK(std::string(buf) == "   42");

  CHECK(PDB::sprintf(buf, "%3D", 0) == 3);
  CHECK(std::string(buf) == "   ");

  CHECK(PDB::sprintf(buf, "%05d", 42) == 5);
  CHECK(std::string(buf) == "00042");

  CHECK(PDB::sprintf(buf, "%-5d", 42) == 5);
  CHECK(std::string(buf) == "42   ");

  CHECK(PDB::sprintf(buf, "%3d", 12345) == 3);
  CHECK(std::string(buf) == "***");

  CHECK(PDB::sprintf(buf, "%3d", -12345) == 3);
  CHECK(std::string(buf) == "***");

  CHECK(PDB::sprintf(buf, "%8.3f", 1.2344) == 8);
  CHECK(std::string(buf) == "   1.234");

  CHECK(PDB::sprintf(buf, "%8.3f", -1.5) == 8);
  CHECK(std::string(buf) == "  -1.500");

  CHECK(PDB::sprintf(buf, "%4.2f", 123.45) == 4);
  CHECK(std::string(buf) == "****");

  CHECK(PDB::sprintf(buf, "%-6.6s", "AB") == 6);
  CHECK(std::string(buf) == "AB    ");

  CHECK(PDB::sprintf(buf, "%2s", "ABC") == 2);
  CHECK(std::string(buf) == "**");

  CHECK(PDB::sprintf(buf, "%3c", 'X') == 3);
  CHECK(std::string(buf) == "  X");
}
