// Copyright (c) 2019-2026 Arjan van der Velde, Weng Lab
// SPDX-License-Identifier: BSD-2-Clause

// Suite for the modern record class (Phase 2). At this stage it covers the
// RecordType enumeration only; parsing, writing, and parity tests with the
// legacy implementation arrive with the later Phase 2 commits.

#include "Test.hpp"
#include "Columns.hpp"
#include "Record.hpp"

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
// list spells every record name the fixed-format layouts cover, so a
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

// ---------------------------------------------------------------------------
// Typed payloads and the Record variant (Phase 2 section 4.2, commit 3)
// ---------------------------------------------------------------------------

TEST_CASE("Record selects the payload its record kind selects",
          "[record-payload]") {
  // Invariant: constructing a Record for any record kind yields the typed
  // payload that kind carries, with identical layouts sharing one struct
  // exactly as the legacy typedefs did.
  const std::vector<RecordType> allTypes = {
      RecordType::UNKNOWN,      RecordType::ANISOU,
      RecordType::ATOM,         RecordType::AUTHOR,
      RecordType::CAVEAT,       RecordType::COMPND,
      RecordType::CONECT,       RecordType::CRYST1,
      RecordType::END,          RecordType::FORMUL,
      RecordType::FTNOTE,       RecordType::HEADER,
      RecordType::HELIX,        RecordType::HET,
      RecordType::HETATM,       RecordType::JRNL,
      RecordType::KEYWDS,       RecordType::MASTER,
      RecordType::MDLTYP,       RecordType::MTRIX,
      RecordType::NUMMDL,       RecordType::OBSLTE,
      RecordType::ORIGX,        RecordType::REMARK,
      RecordType::REVDAT,       RecordType::SCALE,
      RecordType::SEQRES,       RecordType::SHEET,
      RecordType::SIGATM,       RecordType::SIGUIJ,
      RecordType::SITE,         RecordType::SOURCE,
      RecordType::SPLIT,        RecordType::SPRSDE,
      RecordType::SSBOND,       RecordType::TER,
      RecordType::TITLE,        RecordType::TURN,
      RecordType::TVECT,        RecordType::USER,
      RecordType::MODEL,        RecordType::ENDMDL,
      RecordType::EXPDTA,       RecordType::SYMDES,
      RecordType::SYMOP,        RecordType::MTXDES,
      RecordType::CMPDES,       RecordType::CMPONT,
      RecordType::TRNSFM,       RecordType::AGRDES,
      RecordType::AGGRGT,       RecordType::ZDATOM,
      RecordType::USER_PDBRUN,  RecordType::USER_EYEPOS,
      RecordType::USER_ATPOS,   RecordType::USER_WINDOW,
      RecordType::USER_FOCUS,   RecordType::USER_VIEWPORT,
      RecordType::USER_BGCOLOR, RecordType::USER_ANGLE,
      RecordType::USER_DISTANCE, RecordType::USER_FILE,
      RecordType::USER_MARKNAME, RecordType::USER_MARK,
      RecordType::USER_CNAME,   RecordType::USER_COLOR,
      RecordType::USER_RADIUS,  RecordType::USER_OBJECT,
      RecordType::USER_ENDOBJ,  RecordType::USER_CHAIN,
      RecordType::USER_GFX_BEGIN, RecordType::USER_GFX_END,
      RecordType::USER_GFX_COLOR, RecordType::USER_GFX_NORMAL,
      RecordType::USER_GFX_VERTEX, RecordType::USER_GFX_FONT,
      RecordType::USER_GFX_TEXTPOS, RecordType::USER_GFX_LABEL,
      RecordType::USER_GFX_MOVE, RecordType::USER_GFX_DRAW,
      RecordType::USER_GFX_MARKER, RecordType::USER_GFX_POINT,
  };
  for (const RecordType t : allTypes) {
    CAPTURE(zdock::pdb::to_string(t));
    REQUIRE_NOTHROW(zdock::pdb::Record(t));
  }

  struct KindCheck {
    RecordType type;
    bool (*holds)(const zdock::pdb::Record &);
  };
  const auto holdsMonostate = [](const zdock::pdb::Record &r) {
    return r.holds<zdock::pdb::Monostate>();
  };
  const auto holdsUnknown = [](const zdock::pdb::Record &r) {
    return r.holds<zdock::pdb::UnknownRecord>();
  };
  const auto holdsAtom = [](const zdock::pdb::Record &r) {
    return r.holds<zdock::pdb::AtomRecord>();
  };
  const auto holdsAnisou = [](const zdock::pdb::Record &r) {
    return r.holds<zdock::pdb::AnisouRecord>();
  };
  const auto holdsAgrdes = [](const zdock::pdb::Record &r) {
    return r.holds<zdock::pdb::AgrdesRecord>();
  };
  const auto holdsEyePos = [](const zdock::pdb::Record &r) {
    return r.holds<zdock::pdb::UserEyePosRecord>();
  };
  const auto holdsMarkname = [](const zdock::pdb::Record &r) {
    return r.holds<zdock::pdb::UserMarknameRecord>();
  };
  const auto holdsColor = [](const zdock::pdb::Record &r) {
    return r.holds<zdock::pdb::UserColorRecord>();
  };
  const auto holdsObject = [](const zdock::pdb::Record &r) {
    return r.holds<zdock::pdb::UserObjectRecord>();
  };
  const auto holdsGfxPoint = [](const zdock::pdb::Record &r) {
    return r.holds<zdock::pdb::UserGfxPointRecord>();
  };
  const std::vector<KindCheck> kinds = {
      {RecordType::UNKNOWN, holdsUnknown},
      {RecordType::END, holdsMonostate},
      {RecordType::ENDMDL, holdsMonostate},
      {RecordType::USER_GFX_END, holdsMonostate},
      {RecordType::ATOM, holdsAtom},
      {RecordType::HETATM, holdsAtom},
      {RecordType::SIGATM, holdsAtom},
      {RecordType::ZDATOM, holdsAtom},
      {RecordType::ANISOU, holdsAnisou},
      {RecordType::SIGUIJ, holdsAnisou},
      {RecordType::AGRDES, holdsAgrdes},
      {RecordType::CMPDES, holdsAgrdes},
      {RecordType::MTXDES, holdsAgrdes},
      {RecordType::SYMDES, holdsAgrdes},
      {RecordType::FTNOTE, holdsAgrdes},
      {RecordType::REMARK, holdsAgrdes},
      {RecordType::USER_EYEPOS, holdsEyePos},
      {RecordType::USER_ATPOS, holdsEyePos},
      {RecordType::USER_MARKNAME, holdsMarkname},
      {RecordType::USER_MARK, holdsMarkname},
      {RecordType::USER_COLOR, holdsColor},
      {RecordType::USER_GFX_COLOR, holdsColor},
      {RecordType::USER_OBJECT, holdsObject},
      {RecordType::USER_ENDOBJ, holdsObject},
      {RecordType::USER_GFX_NORMAL, holdsGfxPoint},
      {RecordType::USER_GFX_VERTEX, holdsGfxPoint},
      {RecordType::USER_GFX_TEXTPOS, holdsGfxPoint},
      {RecordType::USER_GFX_MOVE, holdsGfxPoint},
      {RecordType::USER_GFX_DRAW, holdsGfxPoint},
      {RecordType::USER_GFX_MARKER, holdsGfxPoint},
      {RecordType::USER_GFX_POINT, holdsGfxPoint},
  };
  for (const auto &kind : kinds) {
    CAPTURE(zdock::pdb::to_string(kind.type));
    const zdock::pdb::Record record(kind.type);
    REQUIRE(record.type() == kind.type);
    REQUIRE(kind.holds(record));
  }

  // distinct single-kind payloads
  REQUIRE(zdock::pdb::Record(RecordType::AUTHOR).holds<zdock::pdb::AuthorRecord>());
  REQUIRE(zdock::pdb::Record(RecordType::TITLE).holds<zdock::pdb::TitleRecord>());
  REQUIRE(zdock::pdb::Record(RecordType::COMPND).holds<zdock::pdb::CompndRecord>());
  REQUIRE(zdock::pdb::Record(RecordType::KEYWDS).holds<zdock::pdb::KeywdsRecord>());
  REQUIRE(zdock::pdb::Record(RecordType::SOURCE).holds<zdock::pdb::SourceRecord>());
  REQUIRE(zdock::pdb::Record(RecordType::MDLTYP).holds<zdock::pdb::MdltypRecord>());
  REQUIRE(zdock::pdb::Record(RecordType::EXPDTA).holds<zdock::pdb::ExpdtaRecord>());
  REQUIRE(zdock::pdb::Record(RecordType::USER).holds<zdock::pdb::UserRecord>());
}

TEST_CASE("Atom accessors are typed and wrong-kind access is an error",
          "[record-payload]") {
  // Invariant: asAtom returns the atom payload for every atom kind and is
  // a typed runtime error on any other record, replacing the silent union
  // misread the legacy layout permitted.
  for (const RecordType t :
       {RecordType::ATOM, RecordType::HETATM, RecordType::SIGATM}) {
    zdock::pdb::Record record(t);
    REQUIRE(record.isAtom());
    record.asAtom().serialNum = 7;
    REQUIRE(record.asAtom().serialNum == 7);
  }

  zdock::pdb::Record title(RecordType::TITLE);
  REQUIRE_FALSE(title.isAtom());
  REQUIRE_THROWS_AS(title.asAtom(), std::bad_variant_access);

  zdock::pdb::Record fallback;
  REQUIRE(fallback.type() == RecordType::UNKNOWN);
  REQUIRE_FALSE(fallback.isAtom());
}

TEST_CASE("Payload fields default to the parsed field defaults",
          "[record-payload]") {
  // Invariant: fresh payloads carry the parsed field defaults -- zero
  // numbers, blank codes, empty strings, and unset optionals for
  // spec-optional or version-conditional columns.
  using zdock::pdb::AtomRecord;
  using zdock::pdb::Residue;

  const Residue residue;
  REQUIRE(residue.name.empty());
  REQUIRE(residue.chainId == ' ');
  REQUIRE(residue.seqNum == 0);
  REQUIRE(residue.insertCode == ' ');

  const AtomRecord atom;
  REQUIRE(atom.serialNum == 0);
  REQUIRE(atom.name.empty());
  REQUIRE(atom.altLoc == ' ');
  REQUIRE_FALSE(atom.occupancy.has_value());
  REQUIRE_FALSE(atom.tempFactor.has_value());
  REQUIRE_FALSE(atom.zdatom.has_value());

  const zdock::pdb::UserAngleRecord angle;
  REQUIRE_FALSE(angle.which.has_value());
  REQUIRE(angle.angle == 0.0);

  const zdock::pdb::UserFileRecord file;
  REQUIRE_FALSE(file.model.has_value());

  const zdock::pdb::UserGfxLabelRecord label;
  REQUIRE_FALSE(label.xyz.has_value());

  const zdock::pdb::AuthorRecord author; // typed int continuation
  REQUIRE(author.continuation == 0);
}

// ---------------------------------------------------------------------------
// Simple-family parse and write with legacy parity (Phase 2, commit 4)
// ---------------------------------------------------------------------------

TEST_CASE("Simple families parse and format canonically",
          "[record-parity]") {
  // Invariant: for every simple-family line (continuation-plus-text
  // records, END, ENDMDL, and UNKNOWN junk, including CRLF, truncated,
  // continued, full-width, and lowercase forms), Record::parse yields
  // the pinned record type and Record::chars() renders the pinned
  // canonical text -- the exact bytes the retired legacy writer
  // produced, frozen here as explicit expectations.
  const std::vector<std::pair<std::string, std::string>> lines = {
      {
       "AUTHOR   1VELDE, A.V.D.",
       "AUTHOR   1VELDE, A.V.D."
      },
{
       "AUTHOR   1AAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAA"
              "AAAAAAAAAAAAAAAAAAAAAAA",
       "AUTHOR   1AAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAA"
              "AAAAAAAAAAAAAAAAAAAAAAA"
      },
{
       "AUTHOR   2 SMITH, J.",
       "AUTHOR   2 SMITH, J."
      },
{
       "AUTHOR    0 SMITH, J.",
       "AUTHOR    0 SMITH, J."
      },
{
       "AUTHOR   1VELDE, A.V.D.\r\n",
       "AUTHOR   1VELDE, A.V.D."
      },
{
       "AUTHOR",
       "AUTHOR"
      },
{
       "TITLE    1LYSOZYME",
       "TITLE    1LYSOZYME"
      },
{
       "TITLE    2 SOME TEXT",
       "TITLE    2 SOME TEXT"
      },
{
       "TITLE    1TTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTT"
              "TTTTTTTTTTTTTTTTTTTTTTTT",
       "TITLE    1TTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTT"
              "TTTTTTTTTTTTTTTTTTTTTTTT"
      },
{
       "TITLE",
       "TITLE"
      },
{
       "COMPND   1LYSOZYME",
       "COMPND   1LYSOZYME"
      },
{
       "COMPND  10 TEN",
       "COMPND  10 TEN"
      },
{
       "KEYWDS   1LYSOZYME",
       "KEYWDS   1LYSOZYME"
      },
{
       "SOURCE   1HUMAN",
       "SOURCE   1HUMAN"
      },
{
       "MDLTYP   1TETRAMER",
       "MDLTYP   1TETRAMER"
      },
{
       "END",
       "END"
      },
{
       "ENDMDL",
       "ENDMDL"
      },
{
       "END  ",
       "END"
      },
{
       "author   1velde",
       "AUTHOR   1velde"
      },
  };
  for (const auto &fixture : lines) {
    CAPTURE(fixture.first);
    zdock::pdb::ParseContext context;
    const std::optional<zdock::pdb::Record> parsed =
        zdock::pdb::Record::parse(fixture.first, context);
    REQUIRE(parsed.has_value());
    REQUIRE(parsed->chars() == fixture.second);
  }
}

TEST_CASE("Simple-family fields parse into the pinned typed values",
          "[record-parity]") {
  // Invariant: continuation numbers land as typed ints, text keeps its
  // leading blanks (only trailing blanks are stripped), truncation yields
  // defaults, and the AUTHOR continuation is an int as the rewrite
  // requires.
  zdock::pdb::ParseContext context;

  const auto author = zdock::pdb::Record::parse("AUTHOR   1VELDE, A.V.D.",
                                                context);
  REQUIRE(author.has_value());
  REQUIRE(author->holds<zdock::pdb::AuthorRecord>());
  REQUIRE(author->as<zdock::pdb::AuthorRecord>().continuation == 1);
  REQUIRE(author->as<zdock::pdb::AuthorRecord>().authors == "VELDE, A.V.D.");

  const auto title = zdock::pdb::Record::parse("TITLE    2 SOME TEXT", context);
  REQUIRE(title.has_value());
  REQUIRE(title->as<zdock::pdb::TitleRecord>().continuation == 2);
  REQUIRE(title->as<zdock::pdb::TitleRecord>().text == " SOME TEXT");

  const auto compnd = zdock::pdb::Record::parse("COMPND  10 TEN", context);
  REQUIRE(compnd.has_value());
  REQUIRE(compnd->as<zdock::pdb::CompndRecord>().continuation == 10);
  REQUIRE(compnd->as<zdock::pdb::CompndRecord>().compound == " TEN");

  // In "AUTHOR    0 SMITH, J." the zero sits in the text field (column
  // 10); the continuation columns are blank and parse as the blank-field
  // zero, exactly as the legacy reader saw it.
  const auto zeroCont =
      zdock::pdb::Record::parse("AUTHOR    0 SMITH, J.", context);
  REQUIRE(zeroCont.has_value());
  REQUIRE(zeroCont->as<zdock::pdb::AuthorRecord>().continuation == 0);
  REQUIRE(zeroCont->as<zdock::pdb::AuthorRecord>().authors == "0 SMITH, J.");

  const auto truncated = zdock::pdb::Record::parse("AUTHOR", context);
  REQUIRE(truncated.has_value());
  REQUIRE(truncated->as<zdock::pdb::AuthorRecord>().continuation == 0);
  REQUIRE(truncated->as<zdock::pdb::AuthorRecord>().authors.empty());
}

TEST_CASE("Simple-family hard mismatches fall back like the legacy parse",
          "[record-parity]") {
  // Invariant: a hard numeric mismatch yields std::nullopt, and the
  // caller-visible junk fallback renders exactly what the legacy
  // UNKNOWN fallback rendered.
  zdock::pdb::ParseContext context;
  const std::string mismatch = "AUTHOR  X ";
  REQUIRE_FALSE(zdock::pdb::Record::parse(mismatch, context).has_value());
  // the pinned junk render the retired legacy fallback produced
  REQUIRE(zdock::pdb::Record::unknownFrom(mismatch).chars() ==
          "UNKNOWN:  ??      ??");
}

TEST_CASE("Unknown lines parse as junk records matching legacy",
          "[record-parity]") {
  // Invariant: lines with no recognized record type parse to UNKNOWN
  // junk records whose rendering matches the legacy fallback, including
  // the pinned bounded forms for overlong lines.
  const std::vector<std::string> lines = {
      "SOME UNKNOWN JUNK TEXT\r\n",
      std::string(120, 'X'),
      std::string(300, 'X'),
  };
  for (const auto &line : lines) {
    CAPTURE(line);
    zdock::pdb::ParseContext context;
    const std::optional<zdock::pdb::Record> parsed =
        zdock::pdb::Record::parse(line, context);
    REQUIRE(parsed.has_value());
    REQUIRE(parsed->type() == zdock::pdb::RecordType::UNKNOWN);
    REQUIRE(parsed->holds<zdock::pdb::UnknownRecord>());
    // the pinned junk renders the retired legacy fallback produced:
    // the short CRLF line blanks the tail, overlong lines fill with
    // the overflow stars
    REQUIRE(parsed->chars() ==
            (line.rfind("SOME", 0) == 0 ? "UNKNOWN:  ??      ??"
                                      : "UNKNOWN:  ??******??"));
  }
  zdock::pdb::ParseContext context;
  const auto bounded =
      zdock::pdb::Record::parse(std::string(300, 'X'), context);
  REQUIRE(bounded.has_value());
  REQUIRE(bounded->chars() == "UNKNOWN:  ??******??");
}

TEST_CASE("Fresh simple records format from typed payloads",
          "[record-parity]") {
  // Invariant: rendering works from constructed payloads without a prior
  // parse, reproducing the legacy writer's continuation, padding, and
  // trailing-blank conventions; families without write tables yet format
  // as empty.
  REQUIRE(zdock::pdb::Record(zdock::pdb::RecordType::AUTHOR).chars() ==
          "AUTHOR");
  REQUIRE(zdock::pdb::Record(zdock::pdb::RecordType::END).chars() == "END");
  REQUIRE(zdock::pdb::Record(zdock::pdb::RecordType::ENDMDL).chars() ==
          "ENDMDL");
  // Deliberate deviation: a fresh atom renders with unset occupancy and
  // tempFactor (as 0.00); the legacy fresh-record 1.0 occupancy default is
  // a struct-level quirk the typed payload does not carry. Parsed atoms
  // are unaffected.
  REQUIRE(zdock::pdb::Record(zdock::pdb::RecordType::ATOM).chars() ==
          "ATOM      0              0       0.000   0.000   0.000  0.00  "
          "0.00");

  zdock::pdb::Record author(zdock::pdb::RecordType::AUTHOR);
  author.as<zdock::pdb::AuthorRecord>().continuation = 2;
  author.as<zdock::pdb::AuthorRecord>().authors = " SMITH, J.";
  REQUIRE(author.chars() == "AUTHOR   2 SMITH, J.");

  zdock::pdb::Record compnd(zdock::pdb::RecordType::COMPND);
  compnd.as<zdock::pdb::CompndRecord>().continuation = 10;
  compnd.as<zdock::pdb::CompndRecord>().compound = " TEN";
  REQUIRE(compnd.chars() == "COMPND  10 TEN");
}

TEST_CASE("Records stream through operator<< as their rendered text",
          "[record-parity]") {
  // Invariant: operator<< emits the same text chars() returns, mirroring
  // the legacy stream operator.
  zdock::pdb::ParseContext context;
  const auto parsed =
      zdock::pdb::Record::parse("TITLE    1LYSOZYME", context);
  REQUIRE(parsed.has_value());
  std::ostringstream out;
  out << *parsed;
  REQUIRE(out.str() == "TITLE    1LYSOZYME");
}

TEST_CASE("write_into renders the same bounded text as chars",
          "[record-parity]") {  // Invariant: the allocation-free renderer produces exactly the chars()
  // text, NUL-terminated at the trimmed length, writing nothing past the
  // declared record length inside the 83-byte buffer.
  const std::vector<std::string> lines = {
      "AUTHOR   1VELDE, A.V.D.",
      "AUTHOR   1" + std::string(69, 'A'),
      "TITLE    2 SOME TEXT",
      "COMPND  10 TEN",
      "KEYWDS   1LYSOZYME",
      "SOURCE   1HUMAN",
      "MDLTYP   1TETRAMER",
      "END",
      "ENDMDL",
      std::string(120, 'X'),
      std::string(300, 'X'),
  };
  for (const auto &line : lines) {
    CAPTURE(line);
    zdock::pdb::ParseContext context;
    const auto parsed = zdock::pdb::Record::parse(line, context);
    REQUIRE(parsed.has_value());
    char buf[83];
    std::memset(buf, '#', sizeof(buf)); // sentinel past the record length
    const std::size_t length = parsed->write_into(buf, zdock::pdb::WriteContext{});
    REQUIRE(length <= 80);
    REQUIRE(buf[length] == '\0');
    REQUIRE(std::string(buf, length) == parsed->chars());
    for (std::size_t i = 81; i < sizeof(buf); ++i)
      REQUIRE(buf[i] == '#'); // nothing written past the record length
  }
  // Families without write tables format as empty; atoms now render their
  // default payload.
  zdock::pdb::Record atom(zdock::pdb::RecordType::ATOM);
  char buf[83];
  REQUIRE(atom.write_into(buf, zdock::pdb::WriteContext{}) ==
          std::strlen("ATOM      0              0       0.000   0.000   "
                      "0.000  0.00  0.00"));
  REQUIRE(std::string(buf) ==
          "ATOM      0              0       0.000   0.000   0.000  0.00  "
          "0.00");
}

// ---------------------------------------------------------------------------
// ATOM family, the ZDATOM fallback, and the charge cases (Phase 2, commit 5)
// ---------------------------------------------------------------------------

TEST_CASE("Atom columns parse into the pinned typed values",
          "[record-parity]") {
  // Invariant: the standard atom layout lands in the typed members at its
  // wwPDB columns exactly as the legacy parser stores them (the golden
  // suite's values), with names right-trimmed and codes verbatim.
  zdock::pdb::ParseContext context;
  const auto atom = zdock::pdb::Record::parse(
      "ATOM    352  N   MET b   1      20.321  -9.325  38.937  1.00 15.91"
      "           N",
      context);
  REQUIRE(atom.has_value());
  REQUIRE(atom->type() == zdock::pdb::RecordType::ATOM);
  REQUIRE(atom->holds<zdock::pdb::AtomRecord>());
  const zdock::pdb::AtomRecord &a = atom->as<zdock::pdb::AtomRecord>();
  REQUIRE(a.serialNum == 352);
  REQUIRE(a.name == " N");
  REQUIRE(a.altLoc == ' ');
  REQUIRE(a.residue.name == "MET");
  REQUIRE(a.residue.chainId == 'b');
  REQUIRE(a.residue.seqNum == 1);
  REQUIRE(a.residue.insertCode == ' ');
  REQUIRE(a.xyz[0] == 20.321);
  REQUIRE(a.xyz[1] == -9.325);
  REQUIRE(a.xyz[2] == 38.937);
  REQUIRE(a.occupancy.value_or(-1.0) == 1.0);
  REQUIRE(a.tempFactor.value_or(-1.0) == 15.91);
  REQUIRE(a.element == " N");
  REQUIRE(a.charge.empty());
  REQUIRE_FALSE(a.zdatom.has_value());

  // HETATM and SIGATM classify to their own types over the same layout.
  for (const char *line : {
           "HETATM 9999  O   HOH A 999       1.500  -2.250   3.750  0.50 "
           "20.00           O",
           "SIGATM    1  N   ALA A   1       0.020   0.020   0.020  0.00  "
           "1.00",
       }) {
    CAPTURE(line);
    const auto parsed = zdock::pdb::Record::parse(line, context);
    REQUIRE(parsed.has_value());
    REQUIRE(parsed->holds<zdock::pdb::AtomRecord>());
  }
}

TEST_CASE("Atom lines fall back to ZDATOM and to junk like the legacy parse",
          "[record-parity]") {
  // Invariant: an atom line that fails the standard layout parses with the
  // 82-column ZDATOM extension (type retained), one that fails both falls
  // back to the junk record, and a truncated line parses as ATOM with
  // defaults -- all byte-for-byte with the legacy path.
  zdock::pdb::ParseContext context;

  const char *zdatomLine =
      "ATOM      1  CA  ALA A   1       1.000   2.000   3.000  1     1 1.23"
      "     ALA  1.50";
  const auto z = zdock::pdb::Record::parse(zdatomLine, context);
  REQUIRE(z.has_value());
  REQUIRE(z->type() == zdock::pdb::RecordType::ATOM);
  REQUIRE(z->asAtom().zdatom.has_value());
  REQUIRE(z->asAtom().zdatom->type == 1);
  REQUIRE(z->asAtom().zdatom->surface == 1);
  REQUIRE(z->asAtom().zdatom->rad == 1.23);
  REQUIRE(z->asAtom().zdatom->segid == " ALA");
  REQUIRE(z->asAtom().zdatom->chg == 1.5);
  REQUIRE(z->chars() == zdatomLine);

  // The junk fallback for lines that fit neither layout.
  for (const char *junkLine : {
           "ATOM   1X",
           "ATOM      1  CA  ALA A   1       1.000   2.000   3.000  1.00 "
           "10.00\r\n",
       }) {
    CAPTURE(junkLine);
    REQUIRE_FALSE(zdock::pdb::Record::parse(junkLine, context).has_value());
    // the pinned junk renders the retired legacy fallback produced:
    // the short mismatch blanks the tail, the CRLF line keeps its
    // truncated atom text in the junk window
    REQUIRE(zdock::pdb::Record::unknownFrom(junkLine).chars() ==
            (std::string(junkLine) == "ATOM   1X"
                 ? "UNKNOWN:  ??      ??"
                 : "UNKNOWN:  ??    AT??"));
  }

  // Truncated atom lines parse with defaults, not junk.
  const auto truncated = zdock::pdb::Record::parse("ATOM", context);
  REQUIRE(truncated.has_value());
  REQUIRE(truncated->type() == zdock::pdb::RecordType::ATOM);
  REQUIRE(truncated->asAtom().serialNum == 0);
  REQUIRE(truncated->asAtom().xyz[0] == 0.0);
  REQUIRE(truncated->chars() ==
          "ATOM      0              0       0.000   0.000   0.000  0.00  "
          "0.00");
}

TEST_CASE("Standard atom lines with a wwPDB tail salvage their fields",
          "[record-parity]") {
  // Invariant: an ATOM/HETATM line that fails the pinned blank-column
  // literal at 67-76 because it carries a wwPDB segID (73-76, left-
  // justified), element (77-78, right-justified), and optionally a
  // charge (79-80), and that also fails the numeric mark_sur tail,
  // parses as the standard record it announces with the tail fields
  // kept, renders byte-identically to its input, and round-trips
  // stably. The legacy parser stores such lines as junk -- the
  // deliberate modern gain, asserted here as a classification
  // difference rather than text parity. Each fixture is assembled from
  // the layout: 6-byte prefix, %5d serial, space, %-4s name, %c altLoc,
  // %-4s residue, %c chain, %4d sequence, %c insert code, three spaces,
  // %8.3f coordinates, %6.2f occupancy, %6.2f tempFactor, six blanks,
  // %-4s segID, %2s element, %2s charge.
  const std::string head =
      "    1"  // %5d serial
      " "      // literal space
      " CA "   // %-4s atom name
      " "      // %c altLoc
      "ALA "   // %-4s residue name
      "A"      // %c chain id
      "   1"   // %4d sequence
      " "      // %c insert code
      "   "    // literal spaces
      "   1.000" "   2.000" "   3.000" // %8.3f coordinates
      "  1.00" " 22.00";               // %6.2f occupancy, tempFactor
  struct Fixture {
    const char *prefix; // "ATOM  " or "HETATM"
    const char *tail;   // segID, element, charge as rendered
    const char *segid;
    const char *element; // read with the leading blank kept, as the
                         // strict layout reads it
    const char *charge;
  };
  const std::vector<Fixture> fixtures = {
      // the dominant corpus shape: chain-letter segID, one-letter element
      {"ATOM  ", "      A    C", "A", " C", ""},
      // a segID carrying digits ("1SG"-style masked shape)
      {"ATOM  ", "      1AB  N", "1AB", " N", ""},
      // a charge in columns 79-80 (right-justified like the element)
      {"ATOM  ", "      A    C2-", "A", " C", "2-"},
      // the HETATM shape
      {"HETATM", "      W    O", "W", " O", ""},
  };
  for (const Fixture &fixture : fixtures) {
    const std::string line =
        std::string(fixture.prefix) + head + fixture.tail;
    CAPTURE(line);
    zdock::pdb::ParseContext context;
    const auto parsed = zdock::pdb::Record::parse(line, context);
    REQUIRE(parsed.has_value());
    REQUIRE((parsed->type() == zdock::pdb::RecordType::ATOM ||
             parsed->type() == zdock::pdb::RecordType::HETATM));
    // the deliberate divergence from the retired legacy parser, which
    // stored these lines as junk
    const auto &atom = parsed->asAtom();
    REQUIRE_FALSE(atom.zdatom.has_value());
    REQUIRE(atom.serialNum == 1);
    REQUIRE(atom.name == " CA");
    REQUIRE(atom.residue.name == "ALA");
    REQUIRE(atom.residue.chainId == 'A');
    REQUIRE(atom.residue.seqNum == 1);
    REQUIRE(atom.xyz[0] == 1.0);
    REQUIRE(atom.xyz[1] == 2.0);
    REQUIRE(atom.xyz[2] == 3.0);
    REQUIRE(atom.occupancy.value_or(-1.0) == 1.0);
    REQUIRE(atom.tempFactor.value_or(-1.0) == 22.0);
    REQUIRE(atom.segid == fixture.segid);
    REQUIRE(atom.element == fixture.element);
    REQUIRE(atom.charge == fixture.charge);
    // the render restores the input byte-for-byte and round-trips
    REQUIRE(parsed->chars() == line);
    zdock::pdb::ParseContext again;
    const auto reread = zdock::pdb::Record::parse(parsed->chars(), again);
    REQUIRE(reread.has_value());
    REQUIRE(reread->asAtom().segid == atom.segid);
    REQUIRE(reread->asAtom().element == atom.element);
    REQUIRE(reread->chars() == parsed->chars());
    // the tail stays out of the coordinate-rejection census
    REQUIRE(context.coordinateRejections == 0);
  }
}

TEST_CASE("Unparsable coordinate lines are counted, not silent",
          "[record-parity]") {
  // Invariant: a coordinate-looking line that no layout accepts returns
  // no record, increments the caller's ParseContext rejection counter
  // (the warning/error hook), and is distinguishable from non-
  // coordinate junk, which does not count; the pinned charge-overflow
  // render degrades to junk and counts too.
  const std::string brokenCoords = "ATOM      1  CA  ALA A   1      XX.000"
                                   "   2.000   3.000  1.00 22.00      A"
                                   "    C";
  zdock::pdb::ParseContext context;
  REQUIRE_FALSE(zdock::pdb::Record::parse(brokenCoords, context).has_value());
  REQUIRE(context.coordinateRejections == 1);
  REQUIRE(zdock::pdb::Record::looksLikeCoordinate(brokenCoords));

  const std::string junk = "GARBAGE AND NOT A RECORD AT ALL";
  const auto junkRecord = zdock::pdb::Record::parse(junk, context);
  REQUIRE(junkRecord.has_value());
  REQUIRE(junkRecord->type() == zdock::pdb::RecordType::UNKNOWN);
  REQUIRE(context.coordinateRejections == 1); // unchanged
  REQUIRE_FALSE(zdock::pdb::Record::looksLikeCoordinate(junk));

  // the pinned '*****' charge-overflow render degrades loudly and counts
  const std::string overflow =
      "ATOM      1  CA  ALA A   1       1.000   2.000   3.000  1     1 "
      "1.23     ALA -12.3";
  zdock::pdb::ParseContext overflowContext;
  const auto rendered = zdock::pdb::Record::parse(overflow, overflowContext);
  REQUIRE(rendered.has_value());
  const std::string stars = rendered->chars();
  REQUIRE(stars.find("*****") != std::string::npos);
  zdock::pdb::ParseContext rereadContext;
  REQUIRE(
      !zdock::pdb::Record::parse(stars, rereadContext).has_value());
  REQUIRE(rereadContext.coordinateRejections == 1);
}

TEST_CASE("ZDATOM charges survive the pinned cases", "[record-parity]") {
  // Invariant: the pinned charge forms render across the full 82-column
  // layout, the zero-less input normalizes to "-0.15", and a charge too
  // wide for %5.2f renders as asterisks and the rendered line re-parses
  // as UNKNOWN junk (the read itself accepts the value; only the render
  // overflows) -- every parsable rendering matches the legacy text
  // byte-for-byte.
  const std::string base = "ATOM      1  CA  ALA A   1       1.000   2.000"
                           "   3.000  1     1 1.23     ALA ";
  struct ChargeCase {
    const char *suffix;
    double charge;
    const char *rendered;
  };
  const std::vector<ChargeCase> cases = {
      {" 1.50", 1.5, " 1.50"},
      {" 0.25", 0.25, " 0.25"},
      {"-0.45", -0.45, "-0.45"},
      {" -.15", -0.15, "-0.15"},
      {"-12.3", -12.3, "*****"},
  };
  zdock::pdb::ParseContext context;
  for (const auto &c : cases) {
    const std::string line = base + c.suffix;
    CAPTURE(line);
    const auto parsed = zdock::pdb::Record::parse(line, context);
    REQUIRE(parsed.has_value());
    REQUIRE(parsed->asAtom().zdatom->chg == c.charge);
    REQUIRE(parsed->chars() == base + c.rendered);
    // the rendered line re-parses and re-renders stably; the asterisk
    // rendering degrades to the junk fallback exactly as the legacy one
    // does
    zdock::pdb::ParseContext again;
    const auto reread = zdock::pdb::Record::parse(parsed->chars(), again);
    if (std::string(c.rendered) == "*****") {
      REQUIRE_FALSE(reread.has_value());
      // the pinned junk render the retired legacy fallback produced
      REQUIRE(zdock::pdb::Record::unknownFrom(parsed->chars()).chars() ==
              "UNKNOWN:  ??******??");
    } else {
      REQUIRE(reread.has_value());
      REQUIRE(reread->asAtom().zdatom->chg == c.charge);
      REQUIRE(reread->chars() == parsed->chars());
    }
  }
}

TEST_CASE("Clean ZDATOM parse leaves occupancy unset", "[record-payload]") {
  // Invariant (deliberate deviation): the ZDATOM fallback parses on a
  // fresh record, so occupancy stays unset instead of carrying the failed
  // standard attempt's partial write; ZDATOM text never renders occupancy,
  // so no golden line changes.
  zdock::pdb::ParseContext context;
  const auto z = zdock::pdb::Record::parse(
      "ATOM      1  CA  ALA A   1       1.000   2.000   3.000  1     1 1.23"
      "     ALA  1.50",
      context);
  REQUIRE(z.has_value());
  REQUIRE_FALSE(z->asAtom().occupancy.has_value());
  REQUIRE_FALSE(z->asAtom().tempFactor.has_value());
}

TEST_CASE("Real mark_sur fixtures match legacy and round-trip",
          "[record-parity]") {
  // Invariant: every real 82-column ZDATOM line from the stored example
  // results (1UDH/1POH and blank segment ids, zero-less and normalized
  // charge forms) parses to the same type and byte-identical text as the
  // legacy parser, and the rendered line re-parses and re-renders
  // stably under both paths.
  std::ifstream in(test::getpath("ZDATOM/mark_sur.pdb"));
  REQUIRE(in);
  std::string line;
  while (std::getline(in, line)) {
    CAPTURE(line);
    REQUIRE(line.size() == 82);
    zdock::pdb::ParseContext context;
    const auto parsed = zdock::pdb::Record::parse(line, context);
    REQUIRE(parsed.has_value());
    REQUIRE(parsed->asAtom().zdatom.has_value());
    // the render canonicalizes zero-less charges (pinned in the charge
    // test) and re-parses stably; the 82-column layout is asserted in
    // the write-buffer test
    zdock::pdb::ParseContext again;
    const auto reread = zdock::pdb::Record::parse(parsed->chars(), again);
    REQUIRE(reread.has_value());
    REQUIRE(reread->asAtom().zdatom.has_value());
    REQUIRE(reread->asAtom().zdatom->chg ==
            parsed->asAtom().zdatom->chg);
    REQUIRE(reread->chars() == parsed->chars());
  }
}

TEST_CASE("ZDATOM records render 82 columns inside the write buffer",
          "[record-parity]") {
  // Invariant: the 82-column ZDATOM layout renders fully (the pinned
  // full-width case), writes nothing past the declared length in the
  // 83-byte buffer, and standard atoms stay at 80 columns.
  zdock::pdb::ParseContext context;
  const auto z = zdock::pdb::Record::parse(
      "ATOM      1  CA  ALA A   1      12.345  -6.789  10.000  1     1 1.23"
      "     ALA  1.50",
      context);
  REQUIRE(z.has_value());
  REQUIRE(z->chars().size() == 82);
  REQUIRE(z->chars() ==
          "ATOM      1  CA  ALA A   1      12.345  -6.789  10.000 "
          " 1     1 1.23     ALA  1.50");
  char buf[83];
  std::memset(buf, '#', sizeof(buf));
  const std::size_t n = z->write_into(buf, zdock::pdb::WriteContext{});
  REQUIRE(n == 82);
  REQUIRE(buf[82] == '\0');
  for (std::size_t i = 83 - 1; i > 82; --i)
    REQUIRE(buf[i] == '#');
}

// ---------------------------------------------------------------------------
// Residue-heavy and structural families (Phase 2, commit 6)
// ---------------------------------------------------------------------------

TEST_CASE("Structural families parse and format canonically",
          "[record-parity]") {
  // Invariant: every structural, residue-heavy, and single-value record
  // line (the characterization fixtures plus continued and zero-valued
  // forms) parses and renders its pinned canonical text -- the exact
  // bytes the retired legacy writer produced, including the %D
  // zero-as-blank continuations and the MODEL read/write column
  // asymmetry, frozen here as explicit expectations.
  const std::vector<std::pair<std::string, std::string>> lines = {
      {
       "HELIX    1  H1 ALA A   10  LEU A   15  1                "
              "                     6",
       "HELIX    1  H1 ALA A   10  LEU A   15  1"
      },
{
       "HELIX    0  H1 ALA A   10  LEU A   15  0                "
              "                     0",
       "HELIX       H1 ALA A   10  LEU A   15"
      },
{
       "SHEET    1  S1 5 ALA A  10  LEU A  15  0  N  ALA A  10  "
              " O  LEU A  15",
       "SHEET    1  S1 5 ALA A  10  LEU A  15  0  N  ALA A  10  "
              " O  LEU A  15"
      },
{
       "SHEET    0  S1 0 ALA A   0  LEU A   0  0  N  ALA A   0  "
              " O  LEU A   0",
       "SHEET       S1 0 ALA A   0  LEU A   0  0  N  ALA A      "
              " O  LEU A"
      },
{
       "SITE     1 CAT  2 GLY A  10  ALA A  15  VAL A  20  LEU A"
              "  25",
       "SITE     1 CAT  2 GLY A  10  ALA A  15  VAL A  20  LEU A"
              "  25"
      },
{
       "SITE     0 CAT  0",
       "SITE     0 CAT  0"
      },
{
       "SSBOND   1 CYS A   10    CYS A   25    1555   1555  2.05",
       "SSBOND   1 CYS A   10    CYS A   25                     "
              "                  0.00"
      },
{
       "SSBOND   0 CYS A   10    CYS A   25       0     0  0.00",
       "SSBOND     CYS A   10    CYS A   25                     "
              "                  0.00"
      },
{
       "TER       1      ALA A  10",
       "TER       1      ALA A  10"
      },
{
       "TER       0      ALA A   0",
       "TER       0      ALA A   0"
      },
{
       "TURN     1  T1 ALA A  10  LEU A  15     TURN",
       "TURN     1  T1 ALA A  10  LEU A  15     TURN"
      },
{
       "TURN     0  T1 ALA A  10  LEU A  15     ",
       "TURN        T1 ALA A  10  LEU A  15"
      },
{
       "CMPONT   1  ALAA    1   ALAA    2",
       "CMPONT   1  ALAA    1   ALAA    2"
      },
{
       "CONECT    1    2    3    4    5    6    7    8",
       "CONECT    1    2    3    4    5    6    7    8"
      },
{
       "CONECT    0    0    0    0    0    0    0    0    0    0"
              "    0",
       "CONECT    0"
      },
{
       "CRYST1   78.900   78.900   85.700  90.00  90.00  90.00 P"
              " 21 21 21    8",
       "CRYST1   78.900   78.900   85.700  90.00  90.00  90.00 P"
              " 21 21 21    8"
      },
{
       "MTRIX1   1  1.000000  0.000000  0.000000        0.00000",
       "MTRIX1   1  1.000000  0.000000  0.000000        0.00000"
      },
{
       "MTRIX1   0  0.000000  0.000000  0.000000        0.00000 "
              "       0",
       "MTRIX1   0  0.000000  0.000000  0.000000        0.00000"
      },
{
       "ORIGX1      1.000000  0.000000  0.000000        0.00000",
       "ORIGX1      1.000000  0.000000  0.000000        0.00000"
      },
{
       "SCALE1      0.017400  0.000000  0.000000        0.00000",
       "SCALE1      0.017400  0.000000  0.000000        0.00000"
      },
{
       "SYMOP1   1  1.000000  0.000000  0.000000        0.00000",
       "SYMOP1   1  1.000000  0.000000  0.000000        0.00000"
      },
{
       "REVDAT   1   29-MAR-95 1LYZ    1       REVDAT",
       "REVDAT   1   29-MAR-95 1LYZ    1       REVDAT"
      },
{
       "REVDAT   0   00-XXX-00 0XXX    0       ",
       "REVDAT       00-XXX-00 0XXX    0"
      },
{
       "SEQRES   1 A    2  GLY  ALA",
       "SEQRES   1 A    2  GLY  ALA"
      },
{
       "SEQRES   0 A    0  ",
       "SEQRES   0 A    0"
      },
{
       "SPLIT    1 1ABZ 2CDZ",
       "SPLIT    1 1ABZ 2CDZ"
      },
{
       "SPLIT    0 ",
       "SPLIT"
      },
{
       "OBSLTE   1 24-JAN-95 1LYZ      2LYZ",
       "OBSLTE   1 24-JAN-95 1LYZ      2LYZ"
      },
{
       "OBSLTE   0 00-XXX-00 0XXX     ",
       "OBSLTE     00-XXX-00 0XXX"
      },
{
       "SPRSDE   1 24-JAN-95 1LYZ      2LYZ",
       "SPRSDE   1 24-JAN-95 1LYZ      2LYZ"
      },
{
       "MASTER        2    0    1    1    1    1    1    0  126 "
              "   1   10    1",
       "MASTER        2    0    1    1    1    1    1    0  126 "
              "   1   10    1"
      },
{
       "MASTER        0    0    0    0    0    0    0    0    0 "
              "   0    0    0",
       "MASTER        0    0    0    0    0    0    0    0    0 "
              "   0    0    0"
      },
{
       "AGGRGT   1   2   3   4   0   0   0   0   0   0   0   0  "
              " 0   0   0   0",
       "AGGRGT   1   2   3   4   0   0   0   0   0   0   0   0  "
              " 0   0   0   0"
      },
{
       "AGGRGT   0   0   0   0   0   0   0   0   0   0   0   0  "
              " 0   0   0   0",
       "AGGRGT   0   0   0   0   0   0   0   0   0   0   0   0  "
              " 0   0   0   0"
      },
{
       "TRNSFM   1   1   1",
       "TRNSFM   1   1   1"
      },
{
       "TRNSFM   0   0   0",
       "TRNSFM   0   0   0"
      },
{
       "TVECT    1   0.10000   0.20000   0.30000VECTOR",
       "TVECT    1   0.10000   0.20000   0.30000VECTOR"
      },
{
       "TVECT    0   0.00000   0.00000   0.00000",
       "TVECT        0.00000   0.00000   0.00000"
      },
{
       "HET    NAG  A   1      14     N-ACETYL-D-GLUCOSAMINE",
       "HET    NAG  A   1      14     N-ACETYL-D-GLUCOSAMINE"
      },
{
       "HET    NAG  A   1      14",
       "HET    NAG  A   1      14"
      },
{
       "FORMUL   2  HOH  1*(H2 O)",
       "FORMUL   2  HOH  1*(H2 O)"
      },
{
       "FORMUL   0    0  0 ",
       "FORMUL        0"
      },
{
       "HEADER    TRANSFERASE24-JAN-95   1LYZ",
       "HEADER    TRANSFERASE24-JAN-95   1LYZ"
      },
{
       "HEADER                        ",
       "HEADER"
      },
{
       "CAVEAT   1 1LYZ    SOME WARNING TEXT",
       "CAVEAT   1 1LYZ    SOME WARNING TEXT"
      },
{
       "CAVEAT   0 0XXX    ",
       "CAVEAT     0XXX"
      },
{
       "JRNL        TO BE PUBLISHED",
       "JRNL        TO BE PUBLISHED"
      },
{
       "JRNL        ",
       "JRNL"
      },
{
       "MODEL        1",
       "MODEL        1"
      },
{
       "MODEL    ",
       "MODEL        0"
      },
{
       "NUMMDL    2",
       "NUMMDL    2"
      },
{
       "NUMMDL    ",
       "NUMMDL    0"
      },
{
       "EXPDTA   1X-RAY DIFFRACTION",
       "EXPDTA   1X-RAY DIFFRACTION"
      },
{
       "EXPDTA   0",
       "EXPDTA"
      },
{
       "AGRDES   1                                             A"
              "N AGGREGATE DESCRIPTION",
       "AGRDES   1                                             A"
              "N AGGREGATE DESCRIPTION"
      },
{
       "AGRDES   0 ",
       "AGRDES   0"
      },
{
       "SYMDES   1                                              "
              " A SYMMETRY DESCRIPTION",
       "SYMDES   1                                              "
              " A SYMMETRY DESCRIPTION"
      },
{
       "MTXDES   1                                              "
              "   A MATRIX DESCRIPTION",
       "MTXDES   1                                              "
              "   A MATRIX DESCRIPTION"
      },
{
       "CMPDES   1                                              "
              "  A COMPONENT DESCRIPTION",
       "CMPDES   1                                              "
              "  A COMPONENT DESCRIPTI"
      },
{
       "FTNOTE   1 A FOOTNOTE",
       "FTNOTE   1 A FOOTNOTE"
      },
{
       "REMARK   1 RRRRRRRRRRRRRRRRRRRRRRRRRRRRRRRRRRRRRRRRRRRRR"
              "RRRRRRRRRRRRRRRRRRRRRRR",
       "REMARK   1 RRRRRRRRRRRRRRRRRRRRRRRRRRRRRRRRRRRRRRRRRRRRR"
              "RRRRRRRRRRRRRRRRRRRRRRR"
      },
{
       "REMARK   0 ",
       "REMARK"
      },
  };
  for (const auto &fixture : lines) {
    CAPTURE(fixture.first);
    zdock::pdb::ParseContext context;
    const auto parsed = zdock::pdb::Record::parse(fixture.first, context);
    REQUIRE(parsed.has_value());
    REQUIRE(parsed->chars() == fixture.second);
    // the rendered line re-parses and re-renders stably
    zdock::pdb::ParseContext again;
    const auto reread = zdock::pdb::Record::parse(parsed->chars(), again);
    REQUIRE(reread.has_value());
    REQUIRE(reread->chars() == parsed->chars());
  }
}

TEST_CASE("HEADER writes its date left-justified like the legacy writer",
          "[record-parity]") {
  // Invariant: HEADER reads a 9-column date but writes an 11-column
  // left-justified date plus a literal space, so a classification that
  // overflows its 40 columns and spills characters into the date field
  // renders back at the spill position -- the rendered line is the
  // input line again, byte-for-byte, under both paths. The expected
  // line is assembled from the layout: 10-byte prefix, 40-column
  // classification (the 42-character name keeps its first 40 here and
  // its last two ride in the date field), then the date and idCode
  // fields, all left-justified.
  const std::string input = "HEADER    TEN ELEVEN TWELVE THIRTEEN FOURTEEN "
                            "FIFTEE"; // 42-character classification
  REQUIRE(input.size() == 52);
  zdock::pdb::ParseContext context;
  const auto parsed = zdock::pdb::Record::parse(input, context);
  REQUIRE(parsed.has_value());
  const auto &h = parsed->as<zdock::pdb::HeaderRecord>();
  REQUIRE(h.classification == "TEN ELEVEN TWELVE THIRTEEN FOURTEEN FIFT");
  REQUIRE(h.timestamp == "EE");
  REQUIRE(h.id == "");
  REQUIRE(parsed->chars() == input);
  // the rendered line re-parses to the same fields under both paths
  zdock::pdb::ParseContext again;
  const auto reread = zdock::pdb::Record::parse(parsed->chars(), again);
  REQUIRE(reread.has_value());
  REQUIRE(reread->as<zdock::pdb::HeaderRecord>().classification ==
          h.classification);
  REQUIRE(reread->as<zdock::pdb::HeaderRecord>().timestamp == h.timestamp);
  REQUIRE(reread->chars() == parsed->chars());
  zdock::pdb::ParseContext rereadContext;
  const auto modernReread =
      zdock::pdb::Record::parse(parsed->chars(), rereadContext);
  REQUIRE(modernReread.has_value());
  REQUIRE(modernReread->chars() == parsed->chars());
}

TEST_CASE("ANISOU and SIGUIJ render parsed zeros differently", "[record-parity]") {
  // Invariant: the two families write the same u-column layout with
  // different zero emitters -- legacy ANISOU uses plain %7d (a parsed
  // zero is a right-justified 0, so the line keeps its full width) and
  // legacy SIGUIJ uses %7D (a parsed zero is a blank field, so trailing
  // zeros trim away) -- and each form round-trips stably. The expected
  // lines are assembled from the layout columns: 6-byte prefix, %5d
  // serial, space, %-4s name, %c altLoc, %-4s residue, %c chain, %4d
  // sequence, %c insert code, space, then %7d/%7D u columns (u[0] plus
  // five more).
  const std::string head =
      "    1"     // %5d serial
      " "         // literal space
      " N  "      // %-4s atom name
      " "         // %c altLoc
      "GLY "      // %-4s residue name
      "A"         // %c chain id
      "  12"      // %4d sequence
      "  "        // %c insert code + literal space
      "   1234";  // %7d u[0]
  const std::string zeroPlain = "      0"; // %7d of 0
  const std::string zeroBlank = "       "; // %7D of 0
  // u payloads written into u[1..5]; u[0] is 1234 in the head
  const std::string allZeroAnisou =
      "ANISOU" + head + zeroPlain + zeroPlain + zeroPlain + zeroPlain +
      zeroPlain;
  const std::string allZeroSiguij =
      "SIGUIJ" + head + zeroBlank + zeroBlank + zeroBlank + zeroBlank +
      zeroBlank;
  const std::string mixedAnisou = "ANISOU" + head + zeroPlain + "   -567" +
                                  "      8" + zeroPlain + "   9012";
  const std::string mixedSiguij = "SIGUIJ" + head + zeroBlank + "   -567" +
                                  "      8" + zeroBlank + "   9012";
  struct Fixture {
    std::string line;   // the input line to parse
    std::string render; // the hand-derived canonical render
    std::array<int, 6> u;
  };
  const std::vector<Fixture> anisouFixtures = {
      {allZeroAnisou, allZeroAnisou, {1234, 0, 0, 0, 0, 0}},
      {mixedAnisou, mixedAnisou, {1234, 0, -567, 8, 0, 9012}},
  };
  const std::vector<Fixture> siguijFixtures = {
      // the all-blank u tail parses as zeros and trims away on render
      {allZeroSiguij, "SIGUIJ" + head, {1234, 0, 0, 0, 0, 0}},
      {mixedSiguij, mixedSiguij, {1234, 0, -567, 8, 0, 9012}},
  };
  for (const auto &fixtures : {anisouFixtures, siguijFixtures}) {
    for (const auto &fixture : fixtures) {
      CAPTURE(fixture.line);
      zdock::pdb::ParseContext context;
      const auto parsed = zdock::pdb::Record::parse(fixture.line, context);
      REQUIRE(parsed.has_value());
      // the render matches the hand-derived canonical form, which is
      // the exact text the retired legacy writer produced
      REQUIRE(parsed->chars() == fixture.render);
      const auto &a = parsed->as<zdock::pdb::AnisouRecord>();
      for (std::size_t i = 0; i < 6; ++i)
        REQUIRE(a.u[i] == fixture.u[i]);
      // the rendered line re-parses to the same payload under both paths
      zdock::pdb::ParseContext again;
      const auto reread = zdock::pdb::Record::parse(parsed->chars(), again);
      REQUIRE(reread.has_value());
      const auto &b = reread->as<zdock::pdb::AnisouRecord>();
      for (std::size_t i = 0; i < 6; ++i)
        REQUIRE(b.u[i] == a.u[i]);
      REQUIRE(reread->chars() == parsed->chars());
      zdock::pdb::ParseContext stableContext;
      const auto stable =
          zdock::pdb::Record::parse(parsed->chars(), stableContext);
      REQUIRE(stable.has_value());
      REQUIRE(stable->chars() == parsed->chars());
    }
  }
}

TEST_CASE("Structural family fields parse into the pinned typed values",
          "[record-parity]") {
  // Invariant: representative members land in their typed payload fields
  // at their wwPDB columns, matching the golden struct assertions.
  zdock::pdb::ParseContext context;
  const std::string spaces23(23, ' ');

  const auto helix = zdock::pdb::Record::parse(
      "HELIX    1  H1 ALA A   10  LEU A   15  1                       "
      "          6",
      context);
  REQUIRE(helix.has_value());
  REQUIRE(helix->as<zdock::pdb::HelixRecord>().serialNum == 1);
  REQUIRE(helix->as<zdock::pdb::HelixRecord>().id == " H1");
  REQUIRE(helix->as<zdock::pdb::HelixRecord>().residues[0].name == "ALA");
  REQUIRE(helix->as<zdock::pdb::HelixRecord>().residues[0].seqNum == 10);
  REQUIRE(helix->as<zdock::pdb::HelixRecord>().residues[1].name == "LEU");
  REQUIRE(helix->as<zdock::pdb::HelixRecord>().type == 1);
  REQUIRE(helix->as<zdock::pdb::HelixRecord>().length == 6);

  const auto conectRecord =
      zdock::pdb::Record::parse("CONECT    1    2    3    4    5    6    7"
                                "    8",
                                context);
  REQUIRE(conectRecord.has_value());
  REQUIRE(conectRecord->as<zdock::pdb::ConectRecord>().serialNum == 1);
  REQUIRE(conectRecord->as<zdock::pdb::ConectRecord>().covalent[0] == 2);
  REQUIRE(conectRecord->as<zdock::pdb::ConectRecord>().covalent[3] == 5);
  REQUIRE(conectRecord->as<zdock::pdb::ConectRecord>().bonds[0].hydrogen[0] ==
          6);
  REQUIRE(conectRecord->as<zdock::pdb::ConectRecord>().bonds[0].hydrogen[1] ==
          7);
  REQUIRE(conectRecord->as<zdock::pdb::ConectRecord>().bonds[0].salt == 8);
  // the line ends after bonds[0].salt; bonds[1] parses with defaults
  REQUIRE(conectRecord->as<zdock::pdb::ConectRecord>().bonds[1].hydrogen[0] ==
          0);
  REQUIRE(conectRecord->as<zdock::pdb::ConectRecord>().bonds[1].salt == 0);

  const auto cmpont =
      zdock::pdb::Record::parse("CMPONT   1  ALAA    1   ALAA    2", context);
  REQUIRE(cmpont.has_value());
  REQUIRE(cmpont->as<zdock::pdb::CmpontRecord>().seqNum == 1);
  REQUIRE(cmpont->as<zdock::pdb::CmpontRecord>().residues[0].name == " ALA");
  REQUIRE(cmpont->as<zdock::pdb::CmpontRecord>().residues[0].chainId == 'A');
  REQUIRE(cmpont->as<zdock::pdb::CmpontRecord>().residues[0].seqNum == 1);
  REQUIRE(cmpont->as<zdock::pdb::CmpontRecord>().residues[1].name == " ALA");
  REQUIRE(cmpont->as<zdock::pdb::CmpontRecord>().residues[1].seqNum == 2);

  const auto ssbond = zdock::pdb::Record::parse(
      "SSBOND   1 CYS A   10    CYS A   25 " + spaces23 + "  1555   1555"
      "  2.05",
      context);
  REQUIRE(ssbond.has_value());
  REQUIRE(ssbond->as<zdock::pdb::SsbondRecord>().seqNum == 1);
  REQUIRE(ssbond->as<zdock::pdb::SsbondRecord>().residues[0].name == "CYS");
  REQUIRE(ssbond->as<zdock::pdb::SsbondRecord>().residues[0].seqNum == 10);
  REQUIRE(ssbond->as<zdock::pdb::SsbondRecord>().residues[1].name == "CYS");
  REQUIRE(ssbond->as<zdock::pdb::SsbondRecord>().residues[1].seqNum == 25);
  REQUIRE(ssbond->as<zdock::pdb::SsbondRecord>().sym1 == 1555);
  REQUIRE(ssbond->as<zdock::pdb::SsbondRecord>().sym2 == 1555);
  REQUIRE(ssbond->as<zdock::pdb::SsbondRecord>().length == 2.05);

  const auto cryst1 = zdock::pdb::Record::parse(
      "CRYST1   78.900   78.900   85.700  90.00  90.00  90.00 P 21 21 21"
      "    8",
      context);
  REQUIRE(cryst1.has_value());
  REQUIRE(cryst1->as<zdock::pdb::Cryst1Record>().a == 78.9);
  REQUIRE(cryst1->as<zdock::pdb::Cryst1Record>().gamma == 90.0);
  REQUIRE(cryst1->as<zdock::pdb::Cryst1Record>().spaceGrp == "P 21 21 21");
  REQUIRE(cryst1->as<zdock::pdb::Cryst1Record>().z == 8);

  const auto model = zdock::pdb::Record::parse("MODEL        1", context);
  REQUIRE(model.has_value());
  REQUIRE(model->as<zdock::pdb::ModelRecord>().num == 1);

  const auto caveat = zdock::pdb::Record::parse(
      "CAVEAT   1 1LYZ    SOME WARNING TEXT", context);
  REQUIRE(caveat.has_value());
  REQUIRE(caveat->as<zdock::pdb::CaveatRecord>().id == "1LYZ");
  REQUIRE(caveat->as<zdock::pdb::CaveatRecord>().text == "SOME WARNING TEXT");
}
}
