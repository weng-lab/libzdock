// Copyright (c) 2019-2026 Arjan van der Velde, Weng Lab
// SPDX-License-Identifier: BSD-2-Clause

#pragma once

#include <array>
#include <cstdint>
#include <iosfwd>
#include <optional>
#include <ostream>
#include <string>
#include <string_view>
#include <variant>

namespace zdock::pdb {

// Record types of the 1992 PDB format and the PDBRUN 5/6 USER extensions.
// The names carry the record spellings the fixed-format
// one-for-one, in the same order, so the two can be correlated during the
// Phase 3 migration.
enum class RecordType : std::uint8_t {
  UNKNOWN,
  ANISOU,
  ATOM,
  AUTHOR,
  CAVEAT,
  COMPND,
  CONECT,
  CRYST1,
  END,
  FORMUL,
  FTNOTE,
  HEADER,
  HELIX,
  HET,
  HETATM,
  JRNL,
  KEYWDS,
  MASTER,
  MDLTYP,
  MTRIX,
  NUMMDL,
  OBSLTE,
  ORIGX,
  REMARK,
  REVDAT,
  SCALE,
  SEQRES,
  SHEET,
  SIGATM,
  SIGUIJ,
  SITE,
  SOURCE,
  SPLIT,
  SPRSDE,
  SSBOND,
  TER,
  TITLE,
  TURN,
  TVECT,
  USER,
  MODEL,
  ENDMDL,
  EXPDTA,
  SYMDES,
  SYMOP,
  MTXDES,
  CMPDES,
  CMPONT,
  TRNSFM,
  AGRDES,
  AGGRGT,
  ZDATOM,
  USER_PDBRUN,
  USER_EYEPOS,
  USER_ATPOS,
  USER_WINDOW,
  USER_FOCUS,
  USER_VIEWPORT,
  USER_BGCOLOR,
  USER_ANGLE,
  USER_DISTANCE,
  USER_FILE,
  USER_MARKNAME,
  USER_MARK,
  USER_CNAME,
  USER_COLOR,
  USER_RADIUS,
  USER_OBJECT,
  USER_ENDOBJ,
  USER_CHAIN,
  USER_GFX_BEGIN,
  USER_GFX_END,
  USER_GFX_COLOR,
  USER_GFX_NORMAL,
  USER_GFX_VERTEX,
  USER_GFX_FONT,
  USER_GFX_TEXTPOS,
  USER_GFX_LABEL,
  USER_GFX_MOVE,
  USER_GFX_DRAW,
  USER_GFX_MARKER,
  USER_GFX_POINT,
};

// Returns the enumerator spelling (for example "ATOM") for diagnostics;
// values outside the enumeration render as "<invalid>".
std::string to_string(RecordType type);

// Graphics primitive names for the USER GFX BEGIN record, mirroring the
// legacy GfxType values.
enum class GfxType : std::uint8_t {
  UNKNOWN,
  POINTS,
  MARKERS,
  LINES,
  LINE_STRIP,
  LINE_LOOP,
  TRIANGLES,
  TRIANGLE_STRIP,
  TRIANGLE_FAN,
  QUADS,
  QUAD_STRIP,
  POLYGON,
};

// Residue identification shared by every residue-bearing record: name as
// parsed (trailing blanks stripped), chain and insertion codes as single
// characters defaulting to blanks.
struct Residue {
  std::string name;
  char chainId = ' ';
  int seqNum = 0;
  char insertCode = ' ';
};

// Payload for records with no fields (END, ENDMDL, USER GFX END).
struct Monostate {};

struct UnknownRecord {
  std::string junk; // the line right-aligned in 72 columns, '*' overflow
};

// The zlab mark_sur extension carried on an atom record; the standard
// ATOM columns are never overloaded with it.
struct Zdatom {
  int type = 0;
  int surface = 0;
  double rad = 0.0;
  std::string segid;
  double chg = 0.0;
};

// Shared by ATOM, HETATM, and SIGATM (identical layouts). Occupancy and
// tempFactor stay unset unless their columns parsed, so a clean ZDATOM
// parse never carries a leaked value from a failed standard attempt.
struct AtomRecord {
  int serialNum = 0;
  std::string name;
  char altLoc = ' ';
  Residue residue;
  std::array<double, 3> xyz{};
  std::optional<double> occupancy;
  std::optional<double> tempFactor;
  std::string element;
  std::string charge;
  // wwPDB segID (columns 73-76), kept when present on a line whose tail
  // breaks the pinned blank-column literal; empty for records parsed
  // through the strict standard layout, which requires those columns
  // blank
  std::string segid;
  std::optional<Zdatom> zdatom;
};

// Continuation numbers are typed ints; the legacy char union member that
// produced the misaligned %2d store is unrepresentable here.
struct AuthorRecord {
  int continuation = 0;
  std::string authors;
};

struct CompndRecord {
  int continuation = 0;
  std::string compound;
};

struct KeywdsRecord {
  int continuation = 0;
  std::string keywds;
};

struct SourceRecord {
  int continuation = 0;
  std::string name;
};

struct MdltypRecord {
  int continuation = 0;
  std::string comment;
};

struct ExpdtaRecord {
  int continuation = 0;
  std::string technique;
};

struct TitleRecord {
  int continuation = 0;
  std::string text;
};

// Shared by AGRDES, CMPDES, MTXDES, SYMDES, FTNOTE, and REMARK.
struct AgrdesRecord {
  int num = 0;
  std::string text;
};

struct CaveatRecord {
  int continuation = 0;
  std::string id;
  std::string text;
};

struct HeaderRecord {
  std::string classification;
  std::string timestamp;
  std::string id;
};

struct AggrgtRecord {
  int serialNum = 0;
  int numComponents = 0;
  std::array<int, 14> cmpontSerialNums{};
};

struct AnisouRecord { // shared by ANISOU and SIGUIJ
  int serialNum = 0;
  std::string name;
  char altLoc = ' ';
  Residue residue;
  std::array<int, 6> u{};
};

struct CmpontRecord {
  int seqNum = 0;
  std::array<Residue, 2> residues;
};

struct ConectRecord {
  struct Bond {
    std::array<int, 2> hydrogen{};
    int salt = 0;
  };
  int serialNum = 0;
  std::array<int, 4> covalent{};
  std::array<Bond, 2> bonds;
};

struct Cryst1Record {
  double a = 0.0;
  double b = 0.0;
  double c = 0.0;
  double alpha = 0.0;
  double beta = 0.0;
  double gamma = 0.0;
  std::string spaceGrp;
  int z = 0;
};

struct FormulRecord {
  int component = 0;
  std::string hetId;
  int continuation = 0;
  char exclude = ' '; // '*' excludes the residue from the formula
  std::string formula;
};

struct HelixRecord {
  int serialNum = 0;
  std::string id;
  std::array<Residue, 2> residues;
  int type = 0;
  std::string comment;
  int length = 0;
};

struct HetRecord {
  Residue hetGrp;
  int numAtoms = 0;
  std::string text;
};

struct JrnlRecord {
  std::string text;
};

struct MasterRecord {
  int numRemark = 0;
  int numFtnote = 0;
  int numHet = 0;
  int numHelix = 0;
  int numSheet = 0;
  int numTurn = 0;
  int numSite = 0;
  int numTransform = 0;
  int numCoordinate = 0;
  int numTer = 0;
  int numConect = 0;
  int numSeqres = 0;
};

struct ModelRecord {
  int num = 0;
};

struct NummdlRecord {
  int number = 0;
};

struct MtrixRecord {
  int rowNum = 0;
  int serialNum = 0;
  double m1 = 0.0;
  double m2 = 0.0;
  double m3 = 0.0;
  double v = 0.0;
  int given = 0;
};

struct OrigxRecord {
  int rowNum = 0;
  double o1 = 0.0;
  double o2 = 0.0;
  double o3 = 0.0;
  double t = 0.0;
};

struct ScaleRecord {
  int rowNum = 0;
  double s1 = 0.0;
  double s2 = 0.0;
  double s3 = 0.0;
  double u = 0.0;
};

struct SymopRecord {
  int rowNum = 0;
  int serialNum = 0;
  double s1 = 0.0;
  double s2 = 0.0;
  double s3 = 0.0;
  double t = 0.0;
};

struct ObslteRecord {
  int continuation = 0;
  std::string timestamp;
  std::string oldId;
  std::array<std::string, 9> idMap;
};

struct RevdatRecord {
  int modification = 0;
  int continuation = 0;
  std::string timestamp;
  std::string id;
  int modType = 0;
  std::array<std::string, 4> record;
};

struct SeqresRecord {
  int serialNum = 0;
  char chainId = ' ';
  int count = 0;
  std::array<std::string, 13> names;
};

struct SheetRecord {
  struct Atom {
    std::string name;
    Residue residue;
  };
  int strandNum = 0;
  std::string id;
  int count = 0;
  std::array<Residue, 2> residues;
  int sense = 0;
  std::array<Atom, 2> atoms;
};

struct SiteRecord {
  int seqNum = 0;
  std::string id;
  int count = 0;
  std::array<Residue, 4> residues;
};

struct SplitRecord {
  int continuation = 0;
  std::array<std::string, 14> idMap;
};

struct SprsdeRecord {
  int continuation = 0;
  std::string timestamp;
  std::string id;
  std::array<std::string, 9> supersede;
};

struct SsbondRecord {
  int seqNum = 0;
  std::array<Residue, 2> residues;
  int sym1 = 0;
  int sym2 = 0;
  double length = 0.0;
};

struct TerRecord {
  int serialNum = 0;
  Residue residue;
};

struct TrnsfmRecord {
  int resultSerialNum = 0;
  int applySerialNum = 0;
  int sourceSerialNum = 0;
};

struct TurnRecord {
  int seqNum = 0;
  std::string id;
  std::array<Residue, 2> residues;
  std::string comment;
};

struct TvectRecord {
  int serialNum = 0;
  double t1 = 0.0;
  double t2 = 0.0;
  double t3 = 0.0;
  std::string comment;
};

struct UserRecord {
  std::string subtype;
  std::string text;
};

struct UserPdbrunRecord {
  int version = 0;
};

struct UserEyePosRecord { // shared by USER EYEPOS and USER ATPOS
  std::array<double, 3> xyz{};
};

struct UserWindowRecord {
  double left = 0.0;
  double right = 0.0;
  double bottom = 0.0;
  double top = 0.0;
  double hither = 0.0;
  double yon = 0.0;
};

struct UserFocusRecord {
  double focus = 0.0;
};

struct UserViewportRecord {
  double xmin = 0.0;
  double xmax = 0.0;
  double ymin = 0.0;
  double ymax = 0.0;
};

struct UserBgColorRecord {
  std::array<double, 3> rgb{};
};

// PDBRUN v5 carries the obsolete leading "which" field; it stays unset
// under v6 instead of shifting the atom numbers.
struct UserAngleRecord {
  std::optional<int> which;
  int atom0 = 0;
  int atom1 = 0;
  int atom2 = 0;
  int atom3 = 0;
  double angle = 0.0;
};

struct UserDistanceRecord {
  std::optional<int> which;
  int atom0 = 0;
  int atom1 = 0;
  double distance = 0.0;
};

// v5 has no model column; it stays unset under v5.
struct UserFileRecord {
  std::optional<int> model;
  std::string filename;
};

struct UserMarknameRecord { // shared by USER MARKNAME and USER MARK
  std::string markname;
};

struct UserCNameRecord {
  std::array<double, 3> rgb{};
  std::string name;
};

struct UserColorRecord { // shared by USER COLOR and USER GFX COLOR
  std::array<double, 3> rgb{};
  std::string spec;
};

struct UserRadiusRecord {
  double radius = 0.0;
};

struct UserObjectRecord { // shared by USER OBJECT and USER ENDOBJ
  int model = 0;
};

struct UserChainRecord {
  int atom0 = 0;
  int atom1 = 0;
};

struct UserGfxBeginRecord {
  GfxType primitive = GfxType::UNKNOWN;
  std::string unknown; // the captured text when the primitive is unknown
};

struct UserGfxFontRecord {
  int size = 0;
  std::string name;
};

// v5 labels carry coordinates before the text; they stay unset under v6.
struct UserGfxLabelRecord {
  std::optional<std::array<double, 3>> xyz;
  std::string text;
};

// Shared by USER GFX NORMAL, VERTEX, TEXTPOS, MOVE, DRAW, MARKER, and
// POINT (identical layouts).
struct UserGfxPointRecord {
  std::array<double, 3> xyz{};
};

// One payload alternative per distinct record layout; record kinds with
// identical columns share a struct (as the legacy typedefs did).
using RecordPayload = std::variant<
    Monostate, UnknownRecord, AtomRecord, AuthorRecord, TitleRecord,
    CompndRecord, KeywdsRecord, SourceRecord, MdltypRecord, ExpdtaRecord,
    AgrdesRecord, CaveatRecord, HeaderRecord, AggrgtRecord, AnisouRecord,
    CmpontRecord, ConectRecord, Cryst1Record, FormulRecord, HelixRecord,
    HetRecord, JrnlRecord, MasterRecord, ModelRecord, MtrixRecord, NummdlRecord,
    ObslteRecord, OrigxRecord, RevdatRecord, ScaleRecord, SeqresRecord,
    SheetRecord, SiteRecord, SplitRecord, SprsdeRecord, SsbondRecord,
    SymopRecord, TerRecord, TrnsfmRecord, TurnRecord, TvectRecord, UserRecord,
    UserPdbrunRecord, UserEyePosRecord, UserWindowRecord, UserFocusRecord,
    UserViewportRecord, UserBgColorRecord, UserAngleRecord, UserDistanceRecord,
    UserFileRecord, UserMarknameRecord, UserCNameRecord, UserColorRecord,
    UserRadiusRecord, UserObjectRecord, UserChainRecord, UserGfxBeginRecord,
    UserGfxFontRecord, UserGfxLabelRecord, UserGfxPointRecord>;

// Caller-owned PDBRUN version state: parsing a USER PDBRUN record updates
// the context owned by the caller instead of process globals, and
// formatting takes the output version explicitly. (The USER records
// themselves land in a later Phase 2 commit.)
struct ParseContext {
  int pdbrunVersion = 6;
  // Input-quality diagnostic: the number of lines this context has seen
  // that announce themselves as coordinate data (ATOM, HETATM, or a
  // mark_sur line) yet parse as no record at all -- the caller's signal
  // that a file is losing atoms, never silently swallowed as UNKNOWN.
  std::size_t coordinateRejections = 0;
};

struct WriteContext {
  int pdbrunVersion = 6;
};

// Default-constructs the payload a record of the given type carries.
RecordPayload payloadFor(RecordType type);

// A typed record: the record kind plus its payload, with wrong-kind field
// access a typed error instead of a silent union misread.
class Record {
public:
  Record(); // UNKNOWN with empty junk
  explicit Record(RecordType type);

  RecordType type() const noexcept { return type_; }

  template <typename T> bool holds() const {
    return std::holds_alternative<T>(payload_);
  }

  template <typename T> const T &as() const { return std::get<T>(payload_); }

  template <typename T> T &as() { return std::get<T>(payload_); }

  bool isAtom() const noexcept {
    return type_ == RecordType::ATOM || type_ == RecordType::HETATM ||
           type_ == RecordType::SIGATM;
  }

  const AtomRecord &asAtom() const { return std::get<AtomRecord>(payload_); }
  AtomRecord &asAtom() { return std::get<AtomRecord>(payload_); }

  // Parsing: classifies the line and walks the layout table for its
  // family. Returns std::nullopt on a hard field mismatch (the caller
  // falls back to unknownFrom, as the legacy parse did) and for record
  // families whose tables have not landed yet. USER PDBRUN parsing
  // updates the context owned by the caller, never process-global state.
  static std::optional<Record> parse(std::string_view line,
                                     ParseContext &context);

  // True when the line announces itself as coordinate data: an ATOM or
  // HETATM record name (mark_sur lines included -- they also start with
  // ATOM). Callers use it together with ParseContext::
  // coordinateRejections to warn or fail loudly when a file loses atoms.
  static bool looksLikeCoordinate(std::string_view line);

  // The UNKNOWN junk fallback: the line right-aligned in 72 columns with
  // '*' overflow, exactly as the legacy setUnknown produced.
  static Record unknownFrom(std::string_view line);

  // Formatting: the record text with trailing blanks trimmed, bounded by
  // the declared record length; thread-safe by construction. Record
  // families whose write tables have not landed yet format as empty.
  std::string chars() const;
  std::string chars(const WriteContext &context) const;

  // Allocation-free formatting into the caller's record buffer (sized for
  // the widest declared layout plus a terminator); renders the same text
  // as chars(), NUL-terminated, returning the trimmed length. The
  // hot-path counterpart of chars() for streaming loops.
  std::size_t write_into(char (&out)[83], const WriteContext &context) const;

private:
  RecordType type_;
  RecordPayload payload_;
};

inline std::ostream &operator<<(std::ostream &s, const Record &record) {
  s << record.chars();
  return s;
}

// The pinned input-line semantics, carried over from the retired
// the retired upstream reader: read one record line through a 4x81-character
// std::istream::getline buffer, strip exactly one trailing CR, and
// fail the stream for lines that do not fit (production read loops
// stop at that point). Returns false when the stream failed.
bool readRecordLine(std::istream &in, std::string &line);

} // namespace zdock::pdb
