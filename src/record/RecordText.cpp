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
  assert(std::strlen(prefix) == cols[0].offset);
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
