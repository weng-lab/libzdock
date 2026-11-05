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
  // The declared record length; the ZDATOM layout extends to 82 columns in
  // a later commit and still fits the buffer.
  constexpr std::size_t kLength = 80;
  static_assert(kLength + 1 <= sizeof(out), "record buffer too small");
  std::memset(out, ' ', kLength);
  ColumnWriter writer(out, kLength);
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
  std::size_t length = kLength;
  while (length > 1 &&
         std::isspace(static_cast<unsigned char>(out[length - 1])))
    --length;
  out[length] = '\0';
  return length;
}

} // namespace zdock::pdb
