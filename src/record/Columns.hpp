// Copyright (c) 2019-2026 Arjan van der Velde, Weng Lab
// SPDX-License-Identifier: BSD-2-Clause

#pragma once

#include <cstddef>
#include <cstdint>
#include <string>
#include <string_view>

namespace zdock::pdb {

// The largest total digit count the fixed-format parser accepts in a Real
// column. Every layout table keeps its Real fields within this bound: the
// exact decimal D / 10^m with D < 10^10 and m <= 10 is exactly
// representable, so one correctly-rounded division by a static power of
// ten reproduces strtod bit-for-bit. A table claiming a wider Real field
// is a table-definition error, not a runtime case.
inline constexpr int kMaxRealDigits = 10;

// The largest total digit count the fixed-format parser accepts in an Int
// column; D < 10^18 always fits a 64-bit accumulator.
inline constexpr int kMaxIntDigits = 18;

// Column kinds of a record layout. Skip columns advance without reading
// (the record-name field, blank gaps); Int and Real columns convert with
// the shared fixed-format parser; Str columns carry trailing-blank-
// stripped text; Char columns capture verbatim, space-padded.
enum class ColumnKind : std::uint8_t { Skip, Int, Real, Str, Char };

// One fixed-width column of a record layout: its start offset in the
// 0-based record line, its width, and its kind. The write-side emitter
// style rides along (ignored while parsing) so one table per record can
// drive both walks where the read and write layouts coincide.
struct ColumnSpec {
  std::uint16_t offset;
  std::uint8_t width;
  ColumnKind kind;
  std::uint8_t places = 0;    // Real: digits after the decimal point
  bool leftJustify = false;   // write: left-justified field
  char fill = ' ';            // write: pad character (' ' or '0')
  bool zeroAsBlank = false;   // write Int: zero renders as blanks (%D)
  std::int8_t precision = -1; // write Str: maximum length, -1 unbounded
};

// Bounded, table-driven field extraction over one record line, with the
// pinned legacy field semantics: fields starting at or past end-of-line
// yield their documented defaults instead of failing (0, 0.0, empty
// string, or spaces), every field started on the line counts toward the
// match count including skipped fields, and a numeric conversion is a
// hard mismatch unless the whole field matches the fixed-format grammar.
class ColumnReader {
public:
  explicit ColumnReader(std::string_view line) noexcept;

  bool atEnd() const noexcept; // at end-of-line (NUL, newline, or past it)
  int nmatch() const noexcept; // number of fields started on the line

  // Consume one literal character; tolerated without consuming at
  // end-of-line. Returns false on a mismatch.
  bool matchLiteral(char expected) noexcept;

  bool readInt(int width, std::int64_t *value);     // default 0
  bool readReal(int width, double *value);          // default 0.0
  std::string readString(int width);                // trailing blanks stripped
  std::string readChars(int width);                 // verbatim, space-padded
  void skipChars(int width);                        // ignored field

private:
  bool beginField(); // counts a started field; false at end-of-line
  std::string_view slice(int width); // consume up to width chars to EOL

  std::string_view line_;
  std::size_t cursor_;
  int nmatch_;
};

// Fixed-width column formatting into a caller-owned buffer sized to the
// declared record length and pre-filled with spaces, with the pinned
// legacy emitter semantics: right or left justification, blank or zero
// fill, '*' overflow markers, and %D integers that render zero as blanks.
// Every field renders exactly its descriptor width at its descriptor
// offset, so the writer is bounded by construction.
class ColumnWriter {
public:
  explicit ColumnWriter(std::string &out) noexcept;
  ColumnWriter(char *out, std::size_t size) noexcept;

  // Writes raw text (the record-name prefix) at an absolute offset.
  void writeLiteral(std::size_t offset, std::string_view text);

  void writeInt(std::int64_t value, const ColumnSpec &spec);
  void writeReal(double value, const ColumnSpec &spec);
  void writeString(std::string_view s, const ColumnSpec &spec);
  void writeChar(char c, const ColumnSpec &spec);

private:
  char *field(const ColumnSpec &spec); // start of the field window

  char *out_;
  std::size_t size_;
};

} // namespace zdock::pdb
