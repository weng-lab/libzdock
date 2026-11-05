// Copyright (c) 2019-2026 Arjan van der Velde, Weng Lab
// SPDX-License-Identifier: BSD-2-Clause

#include "Columns.hpp"

#include <cassert>
#include <cctype>
#include <algorithm>

namespace zdock::pdb {
namespace {

// Exact powers of ten for Real conversion: the layouts are width-only, so
// the scale comes from the field itself, exactly as legacy strtod saw it.
constexpr double kPow10[kMaxRealDigits + 1] = {1e0,  1e1,  1e2,  1e3,
                                               1e4,  1e5,  1e6,  1e7,
                                               1e8,  1e9,  1e10};

// The supported fixed-format grammar (everything else is rejected, never
// guessed): optional leading blanks, an optional sign, and a decimal
// number with at least one digit — digits with an optional point and
// fraction, or a point followed by digits. Inner blanks, tabs, carriage
// returns, exponents, inf/nan, and hex are hard mismatches; legacy
// strtod's acceptance of leading tabs and exponent forms is a deliberate,
// pinned divergence (fixed-format PDB never contains them).
//
// On rejection the accumulated prefix value is still stored (legacy
// strtol/strtod semantics), so callers observe the same partial value.
//
// keepPoint (Real) allows a trailing point with no fraction digits, as
// strtod does; the Int grammar rejects it. The return value reports a
// full match; *digits and *fraction carry the accumulated numerator and
// its scale on both full matches and rejections.
bool scanFixed(std::string_view field, bool keepPoint, std::uint64_t *digits,
               int *fraction, bool *negative) {
  std::size_t i = 0;
  const std::size_t n = field.size();
  while (i < n && field[i] == ' ')
    ++i;
  *negative = i < n && field[i] == '-';
  if (*negative || (i < n && field[i] == '+'))
    ++i;

  std::uint64_t d = 0;
  int m = 0;
  int total = 0;
  bool anyDigit = false;
  while (i < n && std::isdigit(static_cast<unsigned char>(field[i]))) {
    d = d * 10 + static_cast<std::uint64_t>(field[i] - '0');
    ++total;
    anyDigit = true;
    ++i;
  }
  if (keepPoint && i < n && field[i] == '.') {
    ++i;
    while (i < n && std::isdigit(static_cast<unsigned char>(field[i]))) {
      d = d * 10 + static_cast<std::uint64_t>(field[i] - '0');
      ++m;
      ++total;
      anyDigit = true;
      ++i;
    }
  }
  *digits = d;
  *fraction = m;
  if (!anyDigit || i != n)
    return false;
  const int cap = keepPoint ? kMaxRealDigits : kMaxIntDigits;
  return total <= cap;
}

} // namespace

ColumnReader::ColumnReader(std::string_view line) noexcept
    : line_(line), cursor_(0), nmatch_(0) {}

bool ColumnReader::atEnd() const noexcept {
  return cursor_ >= line_.size() || line_[cursor_] == '\0' ||
         line_[cursor_] == '\n';
}

int ColumnReader::nmatch() const noexcept { return nmatch_; }

bool ColumnReader::matchLiteral(char expected) noexcept {
  if (cursor_ < line_.size() && line_[cursor_] == expected) {
    ++cursor_;
    return true;
  }
  return atEnd(); // a literal past end-of-line is tolerated
}

bool ColumnReader::beginField() {
  if (atEnd())
    return false;
  ++nmatch_;
  return true;
}

std::string_view ColumnReader::slice(int width) {
  std::size_t end = cursor_;
  while (end < line_.size() && line_[end] != '\0' && line_[end] != '\n' &&
         end - cursor_ < static_cast<std::size_t>(width))
    ++end;
  std::string_view s = line_.substr(cursor_, end - cursor_);
  cursor_ = end;
  return s;
}

bool ColumnReader::readInt(int width, std::int64_t *value) {
  if (!beginField()) {
    *value = 0;
    return true;
  }
  std::string_view field = slice(width);
  while (!field.empty() &&
         std::isspace(static_cast<unsigned char>(field.back())))
    field.remove_suffix(1);
  if (field.empty()) { // an all-blank field parses as 0, exactly as
    *value = 0;        // legacy strtol("") did, with the field started
    return true;
  }
  std::uint64_t d = 0;
  int m = 0;
  bool negative = false;
  const bool ok = scanFixed(field, false, &d, &m, &negative);
  std::int64_t v = static_cast<std::int64_t>(d);
  *value = negative ? -v : v;
  return ok;
}

bool ColumnReader::readReal(int width, double *value) {
  if (!beginField()) {
    *value = 0.0;
    return true;
  }
  std::string_view field = slice(width);
  while (!field.empty() &&
         std::isspace(static_cast<unsigned char>(field.back())))
    field.remove_suffix(1);
  if (field.empty()) { // an all-blank field parses as 0.0, exactly as
    *value = 0.0;      // legacy strtod("") did, with the field started
    return true;
  }
  std::uint64_t d = 0;
  int m = 0;
  bool negative = false;
  const bool ok = scanFixed(field, true, &d, &m, &negative);
  double v = static_cast<double>(d);
  if (negative)
    v = -v;
  *value = v / kPow10[m];
  return ok;
}

std::string ColumnReader::readString(int width) {
  if (!beginField())
    return {};
  std::string_view field = slice(width);
  while (!field.empty() &&
         std::isspace(static_cast<unsigned char>(field.back())))
    field.remove_suffix(1);
  return std::string(field);
}

std::string ColumnReader::readChars(int width) {
  std::string value(static_cast<std::size_t>(width), ' ');
  if (!beginField())
    return value;
  std::string_view field = slice(width);
  field.copy(value.data(), field.size());
  return value;
}

void ColumnReader::skipChars(int width) {
  if (!beginField())
    return;
  slice(width);
}

ColumnWriter::ColumnWriter(std::string &out) noexcept
    : out_(out.data()), size_(out.size()) {}

ColumnWriter::ColumnWriter(char *out, std::size_t size) noexcept
    : out_(out), size_(size) {}

char *ColumnWriter::field(const ColumnSpec &spec) {
  assert(static_cast<std::size_t>(spec.offset) + spec.width <= size_);
  return out_ + spec.offset;
}

void ColumnWriter::writeLiteral(std::size_t offset, std::string_view text) {
  assert(offset + text.size() <= size_);
  text.copy(out_ + offset, text.size());
}

// The pinned integer emitter: right or left justification within the
// fixed width, ' ' or '0' fill (zero fill moves the sign to the field
// start), %D zero-as-blank, and '*' overflow markers filling the field.
void ColumnWriter::writeInt(std::int64_t value, const ColumnSpec &spec) {
  char *const p = field(spec);
  char *cursor = p;
  int width = static_cast<int>(spec.width);
  const char fill = spec.fill;
  const bool leftJustify = spec.leftJustify;
  const char zero = spec.zeroAsBlank ? ' ' : '0';
  const bool negative = value < 0;
  const std::uint64_t magnitude =
      negative ? 0u - static_cast<std::uint64_t>(value)
               : static_cast<std::uint64_t>(value);

  if (negative)
    width--; // the sign shares the field, as in the legacy emitter
  char digits[24];
  int n = 0;
  std::uint64_t v = magnitude;
  if (v)
    do {
      digits[n++] = static_cast<char>('0' + v % 10);
      v /= 10;
    } while (v);
  else
    digits[n++] = zero;
  const int ndigits = n;
  if (n > width) {
    std::fill_n(p, spec.width, '*');
    return;
  }

  if (negative && fill == '0')
    *cursor++ = '-';
  if (!leftJustify)
    for (int pad = width - n; pad > 0; pad--)
      *cursor++ = fill;
  if (negative && fill == ' ')
    *cursor++ = '-';
  while (n > 0)
    *cursor++ = digits[--n];
  if (leftJustify)
    for (int pad = width - ndigits; pad > 0; pad--)
      *cursor++ = fill;
  assert(cursor == p + spec.width);
}

// The pinned fixed-point emitter, reproducing the legacy rounding
// ((value + 0.5) truncated at the scaled integer), the '*', and the sign
// placement conventions of the legacy %<w>.<p>f writer exactly.
void ColumnWriter::writeReal(double value, const ColumnSpec &spec) {
  const std::size_t offset = spec.offset;
  const std::size_t width = spec.width;
  const int places = spec.places;
  const char fill = spec.fill;
  const bool leftJustify = spec.leftJustify;

  bool negative = value < 0.0;
  if (negative)
    value = -value;
  for (int i = 0; i < places; i++)
    value *= 10.0;
  int intval = static_cast<int>(value + 0.5);

  if (places + (places == 0 ? 1 : 2) > static_cast<int>(width)) {
    std::fill_n(out_ + offset, width, '*');
    return;
  }

  char *const p = out_ + offset;
  char *place = p + width - 1;
  for (; place >= p + width - places; place--) {
    *place = static_cast<char>('0' + intval % 10);
    intval /= 10;
  }
  if (places > 0)
    *place-- = '.';
  if (intval == 0)
    *place-- = '0';
  for (; place >= p; place--) {
    if (intval == 0)
      break;
    *place = static_cast<char>('0' + intval % 10);
    intval /= 10;
  }
  if (intval != 0) {
    std::fill_n(p, width, '*');
    return;
  }
  if (place < p && negative) {
    std::fill_n(p, width, '*');
    return;
  }

  if (leftJustify) {
    char *from = place + 1;
    char *to = negative ? p + 1 : p;
    for (; from < p + width; from++, to++)
      *to = *from;
    for (; to < p + width; to++)
      *to = fill;
    if (negative)
      *p = '-';
  } else {
    for (char *to = place; to >= p; to--)
      *to = fill;
    if (negative) {
      if (fill == ' ')
        *place = '-';
      else
        *p = '-';
    }
  }
}

// The pinned string emitter: an optional precision truncates, an
// overlong value overflows to '*', otherwise the value is justified in
// the field with the fill character.
void ColumnWriter::writeString(std::string_view s, const ColumnSpec &spec) {
  const std::size_t offset = spec.offset;
  const std::size_t width = spec.width;
  std::size_t len = s.size();
  if (spec.precision >= 0 && len > static_cast<std::size_t>(spec.precision))
    len = static_cast<std::size_t>(spec.precision);
  if (len > width) {
    std::fill_n(out_ + offset, width, '*');
    return;
  }

  char *p = out_ + offset;
  char *const end = p + width;
  if (!spec.leftJustify)
    while (end - p > static_cast<std::ptrdiff_t>(len))
      *p++ = spec.fill;
  for (std::size_t i = 0; i < len; i++)
    *p++ = s[i];
  if (spec.leftJustify)
    while (p < end)
      *p++ = spec.fill;
}

void ColumnWriter::writeChar(char c, const ColumnSpec &spec) {
  if (c == '\0')
    c = ' ';
  char *p = out_ + spec.offset;
  char *const end = p + spec.width;
  if (spec.leftJustify)
    *p++ = c;
  while (end - p > 1)
    *p++ = spec.fill;
  if (!spec.leftJustify)
    *p = c;
}

} // namespace zdock::pdb
