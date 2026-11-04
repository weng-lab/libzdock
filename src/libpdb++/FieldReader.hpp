// Copyright (c) 2019-2026 Arjan van der Velde, Weng Lab
// SPDX-License-Identifier: BSD-2-Clause

#pragma once

namespace libpdb {

// Stateful cursor over one PDB record line for fixed-width field
// extraction. Field semantics match the legacy varargs PDB::sscanf:
// fields starting at or past end-of-line yield their documented
// defaults (0, 0.0, empty string, or spaces) instead of failing,
// a value that does not fill its field is a hard mismatch, and
// every field started on the line counts toward the match count,
// including ignored fields.
class FieldReader {
public:
  explicit FieldReader(const char *line) noexcept;

  bool atEnd() const noexcept; // at end-of-line ('\0' or '\n')
  int nmatch() const noexcept; // number of fields started on the line

  // Consume one literal character; tolerated without consuming at
  // end-of-line. Returns false on a mismatch.
  bool matchLiteral(char expected) noexcept;

  bool readInt(int width, int *value);     // '%d'; default 0
  bool readReal(int width, double *value); // '%f'; default 0.0
  void readString(int width, char *value); // '%s'; trailing blanks stripped
  void readChars(int width, char *value);  // '%c'; no stripping, no terminator
  void skipChars(int width);               // '% '; ignored field

private:
  bool beginField(); // counts a started field; false at end-of-line

  const char *cursor_;
  int nmatch_;
};

} // namespace libpdb
