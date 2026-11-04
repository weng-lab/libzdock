// Copyright (c) 2019-2026 Arjan van der Velde, Weng Lab
// SPDX-License-Identifier: BSD-2-Clause

#pragma once

namespace libpdb {

// Fixed-width field formatting for PDB record text, matching the legacy
// varargs PDB::sprintf field semantics: fixed field widths with left or
// right justification, blank or zero fill, and '*' overflow markers.
// Fields are appended to the buffer supplied at construction; as in the
// legacy formatter, the caller owns that buffer and its capacity.
class FieldWriter {
public:
  explicit FieldWriter(char *out) noexcept;

  void writeChar(char c);
  void writePaddedChar(char c, int width, char fill, bool leftJustify);

  void writeInt(int value, int width, int radix, char fill, char hex,
                bool leftJustify, char zero);
  void writeUnsigned(unsigned value, int width, char fill, bool leftJustify);
  void writeString(const char *s, int width, int maxLen, char fill,
                   bool leftJustify);
  void writeFloat(double value, int width, int places, char fill,
                  bool leftJustify);
  void writeExp(double value, int width, int places, char fill,
                bool leftJustify);

  void terminate();           // NUL-terminate at the write cursor
  int count() const noexcept; // characters written so far

private:
  void overflow(int width); // emit '*' markers for an overflowing field

  char *out_;      // start of the caller-owned buffer
  char *cursor_;   // write cursor
  char scratch_[256]; // transient digit buffer, must hold the largest number
};

} // namespace libpdb
