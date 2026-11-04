// Copyright (c) 2019-2026 Arjan van der Velde, Weng Lab
// SPDX-License-Identifier: BSD-2-Clause

#pragma once

#include "pdb++.h"

#include <cstddef>

namespace libpdb {

// Format and field-descriptor tables for the 1992 record format and the
// PDBRUN 5/6 USER extensions. Owns the legacy format string tables and the
// per-record field lists (field kind plus byte offset within the record
// storage) that replace the per-type varargs argument bindings. Field
// widths and literals stay in the format strings; the field lists pair
// with the value conversions in order.
class PdbFormat {
public:
  enum class FieldKind { Int, Real, Str, Char, End };

  struct Field {
    FieldKind kind;
    std::size_t offset; // byte offset of the target within the record storage
  };

  struct ReadSpec {
    const char *fmt;      // nullptr: nothing is parsed
    bool libcBackend;     // v5 USER paths consumed by libc sscanf
    bool fallbackToUser;  // parse failure falls back to a USER record
    const Field *fields;  // End-terminated, pairs with the value conversions
  };

  struct WriteSpec {
    const char *fmt;     // nullptr: the record formats as empty (count 0)
    bool libcBackend;    // v5 USER paths consumed by libc snprintf
    const Field *fields; // End-terminated, pairs with the value conversions
  };

  // Parses buf into record according to the read spec for `type` and
  // returns the final record type (the type itself, USER after the USER
  // fallback, or UNKNOWN for junk records). Applies the pinned side
  // effects: USER PDBRUN updates the process-global input version and a
  // parsed USER GFX BEGIN records its primitive.
  static PDB::RecordType parse(PDB &record, const char *buf,
                               PDB::RecordType type, int pdbrunInputVersion);

  static ReadSpec readSpec(PDB::RecordType type, int pdbrunInputVersion);

  // Formats record into buf according to the write spec for `type` and
  // returns the produced length. Preserves the pinned static-buffer and
  // global-version behaviors of PDB::chars: the ZDATOM layout still
  // overruns BufLen before the caller clamps, and formatting USER PDBRUN
  // updates the process-global output version.
  static int formatRecord(const PDB &record, char *buf, PDB::RecordType type,
                          int pdbrunOutputVersion);

  static WriteSpec writeSpec(PDB::RecordType type, int pdbrunOutputVersion);
};

} // namespace libpdb
