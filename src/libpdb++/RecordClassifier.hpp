// Copyright (c) 2019-2026 Arjan van der Velde, Weng Lab
// SPDX-License-Identifier: BSD-2-Clause

#pragma once

#include "pdb++.h"

namespace libpdb {

// Record-name dispatch for the 1992 format and the PDBRUN 5/6 USER
// extensions: maps the leading characters of a record line to a record
// type and maps graphics primitive names to their enumeration values
// and back. The PDBRUN input version is supplied per call; this class
// holds no state.
class RecordClassifier {
public:
  static PDB::RecordType getType(const char *buf, int pdbrunInputVersion);
  static PDB::GfxType getGfxType(const char *buf);
  static const char *gfxChars(PDB::GfxType gt);

private:
  RecordClassifier() = delete;

  static PDB::RecordType pdbrun5Type(const char *buf);
  static PDB::RecordType pdbrun6Type(const char *buf);
};

} // namespace libpdb
