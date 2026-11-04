// Copyright (c) 2019-2026 Arjan van der Velde, Weng Lab
// SPDX-License-Identifier: BSD-2-Clause

#pragma once

#include <cstdint>
#include <string>

namespace zdock::pdb {

// Record types of the 1992 PDB format and the PDBRUN 5/6 USER extensions.
// The names mirror the legacy libpdb::PDB::RecordType enumerators
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

} // namespace zdock::pdb
