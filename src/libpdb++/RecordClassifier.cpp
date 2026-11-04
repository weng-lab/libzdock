// Copyright (c) 2019-2026 Arjan van der Velde, Weng Lab
// SPDX-License-Identifier: BSD-2-Clause

#include "RecordClassifier.hpp"

extern "C" {
#include <ctype.h>
#include <string.h>
#include <strings.h>
}

#ifndef _toupper
#define _toupper toupper
#endif

namespace libpdb {

PDB::GfxType RecordClassifier::getGfxType(const char *buf) {
  switch (buf[0]) {
  case 'L':
  case 'l':
    if (strcasecmp(buf + 1, "INE-LOOP") == 0)
      return PDB::GFX_LINE_LOOP;
    if (strcasecmp(buf + 1, "INE-STRIP") == 0)
      return PDB::GFX_LINE_STRIP;
    if (strcasecmp(buf + 1, "INES") == 0)
      return PDB::GFX_LINES;
    break;
  case 'M':
  case 'm':
    if (strcasecmp(buf + 1, "ARKERS") == 0)
      return PDB::GFX_MARKERS;
    break;
  case 'P':
  case 'p':
    if (strcasecmp(buf + 1, "OINTS") == 0)
      return PDB::GFX_POINTS;
    if (strcasecmp(buf + 1, "OLYGON") == 0)
      return PDB::GFX_POLYGON;
    break;
  case 'Q':
  case 'q':
    if (strcasecmp(buf + 1, "UAD-STRIP") == 0)
      return PDB::GFX_QUAD_STRIP;
    if (strcasecmp(buf + 1, "UADS") == 0)
      return PDB::GFX_QUADS;
    break;
  case 'T':
  case 't':
    if (strcasecmp(buf + 1, "RIANGLE-FAN") == 0)
      return PDB::GFX_TRIANGLE_FAN;
    if (strcasecmp(buf + 1, "RIANGLE-STRIP") == 0)
      return PDB::GFX_TRIANGLE_STRIP;
    if (strcasecmp(buf + 1, "RIANGLES") == 0)
      return PDB::GFX_TRIANGLES;
    break;
  }
  return PDB::GFX_UNKNOWN;
}

const char *RecordClassifier::gfxChars(PDB::GfxType gt) {
  switch (gt) {
  default:
    return "UNKNOWN";
  case PDB::GFX_POINTS:
    return "POINTS";
  case PDB::GFX_MARKERS:
    return "MARKERS";
  case PDB::GFX_LINES:
    return "LINES";
  case PDB::GFX_LINE_STRIP:
    return "LINE-STRIP";
  case PDB::GFX_LINE_LOOP:
    return "LINE-LOOP";
  case PDB::GFX_TRIANGLES:
    return "TRIANGLES";
  case PDB::GFX_TRIANGLE_STRIP:
    return "TRIANGLE-STRIP";
  case PDB::GFX_TRIANGLE_FAN:
    return "TRIANGLE-FAN";
  case PDB::GFX_QUADS:
    return "QUADS";
  case PDB::GFX_QUAD_STRIP:
    return "QUAD-STRIP";
  case PDB::GFX_POLYGON:
    return "POLYGON";
  }
}

PDB::RecordType RecordClassifier::pdbrun5Type(const char *buf) {
  switch (buf[0]) {
  case 'A':
  case 'a':
    if (strncasecmp(buf + 1, "NGLE ", 5) == 0)
      return PDB::USER_ANGLE;
    if (strncasecmp(buf + 1, "TPOS ", 5) == 0)
      return PDB::USER_ATPOS;
    break;
  case 'B':
  case 'b':
    if (strncasecmp(buf + 1, "GCOLOR ", 7) == 0)
      return PDB::USER_BGCOLOR;
    break;
  case 'C':
  case 'c':
    if (strncasecmp(buf + 1, "HAIN ", 5) == 0)
      return PDB::USER_CHAIN;
    if (strncasecmp(buf + 1, "NAME ", 5) == 0)
      return PDB::USER_CNAME;
    if (strncasecmp(buf + 1, "OLOR ", 5) == 0)
      return PDB::USER_COLOR;
    break;
  case 'D':
  case 'd':
    if (strncasecmp(buf + 1, "ISTANCE ", 8) == 0)
      return PDB::USER_DISTANCE;
    break;
  case 'E':
  case 'e':
    if (strncasecmp(buf + 1, "NDOBJ ", 6) == 0)
      return PDB::USER_ENDOBJ;
    if (strncasecmp(buf + 1, "YEPOS ", 6) == 0)
      return PDB::USER_EYEPOS;
    break;
  case 'F':
  case 'f':
    if (strncasecmp(buf + 1, "ILE ", 4) == 0)
      return PDB::USER_FILE;
    if (strncasecmp(buf + 1, "OCUS ", 5) == 0)
      return PDB::USER_FOCUS;
    break;
  case 'G':
  case 'g':
    if (strncasecmp(buf + 1, "FX ", 3) != 0)
      break;
    if (strncasecmp(buf + 4, "COLOR ", 6) == 0)
      return PDB::USER_GFX_COLOR;
    if (strncasecmp(buf + 4, "DRAW ", 5) == 0)
      return PDB::USER_GFX_DRAW;
    if (strncasecmp(buf + 4, "FONT ", 5) == 0)
      return PDB::USER_GFX_FONT;
    if (strncasecmp(buf + 4, "LABEL ", 6) == 0)
      return PDB::USER_GFX_LABEL;
    if (strncasecmp(buf + 4, "MARKER ", 7) == 0)
      return PDB::USER_GFX_MARKER;
    if (strncasecmp(buf + 4, "MOVE ", 5) == 0)
      return PDB::USER_GFX_MOVE;
    if (strncasecmp(buf + 4, "POINT ", 6) == 0)
      return PDB::USER_GFX_POINT;
    break;
  case 'O':
  case 'o':
    if (strncasecmp(buf + 1, "BJECT ", 6) == 0)
      return PDB::USER_OBJECT;
    break;
  case 'P':
  case 'p':
    if (strncasecmp(buf + 1, "DBRUN ", 6) == 0)
      return PDB::USER_PDBRUN;
    break;
  case 'R':
  case 'r':
    if (strncasecmp(buf + 1, "ADIUS ", 6) == 0)
      return PDB::USER_RADIUS;
    break;
  case 'V':
  case 'v':
    if (strncasecmp(buf + 1, "IEWPORT ", 8) == 0)
      return PDB::USER_VIEWPORT;
    break;
  case 'W':
  case 'w':
    if (strncasecmp(buf + 1, "INDOW ", 6) == 0)
      return PDB::USER_WINDOW;
    break;
  }
  return PDB::USER;
}

PDB::RecordType RecordClassifier::pdbrun6Type(const char *buf) {
  switch (buf[0]) {
  case 'A':
  case 'a':
    if (strncasecmp(buf + 1, "NGLE ", 5) == 0)
      return PDB::USER_ANGLE;
    if (strncasecmp(buf + 1, "TPOS ", 5) == 0)
      return PDB::USER_ATPOS;
    break;
  case 'B':
  case 'b':
    if (strncasecmp(buf + 1, "GCOLOR ", 7) == 0)
      return PDB::USER_BGCOLOR;
    break;
  case 'C':
  case 'c':
    if (strncasecmp(buf + 1, "HAIN ", 5) == 0)
      return PDB::USER_CHAIN;
    if (strncasecmp(buf + 1, "NAME ", 5) == 0)
      return PDB::USER_CNAME;
    if (strncasecmp(buf + 1, "OLOR ", 5) == 0)
      return PDB::USER_COLOR;
    break;
  case 'D':
  case 'd':
    if (strncasecmp(buf + 1, "ISTANCE ", 8) == 0)
      return PDB::USER_DISTANCE;
    break;
  case 'E':
  case 'e':
    if (strncasecmp(buf + 1, "NDOBJ", 5) == 0 &&
        (buf[6] == '\0' || buf[6] == '\n' || buf[6] == ' '))
      return PDB::USER_ENDOBJ;
    if (strncasecmp(buf + 1, "YEPOS ", 6) == 0)
      return PDB::USER_EYEPOS;
    break;
  case 'F':
  case 'f':
    if (strncasecmp(buf + 1, "ILE ", 4) == 0)
      return PDB::USER_FILE;
    if (strncasecmp(buf + 1, "OCUS ", 5) == 0)
      return PDB::USER_FOCUS;
    break;
  case 'G':
  case 'g':
    if (buf[1] != 'F' || buf[2] != 'X' || buf[3] != ' ')
      break;
    switch (buf[4]) {
    case 'B':
    case 'b':
      if (strncasecmp(buf + 5, "EGIN ", 5) == 0)
        return PDB::USER_GFX_BEGIN;
      break;
    case 'C':
    case 'c':
      if (strncasecmp(buf + 5, "OLOR ", 5) == 0)
        return PDB::USER_GFX_COLOR;
      break;
    case 'E':
    case 'e':
      if (buf[5] == 'N' && buf[6] == 'D' &&
          (buf[7] == '\0' || buf[7] == '\n' || buf[7] == ' '))
        return PDB::USER_GFX_END;
      break;
    case 'F':
    case 'f':
      if (strncasecmp(buf + 5, "ONT ", 4) == 0)
        return PDB::USER_GFX_FONT;
      break;
    case 'L':
    case 'l':
      if (strncasecmp(buf + 5, "ABEL ", 5) == 0)
        return PDB::USER_GFX_LABEL;
      break;
    case 'N':
    case 'n':
      if (strncasecmp(buf + 5, "ORMAL ", 6) == 0)
        return PDB::USER_GFX_NORMAL;
      break;
    case 'T':
    case 't':
      if (strncasecmp(buf + 5, "EXTPOS ", 7) == 0)
        return PDB::USER_GFX_TEXTPOS;
      break;
    case 'V':
    case 'v':
      if (strncasecmp(buf + 5, "ERTEX ", 6) == 0)
        return PDB::USER_GFX_VERTEX;
      break;
    }
    break;
  case 'M':
  case 'm':
    if (strncasecmp(buf + 1, "ARK ", 4) == 0)
      return PDB::USER_MARK;
    if (strncasecmp(buf + 1, "ARKNAME ", 6) == 0)
      return PDB::USER_MARKNAME;
    break;
  case 'O':
  case 'o':
    if (strncasecmp(buf + 1, "BJECT", 5) == 0 &&
        (buf[6] == '\0' || buf[6] == '\n' || buf[6] == ' '))
      return PDB::USER_OBJECT;
    break;
  case 'P':
  case 'p':
    if (strncasecmp(buf + 1, "DBRUN ", 6) == 0)
      return PDB::USER_PDBRUN;
    break;
  case 'R':
  case 'r':
    if (strncasecmp(buf + 1, "ADIUS ", 6) == 0)
      return PDB::USER_RADIUS;
    break;
  case 'V':
  case 'v':
    if (strncasecmp(buf + 1, "IEWPORT ", 6) == 0)
      return PDB::USER_VIEWPORT;
    break;
  case 'W':
  case 'w':
    if (strncasecmp(buf + 1, "INDOW ", 6) == 0)
      return PDB::USER_WINDOW;
    break;
  }
  return PDB::USER;
}

PDB::RecordType RecordClassifier::getType(const char *buf,
                                          int pdbrunInputVersion) {
  char rt[4]; // PDB record type
  int i;

  for (i = 0; buf[i] != '\0' && buf[i] != '\n' && i < 4; i += 1) {
    if (islower(static_cast<unsigned char>(buf[i])))
      rt[i] = _toupper(buf[i]);
    else
      rt[i] = buf[i];
  }
  if (i < 4)
    for (; i < 4; i += 1)
      rt[i] = ' ';

  switch (rt[0]) {

  case 'A':
    switch (rt[1]) {
    case 'G':
      if (rt[2] == 'R' && rt[3] == 'D')
        return PDB::AGRDES;
      if (rt[2] == 'G' && rt[3] == 'R')
        return PDB::AGGRGT;
      break;
    case 'N':
      if (rt[2] == 'I' && rt[3] == 'S')
        return PDB::ANISOU;
      break;
    case 'T':
      if (rt[2] == 'O' && rt[3] == 'M')
        return PDB::ATOM;
      break;
    case 'U':
      if (rt[2] == 'T' && rt[3] == 'H')
        return PDB::AUTHOR;
      break;
    }
    break;

  case 'C':
    switch (rt[1]) {
    case 'A':
      if (rt[2] == 'V' && rt[3] == 'E')
        return PDB::CAVEAT;
      break;
    case 'M':
      if (rt[2] == 'P' && rt[3] == 'D')
        return PDB::CMPDES;
      if (rt[2] == 'P' && rt[3] == 'O')
        return PDB::CMPONT;
      break;
    case 'O':
      if (rt[2] == 'M' && rt[3] == 'P')
        return PDB::COMPND;
      if (rt[2] == 'N' && rt[3] == 'E')
        return PDB::CONECT;
      break;
    case 'R':
      if (rt[2] == 'Y' && rt[3] == 'S')
        return PDB::CRYST1;
      break;
    }
    break;

  case 'E':
    switch (rt[1]) {
    case 'N':
      if (rt[2] == 'D' && rt[3] == ' ')
        return PDB::END;
      if (rt[2] == 'D' && rt[3] == 'M')
        return PDB::ENDMDL;
      break;
    case 'X':
      if (rt[2] == 'P' && rt[3] == 'D')
        return PDB::EXPDTA;
      break;
    }
    break;

  case 'F':
    switch (rt[1]) {
    case 'T':
      if (rt[2] == 'N' && rt[3] == 'O')
        return PDB::FTNOTE;
      break;
    case 'O':
      if (rt[2] == 'R' && rt[3] == 'M')
        return PDB::FORMUL;
      break;
    }
    break;

  case 'H':
    if (rt[1] != 'E')
      break;
    if (rt[2] == 'T' && rt[3] == 'A')
      return PDB::HETATM;
    if (rt[2] == 'A' && rt[3] == 'D')
      return PDB::HEADER;
    if (rt[2] == 'T' && rt[3] == ' ')
      return PDB::HET;
    if (rt[2] == 'L' && rt[3] == 'I')
      return PDB::HELIX;
    break;

  case 'J':
    if (rt[1] == 'R' && rt[2] == 'N' && rt[3] == 'L')
      return PDB::JRNL;
    break;

  case 'K':
    if (rt[1] == 'E' && rt[2] == 'Y' && rt[3] == 'W')
      return PDB::KEYWDS;
    break;

  case 'M':
    switch (rt[1]) {
    case 'A':
      if (rt[2] == 'S' && rt[3] == 'T')
        return PDB::MASTER;
      break;
    case 'D':
      if (rt[2] == 'L' && rt[3] == 'T')
        return PDB::MDLTYP;
      break;
    case 'O':
      if (rt[2] == 'D' && rt[3] == 'E')
        return PDB::MODEL;
      break;
    case 'T':
      if (rt[2] == 'R' && rt[3] == 'I')
        return PDB::MTRIX;
      if (rt[2] == 'X' && rt[3] == 'D')
        return PDB::MTXDES;
      break;
    }
    break;

  case 'N':
    if (rt[1] == 'U' && rt[2] == 'M' && rt[3] == 'M')
      return PDB::NUMMDL;
    break;

  case 'O':
    switch (rt[1]) {
    case 'B':
      if (rt[2] == 'S' && rt[3] == 'L')
        return PDB::OBSLTE;
      break;
    case 'R':
      if (rt[2] == 'I' && rt[3] == 'G')
        return PDB::ORIGX;
      break;
    }
    break;

  case 'R':
    if (rt[1] != 'E')
      break;
    if (rt[2] == 'M' && rt[3] == 'A')
      return PDB::REMARK;
    if (rt[2] == 'V' && rt[3] == 'D')
      return PDB::REVDAT;
    break;

  case 'S':
    switch (rt[1]) {

    case 'C':
      if (rt[2] == 'A' && rt[3] == 'L')
        return PDB::SCALE;
      break;

    case 'E':
      if (rt[2] == 'Q' && rt[3] == 'R')
        return PDB::SEQRES;
      break;

    case 'H':
      if (rt[2] == 'E' && rt[3] == 'E')
        return PDB::SHEET;
      break;

    case 'I':
      if (rt[2] == 'T' && rt[3] == 'E')
        return PDB::SITE;
      if (rt[2] == 'G' && rt[3] == 'A')
        return PDB::SIGATM;
      if (rt[2] == 'G' && rt[3] == 'U')
        return PDB::SIGUIJ;
      break;

    case 'O':
      if (rt[2] == 'U' && rt[3] == 'R')
        return PDB::SOURCE;
      break;

    case 'P':
      if (rt[2] == 'L' && rt[3] == 'I')
        return PDB::SPLIT;
      if (rt[2] == 'R' && rt[3] == 'S')
        return PDB::SPRSDE;
      break;

    case 'S':
      if (rt[2] == 'B' && rt[3] == 'O')
        return PDB::SSBOND;
      break;

    case 'Y':
      if (rt[2] == 'M' && rt[3] == 'D')
        return PDB::SYMDES;
      if (rt[2] == 'M' && rt[3] == 'O')
        return PDB::SYMOP;
      break;
    }
    break;

  case 'T':
    switch (rt[1]) {
    case 'E':
      if (rt[2] == 'R' && rt[3] == ' ')
        return PDB::TER;
      break;
    case 'I':
      if (rt[2] == 'T' && rt[3] == 'L')
        return PDB::TITLE;
      break;
    case 'R':
      if (rt[2] == 'N' && rt[3] == 'S')
        return PDB::TRNSFM;
      break;
    case 'U':
      if (rt[2] == 'R' && rt[3] == 'N')
        return PDB::TURN;
      break;
    case 'V':
      if (rt[2] == 'E' && rt[3] == 'C')
        return PDB::TVECT;
      break;
    }
    break;

  case 'U':
    if (rt[1] == 'S' && rt[2] == 'E' && rt[3] == 'R')
      switch (pdbrunInputVersion) {
      case 1:
      case 2:
      case 3:
      case 4:
      case 5:
        return pdbrun5Type(buf + 6);
      case 6:
        return pdbrun6Type(buf + 6);
      default:
        if (strncasecmp(buf + 6, "PDBRUN ", 7) == 0)
          return PDB::USER_PDBRUN;
        return PDB::USER;
      }
    break;
  }
  return PDB::UNKNOWN;
}

} // namespace libpdb
