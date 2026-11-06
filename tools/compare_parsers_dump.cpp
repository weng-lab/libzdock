// Copyright (c) 2019-2026 Arjan van der Velde, Weng Lab
// SPDX-License-Identifier: BSD-2-Clause

// Baseline record dump: parses one PDB file with the pre-rewrite parser
// and prints one line per record. Compiled out-of-tree against a bb0032c
// worktree (the in-tree legacy parser is retired once the snapshot
// exists); see tools/parity/README.md for the build command, which
// force-includes a small generated meta header carrying the commit,
// compiler flags, and harness version that --meta reports.
//
// Modes:
//   baseline_dump <file>            one "<type>\t<text>" line per record
//   baseline_dump --json <file>     one JSON line for the whole file:
//                                   {"records":[{"line":N,"type":T,
//                                   "text":"...","fields":{...}},...]}
//   baseline_dump --meta            one JSON line of build metadata
//
// The JSON field names and value formats mirror the parity harness's
// field comparison exactly, so a snapshot compare against these records
// reproduces the live legacy-vs-modern comparison field-for-field.

// built only against the bb0032c worktree (see README.md);
// the in-tree legacy header is retired
#include "pdb++.h"

#include <cstdio>
#include <cstring>
#include <iostream>
#include <sstream>
#include <string>

#ifndef HARNESS_VERSION
#define HARNESS_VERSION "unknown"
#endif
#ifndef HARNESS_COMMIT
#define HARNESS_COMMIT "unknown"
#endif
#ifndef HARNESS_FLAGS
#define HARNESS_FLAGS "unknown"
#endif

namespace {

std::string jsonEscape(const std::string &s) {
  std::string out;
  out.reserve(s.size() + 8);
  for (unsigned char c : s) {
    switch (c) {
    case '"':
      out += "\\\"";
      break;
    case '\\':
      out += "\\\\";
      break;
    case '\b':
      out += "\\b";
      break;
    case '\f':
      out += "\\f";
      break;
    case '\n':
      out += "\\n";
      break;
    case '\r':
      out += "\\r";
      break;
    case '\t':
      out += "\\t";
      break;
    default:
      if (c < 0x20 || c >= 0x80) {
        // control bytes and non-ASCII bytes escape as \u00XX so the
        // JSON stays valid UTF-8 and round-trips byte-exactly
        char buf[8];
        std::snprintf(buf, sizeof(buf), "\\u%04x", c);
        out += buf;
      } else {
        out += static_cast<char>(c);
      }
    }
  }
  return out;
}

std::string jsonStr(const char *s) {
  return "\"" + jsonEscape(s == nullptr ? "" : s) + "\"";
}

std::string jsonStr(const std::string &s) { return "\"" + jsonEscape(s) + "\""; }

std::string jsonInt(long v) {
  std::ostringstream os;
  os << v;
  return os.str();
}

std::string jsonReal(double v) {
  char buf[40];
  std::snprintf(buf, sizeof(buf), "%.17g", v);
  return buf;
}

std::string jsonChar(char c) {
  const std::string s(1, c);
  return jsonStr(s);
}

std::string jsonCharNull(char c) { return jsonChar(c == '\0' ? ' ' : c); }

// The parsed-field emission for one record: the same field names the
// parity harness compares, in the same value formats.
std::string emitFields(const libpdb::PDB &r) {
  std::ostringstream os;
  os << "{";
  const bool atomish = r.type() == libpdb::PDB::ATOM ||
                       r.type() == libpdb::PDB::HETATM ||
                       r.type() == libpdb::PDB::SIGATM;
  if (atomish && r.atom.iszdatom) {
    os << "\"zdatom.type\":" << jsonInt(r.atom.type) << ","
       << "\"zdatom.surface\":" << jsonInt(r.atom.surface) << ","
       << "\"zdatom.rad\":" << jsonReal(r.atom.rad) << ","
       << "\"zdatom.segid\":" << jsonStr(r.atom.segid) << ","
       << "\"zdatom.chg\":" << jsonReal(r.atom.chg);
  } else if (atomish) {
    os << "\"serialNum\":" << jsonInt(r.atom.serialNum) << ","
       << "\"name\":" << jsonStr(r.atom.name) << ","
       << "\"altLoc\":" << jsonChar(r.atom.altLoc) << ","
       << "\"residue.name\":" << jsonStr(r.atom.residue.name) << ","
       << "\"residue.chainId\":" << jsonChar(r.atom.residue.chainId) << ","
       << "\"residue.seqNum\":" << jsonInt(r.atom.residue.seqNum) << ","
       << "\"residue.insertCode\":" << jsonChar(r.atom.residue.insertCode);
    for (int i = 0; i < 3; ++i)
      os << ",\"xyz[" << i << "]\":" << jsonReal(r.atom.xyz[i]);
    os << ",\"occupancy\":" << jsonReal(r.atom.occupancy)
       << ",\"tempFactor\":" << jsonReal(r.atom.tempFactor)
       << ",\"element\":" << jsonStr(r.atom.element)
       << ",\"charge\":" << jsonStr(r.atom.charge);
  } else
    switch (r.type()) {
    case libpdb::PDB::ANISOU:
    case libpdb::PDB::SIGUIJ: {
      os << "\"serialNum\":" << jsonInt(r.anisou.serialNum)
         << ",\"name\":" << jsonStr(r.anisou.name)
         << ",\"altLoc\":" << jsonChar(r.anisou.altLoc)
         << ",\"residue.name\":" << jsonStr(r.anisou.residue.name)
         << ",\"residue.chainId\":" << jsonChar(r.anisou.residue.chainId)
         << ",\"residue.seqNum\":" << jsonInt(r.anisou.residue.seqNum)
         << ",\"residue.insertCode\":" << jsonChar(r.anisou.residue.insertCode);
      for (int i = 0; i < 6; ++i)
        os << ",\"u[" << i << "]\":" << jsonInt(r.anisou.u[i]);
      break;
    }
    case libpdb::PDB::AUTHOR:
      os << "\"continuation\":" << jsonInt(r.author.continuation)
         << ",\"authors\":" << jsonStr(r.author.authors);
      break;
    case libpdb::PDB::CAVEAT:
      os << "\"continuation\":" << jsonInt(r.caveat.continuation)
         << ",\"id\":" << jsonStr(r.caveat.id)
         << ",\"text\":" << jsonStr(r.caveat.text);
      break;
    case libpdb::PDB::COMPND:
      os << "\"continuation\":" << jsonInt(r.compnd.continuation)
         << ",\"compound\":" << jsonStr(r.compnd.compound);
      break;
    case libpdb::PDB::CONECT: {
      os << "\"serialNum\":" << jsonInt(r.conect.serialNum);
      for (int i = 0; i < 4; ++i)
        os << ",\"covalent[" << i << "]\":" << jsonInt(r.conect.covalent[i]);
      for (int b = 0; b < 2; ++b) {
        for (int i = 0; i < 2; ++i)
          os << ",\"bonds[" << b << "].hydrogen[" << i
             << "]\":" << jsonInt(r.conect.bonds[b].hydrogen[i]);
        os << ",\"bonds[" << b << "].salt\":" << jsonInt(r.conect.bonds[b].salt);
      }
      break;
    }
    case libpdb::PDB::CRYST1:
      os << "\"a\":" << jsonReal(r.cryst1.a) << ",\"b\":" << jsonReal(r.cryst1.b)
         << ",\"c\":" << jsonReal(r.cryst1.c)
         << ",\"alpha\":" << jsonReal(r.cryst1.alpha)
         << ",\"beta\":" << jsonReal(r.cryst1.beta)
         << ",\"gamma\":" << jsonReal(r.cryst1.gamma)
         << ",\"spaceGrp\":" << jsonStr(r.cryst1.spaceGrp)
         << ",\"z\":" << jsonInt(r.cryst1.z);
      break;
    case libpdb::PDB::FORMUL:
      os << "\"component\":" << jsonInt(r.formul.component)
         << ",\"hetId\":" << jsonStr(r.formul.hetId)
         << ",\"continuation\":" << jsonInt(r.formul.continuation)
         << ",\"exclude\":" << jsonCharNull(r.formul.exclude)
         << ",\"formula\":" << jsonStr(r.formul.formula);
      break;
    case libpdb::PDB::FTNOTE:
    case libpdb::PDB::REMARK:
    case libpdb::PDB::SYMDES:
    case libpdb::PDB::MTXDES:
    case libpdb::PDB::CMPDES:
    case libpdb::PDB::AGRDES:
      os << "\"num\":" << jsonInt(r.agrdes.num)
         << ",\"text\":" << jsonStr(r.agrdes.text);
      break;
    case libpdb::PDB::HEADER:
      os << "\"classification\":" << jsonStr(r.header.classification)
         << ",\"timestamp\":" << jsonStr(r.header.timestamp)
         << ",\"id\":" << jsonStr(r.header.id);
      break;
    case libpdb::PDB::HELIX:
      os << "\"serialNum\":" << jsonInt(r.helix.serialNum)
         << ",\"id\":" << jsonStr(r.helix.id);
      for (int i = 0; i < 2; ++i)
        os << ",\"residues[" << i << "].name\":" << jsonStr(r.helix.residues[i].name)
           << ",\"residues[" << i << "].chainId\":" << jsonChar(r.helix.residues[i].chainId)
           << ",\"residues[" << i << "].seqNum\":" << jsonInt(r.helix.residues[i].seqNum)
           << ",\"residues[" << i << "].insertCode\":" << jsonChar(r.helix.residues[i].insertCode);
      os << ",\"type\":" << jsonInt(r.helix.type)
         << ",\"comment\":" << jsonStr(r.helix.comment)
         << ",\"length\":" << jsonInt(r.helix.length);
      break;
    case libpdb::PDB::HET:
      os << "\"hetGrp.name\":" << jsonStr(r.het.hetGrp.name)
         << ",\"hetGrp.chainId\":" << jsonChar(r.het.hetGrp.chainId)
         << ",\"hetGrp.seqNum\":" << jsonInt(r.het.hetGrp.seqNum)
         << ",\"hetGrp.insertCode\":" << jsonChar(r.het.hetGrp.insertCode)
         << ",\"numAtoms\":" << jsonInt(r.het.numAtoms)
         << ",\"text\":" << jsonStr(r.het.text);
      break;
    case libpdb::PDB::JRNL:
      os << "\"text\":" << jsonStr(r.jrnl.text);
      break;
    case libpdb::PDB::KEYWDS:
      os << "\"continuation\":" << jsonInt(r.keywds.continuation)
         << ",\"keywds\":" << jsonStr(r.keywds.keywds);
      break;
    case libpdb::PDB::MASTER: {
      os << "\"numRemark\":" << jsonInt(r.master.numRemark)
         << ",\"numFtnote\":" << jsonInt(r.master.numFtnote)
         << ",\"numHet\":" << jsonInt(r.master.numHet)
         << ",\"numHelix\":" << jsonInt(r.master.numHelix)
         << ",\"numSheet\":" << jsonInt(r.master.numSheet)
         << ",\"numTurn\":" << jsonInt(r.master.numTurn)
         << ",\"numSite\":" << jsonInt(r.master.numSite)
         << ",\"numTransform\":" << jsonInt(r.master.numTransform)
         << ",\"numCoordinate\":" << jsonInt(r.master.numCoordinate)
         << ",\"numTer\":" << jsonInt(r.master.numTer)
         << ",\"numConect\":" << jsonInt(r.master.numConect)
         << ",\"numSeqres\":" << jsonInt(r.master.numSeqres);
      break;
    }
    case libpdb::PDB::MDLTYP:
      os << "\"continuation\":" << jsonInt(r.mdltyp.continuation)
         << ",\"comment\":" << jsonStr(r.mdltyp.comment);
      break;
    case libpdb::PDB::MODEL:
      os << "\"num\":" << jsonInt(r.model.num);
      break;
    case libpdb::PDB::NUMMDL:
      os << "\"number\":" << jsonInt(r.nummdl.number);
      break;
    case libpdb::PDB::MTRIX: {
      os << "\"rowNum\":" << jsonInt(r.mtrix.rowNum)
         << ",\"serialNum\":" << jsonInt(r.mtrix.serialNum)
         << ",\"m1\":" << jsonReal(r.mtrix.m1)
         << ",\"m2\":" << jsonReal(r.mtrix.m2)
         << ",\"m3\":" << jsonReal(r.mtrix.m3)
         << ",\"v\":" << jsonReal(r.mtrix.v)
         << ",\"given\":" << jsonInt(r.mtrix.given);
      break;
    }
    case libpdb::PDB::OBSLTE: {
      os << "\"continuation\":" << jsonInt(r.obslte.continuation)
         << ",\"timestamp\":" << jsonStr(r.obslte.timestamp)
         << ",\"oldId\":" << jsonStr(r.obslte.oldId);
      for (int i = 0; i < 9; ++i)
        os << ",\"idMap[" << i << "]\":" << jsonStr(r.obslte.idMap[i]);
      break;
    }
    case libpdb::PDB::ORIGX:
      os << "\"rowNum\":" << jsonInt(r.origx.rowNum)
         << ",\"o1\":" << jsonReal(r.origx.o1)
         << ",\"o2\":" << jsonReal(r.origx.o2)
         << ",\"o3\":" << jsonReal(r.origx.o3)
         << ",\"t\":" << jsonReal(r.origx.t);
      break;
    case libpdb::PDB::REVDAT: {
      os << "\"modification\":" << jsonInt(r.revdat.modification)
         << ",\"continuation\":" << jsonInt(r.revdat.continuation)
         << ",\"timestamp\":" << jsonStr(r.revdat.timestamp)
         << ",\"id\":" << jsonStr(r.revdat.id)
         << ",\"modType\":" << jsonInt(r.revdat.modType);
      for (int i = 0; i < 4; ++i)
        os << ",\"record[" << i << "]\":" << jsonStr(r.revdat.record[i]);
      break;
    }
    case libpdb::PDB::SCALE:
      os << "\"rowNum\":" << jsonInt(r.scale.rowNum)
         << ",\"s1\":" << jsonReal(r.scale.s1)
         << ",\"s2\":" << jsonReal(r.scale.s2)
         << ",\"s3\":" << jsonReal(r.scale.s3)
         << ",\"u\":" << jsonReal(r.scale.u);
      break;
    case libpdb::PDB::SEQRES: {
      os << "\"serialNum\":" << jsonInt(r.seqres.serialNum)
         << ",\"chainId\":" << jsonChar(r.seqres.chainId)
         << ",\"count\":" << jsonInt(r.seqres.count);
      for (int i = 0; i < 13; ++i)
        os << ",\"names[" << i << "]\":" << jsonStr(r.seqres.names[i]);
      break;
    }
    case libpdb::PDB::SHEET:
      os << "\"strandNum\":" << jsonInt(r.sheet.strandNum)
         << ",\"id\":" << jsonStr(r.sheet.id)
         << ",\"count\":" << jsonInt(r.sheet.count);
      for (int i = 0; i < 2; ++i)
        os << ",\"residues[" << i << "].name\":" << jsonStr(r.sheet.residues[i].name)
           << ",\"residues[" << i << "].chainId\":" << jsonChar(r.sheet.residues[i].chainId)
           << ",\"residues[" << i << "].seqNum\":" << jsonInt(r.sheet.residues[i].seqNum)
           << ",\"residues[" << i << "].insertCode\":" << jsonChar(r.sheet.residues[i].insertCode);
      os << ",\"sense\":" << jsonInt(r.sheet.sense);
      for (int i = 0; i < 2; ++i) {
        os << ",\"atoms[" << i << "].name\":" << jsonStr(r.sheet.atoms[i].name)
           << ",\"atoms[" << i << "].residue.name\":" << jsonStr(r.sheet.atoms[i].residue.name)
           << ",\"atoms[" << i << "].residue.chainId\":" << jsonChar(r.sheet.atoms[i].residue.chainId)
           << ",\"atoms[" << i << "].residue.seqNum\":" << jsonInt(r.sheet.atoms[i].residue.seqNum)
           << ",\"atoms[" << i << "].residue.insertCode\":" << jsonChar(r.sheet.atoms[i].residue.insertCode);
      }
      break;
    case libpdb::PDB::SITE:
      os << "\"seqNum\":" << jsonInt(r.site.seqNum)
         << ",\"id\":" << jsonStr(r.site.id)
         << ",\"count\":" << jsonInt(r.site.count);
      for (int i = 0; i < 4; ++i)
        os << ",\"residues[" << i << "].name\":" << jsonStr(r.site.residues[i].name)
           << ",\"residues[" << i << "].chainId\":" << jsonChar(r.site.residues[i].chainId)
           << ",\"residues[" << i << "].seqNum\":" << jsonInt(r.site.residues[i].seqNum)
           << ",\"residues[" << i << "].insertCode\":" << jsonChar(r.site.residues[i].insertCode);
      break;
    case libpdb::PDB::SOURCE:
      os << "\"continuation\":" << jsonInt(r.source.continuation)
         << ",\"name\":" << jsonStr(r.source.name);
      break;
    case libpdb::PDB::SPLIT: {
      os << "\"continuation\":" << jsonInt(r.split.continuation);
      for (int i = 0; i < 14; ++i)
        os << ",\"idMap[" << i << "]\":" << jsonStr(r.split.idMap[i]);
      break;
    }
    case libpdb::PDB::SPRSDE: {
      os << "\"continuation\":" << jsonInt(r.sprsde.continuation)
         << ",\"timestamp\":" << jsonStr(r.sprsde.timestamp)
         << ",\"id\":" << jsonStr(r.sprsde.id);
      for (int i = 0; i < 9; ++i)
        os << ",\"supersede[" << i << "]\":" << jsonStr(r.sprsde.supersede[i]);
      break;
    }
    case libpdb::PDB::SSBOND:
      os << "\"seqNum\":" << jsonInt(r.ssbond.seqNum);
      for (int i = 0; i < 2; ++i)
        os << ",\"residues[" << i << "].name\":" << jsonStr(r.ssbond.residues[i].name)
           << ",\"residues[" << i << "].chainId\":" << jsonChar(r.ssbond.residues[i].chainId)
           << ",\"residues[" << i << "].seqNum\":" << jsonInt(r.ssbond.residues[i].seqNum)
           << ",\"residues[" << i << "].insertCode\":" << jsonChar(r.ssbond.residues[i].insertCode);
      os << ",\"sym1\":" << jsonInt(r.ssbond.sym1)
         << ",\"sym2\":" << jsonInt(r.ssbond.sym2)
         << ",\"length\":" << jsonReal(r.ssbond.length);
      break;
    case libpdb::PDB::SYMOP:
      os << "\"rowNum\":" << jsonInt(r.symop.rowNum)
         << ",\"serialNum\":" << jsonInt(r.symop.serialNum)
         << ",\"s1\":" << jsonReal(r.symop.s1)
         << ",\"s2\":" << jsonReal(r.symop.s2)
         << ",\"s3\":" << jsonReal(r.symop.s3)
         << ",\"t\":" << jsonReal(r.symop.t);
      break;
    case libpdb::PDB::TER:
      os << "\"serialNum\":" << jsonInt(r.ter.serialNum)
         << ",\"residue.name\":" << jsonStr(r.ter.residue.name)
         << ",\"residue.chainId\":" << jsonChar(r.ter.residue.chainId)
         << ",\"residue.seqNum\":" << jsonInt(r.ter.residue.seqNum)
         << ",\"residue.insertCode\":" << jsonChar(r.ter.residue.insertCode);
      break;
    case libpdb::PDB::TITLE:
      os << "\"continuation\":" << jsonInt(r.title.continuation)
         << ",\"text\":" << jsonStr(r.title.text);
      break;
    case libpdb::PDB::TRNSFM:
      os << "\"resultSerialNum\":" << jsonInt(r.trnsfm.resultSerialNum)
         << ",\"applySerialNum\":" << jsonInt(r.trnsfm.applySerialNum)
         << ",\"sourceSerialNum\":" << jsonInt(r.trnsfm.sourceSerialNum);
      break;
    case libpdb::PDB::TURN:
      os << "\"seqNum\":" << jsonInt(r.turn.seqNum)
         << ",\"id\":" << jsonStr(r.turn.id);
      for (int i = 0; i < 2; ++i)
        os << ",\"residues[" << i << "].name\":" << jsonStr(r.turn.residues[i].name)
           << ",\"residues[" << i << "].chainId\":" << jsonChar(r.turn.residues[i].chainId)
           << ",\"residues[" << i << "].seqNum\":" << jsonInt(r.turn.residues[i].seqNum)
           << ",\"residues[" << i << "].insertCode\":" << jsonChar(r.turn.residues[i].insertCode);
      os << ",\"comment\":" << jsonStr(r.turn.comment);
      break;
    case libpdb::PDB::TVECT:
      os << "\"serialNum\":" << jsonInt(r.tvect.serialNum)
         << ",\"t1\":" << jsonReal(r.tvect.t1)
         << ",\"t2\":" << jsonReal(r.tvect.t2)
         << ",\"t3\":" << jsonReal(r.tvect.t3)
         << ",\"comment\":" << jsonStr(r.tvect.comment);
      break;
    case libpdb::PDB::USER:
      os << "\"subtype\":" << jsonStr(r.user.subtype)
         << ",\"text\":" << jsonStr(r.user.text);
      break;
    case libpdb::PDB::USER_PDBRUN:
      os << "\"version\":" << jsonInt(r.userPdbrun.version);
      break;
    case libpdb::PDB::USER_EYEPOS:
    case libpdb::PDB::USER_ATPOS:
      for (int i = 0; i < 3; ++i)
        os << (i ? "," : "") << "\"xyz[" << i << "]\":" << jsonReal(r.userEyePos.xyz[i]);
      break;
    case libpdb::PDB::USER_WINDOW:
      os << "\"left\":" << jsonReal(r.userWindow.left)
         << ",\"right\":" << jsonReal(r.userWindow.right)
         << ",\"bottom\":" << jsonReal(r.userWindow.bottom)
         << ",\"top\":" << jsonReal(r.userWindow.top)
         << ",\"hither\":" << jsonReal(r.userWindow.hither)
         << ",\"yon\":" << jsonReal(r.userWindow.yon);
      break;
    case libpdb::PDB::USER_FOCUS:
      os << "\"focus\":" << jsonReal(r.userFocus.focus);
      break;
    case libpdb::PDB::USER_VIEWPORT:
      os << "\"xmin\":" << jsonReal(r.userViewport.xmin)
         << ",\"xmax\":" << jsonReal(r.userViewport.xmax)
         << ",\"ymin\":" << jsonReal(r.userViewport.ymin)
         << ",\"ymax\":" << jsonReal(r.userViewport.ymax);
      break;
    case libpdb::PDB::USER_BGCOLOR:
      for (int i = 0; i < 3; ++i)
        os << (i ? "," : "") << "\"rgb[" << i << "]\":" << jsonReal(r.userBgColor.rgb[i]);
      break;
    case libpdb::PDB::USER_ANGLE:
      os << "\"which\":" << jsonInt(r.userAngle.which)
         << ",\"atom0\":" << jsonInt(r.userAngle.atom0)
         << ",\"atom1\":" << jsonInt(r.userAngle.atom1)
         << ",\"atom2\":" << jsonInt(r.userAngle.atom2)
         << ",\"atom3\":" << jsonInt(r.userAngle.atom3)
         << ",\"angle\":" << jsonReal(r.userAngle.angle);
      break;
    case libpdb::PDB::USER_DISTANCE:
      os << "\"which\":" << jsonInt(r.userDistance.which)
         << ",\"atom0\":" << jsonInt(r.userDistance.atom0)
         << ",\"atom1\":" << jsonInt(r.userDistance.atom1)
         << ",\"distance\":" << jsonReal(r.userDistance.distance);
      break;
    case libpdb::PDB::USER_FILE:
      os << "\"filename\":" << jsonStr(r.userFile.filename)
         << ",\"model\":" << jsonInt(r.userFile.model);
      break;
    case libpdb::PDB::USER_MARKNAME:
    case libpdb::PDB::USER_MARK:
      os << "\"markname\":" << jsonStr(r.userMarkname.markname);
      break;
    case libpdb::PDB::USER_CNAME:
      for (int i = 0; i < 3; ++i)
        os << (i ? "," : "") << "\"rgb[" << i << "]\":" << jsonReal(r.userCName.rgb[i]);
      os << ",\"name\":" << jsonStr(r.userCName.name);
      break;
    case libpdb::PDB::USER_COLOR:
    case libpdb::PDB::USER_GFX_COLOR:
      for (int i = 0; i < 3; ++i)
        os << (i ? "," : "") << "\"rgb[" << i << "]\":" << jsonReal(r.userColor.rgb[i]);
      os << ",\"spec\":" << jsonStr(r.userColor.spec);
      break;
    case libpdb::PDB::USER_RADIUS:
      os << "\"radius\":" << jsonReal(r.userRadius.radius);
      break;
    case libpdb::PDB::USER_OBJECT:
    case libpdb::PDB::USER_ENDOBJ:
      os << "\"model\":" << jsonInt(r.userObject.model);
      break;
    case libpdb::PDB::USER_CHAIN:
      os << "\"atom0\":" << jsonInt(r.userChain.atom0)
         << ",\"atom1\":" << jsonInt(r.userChain.atom1);
      break;
    case libpdb::PDB::USER_GFX_BEGIN:
      os << "\"primitive\":" << jsonInt(r.userGfxBegin.primitive);
      break;
    case libpdb::PDB::USER_GFX_NORMAL:
    case libpdb::PDB::USER_GFX_VERTEX:
    case libpdb::PDB::USER_GFX_TEXTPOS:
    case libpdb::PDB::USER_GFX_MOVE:
    case libpdb::PDB::USER_GFX_DRAW:
    case libpdb::PDB::USER_GFX_MARKER:
    case libpdb::PDB::USER_GFX_POINT:
      for (int i = 0; i < 3; ++i)
        os << (i ? "," : "") << "\"xyz[" << i << "]\":" << jsonReal(r.userGfxNormal.xyz[i]);
      break;
    case libpdb::PDB::USER_GFX_FONT:
      os << "\"size\":" << jsonInt(r.userGfxFont.size)
         << ",\"name\":" << jsonStr(r.userGfxFont.name);
      break;
    case libpdb::PDB::USER_GFX_LABEL:
      for (int i = 0; i < 3; ++i)
        os << (i ? "," : "") << "\"xyz[" << i << "]\":" << jsonReal(r.userGfxLabel.xyz[i]);
      os << ",\"text\":" << jsonStr(r.userGfxLabel.text);
      break;
    default:
      break;
    }
  os << "}";
  return os.str();
}

} // namespace

int main(int argc, char **argv) {
  bool json = false;
  const char *path = nullptr;
  for (int i = 1; i < argc; ++i) {
    const std::string arg(argv[i]);
    if (arg == "--json")
      json = true;
    else if (arg == "--meta") {
      std::printf("{\"commit\":\"%s\",\"compiler\":\"%s\",\"flags\":\"%s\","
                  "\"harness\":\"%s\"}\n",
                  HARNESS_COMMIT, __VERSION__, HARNESS_FLAGS, HARNESS_VERSION);
      return 0;
    } else if (!path)
      path = argv[i];
    else {
      std::fprintf(stderr, "usage: %s [--json|--meta] <pdb-file>\n", argv[0]);
      return 2;
    }
  }
  if (!path) {
    std::fprintf(stderr, "usage: %s [--json|--meta] <pdb-file>\n", argv[0]);
    return 2;
  }
  FILE *in = std::fopen(path, "r");
  if (!in) {
    std::fprintf(stderr, "cannot open %s\n", path);
    return 2;
  }
  char buf[libpdb::PDB::BufLen * 4];
  if (json) {
    // one JSON line for the whole file: an array of record objects
    std::printf("{\"records\":[");
    bool first = true;
    std::size_t lineNo = 0;
    while (std::fgets(buf, sizeof(buf), in)) {
      ++lineNo;
      std::string line(buf);
      while (!line.empty() && (line.back() == '\n' || line.back() == '\r'))
        line.pop_back();
      libpdb::PDB record(line.c_str());
      if (!first)
        std::printf(",");
      first = false;
      std::printf("{\"line\":%zu,\"raw\":%s,\"type\":%d,\"text\":%s,"
                  "\"fields\":%s}",
                  lineNo, jsonStr(line).c_str(), (int)record.type(),
                  jsonStr(record.chars()).c_str(),
                  emitFields(record).c_str());
    }
    std::printf("]}\n");
  } else {
    while (std::fgets(buf, sizeof(buf), in)) {
      // strip the trailing newline (operator>> semantics)
      std::string line(buf);
      while (!line.empty() && (line.back() == '\n' || line.back() == '\r'))
        line.pop_back();
      libpdb::PDB record(line.c_str());
      std::printf("%d\t%s\n", (int)record.type(), record.chars());
    }
  }
  std::fclose(in);
  return 0;
}
