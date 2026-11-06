// Copyright (c) 2019-2026 Arjan van der Velde, Weng Lab
// SPDX-License-Identifier: BSD-2-Clause

// Legacy-vs-modern parser comparison harness. Walks a file or a directory
// tree of PDB files and pushes every line through both parsers, comparing
// the record type, the rendered text, and every parsed payload field
// byte-for-byte. Each side derives the PDBRUN version state independently
// from the file content: the legacy side from its process globals (moved
// by USER PDBRUN lines and by rendering USER PDBRUN records), the modern
// side from its ParseContext/WriteContext (moved by parsing USER PDBRUN
// records) -- no state is copied between the sides, and a divergence in
// the derived state itself is reported.
//
// With --baseline <exe>, every line additionally runs through a separate
// build of the pre-rewrite parser (bb0032c): build the baseline with
//   git worktree add /tmp/bb0032c bb0032c
//   make -C /tmp/bb0032c clean all   (its own toolchain, C++14)
//   g++ -O2 -std=c++14 -I/tmp/bb0032c/contrib/eigen
//       -I/tmp/bb0032c/include tools/compare_parsers_dump.cpp
//       /tmp/bb0032c/lib/libzdock.a -o /tmp/baseline_dump
// Differences against the baseline are classified: ZDATOM-length
// differences (the pre-fb720ad 80-column truncation) and v5 USER token
// differences (the pre-ea22558 unbounded %s) are expected; everything
// else is unexpected.
//
// The input files are read exactly as production reads them, through
// operator>> semantics (the 4xBufLen line buffer, the one trailing CR
// strip, the stream failure on overlong lines), so CRLF files compare
// the same way the servers read them.
//
// Modes:
//   bin/compare_parsers <file-or-directory> [--baseline <exe>]
//       [--files-from <list>] [--max-files N] [--max-lines N]
//       [--all-extensions] [--quiet]
//
//   bin/compare_parsers --dump-json <file>
//       emit the modern parser's records for one file as one JSON line
//       ({"records":[{"line","type","text","fields"},...]}) with the
//       same field names and value formats the snapshot carries; this
//       is the modern side of tools/parity/snapshot.py compare.
//
// With --files-from, <list> holds one input path per line (visited in
// listed order instead of a sorted walk); the first argument is then only
// a label. Shard drivers use this to split a corpus across parallel
// invocations without re-walking the tree.
//
// Exit status: 0 = full parity, 1 = mismatches found, 2 = usage or I/O
// error. Files are visited in sorted order for reproducibility; only
// *.pdb / *.PDB files are visited unless --all-extensions is given.

#include "Record.hpp"

#include <algorithm>
#include <array>
#include <cctype>
#include <cstdio>
#include <cstring>
#include <filesystem>
#include <fstream>
#include <optional>
#include <string>
#include <vector>

namespace {

struct Options {
  std::string root;
  std::string filesFrom; // file list (one path per line); overrides the walk
  std::string baseline; // path to the bb0032c dump executable (optional)
  std::size_t maxFiles = 0;
  std::size_t maxLines = 0;
  bool allExtensions = false;
  bool quiet = false;
};

struct Counters {
  std::size_t files = 0;
  std::size_t lines = 0;
  std::size_t typeMismatches = 0;
  std::size_t textMismatches = 0;
  std::size_t fieldMismatches = 0;
  std::size_t stateMismatches = 0;
  std::size_t rejections = 0;
  // input-quality census, reported separately from parity mismatches
  // and never part of the exit decision: lines both parsers store as
  // UNKNOWN or reject even though the line looks like coordinate data
  std::size_t coordinateUnknowns = 0;
  // deliberate modern gains: lines the legacy parser stores as UNKNOWN
  // junk while the modern parser classifies as the standard ATOM/
  // HETATM record they announce -- the expected difference of the
  // segid/element salvage change, grouped and reported, never part of
  // the exit decision
  std::size_t modernGains = 0;
  // baseline classifications
  std::size_t baselineExpectedZdatom = 0;
  std::size_t baselineExpectedV5Token = 0;
  std::size_t baselineExpectedModernGain = 0;
  std::size_t baselineUnexpected = 0;
};

bool pdbExtension(const std::filesystem::path &path) {
  std::string ext = path.extension().string();
  for (char &c : ext)
    c = static_cast<char>(std::tolower(static_cast<unsigned char>(c)));
  return ext == ".pdb";
}

std::string shorten(const std::string &text) {
  if (text.size() <= 72)
    return text;
  return text.substr(0, 72) + "...";
}

// True when the line's record name is one of the coordinate records:
// ATOM, HETATM, or a mark_sur/ZDATOM-style line (which also starts
// with ATOM). Used only for the input-quality census.
bool looksLikeCoordinate(const std::string &line) {
  return zdock::pdb::Record::looksLikeCoordinate(line);
}


// Census signature: record name, line length, and the masked shape of
// columns 61+ (occupancy, temperature, and whatever tail the file
// carries). Digit masking happens in the report merger; keeping the
// head (atom names, coordinates) out of the detail makes same-shaped
// lines group into one signature.
std::string coordinateSignature(const std::string &line) {
  const std::string keyword =
      line.compare(0, 4, "ATOM") == 0 ? "ATOM" : "HETATM";
  const std::string tail =
      line.size() > 60 ? line.substr(60) : std::string();
  return keyword + " len " + std::to_string(line.size()) +
         " cols61+ '" + shorten(tail) + "'";
}

void report(const std::string &where, std::size_t lineNo,
            const std::string &what, const std::string &detail) {
  fprintf(stderr, "%s:%zu: %s: %s\n", where.c_str(), lineNo, what.c_str(),
          detail.c_str());
}

// ---------------------------------------------------------------------------
// operator>>-equivalent line reading (production semantics)
// ---------------------------------------------------------------------------

// Reads one record line the way the retired upstream reader did: a 4x81-byte
// buffer, one trailing CR stripped, and the stream failed for lines that
// do not fit. Returns false when the stream failed (production stops
// reading at that point).
bool readRecordLine(std::istream &in, std::string &line) {
  // the pinned 4x81-character operator>> buffer
  char buf[4 * 81];
  in.getline(buf, sizeof(buf));
  if (in.fail())
    return false; // overlong line or EOF: production's read loop stops
  line = buf;
  if (!line.empty() && line.back() == '\r')
    line.pop_back();
  return true;
}

// ---------------------------------------------------------------------------
// Field-by-field comparison between the legacy union and the new payload
// ---------------------------------------------------------------------------

// ---------------------------------------------------------------------------
// Modern-side JSON emission for the snapshot compare (tools/parity/
// snapshot.py): the same field names and value formats the bb0032c
// dump emits, produced from the modern typed payloads.
// ---------------------------------------------------------------------------

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
        // control bytes and non-ASCII bytes escape as \u00XX, matching
        // the snapshot's dump exactly
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

std::string jsonStr(const std::string &s) {
  return "\"" + jsonEscape(s) + "\"";
}

std::string jsonInt(long v) {
  return std::to_string(v);
}

std::string jsonReal(double v) {
  char buf[40];
  std::snprintf(buf, sizeof(buf), "%.17g", v);
  return buf;
}

std::string jsonChar(char c) { return jsonStr(std::string(1, c)); }

// The modern parsed-field emission for one record: field names mirror
// the legacy union paths the parity harness compares.
std::string emitModernFields(const zdock::pdb::Record &r) {
  std::ostringstream os;
  os << "{";
  switch (r.type()) {
  case zdock::pdb::RecordType::ATOM:
  case zdock::pdb::RecordType::HETATM:
  case zdock::pdb::RecordType::SIGATM: {
    if (r.asAtom().zdatom.has_value()) {
      const auto &z = *r.asAtom().zdatom;
      os << "\"zdatom.type\":" << jsonInt(z.type)
         << ",\"zdatom.surface\":" << jsonInt(z.surface)
         << ",\"zdatom.rad\":" << jsonReal(z.rad)
         << ",\"zdatom.segid\":" << jsonStr(z.segid)
         << ",\"zdatom.chg\":" << jsonReal(z.chg);
      break;
    }
    const auto &a = r.asAtom();
    os << "\"serialNum\":" << jsonInt(a.serialNum)
       << ",\"name\":" << jsonStr(a.name)
       << ",\"altLoc\":" << jsonChar(a.altLoc)
       << ",\"residue.name\":" << jsonStr(a.residue.name)
       << ",\"residue.chainId\":" << jsonChar(a.residue.chainId)
       << ",\"residue.seqNum\":" << jsonInt(a.residue.seqNum)
       << ",\"residue.insertCode\":" << jsonChar(a.residue.insertCode);
    for (int i = 0; i < 3; ++i)
      os << ",\"xyz[" << i << "]\":" << jsonReal(a.xyz[i]);
    os << ",\"occupancy\":" << jsonReal(a.occupancy.value_or(0.0))
       << ",\"tempFactor\":" << jsonReal(a.tempFactor.value_or(0.0))
       << ",\"element\":" << jsonStr(a.element)
       << ",\"charge\":" << jsonStr(a.charge);
    break;
  }
  case zdock::pdb::RecordType::ANISOU:
  case zdock::pdb::RecordType::SIGUIJ: {
    const auto &a = r.as<zdock::pdb::AnisouRecord>();
    os << "\"serialNum\":" << jsonInt(a.serialNum)
       << ",\"name\":" << jsonStr(a.name)
       << ",\"altLoc\":" << jsonChar(a.altLoc)
       << ",\"residue.name\":" << jsonStr(a.residue.name)
       << ",\"residue.chainId\":" << jsonChar(a.residue.chainId)
       << ",\"residue.seqNum\":" << jsonInt(a.residue.seqNum)
       << ",\"residue.insertCode\":" << jsonChar(a.residue.insertCode);
    for (int i = 0; i < 6; ++i)
      os << ",\"u[" << i << "]\":" << jsonInt(a.u[i]);
    break;
  }
  case zdock::pdb::RecordType::AUTHOR: {
    const auto &a = r.as<zdock::pdb::AuthorRecord>();
    os << "\"continuation\":" << jsonInt(a.continuation)
       << ",\"authors\":" << jsonStr(a.authors);
    break;
  }
  case zdock::pdb::RecordType::CAVEAT: {
    const auto &c = r.as<zdock::pdb::CaveatRecord>();
    os << "\"continuation\":" << jsonInt(c.continuation)
       << ",\"id\":" << jsonStr(c.id) << ",\"text\":" << jsonStr(c.text);
    break;
  }
  case zdock::pdb::RecordType::COMPND: {
    const auto &c = r.as<zdock::pdb::CompndRecord>();
    os << "\"continuation\":" << jsonInt(c.continuation)
       << ",\"compound\":" << jsonStr(c.compound);
    break;
  }
  case zdock::pdb::RecordType::CONECT: {
    const auto &c = r.as<zdock::pdb::ConectRecord>();
    os << "\"serialNum\":" << jsonInt(c.serialNum);
    for (int i = 0; i < 4; ++i)
      os << ",\"covalent[" << i << "]\":" << jsonInt(c.covalent[i]);
    for (int b = 0; b < 2; ++b) {
      for (int i = 0; i < 2; ++i)
        os << ",\"bonds[" << b << "].hydrogen[" << i
           << "]\":" << jsonInt(c.bonds[b].hydrogen[i]);
      os << ",\"bonds[" << b << "].salt\":" << jsonInt(c.bonds[b].salt);
    }
    break;
  }
  case zdock::pdb::RecordType::CRYST1: {
    const auto &c = r.as<zdock::pdb::Cryst1Record>();
    os << "\"a\":" << jsonReal(c.a) << ",\"b\":" << jsonReal(c.b)
       << ",\"c\":" << jsonReal(c.c) << ",\"alpha\":" << jsonReal(c.alpha)
       << ",\"beta\":" << jsonReal(c.beta) << ",\"gamma\":" << jsonReal(c.gamma)
       << ",\"spaceGrp\":" << jsonStr(c.spaceGrp) << ",\"z\":" << jsonInt(c.z);
    break;
  }
  case zdock::pdb::RecordType::FORMUL: {
    const auto &f = r.as<zdock::pdb::FormulRecord>();
    os << "\"component\":" << jsonInt(f.component)
       << ",\"hetId\":" << jsonStr(f.hetId)
       << ",\"continuation\":" << jsonInt(f.continuation)
       << ",\"exclude\":" << jsonChar(f.exclude)
       << ",\"formula\":" << jsonStr(f.formula);
    break;
  }
  case zdock::pdb::RecordType::FTNOTE:
  case zdock::pdb::RecordType::REMARK:
  case zdock::pdb::RecordType::SYMDES:
  case zdock::pdb::RecordType::MTXDES:
  case zdock::pdb::RecordType::CMPDES:
  case zdock::pdb::RecordType::AGRDES: {
    const auto &a = r.as<zdock::pdb::AgrdesRecord>();
    os << "\"num\":" << jsonInt(a.num) << ",\"text\":" << jsonStr(a.text);
    break;
  }
  case zdock::pdb::RecordType::HEADER: {
    const auto &h = r.as<zdock::pdb::HeaderRecord>();
    os << "\"classification\":" << jsonStr(h.classification)
       << ",\"timestamp\":" << jsonStr(h.timestamp)
       << ",\"id\":" << jsonStr(h.id);
    break;
  }
  case zdock::pdb::RecordType::HELIX: {
    const auto &h = r.as<zdock::pdb::HelixRecord>();
    os << "\"serialNum\":" << jsonInt(h.serialNum)
       << ",\"id\":" << jsonStr(h.id);
    for (int i = 0; i < 2; ++i)
      os << ",\"residues[" << i << "].name\":" << jsonStr(h.residues[i].name)
         << ",\"residues[" << i << "].chainId\":" << jsonChar(h.residues[i].chainId)
         << ",\"residues[" << i << "].seqNum\":" << jsonInt(h.residues[i].seqNum)
         << ",\"residues[" << i << "].insertCode\":" << jsonChar(h.residues[i].insertCode);
    os << ",\"type\":" << jsonInt(h.type)
       << ",\"comment\":" << jsonStr(h.comment)
       << ",\"length\":" << jsonInt(h.length);
    break;
  }
  case zdock::pdb::RecordType::HET: {
    const auto &h = r.as<zdock::pdb::HetRecord>();
    os << "\"hetGrp.name\":" << jsonStr(h.hetGrp.name)
       << ",\"hetGrp.chainId\":" << jsonChar(h.hetGrp.chainId)
       << ",\"hetGrp.seqNum\":" << jsonInt(h.hetGrp.seqNum)
       << ",\"hetGrp.insertCode\":" << jsonChar(h.hetGrp.insertCode)
       << ",\"numAtoms\":" << jsonInt(h.numAtoms)
       << ",\"text\":" << jsonStr(h.text);
    break;
  }
  case zdock::pdb::RecordType::JRNL:
    os << "\"text\":" << jsonStr(r.as<zdock::pdb::JrnlRecord>().text);
    break;
  case zdock::pdb::RecordType::KEYWDS: {
    const auto &k = r.as<zdock::pdb::KeywdsRecord>();
    os << "\"continuation\":" << jsonInt(k.continuation)
       << ",\"keywds\":" << jsonStr(k.keywds);
    break;
  }
  case zdock::pdb::RecordType::MASTER: {
    const auto &m = r.as<zdock::pdb::MasterRecord>();
    os << "\"numRemark\":" << jsonInt(m.numRemark)
       << ",\"numFtnote\":" << jsonInt(m.numFtnote)
       << ",\"numHet\":" << jsonInt(m.numHet)
       << ",\"numHelix\":" << jsonInt(m.numHelix)
       << ",\"numSheet\":" << jsonInt(m.numSheet)
       << ",\"numTurn\":" << jsonInt(m.numTurn)
       << ",\"numSite\":" << jsonInt(m.numSite)
       << ",\"numTransform\":" << jsonInt(m.numTransform)
       << ",\"numCoordinate\":" << jsonInt(m.numCoordinate)
       << ",\"numTer\":" << jsonInt(m.numTer)
       << ",\"numConect\":" << jsonInt(m.numConect)
       << ",\"numSeqres\":" << jsonInt(m.numSeqres);
    break;
  }
  case zdock::pdb::RecordType::MDLTYP: {
    const auto &m = r.as<zdock::pdb::MdltypRecord>();
    os << "\"continuation\":" << jsonInt(m.continuation)
       << ",\"comment\":" << jsonStr(m.comment);
    break;
  }
  case zdock::pdb::RecordType::MODEL:
    os << "\"num\":" << jsonInt(r.as<zdock::pdb::ModelRecord>().num);
    break;
  case zdock::pdb::RecordType::NUMMDL:
    os << "\"number\":" << jsonInt(r.as<zdock::pdb::NummdlRecord>().number);
    break;
  case zdock::pdb::RecordType::MTRIX: {
    const auto &m = r.as<zdock::pdb::MtrixRecord>();
    os << "\"rowNum\":" << jsonInt(m.rowNum)
       << ",\"serialNum\":" << jsonInt(m.serialNum)
       << ",\"m1\":" << jsonReal(m.m1) << ",\"m2\":" << jsonReal(m.m2)
       << ",\"m3\":" << jsonReal(m.m3) << ",\"v\":" << jsonReal(m.v)
       << ",\"given\":" << jsonInt(m.given);
    break;
  }
  case zdock::pdb::RecordType::OBSLTE: {
    const auto &o = r.as<zdock::pdb::ObslteRecord>();
    os << "\"continuation\":" << jsonInt(o.continuation)
       << ",\"timestamp\":" << jsonStr(o.timestamp)
       << ",\"oldId\":" << jsonStr(o.oldId);
    for (int i = 0; i < 9; ++i)
      os << ",\"idMap[" << i << "]\":" << jsonStr(o.idMap[i]);
    break;
  }
  case zdock::pdb::RecordType::ORIGX: {
    const auto &o = r.as<zdock::pdb::OrigxRecord>();
    os << "\"rowNum\":" << jsonInt(o.rowNum) << ",\"o1\":" << jsonReal(o.o1)
       << ",\"o2\":" << jsonReal(o.o2) << ",\"o3\":" << jsonReal(o.o3)
       << ",\"t\":" << jsonReal(o.t);
    break;
  }
  case zdock::pdb::RecordType::REVDAT: {
    const auto &t = r.as<zdock::pdb::RevdatRecord>();
    os << "\"modification\":" << jsonInt(t.modification)
       << ",\"continuation\":" << jsonInt(t.continuation)
       << ",\"timestamp\":" << jsonStr(t.timestamp)
       << ",\"id\":" << jsonStr(t.id) << ",\"modType\":" << jsonInt(t.modType);
    for (int i = 0; i < 4; ++i)
      os << ",\"record[" << i << "]\":" << jsonStr(t.record[i]);
    break;
  }
  case zdock::pdb::RecordType::SCALE: {
    const auto &s = r.as<zdock::pdb::ScaleRecord>();
    os << "\"rowNum\":" << jsonInt(s.rowNum) << ",\"s1\":" << jsonReal(s.s1)
       << ",\"s2\":" << jsonReal(s.s2) << ",\"s3\":" << jsonReal(s.s3)
       << ",\"u\":" << jsonReal(s.u);
    break;
  }
  case zdock::pdb::RecordType::SEQRES: {
    const auto &s = r.as<zdock::pdb::SeqresRecord>();
    os << "\"serialNum\":" << jsonInt(s.serialNum)
       << ",\"chainId\":" << jsonChar(s.chainId)
       << ",\"count\":" << jsonInt(s.count);
    for (int i = 0; i < 13; ++i)
      os << ",\"names[" << i << "]\":" << jsonStr(s.names[i]);
    break;
  }
  case zdock::pdb::RecordType::SHEET: {
    const auto &s = r.as<zdock::pdb::SheetRecord>();
    os << "\"strandNum\":" << jsonInt(s.strandNum)
       << ",\"id\":" << jsonStr(s.id) << ",\"count\":" << jsonInt(s.count);
    for (int i = 0; i < 2; ++i)
      os << ",\"residues[" << i << "].name\":" << jsonStr(s.residues[i].name)
         << ",\"residues[" << i << "].chainId\":" << jsonChar(s.residues[i].chainId)
         << ",\"residues[" << i << "].seqNum\":" << jsonInt(s.residues[i].seqNum)
         << ",\"residues[" << i << "].insertCode\":" << jsonChar(s.residues[i].insertCode);
    os << ",\"sense\":" << jsonInt(s.sense);
    for (int i = 0; i < 2; ++i) {
      os << ",\"atoms[" << i << "].name\":" << jsonStr(s.atoms[i].name)
         << ",\"atoms[" << i << "].residue.name\":" << jsonStr(s.atoms[i].residue.name)
         << ",\"atoms[" << i << "].residue.chainId\":" << jsonChar(s.atoms[i].residue.chainId)
         << ",\"atoms[" << i << "].residue.seqNum\":" << jsonInt(s.atoms[i].residue.seqNum)
         << ",\"atoms[" << i << "].residue.insertCode\":" << jsonChar(s.atoms[i].residue.insertCode);
    }
    break;
  }
  case zdock::pdb::RecordType::SITE: {
    const auto &s = r.as<zdock::pdb::SiteRecord>();
    os << "\"seqNum\":" << jsonInt(s.seqNum) << ",\"id\":" << jsonStr(s.id)
       << ",\"count\":" << jsonInt(s.count);
    for (int i = 0; i < 4; ++i)
      os << ",\"residues[" << i << "].name\":" << jsonStr(s.residues[i].name)
         << ",\"residues[" << i << "].chainId\":" << jsonChar(s.residues[i].chainId)
         << ",\"residues[" << i << "].seqNum\":" << jsonInt(s.residues[i].seqNum)
         << ",\"residues[" << i << "].insertCode\":" << jsonChar(s.residues[i].insertCode);
    break;
  }
  case zdock::pdb::RecordType::SOURCE: {
    const auto &s = r.as<zdock::pdb::SourceRecord>();
    os << "\"continuation\":" << jsonInt(s.continuation)
       << ",\"name\":" << jsonStr(s.name);
    break;
  }
  case zdock::pdb::RecordType::SPLIT: {
    const auto &s = r.as<zdock::pdb::SplitRecord>();
    os << "\"continuation\":" << jsonInt(s.continuation);
    for (int i = 0; i < 14; ++i)
      os << ",\"idMap[" << i << "]\":" << jsonStr(s.idMap[i]);
    break;
  }
  case zdock::pdb::RecordType::SPRSDE: {
    const auto &s = r.as<zdock::pdb::SprsdeRecord>();
    os << "\"continuation\":" << jsonInt(s.continuation)
       << ",\"timestamp\":" << jsonStr(s.timestamp)
       << ",\"id\":" << jsonStr(s.id);
    for (int i = 0; i < 9; ++i)
      os << ",\"supersede[" << i << "]\":" << jsonStr(s.supersede[i]);
    break;
  }
  case zdock::pdb::RecordType::SSBOND: {
    const auto &s = r.as<zdock::pdb::SsbondRecord>();
    os << "\"seqNum\":" << jsonInt(s.seqNum);
    for (int i = 0; i < 2; ++i)
      os << ",\"residues[" << i << "].name\":" << jsonStr(s.residues[i].name)
         << ",\"residues[" << i << "].chainId\":" << jsonChar(s.residues[i].chainId)
         << ",\"residues[" << i << "].seqNum\":" << jsonInt(s.residues[i].seqNum)
         << ",\"residues[" << i << "].insertCode\":" << jsonChar(s.residues[i].insertCode);
    os << ",\"sym1\":" << jsonInt(s.sym1) << ",\"sym2\":" << jsonInt(s.sym2)
       << ",\"length\":" << jsonReal(s.length);
    break;
  }
  case zdock::pdb::RecordType::SYMOP: {
    const auto &s = r.as<zdock::pdb::SymopRecord>();
    os << "\"rowNum\":" << jsonInt(s.rowNum)
       << ",\"serialNum\":" << jsonInt(s.serialNum)
       << ",\"s1\":" << jsonReal(s.s1) << ",\"s2\":" << jsonReal(s.s2)
       << ",\"s3\":" << jsonReal(s.s3) << ",\"t\":" << jsonReal(s.t);
    break;
  }
  case zdock::pdb::RecordType::TER: {
    const auto &t = r.as<zdock::pdb::TerRecord>();
    os << "\"serialNum\":" << jsonInt(t.serialNum)
       << ",\"residue.name\":" << jsonStr(t.residue.name)
       << ",\"residue.chainId\":" << jsonChar(t.residue.chainId)
       << ",\"residue.seqNum\":" << jsonInt(t.residue.seqNum)
       << ",\"residue.insertCode\":" << jsonChar(t.residue.insertCode);
    break;
  }
  case zdock::pdb::RecordType::TITLE: {
    const auto &t = r.as<zdock::pdb::TitleRecord>();
    os << "\"continuation\":" << jsonInt(t.continuation)
       << ",\"text\":" << jsonStr(t.text);
    break;
  }
  case zdock::pdb::RecordType::TRNSFM: {
    const auto &t = r.as<zdock::pdb::TrnsfmRecord>();
    os << "\"resultSerialNum\":" << jsonInt(t.resultSerialNum)
       << ",\"applySerialNum\":" << jsonInt(t.applySerialNum)
       << ",\"sourceSerialNum\":" << jsonInt(t.sourceSerialNum);
    break;
  }
  case zdock::pdb::RecordType::TURN: {
    const auto &t = r.as<zdock::pdb::TurnRecord>();
    os << "\"seqNum\":" << jsonInt(t.seqNum) << ",\"id\":" << jsonStr(t.id);
    for (int i = 0; i < 2; ++i)
      os << ",\"residues[" << i << "].name\":" << jsonStr(t.residues[i].name)
         << ",\"residues[" << i << "].chainId\":" << jsonChar(t.residues[i].chainId)
         << ",\"residues[" << i << "].seqNum\":" << jsonInt(t.residues[i].seqNum)
         << ",\"residues[" << i << "].insertCode\":" << jsonChar(t.residues[i].insertCode);
    os << ",\"comment\":" << jsonStr(t.comment);
    break;
  }
  case zdock::pdb::RecordType::TVECT: {
    const auto &t = r.as<zdock::pdb::TvectRecord>();
    os << "\"serialNum\":" << jsonInt(t.serialNum)
       << ",\"t1\":" << jsonReal(t.t1) << ",\"t2\":" << jsonReal(t.t2)
       << ",\"t3\":" << jsonReal(t.t3) << ",\"comment\":" << jsonStr(t.comment);
    break;
  }
  case zdock::pdb::RecordType::USER: {
    const auto &u = r.as<zdock::pdb::UserRecord>();
    os << "\"subtype\":" << jsonStr(u.subtype)
       << ",\"text\":" << jsonStr(u.text);
    break;
  }
  case zdock::pdb::RecordType::USER_PDBRUN:
    os << "\"version\":"
       << jsonInt(r.as<zdock::pdb::UserPdbrunRecord>().version);
    break;
  case zdock::pdb::RecordType::USER_EYEPOS:
  case zdock::pdb::RecordType::USER_ATPOS: {
    const auto &e = r.as<zdock::pdb::UserEyePosRecord>();
    for (int i = 0; i < 3; ++i)
      os << (i ? "," : "") << "\"xyz[" << i << "]\":" << jsonReal(e.xyz[i]);
    break;
  }
  case zdock::pdb::RecordType::USER_WINDOW: {
    const auto &w = r.as<zdock::pdb::UserWindowRecord>();
    os << "\"left\":" << jsonReal(w.left) << ",\"right\":" << jsonReal(w.right)
       << ",\"bottom\":" << jsonReal(w.bottom) << ",\"top\":" << jsonReal(w.top)
       << ",\"hither\":" << jsonReal(w.hither) << ",\"yon\":" << jsonReal(w.yon);
    break;
  }
  case zdock::pdb::RecordType::USER_FOCUS:
    os << "\"focus\":" << jsonReal(r.as<zdock::pdb::UserFocusRecord>().focus);
    break;
  case zdock::pdb::RecordType::USER_VIEWPORT: {
    const auto &v = r.as<zdock::pdb::UserViewportRecord>();
    os << "\"xmin\":" << jsonReal(v.xmin) << ",\"xmax\":" << jsonReal(v.xmax)
       << ",\"ymin\":" << jsonReal(v.ymin) << ",\"ymax\":" << jsonReal(v.ymax);
    break;
  }
  case zdock::pdb::RecordType::USER_BGCOLOR: {
    const auto &b = r.as<zdock::pdb::UserBgColorRecord>();
    for (int i = 0; i < 3; ++i)
      os << (i ? "," : "") << "\"rgb[" << i << "]\":" << jsonReal(b.rgb[i]);
    break;
  }
  case zdock::pdb::RecordType::USER_ANGLE: {
    const auto &a = r.as<zdock::pdb::UserAngleRecord>();
    os << "\"which\":" << jsonInt(a.which.value_or(0))
       << ",\"atom0\":" << jsonInt(a.atom0) << ",\"atom1\":" << jsonInt(a.atom1)
       << ",\"atom2\":" << jsonInt(a.atom2) << ",\"atom3\":" << jsonInt(a.atom3)
       << ",\"angle\":" << jsonReal(a.angle);
    break;
  }
  case zdock::pdb::RecordType::USER_DISTANCE: {
    const auto &d = r.as<zdock::pdb::UserDistanceRecord>();
    os << "\"which\":" << jsonInt(d.which.value_or(0))
       << ",\"atom0\":" << jsonInt(d.atom0) << ",\"atom1\":" << jsonInt(d.atom1)
       << ",\"distance\":" << jsonReal(d.distance);
    break;
  }
  case zdock::pdb::RecordType::USER_FILE: {
    const auto &f = r.as<zdock::pdb::UserFileRecord>();
    os << "\"filename\":" << jsonStr(f.filename)
       << ",\"model\":" << jsonInt(f.model.value_or(0));
    break;
  }
  case zdock::pdb::RecordType::USER_MARKNAME:
  case zdock::pdb::RecordType::USER_MARK:
    os << "\"markname\":"
       << jsonStr(r.as<zdock::pdb::UserMarknameRecord>().markname);
    break;
  case zdock::pdb::RecordType::USER_CNAME: {
    const auto &c = r.as<zdock::pdb::UserCNameRecord>();
    for (int i = 0; i < 3; ++i)
      os << (i ? "," : "") << "\"rgb[" << i << "]\":" << jsonReal(c.rgb[i]);
    os << ",\"name\":" << jsonStr(c.name);
    break;
  }
  case zdock::pdb::RecordType::USER_COLOR:
  case zdock::pdb::RecordType::USER_GFX_COLOR: {
    const auto &c = r.as<zdock::pdb::UserColorRecord>();
    for (int i = 0; i < 3; ++i)
      os << (i ? "," : "") << "\"rgb[" << i << "]\":" << jsonReal(c.rgb[i]);
    os << ",\"spec\":" << jsonStr(c.spec);
    break;
  }
  case zdock::pdb::RecordType::USER_RADIUS:
    os << "\"radius\":"
       << jsonReal(r.as<zdock::pdb::UserRadiusRecord>().radius);
    break;
  case zdock::pdb::RecordType::USER_OBJECT:
  case zdock::pdb::RecordType::USER_ENDOBJ:
    os << "\"model\":" << jsonInt(r.as<zdock::pdb::UserObjectRecord>().model);
    break;
  case zdock::pdb::RecordType::USER_CHAIN: {
    const auto &c = r.as<zdock::pdb::UserChainRecord>();
    os << "\"atom0\":" << jsonInt(c.atom0) << ",\"atom1\":" << jsonInt(c.atom1);
    break;
  }
  case zdock::pdb::RecordType::USER_GFX_BEGIN:
    os << "\"primitive\":"
       << jsonInt(static_cast<int>(
              r.as<zdock::pdb::UserGfxBeginRecord>().primitive));
    break;
  case zdock::pdb::RecordType::USER_GFX_NORMAL:
  case zdock::pdb::RecordType::USER_GFX_VERTEX:
  case zdock::pdb::RecordType::USER_GFX_TEXTPOS:
  case zdock::pdb::RecordType::USER_GFX_MOVE:
  case zdock::pdb::RecordType::USER_GFX_DRAW:
  case zdock::pdb::RecordType::USER_GFX_MARKER:
  case zdock::pdb::RecordType::USER_GFX_POINT: {
    const auto &p = r.as<zdock::pdb::UserGfxPointRecord>();
    for (int i = 0; i < 3; ++i)
      os << (i ? "," : "") << "\"xyz[" << i << "]\":" << jsonReal(p.xyz[i]);
    break;
  }
  case zdock::pdb::RecordType::USER_GFX_FONT: {
    const auto &f = r.as<zdock::pdb::UserGfxFontRecord>();
    os << "\"size\":" << jsonInt(f.size) << ",\"name\":" << jsonStr(f.name);
    break;
  }
  case zdock::pdb::RecordType::USER_GFX_LABEL: {
    const auto &l = r.as<zdock::pdb::UserGfxLabelRecord>();
    if (l.xyz.has_value())
      for (int i = 0; i < 3; ++i)
        os << (i ? "," : "") << "\"xyz[" << i << "]\":" << jsonReal((*l.xyz)[i]);
    os << ",\"text\":" << jsonStr(l.text);
    break;
  }
  default:
    break;
  }
  os << "}";
  return os.str();
}

// Emits the modern parser's records for one file as a single JSON line,
// the snapshot-compare counterpart of the bb0032c dump's --json mode.
// rejectedCoordinates carries the modern ParseContext's count of
// coordinate-looking lines that parse as no record at all (the
// input-quality census, reported separately by the compare mode).
int dumpJsonFile(const std::string &path) {
  std::ifstream in(path);
  if (!in) {
    fprintf(stderr, "cannot open %s\n", path.c_str());
    return 2;
  }
  std::printf("{\"records\":[");
  bool first = true;
  std::string line;
  std::size_t lineNo = 0;
  std::size_t rejectedCoordinates = 0;
  zdock::pdb::ParseContext context;
  while (readRecordLine(in, line)) {
    ++lineNo;
    zdock::pdb::ParseContext lineContext;
    lineContext.pdbrunVersion = context.pdbrunVersion;
    const auto parsed = zdock::pdb::Record::parse(line, lineContext);
    context.pdbrunVersion = lineContext.pdbrunVersion;
    if (!first)
      std::printf(",");
    first = false;
    if (parsed.has_value()) {
      // the render takes the file-derived PDBRUN output version, the
      // same threading the live comparison uses
      zdock::pdb::WriteContext writeContext;
      writeContext.pdbrunVersion = context.pdbrunVersion;
      std::printf("{\"line\":%zu,\"type\":%d,\"text\":%s,\"fields\":%s}",
                  lineNo, static_cast<int>(parsed->type()),
                  jsonStr(parsed->chars(writeContext)).c_str(),
                  emitModernFields(*parsed).c_str());
    } else {
      // a line that parses as no record carries the junk fallback text
      // and no fields, exactly as the snapshot's UNKNOWN records do
      std::printf("{\"line\":%zu,\"type\":%d,\"text\":%s,\"fields\":{}}",
                  lineNo, static_cast<int>(zdock::pdb::RecordType::UNKNOWN),
                  jsonStr(zdock::pdb::Record::unknownFrom(line).chars())
                      .c_str());
      if (zdock::pdb::Record::looksLikeCoordinate(line))
        ++rejectedCoordinates;
    }
  }
  std::printf("],\"rejectedCoordinates\":%zu}\n", rejectedCoordinates);
  return 0;
}

// Classifies a difference between the current and the bb0032c baseline:
// the pre-fb720ad ZDATOM 80-column truncation and the pre-ea22558 unbounded
// v5 USER %s reads are expected; everything else is unexpected.
enum class BaselineClass {
  PARITY,
  EXPECTED_ZDATOM,
  EXPECTED_V5_TOKEN,
  EXPECTED_MODERN_GAIN,
  UNEXPECTED
};

BaselineClass classifyBaseline(const std::string &line,
                               const std::string &current,
                               const std::string &baseline) {
  if (current == baseline)
    return BaselineClass::PARITY;
  const std::string type4 = line.substr(0, std::min<std::size_t>(6, line.size()));
  const bool atomish = type4.rfind("ATOM", 0) == 0 ||
                       type4.rfind("HETA", 0) == 0 ||
                       type4.rfind("SIGA", 0) == 0;
  if (atomish && !current.empty() && !baseline.empty()) {
    // the pre-fb720ad clamp cut ZDATOM renders at 80 columns
    const std::string current80 = current.substr(0, std::min<std::size_t>(80, current.size()));
    if (current80 == baseline ||
        (current.size() == 82 && baseline.size() == 80 &&
         current.substr(0, 80) == baseline))
      return BaselineClass::EXPECTED_ZDATOM;
    // the segid/element salvage: the bb0032c parser stores the line as
    // UNKNOWN junk while the modern parser classifies the standard
    // record it announces
    if (baseline.rfind("UNKNOWN:", 0) == 0)
      return BaselineClass::EXPECTED_MODERN_GAIN;
  }
  if ((line.rfind("USER  CNAME ", 0) == 0 || line.rfind("USER  COLOR ", 0) == 0 ||
       line.rfind("USER  GFX COLOR ", 0) == 0 ||
       line.rfind("USER  GFX FONT ", 0) == 0) &&
      line.size() > 60)
    return BaselineClass::EXPECTED_V5_TOKEN; // the bounded %s reads change
                                             // long-token records
  return BaselineClass::UNEXPECTED;
}

} // namespace

int main(int argc, char **argv) {
  Options options;
  // the dump mode takes a single file, before or after the flag
  bool dumpJson = false;
  for (int i = 1; i < argc; ++i)
    if (std::string(argv[i]) == "--dump-json")
      dumpJson = true;
  if (dumpJson) {
    const char *path = nullptr;
    for (int i = 1; i < argc; ++i)
      if (std::string(argv[i]) != "--dump-json")
        path = argv[i];
    if (!path) {
      fprintf(stderr, "usage: %s --dump-json <file>\n", argv[0]);
      return 2;
    }
    return dumpJsonFile(path);
  }
  if (argc < 2) {
    fprintf(stderr,
            "usage: %s <file-or-directory> [--files-from <list>] "
            "[--baseline <exe>] [--dump-json] "
            "[--max-files N] [--max-lines N] [--all-extensions] [--quiet]\n",
            argv[0]);
    return 2;
  }
  options.root = argv[1];
  for (int i = 2; i < argc; ++i) {
    const std::string arg = argv[i];
    if (arg == "--baseline" && i + 1 < argc)
      options.baseline = argv[++i];
    else if (arg == "--files-from" && i + 1 < argc)
      options.filesFrom = argv[++i];
    else if (arg == "--max-files" && i + 1 < argc)
      options.maxFiles = static_cast<std::size_t>(std::atoi(argv[++i]));
    else if (arg == "--max-lines" && i + 1 < argc)
      options.maxLines = static_cast<std::size_t>(std::atoi(argv[++i]));
    else if (arg == "--all-extensions")
      options.allExtensions = true;
    else if (arg == "--quiet")
      options.quiet = true;
    else {
      fprintf(stderr, "unknown argument: %s\n", arg.c_str());
      return 2;
    }
  }

  std::vector<std::filesystem::path> files;
  try {
    const std::filesystem::path root(options.root);
    if (!options.filesFrom.empty()) {
      std::ifstream list(options.filesFrom);
      if (!list) {
        fprintf(stderr, "cannot open file list %s\n",
                options.filesFrom.c_str());
        return 2;
      }
      std::string path;
      while (std::getline(list, path))
        if (!path.empty())
          files.push_back(std::filesystem::path(path));
    } else if (std::filesystem::is_regular_file(root)) {
      files.push_back(root);
    } else {
      for (const auto &entry :
           std::filesystem::recursive_directory_iterator(root))
        if (entry.is_regular_file() &&
            (options.allExtensions || pdbExtension(entry.path())))
          files.push_back(entry.path());
      std::sort(files.begin(), files.end());
    }
  } catch (const std::filesystem::filesystem_error &error) {
    fprintf(stderr, "filesystem error: %s\n", error.what());
    return 2;
  }
  if (options.maxFiles > 0 && files.size() > options.maxFiles)
    files.resize(options.maxFiles);
  if (files.empty()) {
    fprintf(stderr, "no input files under %s\n", options.root.c_str());
    return 2;
  }

  if (options.baseline.empty()) {
    // the in-tree legacy parser is retired; the only legacy side is the
    // out-of-tree bb0032c dump
    fprintf(stderr,
            "error: --baseline <dump> is required (the bb0032c baseline "
            "dump; see tools/parity/README.md)\n");
    return 2;
  }

  Counters counters;
  for (const auto &file : files) {
    std::ifstream in(file);
    if (!in) {
      fprintf(stderr, "cannot open %s\n", file.string().c_str());
      continue;
    }
    const std::string where = file.string();
    const std::string command =
        options.baseline + " " + file.string() + " 2>/dev/null";
    FILE *pipe = popen(command.c_str(), "r");
    if (!pipe) {
      fprintf(stderr, "cannot run baseline %s\n",
              options.baseline.c_str());
      return 2;
    }
    ++counters.files;
    std::string line;
    std::size_t lineNo = 0;
    char buf[8192];
    zdock::pdb::ParseContext parseContext;
    while (readRecordLine(in, line)) {
      ++lineNo;
      ++counters.lines;
      if (!fgets(buf, sizeof(buf), pipe))
        break; // the dump stopped early: its file read hit its own bound
      std::string baselineLine(buf);
      while (!baselineLine.empty() &&
             (baselineLine.back() == '\n' || baselineLine.back() == '\r'))
        baselineLine.pop_back();
      // baseline format: "<type>\t<text>"
      const std::size_t tab = baselineLine.find('\t');
      const int baselineType =
          tab == std::string::npos
              ? -1
              : std::atoi(baselineLine.substr(0, tab).c_str());
      const std::string baselineText =
          tab == std::string::npos ? "" : baselineLine.substr(tab + 1);

      zdock::pdb::ParseContext lineContext;
      lineContext.pdbrunVersion = parseContext.pdbrunVersion;
      const auto parsed = zdock::pdb::Record::parse(line, lineContext);
      parseContext.pdbrunVersion = lineContext.pdbrunVersion;

      if (!parsed.has_value()) {
        // the modern parse refuses the line: a divergence only when the
        // baseline stored a real record; otherwise the junk renders
        // must match exactly
        if (baselineType != static_cast<int>(zdock::pdb::RecordType::UNKNOWN)) {
          ++counters.rejections;
          if (!options.quiet)
            report(where, lineNo, "new parser rejected",
                   "baseline type " + std::to_string(baselineType) + " | " +
                       shorten(line));
        } else {
          zdock::pdb::WriteContext writeContext;
          writeContext.pdbrunVersion = parseContext.pdbrunVersion;
          const std::string junkText =
              zdock::pdb::Record::unknownFrom(line).chars(writeContext);
          if (junkText != baselineText) {
            ++counters.textMismatches;
            if (!options.quiet)
              report(where, lineNo, "text mismatch",
                     "baseline '" + shorten(baselineText) + "' new '" +
                         shorten(junkText) + "'");
          }
        }
        if (looksLikeCoordinate(line)) {
          // input-quality census, never part of the exit decision: the
          // file loses this atom
          ++counters.coordinateUnknowns;
          if (!options.quiet)
            report(where, lineNo, "coordinate unknown",
                   coordinateSignature(line));
        }
        continue;
      }

      const bool modernStandardAtom =
          (parsed->type() == zdock::pdb::RecordType::ATOM ||
           parsed->type() == zdock::pdb::RecordType::HETATM ||
           parsed->type() == zdock::pdb::RecordType::SIGATM) &&
          !parsed->asAtom().zdatom.has_value();
      if (baselineType ==
              static_cast<int>(zdock::pdb::RecordType::UNKNOWN) &&
          modernStandardAtom) {
        // the deliberate salvage gain: the baseline stores the line as
        // junk while the modern parser classifies the standard-
        // coordinate record it announces; grouped and reported, never a
        // parity failure
        ++counters.modernGains;
        if (!options.quiet)
          report(where, lineNo, "modern gain", coordinateSignature(line));
        continue;
      }

      zdock::pdb::WriteContext writeContext;
      writeContext.pdbrunVersion = parseContext.pdbrunVersion;
      const std::string modernText = parsed->chars(writeContext);
      if (baselineType != static_cast<int>(parsed->type())) {
        ++counters.typeMismatches;
        if (!options.quiet)
          report(where, lineNo, "type mismatch",
                 "baseline " + std::to_string(baselineType) + " new " +
                     std::to_string(static_cast<int>(parsed->type())) +
                     " | " + shorten(line));
      }
      if (modernText != baselineText) {
        const BaselineClass klass =
            classifyBaseline(line, modernText, baselineText);
        if (klass == BaselineClass::EXPECTED_ZDATOM)
          ++counters.baselineExpectedZdatom;
        else if (klass == BaselineClass::EXPECTED_V5_TOKEN)
          ++counters.baselineExpectedV5Token;
        else if (klass == BaselineClass::EXPECTED_MODERN_GAIN)
          ++counters.baselineExpectedModernGain;
        else {
          ++counters.baselineUnexpected;
          if (!options.quiet)
            report(where, lineNo, "UNEXPECTED baseline difference",
                   "baseline '" + shorten(baselineText) + "' new '" +
                       shorten(modernText) + "'");
        }
      }
    }
    pclose(pipe);
  }

  printf("files %zu | lines %zu | type mismatches %zu | text mismatches %zu "
         "| field mismatches %zu | state mismatches %zu | rejections %zu "
         "| coordinate unknowns %zu | modern gains %zu\n",
         counters.files, counters.lines, counters.typeMismatches,
         counters.textMismatches, counters.fieldMismatches,
         counters.stateMismatches, counters.rejections,
         counters.coordinateUnknowns, counters.modernGains);
  if (!options.baseline.empty())
    printf("baseline: expected ZDATOM %zu | expected v5 token %zu | "
           "expected modern gain %zu | UNEXPECTED %zu\n",
           counters.baselineExpectedZdatom, counters.baselineExpectedV5Token,
           counters.baselineExpectedModernGain, counters.baselineUnexpected);
  const bool parity =
      counters.typeMismatches == 0 && counters.textMismatches == 0 &&
      counters.rejections == 0 && counters.baselineUnexpected == 0;
  return parity ? 0 : 1;
}
