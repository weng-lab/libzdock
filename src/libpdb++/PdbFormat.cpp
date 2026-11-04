// Copyright (c) 2019-2026 Arjan van der Velde, Weng Lab
// SPDX-License-Identifier: BSD-2-Clause

#include "PdbFormat.hpp"

#include "FieldReader.hpp"
#include "RecordClassifier.hpp"

#include <cstddef>
#include <cstdio>
#include <ctype.h>

#ifndef NULL
#define NULL 0
#endif

//
//	for each pdb record type there is a format reading in the
//	record values and for printing them out.
//
//	The actual format of a line written, is the print format
//	followed by blank padding to 72 characters, followed by
//	8 characters of file and line information.
//

namespace libpdb {
namespace {

using Kind = PdbFormat::FieldKind;
using Field = PdbFormat::Field;

constexpr Field endField() { return {Kind::End, 0}; }

const Field kAnisouFields[] = {
    {Kind::Int, offsetof(PDB::Anisou, serialNum)},
    {Kind::Str, offsetof(PDB::Anisou, name)},
    {Kind::Char, offsetof(PDB::Anisou, altLoc)},
    {Kind::Str, offsetof(PDB::Anisou, residue.name)},
    {Kind::Char, offsetof(PDB::Anisou, residue.chainId)},
    {Kind::Int, offsetof(PDB::Anisou, residue.seqNum)},
    {Kind::Char, offsetof(PDB::Anisou, residue.insertCode)},
    {Kind::Int, offsetof(PDB::Anisou, u[0])},
    {Kind::Int, offsetof(PDB::Anisou, u[1])},
    {Kind::Int, offsetof(PDB::Anisou, u[2])},
    {Kind::Int, offsetof(PDB::Anisou, u[3])},
    {Kind::Int, offsetof(PDB::Anisou, u[4])},
    {Kind::Int, offsetof(PDB::Anisou, u[5])},
    endField(),
};

const Field kAtomFields[] = {
    {Kind::Int, offsetof(PDB::Atom, serialNum)},
    {Kind::Str, offsetof(PDB::Atom, name)},
    {Kind::Char, offsetof(PDB::Atom, altLoc)},
    {Kind::Str, offsetof(PDB::Atom, residue.name)},
    {Kind::Char, offsetof(PDB::Atom, residue.chainId)},
    {Kind::Int, offsetof(PDB::Atom, residue.seqNum)},
    {Kind::Char, offsetof(PDB::Atom, residue.insertCode)},
    {Kind::Real, offsetof(PDB::Atom, xyz[0])},
    {Kind::Real, offsetof(PDB::Atom, xyz[1])},
    {Kind::Real, offsetof(PDB::Atom, xyz[2])},
    {Kind::Real, offsetof(PDB::Atom, occupancy)},
    {Kind::Real, offsetof(PDB::Atom, tempFactor)},
    {Kind::Str, offsetof(PDB::Atom, element)},
    {Kind::Str, offsetof(PDB::Atom, charge)},
    endField(),
};

const Field kZdatomFields[] = {
    {Kind::Int, offsetof(PDB::Atom, serialNum)},
    {Kind::Str, offsetof(PDB::Atom, name)},
    {Kind::Char, offsetof(PDB::Atom, altLoc)},
    {Kind::Str, offsetof(PDB::Atom, residue.name)},
    {Kind::Char, offsetof(PDB::Atom, residue.chainId)},
    {Kind::Int, offsetof(PDB::Atom, residue.seqNum)},
    {Kind::Char, offsetof(PDB::Atom, residue.insertCode)},
    {Kind::Real, offsetof(PDB::Atom, xyz[0])},
    {Kind::Real, offsetof(PDB::Atom, xyz[1])},
    {Kind::Real, offsetof(PDB::Atom, xyz[2])},
    {Kind::Int, offsetof(PDB::Atom, type)},
    {Kind::Int, offsetof(PDB::Atom, surface)},
    {Kind::Real, offsetof(PDB::Atom, rad)},
    {Kind::Str, offsetof(PDB::Atom, segid)},
    {Kind::Real, offsetof(PDB::Atom, chg)},
    endField(),
};

const Field kAgrdesFields[] = {
    {Kind::Int, offsetof(PDB::Agrdes, num)},
    {Kind::Str, offsetof(PDB::Agrdes, text)},
    endField(),
};

const Field kAuthorFields[] = {
    {Kind::Int, offsetof(PDB::Author, continuation)},
    {Kind::Str, offsetof(PDB::Author, authors)},
    endField(),
};

const Field kCompndFields[] = {
    {Kind::Int, offsetof(PDB::Compnd, continuation)},
    {Kind::Str, offsetof(PDB::Compnd, compound)},
    endField(),
};

const Field kExpdtaFields[] = {
    {Kind::Int, offsetof(PDB::Expdta, continuation)},
    {Kind::Str, offsetof(PDB::Expdta, technique)},
    endField(),
};

const Field kJrnlFields[] = {
    {Kind::Str, offsetof(PDB::Jrnl, text)},
    endField(),
};

const Field kKeywdsFields[] = {
    {Kind::Int, offsetof(PDB::Keywds, continuation)},
    {Kind::Str, offsetof(PDB::Keywds, keywds)},
    endField(),
};

const Field kMdltypFields[] = {
    {Kind::Int, offsetof(PDB::Mdltyp, continuation)},
    {Kind::Str, offsetof(PDB::Mdltyp, comment)},
    endField(),
};

const Field kModelFields[] = {
    {Kind::Int, offsetof(PDB::Model, num)},
    endField(),
};

const Field kNummdlFields[] = {
    {Kind::Int, offsetof(PDB::Nummdl, number)},
    endField(),
};

const Field kSourceFields[] = {
    {Kind::Int, offsetof(PDB::Source, continuation)},
    {Kind::Str, offsetof(PDB::Source, name)},
    endField(),
};

const Field kTitleFields[] = {
    {Kind::Int, offsetof(PDB::Title, continuation)},
    {Kind::Str, offsetof(PDB::Title, text)},
    endField(),
};

const Field kUserFields[] = {
    {Kind::Str, offsetof(PDB::User, subtype)},
    {Kind::Str, offsetof(PDB::User, text)},
    endField(),
};

const Field kCaveatFields[] = {
    {Kind::Int, offsetof(PDB::Caveat, continuation)},
    {Kind::Str, offsetof(PDB::Caveat, id)},
    {Kind::Str, offsetof(PDB::Caveat, text)},
    endField(),
};

const Field kFormulFields[] = {
    {Kind::Int, offsetof(PDB::Formul, component)},
    {Kind::Str, offsetof(PDB::Formul, hetId)},
    {Kind::Int, offsetof(PDB::Formul, continuation)},
    {Kind::Char, offsetof(PDB::Formul, exclude)},
    {Kind::Str, offsetof(PDB::Formul, formula)},
    endField(),
};

const Field kHeaderFields[] = {
    {Kind::Str, offsetof(PDB::Header, classification)},
    {Kind::Str, offsetof(PDB::Header, timestamp)},
    {Kind::Str, offsetof(PDB::Header, id)},
    endField(),
};

const Field kHetFields[] = {
    {Kind::Str, offsetof(PDB::Het, hetGrp.name)},
    {Kind::Char, offsetof(PDB::Het, hetGrp.chainId)},
    {Kind::Int, offsetof(PDB::Het, hetGrp.seqNum)},
    {Kind::Char, offsetof(PDB::Het, hetGrp.insertCode)},
    {Kind::Int, offsetof(PDB::Het, numAtoms)},
    {Kind::Str, offsetof(PDB::Het, text)},
    endField(),
};

const Field kConectFields[] = {
    {Kind::Int, offsetof(PDB::Conect, serialNum)},
    {Kind::Int, offsetof(PDB::Conect, covalent[0])},
    {Kind::Int, offsetof(PDB::Conect, covalent[1])},
    {Kind::Int, offsetof(PDB::Conect, covalent[2])},
    {Kind::Int, offsetof(PDB::Conect, covalent[3])},
    {Kind::Int, offsetof(PDB::Conect, bonds[0].hydrogen[0])},
    {Kind::Int, offsetof(PDB::Conect, bonds[0].hydrogen[1])},
    {Kind::Int, offsetof(PDB::Conect, bonds[0].salt)},
    {Kind::Int, offsetof(PDB::Conect, bonds[1].hydrogen[0])},
    {Kind::Int, offsetof(PDB::Conect, bonds[1].hydrogen[1])},
    {Kind::Int, offsetof(PDB::Conect, bonds[1].salt)},
    endField(),
};

const Field kCryst1Fields[] = {
    {Kind::Real, offsetof(PDB::Cryst1, a)},
    {Kind::Real, offsetof(PDB::Cryst1, b)},
    {Kind::Real, offsetof(PDB::Cryst1, c)},
    {Kind::Real, offsetof(PDB::Cryst1, alpha)},
    {Kind::Real, offsetof(PDB::Cryst1, beta)},
    {Kind::Real, offsetof(PDB::Cryst1, gamma)},
    {Kind::Str, offsetof(PDB::Cryst1, spaceGrp)},
    {Kind::Int, offsetof(PDB::Cryst1, z)},
    endField(),
};

const Field kMtrixFields[] = {
    {Kind::Int, offsetof(PDB::Mtrix, rowNum)},
    {Kind::Int, offsetof(PDB::Mtrix, serialNum)},
    {Kind::Real, offsetof(PDB::Mtrix, m1)},
    {Kind::Real, offsetof(PDB::Mtrix, m2)},
    {Kind::Real, offsetof(PDB::Mtrix, m3)},
    {Kind::Real, offsetof(PDB::Mtrix, v)},
    {Kind::Int, offsetof(PDB::Mtrix, given)},
    endField(),
};

const Field kOrigxFields[] = {
    {Kind::Int, offsetof(PDB::Origx, rowNum)},
    {Kind::Real, offsetof(PDB::Origx, o1)},
    {Kind::Real, offsetof(PDB::Origx, o2)},
    {Kind::Real, offsetof(PDB::Origx, o3)},
    {Kind::Real, offsetof(PDB::Origx, t)},
    endField(),
};

const Field kScaleFields[] = {
    {Kind::Int, offsetof(PDB::Scale, rowNum)},
    {Kind::Real, offsetof(PDB::Scale, s1)},
    {Kind::Real, offsetof(PDB::Scale, s2)},
    {Kind::Real, offsetof(PDB::Scale, s3)},
    {Kind::Real, offsetof(PDB::Scale, u)},
    endField(),
};

const Field kSymopFields[] = {
    {Kind::Int, offsetof(PDB::Symop, rowNum)},
    {Kind::Int, offsetof(PDB::Symop, serialNum)},
    {Kind::Real, offsetof(PDB::Symop, s1)},
    {Kind::Real, offsetof(PDB::Symop, s2)},
    {Kind::Real, offsetof(PDB::Symop, s3)},
    {Kind::Real, offsetof(PDB::Symop, t)},
    endField(),
};

const Field kTrnsfmFields[] = {
    {Kind::Int, offsetof(PDB::Trnsfm, resultSerialNum)},
    {Kind::Int, offsetof(PDB::Trnsfm, applySerialNum)},
    {Kind::Int, offsetof(PDB::Trnsfm, sourceSerialNum)},
    endField(),
};

const Field kAggrgtFields[] = {
    {Kind::Int, offsetof(PDB::Aggrgt, serialNum)},
    {Kind::Int, offsetof(PDB::Aggrgt, numComponents)},
    {Kind::Int, offsetof(PDB::Aggrgt, cmpontSerialNums[0])},
    {Kind::Int, offsetof(PDB::Aggrgt, cmpontSerialNums[1])},
    {Kind::Int, offsetof(PDB::Aggrgt, cmpontSerialNums[2])},
    {Kind::Int, offsetof(PDB::Aggrgt, cmpontSerialNums[3])},
    {Kind::Int, offsetof(PDB::Aggrgt, cmpontSerialNums[4])},
    {Kind::Int, offsetof(PDB::Aggrgt, cmpontSerialNums[5])},
    {Kind::Int, offsetof(PDB::Aggrgt, cmpontSerialNums[6])},
    {Kind::Int, offsetof(PDB::Aggrgt, cmpontSerialNums[7])},
    {Kind::Int, offsetof(PDB::Aggrgt, cmpontSerialNums[8])},
    {Kind::Int, offsetof(PDB::Aggrgt, cmpontSerialNums[9])},
    {Kind::Int, offsetof(PDB::Aggrgt, cmpontSerialNums[10])},
    {Kind::Int, offsetof(PDB::Aggrgt, cmpontSerialNums[11])},
    {Kind::Int, offsetof(PDB::Aggrgt, cmpontSerialNums[12])},
    {Kind::Int, offsetof(PDB::Aggrgt, cmpontSerialNums[13])},
    endField(),
};

const Field kHelixFields[] = {
    {Kind::Int, offsetof(PDB::Helix, serialNum)},
    {Kind::Str, offsetof(PDB::Helix, id)},
    {Kind::Str, offsetof(PDB::Helix, residues[0].name)},
    {Kind::Char, offsetof(PDB::Helix, residues[0].chainId)},
    {Kind::Int, offsetof(PDB::Helix, residues[0].seqNum)},
    {Kind::Char, offsetof(PDB::Helix, residues[0].insertCode)},
    {Kind::Str, offsetof(PDB::Helix, residues[1].name)},
    {Kind::Char, offsetof(PDB::Helix, residues[1].chainId)},
    {Kind::Int, offsetof(PDB::Helix, residues[1].seqNum)},
    {Kind::Char, offsetof(PDB::Helix, residues[1].insertCode)},
    {Kind::Int, offsetof(PDB::Helix, type)},
    {Kind::Str, offsetof(PDB::Helix, comment)},
    {Kind::Int, offsetof(PDB::Helix, length)},
    endField(),
};

const Field kSheetFields[] = {
    {Kind::Int, offsetof(PDB::Sheet, strandNum)},
    {Kind::Str, offsetof(PDB::Sheet, id)},
    {Kind::Int, offsetof(PDB::Sheet, count)},
    {Kind::Str, offsetof(PDB::Sheet, residues[0].name)},
    {Kind::Char, offsetof(PDB::Sheet, residues[0].chainId)},
    {Kind::Int, offsetof(PDB::Sheet, residues[0].seqNum)},
    {Kind::Char, offsetof(PDB::Sheet, residues[0].insertCode)},
    {Kind::Str, offsetof(PDB::Sheet, residues[1].name)},
    {Kind::Char, offsetof(PDB::Sheet, residues[1].chainId)},
    {Kind::Int, offsetof(PDB::Sheet, residues[1].seqNum)},
    {Kind::Char, offsetof(PDB::Sheet, residues[1].insertCode)},
    {Kind::Int, offsetof(PDB::Sheet, sense)},
    {Kind::Str, offsetof(PDB::Sheet, atoms[0].name)},
    {Kind::Str, offsetof(PDB::Sheet, atoms[0].residue.name)},
    {Kind::Char, offsetof(PDB::Sheet, atoms[0].residue.chainId)},
    {Kind::Int, offsetof(PDB::Sheet, atoms[0].residue.seqNum)},
    {Kind::Char, offsetof(PDB::Sheet, atoms[0].residue.insertCode)},
    {Kind::Str, offsetof(PDB::Sheet, atoms[1].name)},
    {Kind::Str, offsetof(PDB::Sheet, atoms[1].residue.name)},
    {Kind::Char, offsetof(PDB::Sheet, atoms[1].residue.chainId)},
    {Kind::Int, offsetof(PDB::Sheet, atoms[1].residue.seqNum)},
    {Kind::Char, offsetof(PDB::Sheet, atoms[1].residue.insertCode)},
    endField(),
};

const Field kSiteFields[] = {
    {Kind::Int, offsetof(PDB::Site, seqNum)},
    {Kind::Str, offsetof(PDB::Site, id)},
    {Kind::Int, offsetof(PDB::Site, count)},
    {Kind::Str, offsetof(PDB::Site, residues[0].name)},
    {Kind::Char, offsetof(PDB::Site, residues[0].chainId)},
    {Kind::Int, offsetof(PDB::Site, residues[0].seqNum)},
    {Kind::Char, offsetof(PDB::Site, residues[0].insertCode)},
    {Kind::Str, offsetof(PDB::Site, residues[1].name)},
    {Kind::Char, offsetof(PDB::Site, residues[1].chainId)},
    {Kind::Int, offsetof(PDB::Site, residues[1].seqNum)},
    {Kind::Char, offsetof(PDB::Site, residues[1].insertCode)},
    {Kind::Str, offsetof(PDB::Site, residues[2].name)},
    {Kind::Char, offsetof(PDB::Site, residues[2].chainId)},
    {Kind::Int, offsetof(PDB::Site, residues[2].seqNum)},
    {Kind::Char, offsetof(PDB::Site, residues[2].insertCode)},
    {Kind::Str, offsetof(PDB::Site, residues[3].name)},
    {Kind::Char, offsetof(PDB::Site, residues[3].chainId)},
    {Kind::Int, offsetof(PDB::Site, residues[3].seqNum)},
    {Kind::Char, offsetof(PDB::Site, residues[3].insertCode)},
    endField(),
};

const Field kSsbondFields[] = {
    {Kind::Int, offsetof(PDB::Ssbond, seqNum)},
    {Kind::Str, offsetof(PDB::Ssbond, residues[0].name)},
    {Kind::Char, offsetof(PDB::Ssbond, residues[0].chainId)},
    {Kind::Int, offsetof(PDB::Ssbond, residues[0].seqNum)},
    {Kind::Char, offsetof(PDB::Ssbond, residues[0].insertCode)},
    {Kind::Str, offsetof(PDB::Ssbond, residues[1].name)},
    {Kind::Char, offsetof(PDB::Ssbond, residues[1].chainId)},
    {Kind::Int, offsetof(PDB::Ssbond, residues[1].seqNum)},
    {Kind::Char, offsetof(PDB::Ssbond, residues[1].insertCode)},
    {Kind::Int, offsetof(PDB::Ssbond, sym1)},
    {Kind::Int, offsetof(PDB::Ssbond, sym2)},
    {Kind::Real, offsetof(PDB::Ssbond, length)},
    endField(),
};

const Field kTerFields[] = {
    {Kind::Int, offsetof(PDB::Ter, serialNum)},
    {Kind::Str, offsetof(PDB::Ter, residue.name)},
    {Kind::Char, offsetof(PDB::Ter, residue.chainId)},
    {Kind::Int, offsetof(PDB::Ter, residue.seqNum)},
    {Kind::Char, offsetof(PDB::Ter, residue.insertCode)},
    endField(),
};

const Field kTurnFields[] = {
    {Kind::Int, offsetof(PDB::Turn, seqNum)},
    {Kind::Str, offsetof(PDB::Turn, id)},
    {Kind::Str, offsetof(PDB::Turn, residues[0].name)},
    {Kind::Char, offsetof(PDB::Turn, residues[0].chainId)},
    {Kind::Int, offsetof(PDB::Turn, residues[0].seqNum)},
    {Kind::Char, offsetof(PDB::Turn, residues[0].insertCode)},
    {Kind::Str, offsetof(PDB::Turn, residues[1].name)},
    {Kind::Char, offsetof(PDB::Turn, residues[1].chainId)},
    {Kind::Int, offsetof(PDB::Turn, residues[1].seqNum)},
    {Kind::Char, offsetof(PDB::Turn, residues[1].insertCode)},
    {Kind::Str, offsetof(PDB::Turn, comment)},
    endField(),
};

const Field kCmpontFields[] = {
    {Kind::Int, offsetof(PDB::Cmpont, seqNum)},
    {Kind::Str, offsetof(PDB::Cmpont, residues[0].name)},
    {Kind::Char, offsetof(PDB::Cmpont, residues[0].chainId)},
    {Kind::Int, offsetof(PDB::Cmpont, residues[0].seqNum)},
    {Kind::Char, offsetof(PDB::Cmpont, residues[0].insertCode)},
    {Kind::Str, offsetof(PDB::Cmpont, residues[1].name)},
    {Kind::Char, offsetof(PDB::Cmpont, residues[1].chainId)},
    {Kind::Int, offsetof(PDB::Cmpont, residues[1].seqNum)},
    {Kind::Char, offsetof(PDB::Cmpont, residues[1].insertCode)},
    endField(),
};

const Field kMasterFields[] = {
    {Kind::Int, offsetof(PDB::Master, numRemark)},
    {Kind::Int, offsetof(PDB::Master, numFtnote)},
    {Kind::Int, offsetof(PDB::Master, numHet)},
    {Kind::Int, offsetof(PDB::Master, numHelix)},
    {Kind::Int, offsetof(PDB::Master, numSheet)},
    {Kind::Int, offsetof(PDB::Master, numTurn)},
    {Kind::Int, offsetof(PDB::Master, numSite)},
    {Kind::Int, offsetof(PDB::Master, numTransform)},
    {Kind::Int, offsetof(PDB::Master, numCoordinate)},
    {Kind::Int, offsetof(PDB::Master, numTer)},
    {Kind::Int, offsetof(PDB::Master, numConect)},
    {Kind::Int, offsetof(PDB::Master, numSeqres)},
    endField(),
};

const Field kObslteFields[] = {
    {Kind::Int, offsetof(PDB::Obslte, continuation)},
    {Kind::Str, offsetof(PDB::Obslte, timestamp)},
    {Kind::Str, offsetof(PDB::Obslte, oldId)},
    {Kind::Str, offsetof(PDB::Obslte, idMap[0])},
    {Kind::Str, offsetof(PDB::Obslte, idMap[1])},
    {Kind::Str, offsetof(PDB::Obslte, idMap[2])},
    {Kind::Str, offsetof(PDB::Obslte, idMap[3])},
    {Kind::Str, offsetof(PDB::Obslte, idMap[4])},
    {Kind::Str, offsetof(PDB::Obslte, idMap[5])},
    {Kind::Str, offsetof(PDB::Obslte, idMap[6])},
    {Kind::Str, offsetof(PDB::Obslte, idMap[7])},
    {Kind::Str, offsetof(PDB::Obslte, idMap[8])},
    endField(),
};

const Field kRevdatFields[] = {
    {Kind::Int, offsetof(PDB::Revdat, modification)},
    {Kind::Int, offsetof(PDB::Revdat, continuation)},
    {Kind::Str, offsetof(PDB::Revdat, timestamp)},
    {Kind::Str, offsetof(PDB::Revdat, id)},
    {Kind::Int, offsetof(PDB::Revdat, modType)},
    {Kind::Str, offsetof(PDB::Revdat, record[0])},
    {Kind::Str, offsetof(PDB::Revdat, record[1])},
    {Kind::Str, offsetof(PDB::Revdat, record[2])},
    {Kind::Str, offsetof(PDB::Revdat, record[3])},
    endField(),
};

const Field kSeqresFields[] = {
    {Kind::Int, offsetof(PDB::Seqres, serialNum)},
    {Kind::Char, offsetof(PDB::Seqres, chainId)},
    {Kind::Int, offsetof(PDB::Seqres, count)},
    {Kind::Str, offsetof(PDB::Seqres, names[0])},
    {Kind::Str, offsetof(PDB::Seqres, names[1])},
    {Kind::Str, offsetof(PDB::Seqres, names[2])},
    {Kind::Str, offsetof(PDB::Seqres, names[3])},
    {Kind::Str, offsetof(PDB::Seqres, names[4])},
    {Kind::Str, offsetof(PDB::Seqres, names[5])},
    {Kind::Str, offsetof(PDB::Seqres, names[6])},
    {Kind::Str, offsetof(PDB::Seqres, names[7])},
    {Kind::Str, offsetof(PDB::Seqres, names[8])},
    {Kind::Str, offsetof(PDB::Seqres, names[9])},
    {Kind::Str, offsetof(PDB::Seqres, names[10])},
    {Kind::Str, offsetof(PDB::Seqres, names[11])},
    {Kind::Str, offsetof(PDB::Seqres, names[12])},
    endField(),
};

const Field kSplitFields[] = {
    {Kind::Int, offsetof(PDB::Split, continuation)},
    {Kind::Str, offsetof(PDB::Split, idMap[0])},
    {Kind::Str, offsetof(PDB::Split, idMap[1])},
    {Kind::Str, offsetof(PDB::Split, idMap[2])},
    {Kind::Str, offsetof(PDB::Split, idMap[3])},
    {Kind::Str, offsetof(PDB::Split, idMap[4])},
    {Kind::Str, offsetof(PDB::Split, idMap[5])},
    {Kind::Str, offsetof(PDB::Split, idMap[6])},
    {Kind::Str, offsetof(PDB::Split, idMap[7])},
    {Kind::Str, offsetof(PDB::Split, idMap[8])},
    {Kind::Str, offsetof(PDB::Split, idMap[9])},
    {Kind::Str, offsetof(PDB::Split, idMap[10])},
    {Kind::Str, offsetof(PDB::Split, idMap[11])},
    {Kind::Str, offsetof(PDB::Split, idMap[12])},
    {Kind::Str, offsetof(PDB::Split, idMap[13])},
    endField(),
};

const Field kSprsdeFields[] = {
    {Kind::Int, offsetof(PDB::Sprsde, continuation)},
    {Kind::Str, offsetof(PDB::Sprsde, timestamp)},
    {Kind::Str, offsetof(PDB::Sprsde, id)},
    {Kind::Str, offsetof(PDB::Sprsde, supersede[0])},
    {Kind::Str, offsetof(PDB::Sprsde, supersede[1])},
    {Kind::Str, offsetof(PDB::Sprsde, supersede[2])},
    {Kind::Str, offsetof(PDB::Sprsde, supersede[3])},
    {Kind::Str, offsetof(PDB::Sprsde, supersede[4])},
    {Kind::Str, offsetof(PDB::Sprsde, supersede[5])},
    {Kind::Str, offsetof(PDB::Sprsde, supersede[6])},
    {Kind::Str, offsetof(PDB::Sprsde, supersede[7])},
    {Kind::Str, offsetof(PDB::Sprsde, supersede[8])},
    endField(),
};

const Field kTvectFields[] = {
    {Kind::Int, offsetof(PDB::Tvect, serialNum)},
    {Kind::Real, offsetof(PDB::Tvect, t1)},
    {Kind::Real, offsetof(PDB::Tvect, t2)},
    {Kind::Real, offsetof(PDB::Tvect, t3)},
    {Kind::Str, offsetof(PDB::Tvect, comment)},
    endField(),
};

const Field kUserPdbrunFields[] = {
    {Kind::Int, offsetof(PDB::UserPdbrun, version)},
    endField(),
};

const Field kUserEyePosFields[] = {
    {Kind::Real, offsetof(PDB::UserEyePos, xyz[0])},
    {Kind::Real, offsetof(PDB::UserEyePos, xyz[1])},
    {Kind::Real, offsetof(PDB::UserEyePos, xyz[2])},
    endField(),
};

const Field kUserAtPosFields[] = {
    {Kind::Real, offsetof(PDB::UserAtPos, xyz[0])},
    {Kind::Real, offsetof(PDB::UserAtPos, xyz[1])},
    {Kind::Real, offsetof(PDB::UserAtPos, xyz[2])},
    endField(),
};

const Field kUserWindowFields[] = {
    {Kind::Real, offsetof(PDB::UserWindow, left)},
    {Kind::Real, offsetof(PDB::UserWindow, right)},
    {Kind::Real, offsetof(PDB::UserWindow, bottom)},
    {Kind::Real, offsetof(PDB::UserWindow, top)},
    {Kind::Real, offsetof(PDB::UserWindow, hither)},
    {Kind::Real, offsetof(PDB::UserWindow, yon)},
    endField(),
};

const Field kUserFocusFields[] = {
    {Kind::Real, offsetof(PDB::UserFocus, focus)},
    endField(),
};

const Field kUserViewportFields[] = {
    {Kind::Real, offsetof(PDB::UserViewport, xmin)},
    {Kind::Real, offsetof(PDB::UserViewport, xmax)},
    {Kind::Real, offsetof(PDB::UserViewport, ymin)},
    {Kind::Real, offsetof(PDB::UserViewport, ymax)},
    endField(),
};

const Field kUserBgColorFields[] = {
    {Kind::Real, offsetof(PDB::UserBgColor, rgb[0])},
    {Kind::Real, offsetof(PDB::UserBgColor, rgb[1])},
    {Kind::Real, offsetof(PDB::UserBgColor, rgb[2])},
    endField(),
};

const Field kUserAngle5Fields[] = {
    {Kind::Int, offsetof(PDB::UserAngle, which)},
    {Kind::Int, offsetof(PDB::UserAngle, atom0)},
    {Kind::Int, offsetof(PDB::UserAngle, atom1)},
    {Kind::Int, offsetof(PDB::UserAngle, atom2)},
    {Kind::Int, offsetof(PDB::UserAngle, atom3)},
    {Kind::Real, offsetof(PDB::UserAngle, angle)},
    endField(),
};

const Field kUserAngle6Fields[] = {
    {Kind::Int, offsetof(PDB::UserAngle, atom0)},
    {Kind::Int, offsetof(PDB::UserAngle, atom1)},
    {Kind::Int, offsetof(PDB::UserAngle, atom2)},
    {Kind::Int, offsetof(PDB::UserAngle, atom3)},
    {Kind::Real, offsetof(PDB::UserAngle, angle)},
    endField(),
};

const Field kUserDistance5Fields[] = {
    {Kind::Int, offsetof(PDB::UserDistance, which)},
    {Kind::Int, offsetof(PDB::UserDistance, atom0)},
    {Kind::Int, offsetof(PDB::UserDistance, atom1)},
    {Kind::Real, offsetof(PDB::UserDistance, distance)},
    endField(),
};

const Field kUserDistance6Fields[] = {
    {Kind::Int, offsetof(PDB::UserDistance, atom0)},
    {Kind::Int, offsetof(PDB::UserDistance, atom1)},
    {Kind::Real, offsetof(PDB::UserDistance, distance)},
    endField(),
};

const Field kUserFile5Fields[] = {
    {Kind::Str, offsetof(PDB::UserFile, filename)},
    endField(),
};

const Field kUserFile6Fields[] = {
    {Kind::Int, offsetof(PDB::UserFile, model)},
    {Kind::Str, offsetof(PDB::UserFile, filename)},
    endField(),
};

const Field kUserMarknameFields[] = {
    {Kind::Str, offsetof(PDB::UserMarkname, markname)},
    endField(),
};

const Field kUserMarkFields[] = {
    {Kind::Str, offsetof(PDB::UserMark, markname)},
    endField(),
};

const Field kUserCName6Fields[] = {
    {Kind::Real, offsetof(PDB::UserCName, rgb[0])},
    {Kind::Real, offsetof(PDB::UserCName, rgb[1])},
    {Kind::Real, offsetof(PDB::UserCName, rgb[2])},
    {Kind::Str, offsetof(PDB::UserCName, name)},
    endField(),
};

const Field kUserColor6Fields[] = {
    {Kind::Real, offsetof(PDB::UserColor, rgb[0])},
    {Kind::Real, offsetof(PDB::UserColor, rgb[1])},
    {Kind::Real, offsetof(PDB::UserColor, rgb[2])},
    {Kind::Str, offsetof(PDB::UserColor, spec)},
    endField(),
};

const Field kUserRadiusFields[] = {
    {Kind::Real, offsetof(PDB::UserRadius, radius)},
    endField(),
};

const Field kUserObjectFields[] = {
    {Kind::Int, offsetof(PDB::UserObject, model)},
    endField(),
};

const Field kUserEndObjFields[] = {
    {Kind::Int, offsetof(PDB::UserEndObj, model)},
    endField(),
};

const Field kUserChainFields[] = {
    {Kind::Int, offsetof(PDB::UserChain, atom0)},
    {Kind::Int, offsetof(PDB::UserChain, atom1)},
    endField(),
};

const Field kUserGfxBeginFields[] = {
    {Kind::Str, offsetof(PDB::UserGfxBegin, unknown)},
    endField(),
};

const Field kUserGfxColor6Fields[] = {
    {Kind::Real, offsetof(PDB::UserGfxColor, rgb[0])},
    {Kind::Real, offsetof(PDB::UserGfxColor, rgb[1])},
    {Kind::Real, offsetof(PDB::UserGfxColor, rgb[2])},
    {Kind::Str, offsetof(PDB::UserGfxColor, spec)},
    endField(),
};

const Field kUserGfxNormalFields[] = {
    {Kind::Real, offsetof(PDB::UserGfxNormal, xyz[0])},
    {Kind::Real, offsetof(PDB::UserGfxNormal, xyz[1])},
    {Kind::Real, offsetof(PDB::UserGfxNormal, xyz[2])},
    endField(),
};

const Field kUserGfxVertexFields[] = {
    {Kind::Real, offsetof(PDB::UserGfxVertex, xyz[0])},
    {Kind::Real, offsetof(PDB::UserGfxVertex, xyz[1])},
    {Kind::Real, offsetof(PDB::UserGfxVertex, xyz[2])},
    endField(),
};

const Field kUserGfxFont6Fields[] = {
    {Kind::Int, offsetof(PDB::UserGfxFont, size)},
    {Kind::Str, offsetof(PDB::UserGfxFont, name)},
    endField(),
};

const Field kUserGfxTextPosFields[] = {
    {Kind::Real, offsetof(PDB::UserGfxTextPos, xyz[0])},
    {Kind::Real, offsetof(PDB::UserGfxTextPos, xyz[1])},
    {Kind::Real, offsetof(PDB::UserGfxTextPos, xyz[2])},
    endField(),
};

const Field kUserGfxLabel6Fields[] = {
    {Kind::Str, offsetof(PDB::UserGfxLabel, text)},
    endField(),
};

const Field kUserGfxMoveFields[] = {
    {Kind::Real, offsetof(PDB::UserGfxMove, xyz[0])},
    {Kind::Real, offsetof(PDB::UserGfxMove, xyz[1])},
    {Kind::Real, offsetof(PDB::UserGfxMove, xyz[2])},
    endField(),
};

const Field kUserGfxDrawFields[] = {
    {Kind::Real, offsetof(PDB::UserGfxDraw, xyz[0])},
    {Kind::Real, offsetof(PDB::UserGfxDraw, xyz[1])},
    {Kind::Real, offsetof(PDB::UserGfxDraw, xyz[2])},
    endField(),
};

const Field kUserGfxMarkerFields[] = {
    {Kind::Real, offsetof(PDB::UserGfxMarker, xyz[0])},
    {Kind::Real, offsetof(PDB::UserGfxMarker, xyz[1])},
    {Kind::Real, offsetof(PDB::UserGfxMarker, xyz[2])},
    endField(),
};

const Field kUserGfxPointFields[] = {
    {Kind::Real, offsetof(PDB::UserGfxPoint, xyz[0])},
    {Kind::Real, offsetof(PDB::UserGfxPoint, xyz[1])},
    {Kind::Real, offsetof(PDB::UserGfxPoint, xyz[2])},
    endField(),
};

// Per-record field lists. Offsets are byte offsets within the record
// storage (every parsed record type keeps its fields inside one union
// member struct, so struct member offsets start at the storage base).
// Each list pairs one-to-one with the value conversions of its record's
// read format and ends with an End sentinel.

// Format string tables (upstream .i data, verbatim).
char const *const pdbRecordFormat[PDB::NUM_TYPES] = {
#include "read_format.i"
};

char const *const pdbrun5Read[PDB::NUM_USER] = {
#include "pdbrun5_read.i"
};

char const *const pdbrun6Read[PDB::NUM_USER] = {
#include "pdbrun6_read.i"
};

// Field lists per standard record type; records without a reader carry
// nullptr (UNKNOWN junk, END, and ENDMDL are handled without parsing).
const Field *const kStandardFields[PDB::NUM_TYPES] = {
    /* UNKNOWN */ nullptr,
    /* ANISOU */ kAnisouFields,
    /* ATOM */ kAtomFields,
    /* AUTHOR */ kAuthorFields,
    /* CAVEAT */ kCaveatFields,
    /* COMPND */ kCompndFields,
    /* CONECT */ kConectFields,
    /* CRYST1 */ kCryst1Fields,
    /* END */ nullptr,
    /* FORMUL */ kFormulFields,
    /* FTNOTE */ kAgrdesFields,
    /* HEADER */ kHeaderFields,
    /* HELIX */ kHelixFields,
    /* HET */ kHetFields,
    /* HETATM */ kAtomFields,
    /* JRNL */ kJrnlFields,
    /* KEYWDS */ kKeywdsFields,
    /* MASTER */ kMasterFields,
    /* MDLTYP */ kMdltypFields,
    /* MTRIX */ kMtrixFields,
    /* NUMMDL */ kNummdlFields,
    /* OBSLTE */ kObslteFields,
    /* ORIGX */ kOrigxFields,
    /* REMARK */ kAgrdesFields,
    /* REVDAT */ kRevdatFields,
    /* SCALE */ kScaleFields,
    /* SEQRES */ kSeqresFields,
    /* SHEET */ kSheetFields,
    /* SIGATM */ kAtomFields,
    /* SIGUIJ */ kAnisouFields,
    /* SITE */ kSiteFields,
    /* SOURCE */ kSourceFields,
    /* SPLIT */ kSplitFields,
    /* SPRSDE */ kSprsdeFields,
    /* SSBOND */ kSsbondFields,
    /* TER */ kTerFields,
    /* TITLE */ kTitleFields,
    /* TURN */ kTurnFields,
    /* TVECT */ kTvectFields,
    /* USER */ kUserFields,
    /* MODEL */ kModelFields,
    /* ENDMDL */ nullptr,
    /* EXPDTA */ kExpdtaFields,
    /* SYMDES */ kAgrdesFields,
    /* SYMOP */ kSymopFields,
    /* MTXDES */ kAgrdesFields,
    /* CMPDES */ kAgrdesFields,
    /* CMPONT */ kCmpontFields,
    /* TRNSFM */ kTrnsfmFields,
    /* AGRDES */ kAgrdesFields,
    /* AGGRGT */ kAggrgtFields,
    /* ZDATOM */ kZdatomFields,
};

// PDBRUN 5/6 USER paths. fallbackToUser mirrors the legacy goto user;
// entries without a reader carry no fields and fall back directly.
struct UserReadSpec {
  bool libcBackend;
  bool fallbackToUser;
  const Field *fields;
};

const UserReadSpec kPdbrun5Read[PDB::NUM_USER] = {
    /* PDBRUN */ {false, true, kUserPdbrunFields},
    /* EYEPOS */ {false, true, kUserEyePosFields},
    /* ATPOS */ {false, true, kUserAtPosFields},
    /* WINDOW */ {false, true, kUserWindowFields},
    /* FOCUS */ {false, true, kUserFocusFields},
    /* VIEWPORT */ {false, true, kUserViewportFields},
    /* BGCOLOR */ {true, true, nullptr},
    /* ANGLE */ {false, true, kUserAngle5Fields},
    /* DISTANCE */ {false, true, kUserDistance5Fields},
    /* FILE */ {false, true, kUserFile5Fields},
    /* MARKNAME */ {false, true, nullptr},
    /* MARK */ {false, true, nullptr},
    /* CNAME */ {true, true, nullptr},
    /* COLOR */ {true, true, nullptr},
    /* RADIUS */ {false, true, kUserRadiusFields},
    /* OBJECT */ {false, true, kUserObjectFields},
    /* ENDOBJ */ {false, true, kUserEndObjFields},
    /* CHAIN */ {true, true, nullptr},
    /* GFX BEGIN */ {false, true, nullptr},
    /* GFX END */ {false, true, nullptr},
    /* GFX COLOR */ {true, true, nullptr},
    /* GFX NORMAL */ {false, true, nullptr},
    /* GFX VERTEX */ {false, true, nullptr},
    /* GFX FONT */ {true, true, nullptr},
    /* GFX TEXTPOS */ {false, true, nullptr},
    /* GFX LABEL */ {true, true, nullptr},
    /* GFX MOVE */ {false, true, kUserGfxMoveFields},
    /* GFX DRAW */ {false, true, kUserGfxDrawFields},
    /* GFX MARKER */ {false, true, kUserGfxMarkerFields},
    /* GFX POINT */ {false, true, kUserGfxPointFields},
};

const UserReadSpec kPdbrun6Read[PDB::NUM_USER] = {
    /* PDBRUN */ {false, true, kUserPdbrunFields},
    /* EYEPOS */ {false, true, kUserEyePosFields},
    /* ATPOS */ {false, true, kUserAtPosFields},
    /* WINDOW */ {false, true, kUserWindowFields},
    /* FOCUS */ {false, true, kUserFocusFields},
    /* VIEWPORT */ {false, true, kUserViewportFields},
    /* BGCOLOR */ {false, true, kUserBgColorFields},
    /* ANGLE */ {false, true, kUserAngle6Fields},
    /* DISTANCE */ {false, true, kUserDistance6Fields},
    /* FILE */ {false, true, kUserFile6Fields},
    /* MARKNAME */ {false, true, kUserMarknameFields},
    /* MARK */ {false, true, kUserMarkFields},
    /* CNAME */ {false, true, kUserCName6Fields},
    /* COLOR */ {false, true, kUserColor6Fields},
    /* RADIUS */ {false, true, kUserRadiusFields},
    /* OBJECT */ {false, false, nullptr},
    /* ENDOBJ */ {false, false, nullptr},
    /* CHAIN */ {false, true, kUserChainFields},
    /* GFX BEGIN */ {false, true, kUserGfxBeginFields},
    /* GFX END */ {false, false, nullptr},
    /* GFX COLOR */ {false, true, kUserGfxColor6Fields},
    /* GFX NORMAL */ {false, true, kUserGfxNormalFields},
    /* GFX VERTEX */ {false, true, kUserGfxVertexFields},
    /* GFX FONT */ {false, true, kUserGfxFont6Fields},
    /* GFX TEXTPOS */ {false, true, kUserGfxTextPosFields},
    /* GFX LABEL */ {false, true, kUserGfxLabel6Fields},
    /* GFX MOVE */ {false, true, nullptr},
    /* GFX DRAW */ {false, true, nullptr},
    /* GFX MARKER */ {false, true, nullptr},
    /* GFX POINT */ {false, true, nullptr},
};

enum class RunResult { Parsed, Skipped, Failed };

template <typename T>
T *fieldPtr(char *base, const Field &field) {
  return reinterpret_cast<T *>(base + field.offset);
}

// Custom fixed-width backend: pairs the format's value conversions with
// the spec's field list, with the PDB::sscanf field semantics.
RunResult runCustom(PDB &record, const char *buf, const char *fmt,
                    const Field *fields) {
  char *const base = reinterpret_cast<char *>(&record.unknown);
  FieldReader reader(buf);
  for (; *fmt != '\0'; fmt++) {
    if (*fmt != '%') {
      if (!reader.matchLiteral(*fmt))
        return RunResult::Failed;
      continue;
    }

    int width = 0;
    for (fmt++; isdigit(static_cast<unsigned char>(*fmt)); fmt++)
      width = width * 10 + *fmt - '0';
    if (width == 0)
      width = 1;

    switch (*fmt) {
    case 'd':
      if (fields->kind == Kind::End)
        return RunResult::Failed;
      if (!reader.readInt(width, fieldPtr<int>(base, *fields++)))
        return RunResult::Failed;
      break;
    case 'f':
      if (fields->kind == Kind::End)
        return RunResult::Failed;
      if (!reader.readReal(width, fieldPtr<double>(base, *fields++)))
        return RunResult::Failed;
      break;
    case 's':
      if (fields->kind == Kind::End)
        return RunResult::Failed;
      reader.readString(width, fieldPtr<char>(base, *fields++));
      break;
    case 'c':
      if (fields->kind == Kind::End)
        return RunResult::Failed;
      reader.readChars(width, fieldPtr<char>(base, *fields++));
      break;
    case ' ':
      reader.skipChars(width);
      break;
    default:
      return RunResult::Failed;
    }
  }
  return RunResult::Parsed;
}

// v5 USER paths whose format strings use libc scanf conversions.
RunResult runLibc(PDB &record, const char *buf, const char *fmt,
                  PDB::RecordType type) {
  switch (type) {
  case PDB::USER_BGCOLOR:
    if (std::sscanf(buf, fmt, &record.userBgColor.rgb[0],
                    &record.userBgColor.rgb[1],
                    &record.userBgColor.rgb[2]) < 0)
      return RunResult::Failed;
    return RunResult::Parsed;
  case PDB::USER_CNAME:
    if (std::sscanf(buf, fmt, record.userCName.name, &record.userCName.rgb[0],
                    &record.userCName.rgb[1], &record.userCName.rgb[2]) < 0)
      return RunResult::Failed;
    return RunResult::Parsed;
  case PDB::USER_COLOR:
    if (std::sscanf(buf, fmt, record.userColor.spec, &record.userColor.rgb[0],
                    &record.userColor.rgb[1], &record.userColor.rgb[2]) < 0)
      return RunResult::Failed;
    return RunResult::Parsed;
  case PDB::USER_CHAIN:
    if (std::sscanf(buf, fmt, &record.userChain.atom0,
                    &record.userChain.atom1) < 0)
      return RunResult::Failed;
    return RunResult::Parsed;
  case PDB::USER_GFX_COLOR:
    if (std::sscanf(buf, fmt, record.userGfxColor.spec,
                    &record.userGfxColor.rgb[0], &record.userGfxColor.rgb[1],
                    &record.userGfxColor.rgb[2]) < 0)
      return RunResult::Failed;
    return RunResult::Parsed;
  case PDB::USER_GFX_FONT:
    if (std::sscanf(buf, fmt, record.userGfxFont.name,
                    &record.userGfxFont.size) < 0)
      return RunResult::Failed;
    return RunResult::Parsed;
  case PDB::USER_GFX_LABEL:
    if (std::sscanf(buf, fmt, &record.userGfxLabel.xyz[0],
                    &record.userGfxLabel.xyz[1], &record.userGfxLabel.xyz[2],
                    record.userGfxLabel.text) < 0)
      return RunResult::Failed;
    return RunResult::Parsed;
  default:
    return RunResult::Failed;
  }
}

RunResult run(PDB &record, const char *buf, const PdbFormat::ReadSpec &spec,
              PDB::RecordType type) {
  if (spec.fmt == nullptr)
    return RunResult::Skipped;
  if (spec.libcBackend)
    return runLibc(record, buf, spec.fmt, type);
  return runCustom(record, buf, spec.fmt, spec.fields);
}

void setUnknown(PDB &record, const char *buf) {
  (void)PDB::sprintf(record.unknown.junk, "%72s", buf);
}

} // namespace

PdbFormat::ReadSpec PdbFormat::readSpec(PDB::RecordType type,
                                        int pdbrunInputVersion) {
  if (type < PDB::USER_PDBRUN)
    return {pdbRecordFormat[type], false, false, kStandardFields[type]};
  const int index = type - PDB::USER_PDBRUN;
  const bool v5 = pdbrunInputVersion < 6;
  const UserReadSpec &spec = (v5 ? kPdbrun5Read : kPdbrun6Read)[index];
  return {(v5 ? pdbrun5Read : pdbrun6Read)[index], spec.libcBackend,
          spec.fallbackToUser, spec.fields};
}

PDB::RecordType PdbFormat::parse(PDB &record, const char *buf,
                                 PDB::RecordType type,
                                 int pdbrunInputVersion) {
  if (type == PDB::UNKNOWN) {
    setUnknown(record, buf);
    return PDB::UNKNOWN;
  }

  if (type == PDB::ATOM || type == PDB::HETATM || type == PDB::SIGATM) {
    if (run(record, buf, readSpec(type, pdbrunInputVersion), type) !=
        RunResult::Failed)
      return type;
    if (run(record, buf, readSpec(PDB::ZDATOM, pdbrunInputVersion),
            PDB::ZDATOM) == RunResult::Parsed) {
      record.atom.iszdatom = true;
      return type;
    }
    setUnknown(record, buf);
    return PDB::UNKNOWN;
  }

  const ReadSpec spec = readSpec(type, pdbrunInputVersion);
  const RunResult result = run(record, buf, spec, type);
  if (result == RunResult::Failed) {
    if (spec.fallbackToUser &&
        run(record, buf, readSpec(PDB::USER, pdbrunInputVersion),
            PDB::USER) != RunResult::Failed)
      return PDB::USER;
    setUnknown(record, buf);
    return PDB::UNKNOWN;
  }
  if (result == RunResult::Parsed) {
    if (type == PDB::USER_PDBRUN)
      PDB::PdbrunInputVersion(record.userPdbrun.version);
    else if (type == PDB::USER_GFX_BEGIN)
      record.userGfxBegin.primitive =
          RecordClassifier::getGfxType(record.userGfxBegin.unknown);
  }
  return type;
}

} // namespace libpdb
