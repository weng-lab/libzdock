// Copyright (c) 2019-2026 Arjan van der Velde, Weng Lab
// SPDX-License-Identifier: BSD-2-Clause

#include "PDB.hpp"
#include "Exception.hpp"
#include "Utils.hpp"
#include <algorithm>
#include <fstream>

namespace p = ::zdock::pdb;
namespace e = ::Eigen;

namespace zdock {

PDB::PDB() : filter_([](const zdock::pdb::Record &) { return true; }) {}

PDB::PDB(const PDB &other) {
  RecordMap records;
  copyFrom_(other, records);
}

PDB::PDB(const std::string &filename,
         std::function<bool(const zdock::pdb::Record &)> filter)
    : filter_(filter) {
  read_(filename);
}

PDB &PDB::operator=(const PDB &other) {
  if (this == &other) {
    return *this;
  }

  RecordMap records;
  copyFrom_(other, records);
  return *this;
}

void PDB::copyFrom_(const PDB &other, RecordMap &clonedRecords) {
  const auto cloneRecord = [&clonedRecords](const Record &record) {
    if (!record) {
      return Record();
    }

    const auto existing = clonedRecords.find(record.get());
    if (existing != clonedRecords.end()) {
      return existing->second;
    }

    const Record clone = std::make_shared<zdock::pdb::Record>(*record);
    clonedRecords.emplace(record.get(), clone);
    return clone;
  };

  filter_ = other.filter_;
  models_.clear();
  records_.clear();
  atoms_.clear();
  records_.reserve(other.records_.size());
  for (const auto &record : other.records_) {
    records_.push_back(cloneRecord(record));
  }

  atoms_.reserve(other.atoms_.size());
  for (const auto &atom : other.atoms_) {
    atoms_.push_back(cloneRecord(atom));
  }
  matrix_ = other.matrix_;

  models_.reserve(other.models_.size());
  for (const auto &sourceModel : other.models_) {
    const auto targetModel = std::make_shared<zdock::Model>();
    targetModel->copyFrom_(*sourceModel, clonedRecords);
    targetModel->modelNum_ = sourceModel->modelNum_;
    models_.push_back(targetModel);
  }
}

void PDB::read_(const std::string &fn) {
  size_t m = 0;
  std::ifstream infile(fn);
  if (infile.is_open()) {
    records_.clear();
    std::string line;
    p::ParseContext context;
    // each side reads with the pinned operator>> line semantics; lines
    // that parse as no record at all were stored as UNKNOWN by the
    // legacy reader and are dropped here, and an overlong line fails
    // the stream and stops the loop exactly as it always has
    while (p::readRecordLine(infile, line)) {
      const auto record = p::Record::parse(line, context);
      if (!record.has_value()) {
        continue; // unparsable: dropped, counted on the context
      }
      switch (record->type()) {
      case p::RecordType::MODEL:
        m = record->as<p::ModelRecord>().num;
        while (models_.size() < m) {
          // TODO: probably make this a map instead...
          models_.push_back(std::make_shared<zdock::Model>());
        }
        break;
      case p::RecordType::ENDMDL:
        m = 0; // ground level
        break;
      default:
        break;
      }
      if (p::RecordType::UNKNOWN != record->type()) {
        append(record.value(), m);
      }
    }
  } else {
    throw PDBOpenException(fn);
  }
}

void PDB::append(const p::Record &record, const int model) {
  append(std::make_shared<p::Record>(record), model);
}

void PDB::append(const Record &r, const int model) {
  std::lock_guard<std::mutex> lock(lock_);
  if (model < 0 || static_cast<size_t>(model) > models_.size()) {
    throw Exception("PDB model index is out of range");
  }
  switch (r->type()) {
  case p::RecordType::UNKNOWN:
    break; // silently drop 'UNKNOWN' type records
  case p::RecordType::ATOM:
  case p::RecordType::HETATM:
    if (filter_(*r)) {
      if (0 == model) {
        atoms_.push_back(r);
        matrix_.conservativeResize(matrix_.rows(), matrix_.cols() + 1);
        const auto &atom = r->asAtom();
        matrix_.col(matrix_.cols() - 1) =
            e::Vector3d(atom.xyz[0], atom.xyz[1], atom.xyz[2]);
      } else {
        models_[model - 1]->append(r);
        models_[model - 1]->modelNum_ = model;
      }
    }
    break;
  default:
    break;
  }
  if (p::RecordType::UNKNOWN != r->type()) {
    records_.push_back(r);
  }
}

const PDB::Matrix &PDB::matrix() const {
  if (models_.size() > 0) {
    return models_[0]->matrix(); // first model
  }
  return matrix_; // only model
}

const PDB::Matrix &PDB::setMatrix(const Matrix &m) {
  if (models_.size() > 0) {
    return models_[0]->setMatrix(m); // first model
  } else {                           // only model
    std::lock_guard<std::mutex> lock(lock_);
    if (m.rows() != 3 || m.cols() != static_cast<long>(atoms_.size())) {
      throw Exception("PDB coordinate matrix dimensions do not match its atoms");
    }
    for (size_t i = 0; i < atoms_.size(); ++i) {
      auto &atom = atoms_[i]->asAtom();
      atom.xyz[0] = m(0, i);
      atom.xyz[1] = m(1, i);
      atom.xyz[2] = m(2, i);
    }
    matrix_ = m;
    return matrix_;
  }
}

const std::vector<PDB::Model> &PDB::models() const { return models_; }

size_t PDB::nmodels() const { return models_.size(); }

const std::vector<PDB::Record> &PDB::records() const { return records_; }

const std::vector<PDB::Record> &PDB::atoms() const { return atoms_; }

const PDB::Record &PDB::operator[](const int serial) const {
  for (const auto &a : atoms_) {
    if (serial == a->asAtom().serialNum) {
      return a;
    }
  }
  throw AtomNotFoundException("Atom not found.");
}

const PDB::Record &PDB::operator[](const RecordCoord &coord) const {
  const Record &x = (*this)[coord.serialNum];
  if (Utils::trim_copy(x->asAtom().name) == coord.atomName &&
      Utils::trim_copy(x->asAtom().residue.name) == coord.resName &&
      x->asAtom().residue.chainId == coord.chain &&
      x->asAtom().residue.seqNum == coord.resNum) {
    return x;
  } else {
    throw AtomNotFoundException("Atom not found.");
  }
}

PDB::Coord PDB::centroid() const  {
  if (matrix().cols() == 0) {
    throw Exception("Cannot calculate the centroid of an empty PDB");
  }
  return matrix().rowwise().mean();
}

} // namespace zdock
