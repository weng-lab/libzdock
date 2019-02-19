/**
 * Copyright 2019 Arjan van der Velde, vandervelde.ag [at] gmail
 * Licensed under the Apache License, Version 2.0 (the "License");
 * you may not use this file except in compliance with the License.
 * You may obtain a copy of the License at
 * http://www.apache.org/licenses/LICENSE-2.0
 *
 * Unless required by applicable law or agreed to in writing, software
 * distributed under the License is distributed on an "AS IS" BASIS,
 * WITHOUT WARRANTIES OR CONDITIONS OF ANY KIND, either express or implied.
 * See the License for the specific language governing permissions and
 * limitations under the License.
 */

#include "PDB.hpp"
#include "Exception.hpp"
#include "Utils.hpp"
#include <algorithm>
#include <fstream>

namespace p = ::libpdb;
namespace e = ::Eigen;

namespace zdock {

PDB::PDB(const PDB &p) {
  models_ = p.models_;
  records_ = p.records_;
  atoms_ = p.atoms_;
  matrix_ = p.matrix_;
};

PDB::PDB(const std::string &filename) { read_(filename); }

PDB &PDB::operator=(const PDB &p) {
  models_ = p.models_;
  records_ = p.records_;
  atoms_ = p.atoms_;
  matrix_ = p.matrix_;
  return *this;
}

void PDB::read_(const std::string &fn) {
  p::PDB record;
  size_t m = 0;
  std::ifstream infile(fn);
  if (infile.is_open()) {
    records_.clear();
    while (infile >> record) {
      switch (record.type()) {
      case p::PDB::UNKNOWN:
        break; // ignore unknown
      case p::PDB::MODEL:
        m = record.model.num;
        while (models_.size() < m) {
          // TODO: probably make this a map instead...
          models_.push_back(std::make_shared<zdock::Model>());
        }
        break;
      case p::PDB::ENDMDL:
        m = 0; // ground level
        break;
      default:
        break;
      }
      if (p::PDB::UNKNOWN != record.type()) {
        append(record, m);
      }
    }
  } else {
    throw PDBOpenException(fn);
  }
}

void PDB::append(const libpdb::PDB &record, const int model) {
  append(std::make_shared<libpdb::PDB>(record), model);
}

void PDB::append(const Record &r, const int model) {
  std::lock_guard<std::mutex> lock(lock_);
  switch (r->type()) {
  case p::PDB::UNKNOWN:
    break; // silently drop 'UNKNOWN' type records
  case p::PDB::ATOM:
  case p::PDB::HETATM:
    if (0 == model) {
      atoms_.push_back(r);
      matrix_.conservativeResize(matrix_.rows(), matrix_.cols() + 1);
      matrix_.col(matrix_.cols() - 1) = e::Vector3d(r->atom.xyz);
    } else {
      models_[model - 1]->append(r);
      models_[model - 1]->modelNum_ = model;
    }
    break;
  default:
    break;
  }
  if (p::PDB::UNKNOWN != r->type()) {
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
    assert(m.cols() == static_cast<long>(atoms_.size()));
    for (size_t i = 0; i < atoms_.size(); ++i) {
      atoms_[i]->atom.xyz[0] = m(0, i);
      atoms_[i]->atom.xyz[1] = m(1, i);
      atoms_[i]->atom.xyz[2] = m(2, i);
    }
    return matrix_;
  }
}

const std::vector<PDB::Model> &PDB::models() const { return models_; }

size_t PDB::nmodels() const { return models_.size(); }

const std::vector<PDB::Record> &PDB::records() const { return records_; }

const std::vector<PDB::Record> &PDB::atoms() const { return atoms_; }

const PDB::Record &PDB::operator[](const int serial) const {
  for (const auto &a : atoms_) {
    if (serial == a->atom.serialNum) {
      return a;
    }
  }
  throw AtomNotFoundException("Atom not found.");
}

const PDB::Record &PDB::operator[](const Coord &coord) const {
  const Record &x = (*this)[coord.serialNum];
  if (Utils::trim_copy(x->atom.name) == coord.atomName &&
      Utils::trim_copy(x->atom.residue.name) == coord.resName &&
      x->atom.residue.chainId == coord.chain &&
      x->atom.residue.seqNum == coord.resNum) {
    return x;
  } else {
    throw AtomNotFoundException("Atom not found.");
  }
}

} // namespace zdock
