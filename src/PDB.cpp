#include "PDB.hpp"
#include "Exception.hpp"
#include "Utils.hpp"
#include <fstream>

namespace p = ::libpdb;
namespace e = ::Eigen;

namespace zdock {

PDB::PDB() {
  m_.resize(3, 0);
}

PDB::PDB(const std::string &fn, const int model,
         std::function<bool(const libpdb::PDB &)> filter) {
  read_(fn, model, filter);
}

const PDB::Matrix &PDB::matrix() const { return m_; }

const std::vector<p::PDB> &PDB::records() { return records_; }

const libpdb::PDB &PDB::operator[](const Coord &coord) const {
  const libpdb::PDB &x = (*this)[coord.serialNum];
  if (Utils::trim_copy(x.atom.name) == coord.atomName &&
      Utils::trim_copy(x.atom.residue.name) == coord.resName &&
      x.atom.residue.chainId == coord.chain &&
      x.atom.residue.seqNum == coord.resNum) {
    return x;
  } else {
    throw AtomNotFoundException("Atom not found.");
  }
}

const libpdb::PDB &PDB::operator[](const int serial) const {
  for (const auto &r : records_) {
    switch (r.type()) {
    case libpdb::PDB::ATOM:
    case libpdb::PDB::HETATM:
      if (serial == r.atom.serialNum) {
        return r;
      }
      break;
    default:
      break;
    }
  }
  throw AtomNotFoundException("Atom not found.");
}

const PDB::Matrix &PDB::transform(const PDB::Transform &t) {
  m_ = t * m_;
  // reconstitute atom records from matrix
  size_t i = 0;
  for (const auto j : atoms_) {
    records_[j].atom.xyz[0] = m_(0, i);
    records_[j].atom.xyz[1] = m_(1, i);
    records_[j].atom.xyz[2] = m_(2, i);
    i++;
  }
  return m_;
}

const PDB::Matrix &PDB::setMatrix(const PDB::Matrix &m) {
  // reconstitute atom records from matrix
  m_ = m;
  size_t i = 0;
  for (const auto j : atoms_) {
    records_[j].atom.xyz[0] = m_(0, i);
    records_[j].atom.xyz[1] = m_(1, i);
    records_[j].atom.xyz[2] = m_(2, i);
    i++;
  }
  return m_;
}

void PDB::read_(const std::string &fn, const int model,
                std::function<bool(const libpdb::PDB &)> filter) {
  p::PDB record;
  int firstmodel = 0;
  int m = 0;
  size_t count = 0;
  std::ifstream infile(fn);
  if (infile.is_open()) {
    std::lock_guard<std::mutex> lock(insertmtx_);
    // read pdb lines
    atoms_.clear();
    records_.clear();
    m_.resize(3, 0);
    while (infile >> record) {
      if (!filter(record)) {
        continue;
      }
      switch (record.type()) {

      // keep track of model number
      case p::PDB::MODEL:
        m = record.model.num;
        if (!firstmodel) {
          if (MODEL_ALL != model) {
            atoms_.clear();
          }
          firstmodel = m;
        }
        break;

      // model dependent capturing of atom/hetatm
      case p::PDB::ATOM:
      case p::PDB::HETATM:
        if (model == m || MODEL_ALL == model ||
            (MODEL_FIRST == model && firstmodel == m)) {
          atoms_.push_back(count);
        }
        break;

      default:
        break;
      }
      if (p::PDB::UNKNOWN != record.type()) {
        records_.push_back(record);
        count++;
      }
    }
    // turn atoms into matrix
    m_.resize(3, atoms_.size());
    size_t col = 0;
    for (const auto x : atoms_) {
      m_.col(col) = e::Vector3d(records_[x].atom.xyz);
      col++;
    }
  } else {
    throw PDBOpenException(fn);
  }
}

void PDB::append(const p::PDB &r) {
  if (p::PDB::UNKNOWN != r.type()) { // silently drop 'UNKNOWN' type records
    std::lock_guard<std::mutex> lock(insertmtx_);
    if (p::PDB::ATOM == r.type() || p::PDB::HETATM == r.type()) {
      m_.conservativeResize(m_.rows(), m_.cols() + 1);
      m_.col(m_.cols() - 1) = e::Vector3d(r.atom.xyz);
      atoms_.push_back(records_.size());
    }
    records_.push_back(r);
  }
}

const e::Vector3d PDB::centroid() const { return m_.rowwise().mean(); }

} // namespace zdock
