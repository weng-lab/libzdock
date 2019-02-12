#include "PDB.hpp"
#include "Exception.hpp"
#include <fstream>

namespace p = ::libpdb;
namespace e = ::Eigen;

namespace zdock {

PDB::PDB(const std::string &fn, const int model, const bool alpha) {
  read_(fn, model, alpha);
}

const PDB::Matrix &PDB::matrix() const { return m_; }

const std::vector<p::PDB> &PDB::records() { return records_; }

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

void PDB::read_(const std::string &fn, const int model, const bool alpha) {
  p::PDB record;
  int firstmodel = 0;
  int m = 0;
  size_t count = 0;
  std::string title;
  std::ifstream infile(fn);
  if (infile.is_open()) {
    // read pdb lines
    atoms_.clear();
    records_.clear();
    m_.resize(3, 0);
    while (infile >> record) {
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
        if (alpha && !record.isalpha()) {
          continue; // skip non-CA if alpha is set
        }
        if (model == m || MODEL_ALL == model ||
            (MODEL_FIRST == model && firstmodel == m)) {
          atoms_.push_back(count);
        }
        break;

      // parse title
      case p::PDB::TITLE: // capture title
        if (record.title.continuation) {
          title += std::string(record.title.text).substr(1);
        } else {
          title = record.title.text;
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

} // namespace zdock
