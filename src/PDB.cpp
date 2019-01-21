#include "PDB.hpp"
#include "Exception.hpp"
#include <fstream>
#include <iostream>

namespace p = ::libpdb;
namespace e = ::Eigen;

namespace zlab {

PDB::PDB(const std::string &fn) { read_(fn); }

const e::Matrix<double, 3, e::Dynamic> &PDB::matrix() { return m_; }

const std::vector<p::PDB> &PDB::records() { return records_; }

const Eigen::Matrix<double, 3, Eigen::Dynamic> &
PDB::transform(Eigen::Transform<double, 3, Eigen::Affine> t) {
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

void PDB::read_(const std::string &fn) {
  p::PDB record;
  int model = 0;
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
      case p::PDB::MODEL:
        model = record.model.num;
        if (model < 2) { // read first model (clear upon 'model entry')
          atoms_.clear();
        }
        break;
      case p::PDB::ATOM:
      case p::PDB::HETATM:
        if (model < 2) { // read first model (or 0 if none defined)
          atoms_.push_back(count);
        }
        break;
      case p::PDB::TITLE:
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

} // namespace zlab
