#include "pdb++.h"
#include <iostream>
#include <string>
#include <vector>

#define EIGEN_USE_LAPACKE

#include <Eigen/Dense>

namespace eg = ::Eigen;

int main() {
  eg::Matrix<double, 3, eg::Dynamic> m;
  PDB record;
  int model = 0;
  int count = 0;
  std::string title;
  std::vector<size_t> atoms; // record indexes for atoms
  std::vector<PDB> records;
  // read pdb lines
  while (std::cin >> record && model < 2) { // read first model only
    switch (record.type()) {
    case PDB::MODEL:
      if (record.model.num < 2) { // read first model (clear upon 'model entry')
        model = record.model.num;
        atoms.clear();
      }
      break;
    case PDB::ATOM:
    case PDB::HETATM:
      atoms.push_back(count);
      break;
    case PDB::TITLE:
      if (record.title.continuation) {
        title += std::string(record.title.text).substr(1);
      } else {
        title = record.title.text;
      }
      break;
    default:
      break;
    }
    if (PDB::UNKNOWN != record.type()) {
      records.push_back(record);
      count++;
    }
  }

  // turn atoms into matrix
  m.resize(3, atoms.size());
  size_t col = 0;
  for (const auto x : atoms) {
    m.col(col) = eg::Vector3d(records[x].atom.xyz);
    col++;
  }

  // perform a translation
  eg::Vector3d x = -m.rowwise().mean();
  eg::Transform<double, 3, eg::Affine> t =
      eg::Translation3d(-x) *
      (eg::AngleAxisd(0.0 * M_PI, eg::Vector3d::UnitX()) *
       eg::AngleAxisd(0.25 * M_PI, eg::Vector3d::UnitY()) *
       eg::AngleAxisd(0.25 * M_PI, eg::Vector3d::UnitZ())) *
      eg::Translation3d(x);
  m = t * m;

  // reconstitute atom records from matrix
  count = 0;
  for (auto i : atoms) {
    records[i].atom.xyz[0] = m(0, count);
    records[i].atom.xyz[1] = m(1, count);
    records[i].atom.xyz[2] = m(2, count);
    count++;
  }

  // output pdb
  for (const auto &r : records) {
    std::cout << r << std::endl;
  }
}

// clang++ -isystem/opt/local/include/eigen3 -framework Accelerate /opt/local/lib/lapack/liblapacke.dylib  -O3 -L. -lpdb++ -o bla main.cpp
