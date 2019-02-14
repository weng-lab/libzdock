#pragma once

#include "pdb++.h"
#include <Eigen/Dense>
#include <iostream>
#include <vector>
#include <mutex>

namespace zdock {

class Coord {
public:
  int serialNum;
  std::string atomName;
  std::string resName;
  char chain;
  int resNum;
  Coord() : serialNum(0), atomName(""), resName(""), chain('\0'), resNum(0) {}
};

class PDB {
private:
  typedef Eigen::Matrix<double, 3, Eigen::Dynamic> Matrix;
  typedef Eigen::Transform<double, 3, Eigen::Affine> Transform;

  std::vector<size_t> atoms_;
  std::vector<libpdb::PDB> records_;
  Matrix m_;
  void read_(const std::string &fn, const int model = MODEL_FIRST,
             std::function<bool(const libpdb::PDB &)> filter =
                 [](const libpdb::PDB &) { return true; });
  std::mutex insertmtx_;
public:
  static const int MODEL_ALL = -1;
  static const int MODEL_FIRST = 0;
  PDB();
  PDB(const std::string &fn, const int model = MODEL_FIRST,
      std::function<bool(const libpdb::PDB &)> filter =
          [](const libpdb::PDB &) { return true; });
  const Matrix &matrix() const;
  const Matrix &transform(const Transform &t);
  const Matrix &setMatrix(const Matrix &m);
  const std::vector<libpdb::PDB> &records();
  const Eigen::Vector3d centroid() const;
  const libpdb::PDB& operator[](const int serial) const;
  const libpdb::PDB& operator[](const Coord& coord) const;
  void append(const libpdb::PDB&);
};

inline std::ostream &operator<<(std::ostream &s, const Coord &c) {
  std::ostringstream os;
  os << c.serialNum << '\t';
  os << c.atomName << '\t';
  os << c.resName << '\t';
  os << c.chain << '\t';
  os << c.resNum;
  s << os.str();
  return s;
}

} // namespace zdock
