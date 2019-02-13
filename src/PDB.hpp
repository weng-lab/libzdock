#pragma once

#include "pdb++.h"
#include <Eigen/Dense>
#include <iostream>
#include <vector>

namespace zdock {

class PDB {
private:
  typedef Eigen::Matrix<double, 3, Eigen::Dynamic> Matrix;
  typedef Eigen::Transform<double, 3, Eigen::Affine> Transform;

  std::vector<size_t> atoms_;
  std::vector<libpdb::PDB> records_;
  Matrix m_;
  void read_(const std::string &fn, const int model = MODEL_FIRST,
             const bool alpha = false);

public:
  static const int MODEL_ALL = -1;
  static const int MODEL_FIRST = 0;
  PDB(const std::string &fn, const int model = MODEL_FIRST,
      const bool alpha = false);
  const Matrix &matrix() const;
  const Matrix &transform(const Transform &t);
  const Matrix &setMatrix(const Matrix &m);
  const std::vector<libpdb::PDB> &records();
  const Eigen::Vector3d centroid() const;
};

} // namespace zdock
