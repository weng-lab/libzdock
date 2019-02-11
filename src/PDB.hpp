#pragma once

#include "pdb++.h"
#include <Eigen/Dense>
#include <iostream>
#include <vector>

namespace zdock {
class PDB {
private:
  std::vector<size_t> atoms_;
  std::vector<libpdb::PDB> records_;
  Eigen::Matrix<double, 3, Eigen::Dynamic> m_;
  void read_(const std::string &fn, const int model = MODEL_FIRST,
             const bool alpha = false);

public:
  static const int MODEL_ALL = -1;
  static const int MODEL_FIRST = 0;
  PDB(const std::string &fn, const int model = MODEL_FIRST,
      const bool alpha = false);
  const Eigen::Matrix<double, 3, Eigen::Dynamic> &matrix();
  const Eigen::Matrix<double, 3, Eigen::Dynamic> &
  transform(const Eigen::Transform<double, 3, Eigen::Affine> &t);
  const Eigen::Matrix<double, 3, Eigen::Dynamic> &
  setMatrix(const Eigen::Matrix<double, 3, Eigen::Dynamic> &m);
  const std::vector<libpdb::PDB> &records();
};

} // namespace zdock
