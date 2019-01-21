#include "pdb++.h"
#include <iostream>
#include <vector>

#define EIGEN_USE_LAPACKE

#include <Eigen/Dense>

namespace zlab {
class PDB {
private:
  std::vector<size_t> atoms_;
  std::vector<libpdb::PDB> records_;
  Eigen::Matrix<double, 3, Eigen::Dynamic> m_;
  std::string fn_;
  void read_();

public:
  PDB(const std::string &fn);
  const Eigen::Matrix<double, 3, Eigen::Dynamic> &matrix();
  const Eigen::Matrix<double, 3, Eigen::Dynamic> &
  transform(Eigen::Transform<double, 3, Eigen::Affine>);
  const std::vector<libpdb::PDB> &records();
};

} // namespace zlab
