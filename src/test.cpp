#include "PDB.hpp"
#include <exception>

namespace eg = ::Eigen;
namespace p = ::libpdb;

int main() {
  try {
    zlab::PDB pdb("1f6g.pdb");

    // perform a translation (i.e. rotate around center of mass)
    eg::Vector3d x = -pdb.matrix().rowwise().mean();
    eg::Transform<double, 3, eg::Affine> t =
        eg::Translation3d(-x) *
        (eg::AngleAxisd(0.0 * M_PI, eg::Vector3d::UnitX()) *
         eg::AngleAxisd(0.25 * M_PI, eg::Vector3d::UnitY()) *
         eg::AngleAxisd(0.25 * M_PI, eg::Vector3d::UnitZ())) *
        eg::Translation3d(x);
    pdb.transform(t);

    for (const auto &x : pdb.records()) {
      std::cout << x << std::endl;
    }
  } catch (const std::exception &e) {
    std::cerr << e.what() << std::endl;
  }
}

// clang++ -isystem/opt/local/include/eigen3 -framework Accelerate
// /opt/local/lib/lapack/liblapacke.dylib  -O3 -L. -lpdb++ -o bla main.cpp
