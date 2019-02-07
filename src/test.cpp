#include "ZDOCK.hpp"
#include <iostream>
#include <exception>

//namespace eg = ::Eigen;
//namespace p = ::libpdb;

int main() {
  /*
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
  */

  std::vector<std::string> files = {
    "/Users/vanderva/git/zdockserver/webservice/test/irad/2MTA.zd.out",
    "/Users/vanderva/Desktop/c0d92b1b-0888-4fa3-83a2-0ccc2f7e60af/zdock.out.pruned",
    "/Users/vanderva/Desktop/0116d0ab47/job.154074.mzdock_24.out"
  };

  for (const auto& x : files) {
    zlab::ZDOCK z(x);
    std::cout << "filename: " << z.filename() << std::endl;
    std::cout << "ismzdock: " << z.ismzdock() << std::endl;
    std::cout << "isswitched: " << z.isswitched() << std::endl;
    std::cout << "isfixed: " << z.isfixed() << std::endl;
    std::cout << "version: " << z.version() << std::endl;
    std::cout << "npred: " << z.npredictions() << std::endl;
    std::cout << "receptor: " << z.receptor() << std::endl;
    if (z.ismzdock()) {
      std::cout << "symmetry: " << z.symmetry() << std::endl;
    }
    if (!z.ismzdock()) {
      std::cout << "ligand: " << z.ligand() << std::endl;
    }
    for (const auto& x : z.predictions()) {
      std::cout << x << '\n';
      break;
    }
    std::cout << std::endl;
  }

}

// clang++ -isystem/opt/local/include/eigen3 -framework Accelerate
// /opt/local/lib/lapack/liblapacke.dylib  -O3 -L. -lpdb++ -o bla main.cpp
