#include "PDB.hpp"
#include "ZDOCK.hpp"
#include <exception>
#include <iostream>

namespace e = Eigen;

const e::Transform<double, 3, e::Affine> eulerRotation(const double (&r)[3],
                                                       bool rev = false) {
  e::Transform<double, 3, e::Affine> t;
  t = e::AngleAxisd(r[0], e::Vector3d::UnitZ()) *
      e::AngleAxisd(r[1], e::Vector3d::UnitX()) *
      e::AngleAxisd(r[2], e::Vector3d::UnitZ());
  return (rev ? t.inverse() : t);
}

const e::Transform<double, 3, e::Affine>
boxTranslation(const int (&t)[3], const double spacing, const int boxsize) {
  e::Transform<double, 3, e::Affine> ret;
  e::Vector3d d;
  d << (t[0] >= boxsize / 2 ? t[0] - boxsize : t[0]),
      (t[1] >= boxsize / 2 ? t[1] - boxsize : t[1]),
      (t[2] >= boxsize / 2 ? t[2] - boxsize : t[2]);
  ret = e::Translation3d(spacing * d);
  return ret;
}

int main() {
  const std::string zfile = "/Users/vanderva/git/libpdb++/572b42d9-8406-41df-8200-fabd347a9548/zdock.out.pruned";
  const std::string receptor = "/Users/vanderva/git/libpdb++/572b42d9-8406-41df-8200-fabd347a9548/receptor.pdb";
  const std::string ligand = "/Users/vanderva/git/libpdb++/572b42d9-8406-41df-8200-fabd347a9548/ligand.pdb";

  e::Transform<double, 3, e::Affine> t0, t1, tp;

  try {
    zdock::PDB rec(receptor); //, zdock::PDB::MODEL_FIRST, true);
    zdock::PDB lig(ligand); //, zdock::PDB::MODEL_FIRST, true);
    zdock::ZDOCK z(zfile);

    const auto p = z.predictions()[0];
    const auto &l = z.ligand();
    const auto &r = z.receptor();

    t0 = e::Translation3d(-e::Vector3d(l.translation)) *
         eulerRotation(r.rotation);
    t1 = e::Translation3d(e::Vector3d(r.translation)) *
         eulerRotation(l.rotation, true);
    tp = eulerRotation(p.rotation, true) *
         boxTranslation(p.translation, z.spacing(), z.boxsize());
    lig.transform(t1 * tp * t0);

    for (const auto &x : rec.records()) {
      std::cout << x << '\n';
    }
    for (const auto &x : lig.records()) {
      std::cout << x << '\n';
    }

  } catch (const std::exception &e) {
    std::cerr << "Exception: " << e.what() << std::endl;
    exit(1);
  }
}

