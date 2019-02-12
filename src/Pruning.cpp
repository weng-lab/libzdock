#include "Pruning.hpp"
#include "Utils.hpp"
#include <cmath>
#include <iostream>
#include <limits>

namespace e = Eigen;

namespace zdock {

Pruning::Pruning(const std::string &zdockouput, const std::string &receptorpdb,
                 const std::string &ligandpdb)
    : zdock_(zdockouput) {

  // read receptor pdb
  if ("" == receptorpdb) {
    std::string fn = zdock_.receptor().filename;
    if ('/' != fn[0]) { // relative
      fn = Utils::copath(zdockouput, fn);
    }
    recpdb_ = std::make_unique<PDB>(fn, PDB::MODEL_ALL, true);
  } else {
    recpdb_ = std::make_unique<PDB>(receptorpdb, PDB::MODEL_ALL, true);
  }

  // read ligand pdb
  if ("" == ligandpdb) {
    std::string fn = zdock_.ligand().filename;
    if ('/' != fn[0]) { // relative
      fn = Utils::copath(zdockouput, fn);
    }
    ligpdb_ = std::make_unique<PDB>(fn, PDB::MODEL_ALL, true);
  } else {
    ligpdb_ = std::make_unique<PDB>(ligandpdb, PDB::MODEL_ALL, true);
  }

  // copy relevant info from zdock file
  receptor_ = zdock_.receptor();
  ligand_ = zdock_.ligand();
  rev_ = zdock_.isswitched();
  fixed_ = zdock_.isfixed();
  spacing_ = zdock_.spacing();
  boxsize_ = zdock_.boxsize();

  // precalculate some transformation matrices
  t0_ = e::Translation3d(-e::Vector3d(ligand_.translation)) *
        eulerRotation(receptor_.rotation);
  t1_ = e::Translation3d(e::Vector3d(receptor_.translation)) *
        eulerRotation(ligand_.rotation, true);
  t2_ = eulerRotation(ligand_.rotation) *
        e::Translation3d(-e::Vector3d(ligand_.translation));
}

void Pruning::prune(const double cutoff) {
  const auto v = zdock_.predictions(); // our copy
  const auto n = zdock_.npredictions();
  const double ligsize = ligpdb_->matrix().cols();
  auto &preds = zdock_.predictions(); // our ref
  std::vector<bool> l(n, true);
  double min = std::numeric_limits<double>::max(); // big number
  int clusters = 0;
  zdock_.predictions().clear();
  for (size_t i = 0; i < n; ++i) {
    if (l[i]) {
      l[i] = false;
      preds.push_back(v[i]);
      const auto x = txLigand(*ligpdb_, v[i]);
      for (size_t j = i + 1; j < n; ++j) {
        if (l[j]) {
          const auto y = txLigand(*ligpdb_, v[j]);
          const double rmsd = std::sqrt((x - y).array().pow(2).sum() / ligsize);
          min = std::min(min, rmsd);
          if (rmsd < cutoff) {
            l[j] = false;
          }
        }
      }
      clusters++;
    }
  }
  // print zdock file
  std::cout << zdock_ << std::endl;
  std::cerr << "cutoff: " << cutoff << ", min: " << min
            << ", ligsize: " << ligsize << ", clusters: " << clusters
            << std::endl;
}

void Pruning::makeComplex(const size_t n) {
  const auto &v = zdock_.predictions();
  const auto &p = v[n];
  ligpdb_->setMatrix(txLigand(*ligpdb_, p));
  for (const auto &x : recpdb_->records()) {
    std::cout << x << '\n';
  }
  for (const auto &x : ligpdb_->records()) {
    std::cout << x << '\n';
  }
}

} // namespace zdock
