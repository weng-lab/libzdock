#include "Pruning.hpp"
#include "Utils.hpp"

#include <cmath>
#include <iostream>
#include <limits>
#include <map>

namespace e = Eigen;
namespace p = libpdb;

namespace zdock {

Pruning::Pruning(const std::string &zdockoutput, const std::string &receptorpdb,
                 const std::string &ligandpdb)
    : zdock_(zdockoutput), txl_(zdockoutput) {

  // receptor file name
  if ("" == receptorpdb) {
    recfn_ = zdock_.receptor().filename;
    if ('/' != recfn_[0]) { // relative
      recfn_ = Utils::copath(zdockoutput, recfn_);
    }
  } else {
    recfn_ = receptorpdb;
  }

  // ligand file name
  if (!zdock_.ismzdock()) {
    if ("" == ligandpdb) {
      ligfn_ = zdock_.ligand().filename;
      if ('/' != ligfn_[0]) { // relative
        ligfn_ = Utils::copath(zdockoutput, ligfn_);
      }
    } else {
      ligfn_ = ligandpdb;
    }
  }
}

void Pruning::prune(const double cutoff) {
  const auto v = zdock_.predictions(); // our copy
  const auto n = zdock_.npredictions();
  auto &preds = zdock_.predictions();              // our ref
  double min = std::numeric_limits<double>::max(); // big number
  std::vector<int> l(n, 0);
  int clusters = 0;

  // read pdb files
  PDB recpdb(recfn_, PDB::MODEL_FIRST, [](const p::PDB &r) {
    return p::PDB::ATOM == r.type() && r.isalpha(); // CA only
  });
  PDB ligpdb(ligfn_, PDB::MODEL_FIRST, [](const p::PDB &r) {
    return p::PDB::ATOM == r.type() && r.isalpha(); // CA only
  });
  const double ligsize = ligpdb.matrix().cols();

  // pre-compute all poses
  std::vector<Pruning::Matrix> poses;
  for (size_t i = 0; i < n; ++i) {
    poses.push_back(txl_.txLigand(ligpdb, v[i]));
  }

  // find clusters
  zdock_.predictions().clear();
  for (size_t i = 0; i < n; ++i) {
    if (!(i % 100)) {
      std::cerr << i << '\t' << clusters << std::endl;
    }
    if (!l.at(i)) {
      l[i] = clusters + 1;
      preds.push_back(v[i]);
      for (size_t j = i + 1; j < n; ++j) {
        if (!l.at(j)) {
          const double rmsd =
              std::sqrt((poses.at(i) - poses.at(j)).squaredNorm() / ligsize);
          min = std::min(min, rmsd);
          if (rmsd < cutoff) {
            l[j] = clusters + 1;
          }
        }
      }
      clusters++;
    }
  }

  // print zdock file
  std::cout << zdock_ << std::endl;

  // print some stats
  std::cerr << std::endl
            << "cutoff: " << cutoff << ", min: " << min
            << ", ligsize: " << ligsize << ", clusters: " << clusters
            << std::endl;

  // count cluster members
  std::map<int, int> counts;
  for (auto x : l) {
    counts[x]++;
  }
  for (const auto &p : counts) {
    std::cerr << "(" << p.first << ", " << p.second << ")" << '\n';
  }
}

void Pruning::makeZDOCKComplex(const size_t n) {
  const auto &v = zdock_.predictions();
  const auto &p = v[n];

  // read pdb files
  PDB recpdb(recfn_, PDB::MODEL_FIRST);
  PDB ligpdb(ligfn_, PDB::MODEL_FIRST);


  ligpdb.setMatrix(txl_.txLigand(ligpdb, p));
  for (const auto &x : recpdb.records()) {
    std::cout << x << '\n';
  }
  for (const auto &x : ligpdb.records()) {
    std::cout << x << '\n';
  }
}

void Pruning::makeMZDOCKComplex(const size_t n) {
  std::cerr << "lala" << std::endl;
}

void Pruning::makeComplex(const size_t n) {
  if (zdock_.ismzdock()) {
    makeMZDOCKComplex(n);
  } else {
    makeZDOCKComplex(n);
  }
}

void Pruning::filterConstraints(const std::string &fn) {
  // load constraints file
  Constraints ccc(fn);

  // load full PDB files (first model)
  PDB receptor(recfn_, PDB::MODEL_FIRST,
               [](const p::PDB &r) { return p::PDB::ATOM == r.type(); });
  PDB ligand(ligfn_, PDB::MODEL_FIRST,
             [](const p::PDB &r) { return p::PDB::ATOM == r.type(); });

  // grab atoms for valid constraints
  PDB ligatoms, recatoms;
  std::vector<double> mindist;
  std::vector<double> maxdist;
  if (!ccc.constraints().size()) {
    throw ConstraintException("No constraints specified");
  }
  for (const auto &x : ccc.constraints()) {
    try {
      // these throw exceptions for bad constraints; we want to append
      // both_ l and r or _none_.
      const p::PDB &r = receptor[x.recCoord];
      const p::PDB &l = ligand[x.ligCoord];
      recatoms.append(r);
      ligatoms.append(l);
      if (Constraint::MAX == x.constraintType) {
        maxdist.push_back(x.distance);
        mindist.push_back(-1.0); // negative distance
      } else {
        maxdist.push_back(std::numeric_limits<double>::max()); // big number
        mindist.push_back(x.distance);
      }
    } catch (const AtomNotFoundException e) {
      throw ConstraintException("Constraint Error: " + std::string(e.what()));
    }
  }

  // assess all poses
  const size_t n = maxdist.size();
  const auto v = zdock_.predictions(); // our copy
  auto &preds = zdock_.predictions();  // our ref
  preds.clear();
  for (const Prediction &p : v) {
    Pruning::Matrix pose = txl_.txLigand(ligatoms, p);
    e::Matrix<double, 1, e::Dynamic> m =
        (pose - recatoms.matrix()).colwise().squaredNorm().array().sqrt();
    bool accepted = true;
    for (size_t i = 0; i < n; ++i) {
      if (m(0, i) > maxdist[i] || m(0, i) < mindist[i]) {
        accepted = false;
        break;
      }
    }
    if (accepted) {
      preds.push_back(p);
    }
  }
  std::cout << zdock_ << std::endl;
}

} // namespace zdock
