#include "PDB.hpp"
#include "ZDOCK.hpp"
#include <memory>

namespace zdock {

class Pruning {
private:
  zdock::Structure receptor_, ligand_;
  std::unique_ptr<PDB> recpdb_, ligpdb_;
  ZDOCK zdock_;
  double spacing_;
  int boxsize_;
  bool rev_, fixed_;
  Eigen::Transform<double, 3, Eigen::Affine> t0_, t1_, t2_, t3_;
  const Eigen::Transform<double, 3, Eigen::Affine>
  eulerRotation(const double (&r)[3], bool rev = false) const;
  const Eigen::Transform<double, 3, Eigen::Affine>
  boxTranslation(const int (&t)[3], bool rev = false) const;

public:
  Pruning(const std::string &zdockouput, const std::string &receptorpdb = "",
          const std::string &ligandpdb = "");
  void prune(const double cutoff);
  double rmsd(const Prediction p1, const Prediction p2,
              const double ligsize = 1.0) const;
  void makeComplex(const size_t n);

  inline const Eigen::Matrix<double, 3, Eigen::Dynamic>
  txLigand(const PDB &pdb, const Prediction &pred) const {
    Eigen::Transform<double, 3, Eigen::Affine> t;
    if (rev_) {
      t = eulerRotation(pred.rotation, true) * boxTranslation(pred.translation);
      return t1_ * t * t0_ * pdb.matrix();
    } else {
      t = Eigen::Translation3d(Eigen::Vector3d(receptor_.translation)) *
          boxTranslation(pred.translation, true) *
          eulerRotation(pred.rotation) * t2_;
      if (!fixed_) {
        t = eulerRotation(receptor_.rotation, true) * t;
      }
      return t * pdb.matrix();
    }
  }
};

} // namespace zdock
