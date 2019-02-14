#include "Constraints.hpp"
#include "PDB.hpp"
#include "ZDOCK.hpp"
#include <memory>

namespace zdock {

class Pruning {
private:
  typedef Eigen::Transform<double, 3, Eigen::Affine> Transform;
  typedef Eigen::Matrix<double, 3, Eigen::Dynamic> Matrix;

  zdock::Structure receptor_, ligand_;   // zdock metadata
  std::unique_ptr<PDB> recpdb_, ligpdb_; // pdb file data
  ZDOCK zdock_;                          // zdock output
  double spacing_;                       // grid spacing
  int boxsize_;                          // grid size
  bool rev_, fixed_;                     // fixed / switched flags
  std::string recfn_, ligfn_;            // receptor and ligand filenames

  // precomputed transformation matrices
  Transform t0_, t1_, t2_;

  // Euler angles to Z-X-Z transformation matrix
  inline const Transform eulerRotation(const double (&r)[3],
                                       bool rev = false) const {
    Transform t;

    using Eigen::AngleAxisd;
    using Eigen::Vector3d;

    t = AngleAxisd(r[0], Vector3d::UnitZ()) *
        AngleAxisd(r[1], Vector3d::UnitX()) *
        AngleAxisd(r[2], Vector3d::UnitZ());
    return (rev ? t.inverse() : t);
  }

  // grid to actual translation ('circularized')
  inline const Transform boxTranslation(const int (&t)[3],
                                        bool rev = false) const {
    Transform ret;

    using Eigen::Translation3d;
    using Eigen::Vector3d;

    Vector3d d;
    d << (t[0] >= boxsize_ / 2 ? t[0] - boxsize_ : t[0]),
        (t[1] >= boxsize_ / 2 ? t[1] - boxsize_ : t[1]),
        (t[2] >= boxsize_ / 2 ? t[2] - boxsize_ : t[2]);
    if (rev) {
      ret = Translation3d(-spacing_ * d);
    } else {
      ret = Translation3d(spacing_ * d);
    }
    return ret;
  }

public:
  Pruning(const std::string &zdockouput,
          const std::string &receptorpdb = "", // or grab from zdock.out
          const std::string &ligandpdb = ""    // or grab from zdock.out
          );

  // perform pruning
  void prune(const double cutoff);

  // ligand pdb record to stdout
  void makeComplex(const size_t n);

  // filter predictions based on constraints
  void filterConstraints(const std::string& fn);

  // perform actual ligand transformation
  inline const Matrix txLigand(const PDB &pdb, const Prediction &pred) const {
    Transform t;

    using Eigen::Translation3d;
    using Eigen::Vector3d;

    if (rev_) {

      /* Transformation; reverse (receptor was rotated, rather than ligand)
       *
       *     X(rec trans) * X(lig rot, reverse) *
       *       X(pred rot, reverse) * X(pred trans) *
       *         X(-lig trans) * X(rec rot) * M
       *
       */

      t = eulerRotation(pred.rotation, true) * boxTranslation(pred.translation);
      return t1_ * t * t0_ * pdb.matrix();
    } else {

      /* Transformation; normal (ligand was rotated)
       *
       * fixed case:
       *     X(rec trans) * X(-pred trans) * X(pred rot) *
       *       X(lig rot) * X(-lig trans) * M
       *
       * not fixed case:
       *     X(rec rot, reverse) * [fixed case]
       *
       */

      t = Translation3d(Vector3d(receptor_.translation)) *
          boxTranslation(pred.translation, true) *
          eulerRotation(pred.rotation) * t2_;
      if (!fixed_) {
        // !fixed means initial random rotation of receptor
        // so we need to rotate to rec frame
        t = eulerRotation(receptor_.rotation, true) * t;
      }
      return t * pdb.matrix();
    }
  }
};

} // namespace zdock
