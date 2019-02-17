#include "PDB.hpp"
#include "TransformUtil.hpp"
#include "ZDOCK.hpp"

namespace zdock {

class TransformMultimer {
public:
  TransformMultimer(const std::string &zdock);
  TransformMultimer(const ZDOCK &zdock);

private:
  typedef Eigen::Transform<double, 3, Eigen::Affine> Transform;
  typedef Eigen::Matrix<double, 3, Eigen::Dynamic> Matrix;
  typedef TransformUtil u;

  zdock::Structure structure_; // zdock metadata
  ZDOCK zdock_;                // zdock output
  double spacing_;             // grid spacing
  int boxsize_;                // grid size
  int symmetry_;               // zdock symmetry number
  double alpha_, beta_;        // alpha and beta angles
  double factor_;              // scaling factor for distance d

  // precomputed transformation matrices
  Transform t0_;

  // grid to actual translation ('circularized')
  inline const Transform initialTranslation(const int (&t)[3],
                                            const double alpha) const {
    using Eigen::Translation3d;
    using Eigen::Vector3d;

    Transform ret;
    Vector3d d;

    // circular translation
    d << (t[0] >= boxsize_ / 2 ? t[0] - boxsize_ : t[0]),
        (t[1] >= boxsize_ / 2 ? t[1] - boxsize_ : t[1]),
        (t[2] >= boxsize_ / 2 ? t[2] - boxsize_ : t[2]);

    // rotate vector by angle alpha
    d = u::eulerRotation({0.0, 0.0, alpha}, true) * (spacing_ * factor_ * d);
    ret = Translation3d(d);
    return ret;
  }

public:
  static const std::string CHAINS[52];

  // perform actual structure transformation
  inline const Matrix txMultimer(const PDB &pdb, const Prediction &pred,
                                 int n) const {
    Transform t;
    using Eigen::Translation3d;
    using Eigen::Vector3d;

    assert(n >= 0 && n < symmetry_);

    /* Transformation; M-ZDOCK (one of n-mer)
     *
     * (create_multimer_num somehow has all rotations reversed;
     *  for compatibility we do the same here.)
     *
     *     X(rotate along z, reverse) *
     *      X(translate d, rotated by alpha) *
     *       X(phi/theta rot, reverse) *
     *        X(random rot, reverse) *r
     *         X(-offset trans) * M
     *
     */

    t = u::eulerRotation({0.0, 0.0, beta_ * n}, true) *
        initialTranslation(pred.translation, alpha_) *
        u::eulerRotation({pred.rotation[0], pred.rotation[1]}, true) * t0_;
    return t * pdb.matrix();
  }
};

} // namespace zdock
