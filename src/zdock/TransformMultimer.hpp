/**
 * Copyright 2019 Arjan van der Velde, vandervelde.ag [at] gmail
 * Licensed under the Apache License, Version 2.0 (the "License");
 * you may not use this file except in compliance with the License.
 * You may obtain a copy of the License at
 * http://www.apache.org/licenses/LICENSE-2.0
 *
 * Unless required by applicable law or agreed to in writing, software
 * distributed under the License is distributed on an "AS IS" BASIS,
 * WITHOUT WARRANTIES OR CONDITIONS OF ANY KIND, either express or implied.
 * See the License for the specific language governing permissions and
 * limitations under the License.
 */

#pragma once

#include "PDB.hpp"
#include "TransformUtil.hpp"
#include "ZDOCK.hpp"

namespace zdock {

class TransformMultimer {

  /**
   * Implements M-ZDOCK prediction transformations on PDB structures.
   *
   * see:
   *
   * M-ZDOCK: a grid-based approach for C n symmetric multimer docking
   * Brian Pierce  Weiwei Tong  Zhiping Weng
   *
   * Bioinformatics, Volume 21, Issue 8, 15 April 2005, Pages 1472–1478
   * https://doi.org/10.1093/bioinformatics/bti229
   *
   */

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
  bool isvalid_;               // indicates successfull init

  // precomputed transformation matrices
  Transform t0_;

  // grid to actual translation ('circularized')
  inline const Transform initialTranslation(const int (&t)[3],
                                            const double alpha) const {
    using Eigen::Translation3d;
    using Eigen::Vector3d;

    Transform ret;

    // rotate vector by angle alpha
    Vector3d d = u::boxedGridCoord(t, boxsize_);
    d = u::eulerRotation({0.0, 0.0, alpha}, true) * (spacing_ * factor_ * d);
    ret = Translation3d(d);
    return ret;
  }

public:
  static const std::string CHAINS[52];

  // perform actual structure transformation
  inline const Matrix txMultimer(const Matrix &matrix, const Prediction &pred,
                                 int n) const {
    Transform t;
    using Eigen::Translation3d;
    using Eigen::Vector3d;

    assert(isvalid_); // did we successfully load m-zdock data?
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
    return t * matrix;
  }
};

} // namespace zdock
