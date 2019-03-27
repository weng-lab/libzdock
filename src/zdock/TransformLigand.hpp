/**
 * Copyright (c) 2019, Arjan van der Velde, Weng Lab
 * All rights reserved.
 *
 * Redistribution and use in source and binary forms, with or without
 * modification, are permitted provided that the following conditions are met:
 *
 * * Redistributions of source code must retain the above copyright notice, this
 *   list of conditions and the following disclaimer.
 *
 * * Redistributions in binary form must reproduce the above copyright notice,
 *   this list of conditions and the following disclaimer in the documentation
 *   and/or other materials provided with the distribution.
 *
 * THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND CONTRIBUTORS "AS IS"
 * AND ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT LIMITED TO, THE
 * IMPLIED WARRANTIES OF MERCHANTABILITY AND FITNESS FOR A PARTICULAR PURPOSE ARE
 * DISCLAIMED. IN NO EVENT SHALL THE COPYRIGHT HOLDER OR CONTRIBUTORS BE LIABLE
 * FOR ANY DIRECT, INDIRECT, INCIDENTAL, SPECIAL, EXEMPLARY, OR CONSEQUENTIAL
 * DAMAGES (INCLUDING, BUT NOT LIMITED TO, PROCUREMENT OF SUBSTITUTE GOODS OR
 * SERVICES; LOSS OF USE, DATA, OR PROFITS; OR BUSINESS INTERRUPTION) HOWEVER
 * CAUSED AND ON ANY THEORY OF LIABILITY, WHETHER IN CONTRACT, STRICT LIABILITY,
 * OR TORT (INCLUDING NEGLIGENCE OR OTHERWISE) ARISING IN ANY WAY OUT OF THE USE
 * OF THIS SOFTWARE, EVEN IF ADVISED OF THE POSSIBILITY OF SUCH DAMAGE.
 */

#pragma once

#include "PDB.hpp"
#include "TransformUtil.hpp"
#include "ZDOCK.hpp"

namespace zdock {

class TransformLigand {

  /**
   * Performs ZDOCK prediction transformations on PDB structures.
   */

public:
  TransformLigand(const std::string &zdock);
  TransformLigand(const ZDOCK &zdock);

private:
  typedef Eigen::Transform<double, 3, Eigen::Affine> Transform;
  typedef Eigen::Matrix<double, 3, Eigen::Dynamic> Matrix;
  typedef TransformUtil u;

  zdock::Structure receptor_, ligand_; // zdock metadata
  ZDOCK zdock_;                        // zdock output
  double spacing_;                     // grid spacing
  int boxsize_;                        // grid size
  bool rev_, fixed_;                   // fixed / switched flags
  bool isvalid_;                       // successfull init

  // precomputed transformation matrices
  Transform t0_, t1_, t2_;

  // grid to actual translation ('circularized')
  inline const Transform boxTranslation(const int (&t)[3],
                                        bool rev = false) const {
    Transform ret;

    using Eigen::Translation3d;
    using Eigen::Vector3d;

    Vector3d d = u::boxedGridCoord(t, boxsize_);
    if (rev) {
      ret = Translation3d(-spacing_ * d);
    } else {
      ret = Translation3d(spacing_ * d);
    }
    return ret;
  }

public:
  // perform actual ligand transformation
  inline const Matrix txLigand(const Matrix &matrix,
                               const Prediction &pred) const {
    Transform t;

    using Eigen::Translation3d;
    using Eigen::Vector3d;

    assert(isvalid_); // did we successfully load zdock data?

    if (rev_) {

      /* Transformation; reverse (receptor was rotated, rather than ligand)
       *
       *     X(rec trans) * X(lig rot, reverse) *
       *       X(pred rot, reverse) * X(pred trans) *
       *         X(-lig trans) * X(rec rot) * M
       *
       */

      t = u::eulerRotation(pred.rotation, true) *
          boxTranslation(pred.translation);
      return t1_ * t * t0_ * matrix;
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
          u::eulerRotation(pred.rotation) * t2_;
      if (!fixed_) {
        // !fixed means initial random rotation of receptor
        // so we need to rotate to rec frame
        t = u::eulerRotation(receptor_.rotation, true) * t;
      }
      return t * matrix;
    }
  }
};

} // namespace zdock
