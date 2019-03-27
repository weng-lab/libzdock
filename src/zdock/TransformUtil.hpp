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

#include <Eigen/Dense>
#include <cmath>

namespace zdock {

class TransformUtil {
private:
  typedef Eigen::Transform<double, 3, Eigen::Affine> Transform;
  typedef Eigen::Matrix<double, 3, Eigen::Dynamic> Matrix;

public:
  static const double PI;

  // Euler angles to Z-X-Z transformation matrix
  static inline const Transform eulerRotation(const double (&r)[3],
                                              bool rev = false) {
    Transform t;

    using Eigen::AngleAxisd;
    using Eigen::Vector3d;

    t = AngleAxisd(r[0], Vector3d::UnitZ()) *
        AngleAxisd(r[1], Vector3d::UnitX()) *
        AngleAxisd(r[2], Vector3d::UnitZ());
    return (rev ? t.inverse() : t);
  }

  // Adjust grid coordinates to fall inside box of size 'boxsize',
  // NOTE: returns Vector3d (double).
  static inline const Eigen::Vector3d boxedGridCoord(const int (&v)[3],
                                                     const int boxsize) {
    Eigen::Vector3d d;
    d << (v[0] >= boxsize / 2 ? v[0] - boxsize : v[0]),
        (v[1] >= boxsize / 2 ? v[1] - boxsize : v[1]),
        (v[2] >= boxsize / 2 ? v[2] - boxsize : v[2]);
    return d;
  }
};

} // namespace zdock
