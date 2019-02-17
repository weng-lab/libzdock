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
