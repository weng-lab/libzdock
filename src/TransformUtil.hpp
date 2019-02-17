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

};

} // namespace zdock
