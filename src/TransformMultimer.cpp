#include "TransformMultimer.hpp"

namespace e = Eigen;
namespace p = libpdb;

namespace zdock {

TransformMultimer::TransformMultimer(const std::string &zdock)
    : TransformMultimer(ZDOCK(zdock)) {}

TransformMultimer::TransformMultimer(const ZDOCK &zdock) : zdock_(zdock) {

  using e::Translation3d;
  using e::Vector3d;
  using std::cos;

  // copy relevant info from zdock file
  structure_ = zdock_.receptor();
  spacing_ = zdock_.spacing();
  boxsize_ = zdock_.boxsize();
  symmetry_ = (zdock_.symmetry() ? zdock_.symmetry() : 3);
  beta_ = 2 * u::PI / symmetry_;
  alpha_ = 0.5 * (u::PI - beta_);
  factor_ = 1.0 / (2.0 * cos(alpha_));

  // precalculate some transformation matrices
  t0_ = TransformUtil::eulerRotation(structure_.rotation, true) *
        Translation3d(-Vector3d(structure_.translation));
}

// alphabeth for chains
const std::string TransformMultimer::CHAINS[52] = {
    "A", "B", "C", "D", "E", "F", "G", "H", "I", "J", "K", "L", "M",
    "N", "O", "P", "Q", "R", "S", "T", "U", "V", "W", "X", "Y", "Z",
    "a", "b", "c", "d", "e", "f", "g", "h", "i", "j", "k", "l", "m",
    "n", "o", "p", "q", "r", "s", "r", "u", "v", "w", "x", "y", "z"};

} // namespace zdock
