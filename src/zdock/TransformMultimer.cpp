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
  if (zdock_.ismzdock()) {
    structure_ = zdock_.receptor();
    spacing_ = zdock_.spacing();
    boxsize_ = zdock_.boxsize();
    symmetry_ = (zdock_.symmetry() ? zdock_.symmetry() : 3);
    beta_ = 2 * u::PI / symmetry_;
    alpha_ = 0.5 * (u::PI - beta_);
    factor_ = 1.0 / (2.0 * cos(alpha_));
    isvalid_ = true;

    // precalculate some transformation matrices
    t0_ = TransformUtil::eulerRotation(structure_.rotation, true) *
          Translation3d(-Vector3d(structure_.translation));
  }
}

// alphabeth for chains
const std::string TransformMultimer::CHAINS[52] = {
    "A", "B", "C", "D", "E", "F", "G", "H", "I", "J", "K", "L", "M",
    "N", "O", "P", "Q", "R", "S", "T", "U", "V", "W", "X", "Y", "Z",
    "a", "b", "c", "d", "e", "f", "g", "h", "i", "j", "k", "l", "m",
    "n", "o", "p", "q", "r", "s", "r", "u", "v", "w", "x", "y", "z"};

} // namespace zdock
