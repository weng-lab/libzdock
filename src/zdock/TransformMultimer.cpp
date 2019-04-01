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
    structure_ = zdock_.structure();
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
