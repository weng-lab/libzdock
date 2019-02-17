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

#include "TransformLigand.hpp"

namespace e = Eigen;
namespace p = libpdb;

namespace zdock {

TransformLigand::TransformLigand(const std::string &zdock)
    : TransformLigand(ZDOCK(zdock)) {}

TransformLigand::TransformLigand(const ZDOCK &zdock) : zdock_(zdock) {

  using e::Translation3d;
  using e::Vector3d;

  // copy relevant info from zdock file
  receptor_ = zdock_.receptor();
  ligand_ = zdock_.ligand();
  rev_ = zdock_.isswitched();
  fixed_ = zdock_.isfixed();
  spacing_ = zdock_.spacing();
  boxsize_ = zdock_.boxsize();

  // precalculate some transformation matrices
  t0_ = Translation3d(-Vector3d(ligand_.translation)) *
        u::eulerRotation(receptor_.rotation);
  t1_ = Translation3d(Vector3d(receptor_.translation)) *
        u::eulerRotation(ligand_.rotation, true);
  t2_ = u::eulerRotation(ligand_.rotation) *
        Translation3d(-Vector3d(ligand_.translation));
}

} // namespace zdock
