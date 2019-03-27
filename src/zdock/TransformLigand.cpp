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
  if (zdock_.iszdock()) {
    receptor_ = zdock_.receptor();
    ligand_ = zdock_.ligand();
    rev_ = zdock_.isswitched();
    fixed_ = zdock_.isfixed();
    spacing_ = zdock_.spacing();
    boxsize_ = zdock_.boxsize();
    isvalid_ = true;

    // precalculate some transformation matrices
    t0_ = Translation3d(-Vector3d(ligand_.translation)) *
          u::eulerRotation(receptor_.rotation);
    t1_ = Translation3d(Vector3d(receptor_.translation)) *
          u::eulerRotation(ligand_.rotation, true);
    t2_ = u::eulerRotation(ligand_.rotation) *
          Translation3d(-Vector3d(ligand_.translation));
  }
}

} // namespace zdock
