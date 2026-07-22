// Copyright (c) 2019-2026 Arjan van der Velde, Weng Lab
// SPDX-License-Identifier: BSD-2-Clause

#include "TransformLigand.hpp"

namespace e = Eigen;
namespace p = libpdb;

namespace zdock {

TransformLigand::TransformLigand(const std::string &zdock)
    : TransformLigand(ZDOCK(zdock)) {}

TransformLigand::TransformLigand(const ZDOCK &zdock)
    : zdock_(zdock), spacing_(0.0), boxsize_(0), rev_(false), fixed_(false),
      isvalid_(false) {

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
