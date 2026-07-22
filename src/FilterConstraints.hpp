// Copyright (c) 2019-2026 Arjan van der Velde, Weng Lab
// SPDX-License-Identifier: BSD-2-Clause

#pragma once

#include "Constraints.hpp"
#include "PDB.hpp"
#include "TransformLigand.hpp"
#include "TransformMultimer.hpp"
#include "ZDOCK.hpp"

namespace zdock {

/**
 * @brief Filter a ZDOCK or M-ZDOCK output file using distance constraints
 */
class FilterConstraints {
private:
  //! short hand for transformation
  typedef Eigen::Transform<double, 3, Eigen::Affine> Transform;
  //! short hand for coordinate matrix
  typedef Eigen::Matrix<double, 3, Eigen::Dynamic> Matrix;

  ZDOCK zdock_;                       //!< zdock output
  const TransformLigand txl_;         //!< ligand tranfomation  (zdock)
  const TransformMultimer txm_;       //!< structure tranfomation class (m-zdock)
  std::string confn_; //!< constraint file name
  std::string recfn_; //!< receptor filenames
  std::string ligfn_; //!< ligand filenames

  //! constraints filtering for ZDOCK
  void filterZDOCKConstraints_();
  //! constraints filtering for M-ZDOCK
  void filterMZDOCKConstraints_();

public:
  FilterConstraints(
      const std::string &zdockoutput,      //!< zdock.out file
      const std::string &constraints,      //!< constraints file
      const std::string &receptorpdb = "", //!< or grab from zdock.out
      const std::string &ligandpdb = ""    //!< or grab from zdock.out
  );

  //! filter predictions based on constraints
  void filter();

  //! get zdock file
  const ZDOCK &zdock() const { return zdock_; }
};

} // namespace zdock
