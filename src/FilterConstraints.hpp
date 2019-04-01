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
