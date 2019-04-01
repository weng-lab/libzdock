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

#include "Exception.hpp"
#include "TransformLigand.hpp"
#include "TransformMultimer.hpp"
#include "ZDOCK.hpp"
#include <Eigen/Dense>
#include <string>

namespace zdock {

/**
 * @brief Perform RMSD based pruning on (M-)ZDOCK output
 */
class Pruning {
private:
  typedef Eigen::Transform<double, 3, Eigen::Affine> Transform;
  typedef Eigen::Matrix<double, 3, Eigen::Dynamic> Matrix;

  ZDOCK zdock_;                 // zdock output
  const double cutoff_;         // cutoff
  const TransformLigand txl_;   // ligand tranfomation class
  const TransformMultimer txm_; // multimertranfomation class
  std::string strucfn_;         // receptor and ligand filenames
  const bool getclusters_;      // return all w/ cluster number in score

  // results
  std::vector<int> clusters_; // cluster assignments
  size_t strucsize_;          // structure size
  int nclusters_;             // number of clusters

public:
  /**
   * @brief Constructor
   *
   * @param zdockoutput ZDOCK or M-ZDOCK output file name
   * @param cutoff RMSD cutoff
   * @param structurefn Structure PDB file name
   * @param getclusters Toggle return for full (M-)ZDOCK output with cluster numbers for scores
   */
  Pruning(
      const std::string &zdockoutput, const double cutoff,
      const std::string &structurefn = "", // or grab from zdock.out
      const bool getclusters = false // return all w/ cluster number in score
  );

  /**
   * @brief Actually perform pruning
   */
  void prune();
  /**
   * @brief Get cluster assignments
   *
   * @return vector of cluster numbers, one for each prediction
   */
  const std::vector<int> &clusters() const { return clusters_; }
  /**
   * @brief Get number of clusters
   *
   * @return number of clusters found
   */
  int nclusters() const { return nclusters_; }
  /**
   * @brief Get ZDOCK output with cluster numbers for scores
   *
   * @return ZDOCK output with cluster numbers for scores
   */
  const ZDOCK &zdock() const { return zdock_; }
};

/**
 * @brief General exception during pruning
 */
class PruningException : public Exception {
public:
  PruningException(const std::string &msg) : Exception(msg) {}
};

} // namespace zdock
