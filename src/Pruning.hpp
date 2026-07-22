// Copyright (c) 2019-2026 Arjan van der Velde, Weng Lab
// SPDX-License-Identifier: BSD-2-Clause

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
