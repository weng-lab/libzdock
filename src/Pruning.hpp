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

#include "Exception.hpp"
#include "TransformLigand.hpp"
#include "ZDOCK.hpp"
#include <Eigen/Dense>
#include <string>

namespace zdock {

class ZdockPruning {
private:
  typedef Eigen::Transform<double, 3, Eigen::Affine> Transform;
  typedef Eigen::Matrix<double, 3, Eigen::Dynamic> Matrix;

  ZDOCK zdock_;               // zdock output
  const double cutoff_;       // cutoff
  const TransformLigand txl_; // ligand tranfomation class
  std::string ligfn_;         // receptor and ligand filenames

  // results
  std::vector<int> clusters_; // cluster assignments
  size_t ligsize_;            // ligand size
  int nclusters_;             // number of clusters

public:
  ZdockPruning(const std::string &zdockoutput, const double cutoff,
               const std::string &ligandpdb = "" // or grab from zdock.out
  );

  // perform pruning
  void prune();
  const std::vector<int> &clusters() const { return clusters_; }
  int nclusters() const { return nclusters_; }
  const ZDOCK &zdock() const { return zdock_; }
};

class ZdockPruningException : public Exception {
public:
  ZdockPruningException(const std::string &msg) : Exception(msg) {}
};

} // namespace zdock
