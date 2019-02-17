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

#include "Constraints.hpp"
#include "PDB.hpp"
#include "TransformLigand.hpp"
#include "ZDOCK.hpp"

namespace zdock {

class Pruning {
private:
  typedef Eigen::Transform<double, 3, Eigen::Affine> Transform;
  typedef Eigen::Matrix<double, 3, Eigen::Dynamic> Matrix;

  ZDOCK zdock_;               // zdock output
  const TransformLigand txl_; // ligand tranfomation class
  std::string recfn_, ligfn_; // receptor and ligand filenames

public:
  Pruning(const std::string &zdockoutput,
          const std::string &receptorpdb = "", // or grab from zdock.out
          const std::string &ligandpdb = ""    // or grab from zdock.out
  );

  // perform pruning
  void prune(const double cutoff);

  // ligand pdb record to stdout
  void makeComplex(const size_t n);
  void makeZDOCKComplex(const size_t n);
  void makeMZDOCKComplex(const size_t n);

  // filter predictions based on constraints
  void filterConstraints(const std::string &fn);
};

} // namespace zdock
