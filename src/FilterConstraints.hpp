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
#include "TransformMultimer.hpp"
#include "ZDOCK.hpp"

namespace zdock {

class FilterConstraints {
private:
  typedef Eigen::Transform<double, 3, Eigen::Affine> Transform;
  typedef Eigen::Matrix<double, 3, Eigen::Dynamic> Matrix;

  ZDOCK zdock_;                       // zdock output
  const TransformLigand txl_;         // ligand tranfomation  (zdock)
  const TransformMultimer txm_;       // structure tranfomation class (m-zdock)
  std::string confn_, recfn_, ligfn_; // receptor and ligand filenames

  void filterZDOCKConstraints_();
  void filterMZDOCKConstraints_();

public:
  FilterConstraints(
      const std::string &zdockoutput,      // zdock.out file
      const std::string &constraints,      // constraints file
      const std::string &receptorpdb = "", // or grab from zdock.out
      const std::string &ligandpdb = ""    // or grab from zdock.out
  );

  // filter predictions based on constraints
  void filter();

  // get zdock file
  const ZDOCK &zdock() const { return zdock_; }
};

} // namespace zdock
