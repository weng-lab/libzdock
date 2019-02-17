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
#include "Pruning.hpp"
#include "TransformMultimer.hpp"
#include "Utils.hpp"
#include "ZDOCK.hpp"
#include <iostream>
#include <regex>

int main(int argc, char **argv) {
  try {
    if (argc > 1) {
      const std::string zfile = std::string(argv[1]);
      const zdock::ZDOCK z(zfile);
      const zdock::TransformMultimer txm(z);

      zdock::PDB s(zdock::Utils::copath(zfile, z.receptor().filename));
      zdock::Prediction p = z.predictions()[0];

      Eigen::Matrix<double, 3, Eigen::Dynamic> m;
      int serial = 0;
      for (int n = 0; n < 24; ++n) {
        zdock::PDB pdb = s;
        m = txm.txMultimer(pdb, p, n);
        pdb.setMatrix(m);
        for (auto x : pdb.records()) {
          x.atom.serialNum = ++serial;
          x.atom.residue.chainId =
              zdock::TransformMultimer::CHAINS[n].c_str()[0];
          std::cout << x << '\n';
        }
      }

      // const double cutoff = std::stod(argv[2]);
      // const std::string cfile = std::string(argv[3]);
      // zdock::Pruning p(zfile);
      // p.prune(cutoff);
      // p.filterConstraints(cfile);
      //     std::cout << cutoff << std::endl;
      // p.makeComplex(0);
    } else {
      std::cerr << "Usage: " << argv[0]
                << " <zdock-output-file> <cut-off> <constrains-file>"
                << std::endl;
      exit(1);
    }
  } catch (const std::exception &e) {
    std::cerr << "Exception: " << e.what() << std::endl;
    exit(1);
  }
}
