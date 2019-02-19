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

#include "ZdockPruning.hpp"
#include <string>

namespace zdock {
class Pruning {
private:
  const std::string zdockfn_;
  const std::string ligfn_;
  const double cutoff_;

public:
  Pruning(const std::string &zdockfn, const std::string ligandfn,
          const double cutoff);
  void doPrune();
};
} // namespace zdock
