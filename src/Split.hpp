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
#include "ZDOCK.hpp"

namespace zdock {

class Split {
private:
  const std::string zdockfn_;
  const int chunksize_;
  const std::string prefix_;
  std::string suffix(const uint64_t i) const;
public:
  Split(const std::string &zdockfn, const int chunksize,
        const std::string &prefix);
  void split();
};

class SplitException : public Exception {
public:
  SplitException(const std::string &msg) : Exception(msg) {}
};

} // namespace zdock
