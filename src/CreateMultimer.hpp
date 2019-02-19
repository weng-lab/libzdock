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

#pragma once

#include "Exception.hpp"
#include "ZDOCK.hpp"

namespace zdock {

class CreateMultimer {
private:
  const std::string zdockfn_;
  const std::string receptorfn_;
  const size_t n_;
  const int mer_;

public:
  CreateMultimer(const std::string &zdockoutput, const std::string &receptor,
                 const size_t n, const int mer);
  void doCreate();
};

class CreateMultimerException : public Exception {
private:
  const std::string what_;

public:
  CreateMultimerException(const std::string &msg) : Exception(msg) {}
};

} // namespace zdock
