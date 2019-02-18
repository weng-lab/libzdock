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

#include <exception>
#include <string>

namespace zdock {

class Exception : public std::exception {
protected:
  const std::string what_;

public:
  Exception(const std::string &msg) : what_(msg) {}
  const char *what() const throw() { return what_.c_str(); }
};



class AtomNotFoundException : public Exception {
public:
  AtomNotFoundException(const std::string &msg) : Exception(msg) {}
};

class ConstraintException : public Exception {
public:
  ConstraintException(const std::string &msg) : Exception(msg) {}
};

class PathException : public Exception {
public:
  PathException(const std::string &msg) : Exception(msg) {}
};

class PDBOpenException : public Exception {
public:
  PDBOpenException(const std::string &fn)
      : Exception("Error opening PDB file '" + fn + "'") {}
};

class ZDOCKInvalidFormat : public Exception {
public:
  ZDOCKInvalidFormat(const std::string &fn, const std::string &msg = "")
      : Exception("Error opening ZDOCK file '" + fn + "'" + ("" != msg ? ": " + msg : "")) {}
};

class ZDOCKUnsupported : public Exception {
public:
  ZDOCKUnsupported(const std::string &msg = "") : Exception(msg) {}
};

} // namespace zdock
