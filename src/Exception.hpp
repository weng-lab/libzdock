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

class AtomNotFoundException : public std::exception {
private:
  const std::string what_;

public:
  AtomNotFoundException(const std::string &msg) : what_(msg) {}
  const char *what() const throw() { return what_.c_str(); }
};

class ConstraintException : public std::exception {
private:
  const std::string what_;

public:
  ConstraintException(const std::string &msg) : what_(msg) {}
  const char *what() const throw() { return what_.c_str(); }
};

class PathException : public std::exception {
private:
  const std::string what_;

public:
  PathException(const std::string &msg) : what_(msg) {}
  const char *what() const throw() { return what_.c_str(); }
};

class PDBOpenException : public std::exception {
private:
  const std::string what_;

public:
  PDBOpenException(const std::string &fn)
      : what_("Error opening '" + fn + "'") {}
  const char *what() const throw() { return what_.c_str(); }
};

class ZDOCKInvalidFormat : public std::exception {
private:
  const std::string what_;

public:
  ZDOCKInvalidFormat(const std::string &fn, const std::string &msg = "")
      : what_("Error opening '" + fn + "'" + ("" != msg ? ": " + msg : "")) {}
  const char *what() const throw() { return what_.c_str(); }
};

class ZDOCKUnsupported : public std::exception {
private:
  const std::string what_;

public:
  ZDOCKUnsupported(const std::string &msg = "") : what_(msg) {}
  const char *what() const throw() { return what_.c_str(); }
};

} // namespace zdock
