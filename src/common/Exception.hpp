/**
 * Copyright (c) 2019, Arjan van der Velde, Weng Lab
 * All rights reserved.
 *
 * Redistribution and use in source and binary forms, with or without
 * modification, are permitted provided that the following conditions are met:
 *
 * * Redistributions of source code must retain the above copyright notice, this
 *   list of conditions and the following disclaimer.
 *
 * * Redistributions in binary form must reproduce the above copyright notice,
 *   this list of conditions and the following disclaimer in the documentation
 *   and/or other materials provided with the distribution.
 *
 * THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND CONTRIBUTORS "AS IS"
 * AND ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT LIMITED TO, THE
 * IMPLIED WARRANTIES OF MERCHANTABILITY AND FITNESS FOR A PARTICULAR PURPOSE ARE
 * DISCLAIMED. IN NO EVENT SHALL THE COPYRIGHT HOLDER OR CONTRIBUTORS BE LIABLE
 * FOR ANY DIRECT, INDIRECT, INCIDENTAL, SPECIAL, EXEMPLARY, OR CONSEQUENTIAL
 * DAMAGES (INCLUDING, BUT NOT LIMITED TO, PROCUREMENT OF SUBSTITUTE GOODS OR
 * SERVICES; LOSS OF USE, DATA, OR PROFITS; OR BUSINESS INTERRUPTION) HOWEVER
 * CAUSED AND ON ANY THEORY OF LIABILITY, WHETHER IN CONTRACT, STRICT LIABILITY,
 * OR TORT (INCLUDING NEGLIGENCE OR OTHERWISE) ARISING IN ANY WAY OUT OF THE USE
 * OF THIS SOFTWARE, EVEN IF ADVISED OF THE POSSIBILITY OF SUCH DAMAGE.
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
