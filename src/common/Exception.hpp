// Copyright (c) 2019-2026 Arjan van der Velde, Weng Lab
// SPDX-License-Identifier: BSD-2-Clause

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
