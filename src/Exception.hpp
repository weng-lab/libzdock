#pragma once

#include <exception>
#include <string>

namespace zlab {

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

class ZDOCKUnsupported: public std::exception {
private:
  const std::string what_;

public:
  ZDOCKUnsupported(const std::string &msg = "")
      : what_(msg) {}
  const char *what() const throw() { return what_.c_str(); }
};


} // namespace zlab
