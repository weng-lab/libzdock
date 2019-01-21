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

} // namespace zlab
