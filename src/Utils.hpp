#pragma once

#include "Exception.hpp"
#include <cstdlib>
#include <string>

namespace zdock {
class Utils {
public:
  static std::string realpath(const std::string &file) {
    char *path = ::realpath(file.c_str(), NULL);
    if (path) {
      std::string ret(path);
      free(path);
      return ret;
    } else {
      throw PathException("realpath not found for '" + file + "'");
    }
  }
  static std::string copath(const std::string &path, const std::string &file) {
    return realpath(path + "/../" + file);
  }
};

} // namespace zdock
