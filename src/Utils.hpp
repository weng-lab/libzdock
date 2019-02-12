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
  static std::string dirname(const std::string &path) {
    size_t p = path.find_last_of('/');
    if (std::string::npos != p) {
      return path.substr(0, std::min(p + 1, path.size() - 1));
    } else {
      return "";
    }
  }
  static std::string copath(const std::string &path, const std::string &file) {
    return realpath(dirname(path) + file);
  }
};

} // namespace zdock
