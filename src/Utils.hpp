#pragma once

#include "Exception.hpp"
#include <algorithm>
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

  // https://stackoverflow.com/questions/216823/whats-the-best-way-to-trim-stdstring/217605#217605

  // trim from start (in place)
  static inline void ltrim(std::string &s) {
    s.erase(s.begin(), std::find_if(s.begin(), s.end(),
                                    [](int ch) { return !std::isspace(ch); }));
  }

  // trim from end (in place)
  static inline void rtrim(std::string &s) {
    s.erase(std::find_if(s.rbegin(), s.rend(),
                         [](int ch) { return !std::isspace(ch); })
                .base(),
            s.end());
  }

  // trim from both ends (in place)
  static inline void trim(std::string &s) {
    ltrim(s);
    rtrim(s);
  }

  // trim from start (copying)
  static inline std::string ltrim_copy(std::string s) {
    ltrim(s);
    return s;
  }

  // trim from end (copying)
  static inline std::string rtrim_copy(std::string s) {
    rtrim(s);
    return s;
  }

  // trim from both ends (copying)
  static inline std::string trim_copy(std::string s) {
    trim(s);
    return s;
  }
};

} // namespace zdock
