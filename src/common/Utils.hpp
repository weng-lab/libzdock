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
#include <algorithm>
#include <chrono>
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

  // return now
  static inline auto tic() { return std::chrono::high_resolution_clock::now(); }

  // return now - t1
  static inline auto toc(std::chrono::high_resolution_clock::time_point t1) {
    auto t2 = tic();
    return std::chrono::duration_cast<std::chrono::microseconds>(t2 - t1)
               .count() /
           1000000.0;
  }
};

} // namespace zdock
