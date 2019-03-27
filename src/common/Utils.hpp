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
    if ('/' == file[0]) {
      return realpath(file); // file is absuolute
    }
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
