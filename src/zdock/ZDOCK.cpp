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

#include "ZDOCK.hpp"
#include "Exception.hpp"
#include <fstream>
#include <iostream>
#include <sstream>

namespace zdock {

ZDOCK::ZDOCK(const std::string &fn)
    : boxsize_(0), spacing_(0.0), isswitched_(false), ismzdock_(false),
      isfixed_(false), version_(0), symmetry_(0), filename_(fn) {
  read_();
}

void ZDOCK::read_() {
  std::vector<std::string> header;
  std::ifstream infile(filename_);
  std::string line;
  if (infile.is_open()) {
    predictions_.clear();
    bool headerdone = false;
    int linenum = 1;
    while (std::getline(infile, line)) {
      Prediction p;
      p.ismzdock = false;
      if (7 == std::sscanf(line.c_str(), "%lf\t%lf\t%lf\t%d\t%d\t%d\t%lf",
                           &p.rotation[0], &p.rotation[1], &p.rotation[2],
                           &p.translation[0], &p.translation[1],
                           &p.translation[2], &p.score)) {
        if (ismzdock_) {
          // M-ZDOCK established but 7-column prediction encountered
          throw ZDOCKInvalidFormat(filename_,
                                   "Invalid M-ZDOCK prediction (line " +
                                       std::to_string(linenum) + ")");
        }
        headerdone = true;
        predictions_.push_back(p);
      } else if (5 == std::sscanf(line.c_str(), "%lf\t%lf\t%d\t%d\t%lf",
                                  &p.rotation[0], &p.rotation[1],
                                  &p.translation[0], &p.translation[1],
                                  &p.score)) {
        if (!ismzdock_ && predictions_.size() > 0) {
          // ZDOCK established but 5-column prediction encountered
          throw ZDOCKInvalidFormat(filename_,
                                   "Invalid ZDOCK prediction (line " +
                                       std::to_string(linenum) + ")");
        }
        p.rotation[2] = 0.0;
        p.translation[2] = 0.0;
        ismzdock_ = true;
        p.ismzdock = true;
        headerdone = true;
        predictions_.push_back(p);
      } else if (!headerdone) {
        header.push_back(line);
      } else {
        throw ZDOCKInvalidFormat(filename_, "Invalid prediction (line " +
                                                std::to_string(linenum) + ")");
      }
      linenum++;
    }
  }

  // figure out header
  if (!ismzdock_) {
    if (5 == header.size()) {
      int sw;
      version_ = 1; // new style has five header rows
      if (3 != std::sscanf(header[0].c_str(), "%d\t%lf\t%d", &boxsize_,
                           &spacing_, &sw)) {
        throw ZDOCKInvalidFormat(filename_, "ZDOCK header error");
      }
      isswitched_ = !!sw;
      isfixed_ = false;
    } else if (4 == header.size()) {
      version_ = 0; // old style has 4 header rows
      isfixed_ = true;
      isswitched_ = false;
      if (2 !=
          std::sscanf(header[0].c_str(), "%d\t%lf", &boxsize_, &spacing_)) {
        throw ZDOCKInvalidFormat(filename_, "ZDOCK header error");
      }
    } else {
      // ZDOCK _must_ have 4 or 5 header rows
      throw ZDOCKInvalidFormat(filename_, "ZDOCK header error");
    }
  } else {
    if (3 != header.size()) {
      // M-ZDOCK _must_ have 3 header rows
      throw ZDOCKInvalidFormat(filename_, "M-ZDOCK header error");
    }
    if (3 != std::sscanf(header[0].c_str(), "%d\t%lf\t%d", &boxsize_, &spacing_,
                         &symmetry_)) {
      throw ZDOCKInvalidFormat(filename_, "M-ZDOCK header error");
    }
    if (3 > symmetry_) {
      throw ZDOCKInvalidFormat(filename_,
                               "M-ZDOCK symmetry cannot be less than 3");
    }
  }

  // receptor
  if (!isfixed_) {
    receptor_.rotation[2] = 0.0; // be sure to set this
    int count =
        std::sscanf(header[1].c_str(), "%lf\t%lf\t%lf", &receptor_.rotation[0],
                    &receptor_.rotation[1], &receptor_.rotation[2]);
    if ((3 != count && !ismzdock_) || (count < 2 && ismzdock_)) {
      throw ZDOCKInvalidFormat(filename_, "Unable to obtain receptor initial rotation");
    }
  }
  try {
    std::stringstream ss(header[(isswitched_ ? 4 : 3) - (!version_)]);
    ss >> receptor_.filename;
    ss >> receptor_.translation[0];
    ss >> receptor_.translation[1];
    ss >> receptor_.translation[2];
  } catch (const std::exception &e) {
    throw ZDOCKInvalidFormat(filename_,
                             "Unable to obtain receptor initial translation");
  }

  // ligand
  if (!ismzdock_) {
    if (3 != std::sscanf(header[(isfixed_ ? 1 : 2)].c_str(), "%lf\t%lf\t%lf",
                         &ligand_.rotation[0], &ligand_.rotation[1],
                         &ligand_.rotation[2])) {
      throw ZDOCKInvalidFormat(filename_,
                               "Unable to obtain receptor initial rotation");
    }
    try {
      std::stringstream ss(header[(isswitched_ ? 3 : 4) - (!version_)]);
      ss >> ligand_.filename;
      ss >> ligand_.translation[0];
      ss >> ligand_.translation[1];
      ss >> ligand_.translation[2];
    } catch (const std::exception &e) {
      throw ZDOCKInvalidFormat(filename_,
                               "Unable to obtain receptor initial translation");
    }
  }
}

// get m-zdock symmetry
int ZDOCK::symmetry() const {
  if (!ismzdock_) {
    throw ZDOCKUnsupported("symmetry() not supported for ZDOCK output");
  }
  return symmetry_;
}

// get ligand (unsupported for M-ZDOCK)
Structure &ZDOCK::ligand() {
  if (ismzdock_) {
    throw ZDOCKUnsupported("ligand() not supported for M-ZDOCK output");
  }
  return ligand_;
}

// get ligand (unsupported for M-ZDOCK)
const Structure &ZDOCK::ligand() const {
  if (ismzdock_) {
    throw ZDOCKUnsupported("ligand() not supported for M-ZDOCK output");
  }
  return ligand_;
}

// text representation of structure (in (m-)zdock format)
std::ostream &operator<<(std::ostream &os, const Structure &obj) {
  char s[256];
  snprintf(s, sizeof(s), "%s\t%.3f\t%.3f\t%.3f", obj.filename.c_str(),
           obj.translation[0], obj.translation[1], obj.translation[2]);
  os << s;
  return os;
}

// text representation of prediction (actual (m-)zdock.out format)
std::ostream &operator<<(std::ostream &os, const Prediction &obj) {
  char s[256];
  if (obj.ismzdock) {
    snprintf(s, sizeof(s), "%.6f\t%.6f\t%d\t%d\t%.2f", obj.rotation[0],
             obj.rotation[1], obj.translation[0], obj.translation[1],
             obj.score);
  } else {
    snprintf(s, sizeof(s), "%.6f\t%.6f\t%.6f\t%d\t%d\t%d\t%.3f",
             obj.rotation[0], obj.rotation[1], obj.rotation[2],
             obj.translation[0], obj.translation[1], obj.translation[2],
             obj.score);
  }
  os << s;
  return os;
}

// text representation of m-)zdock
std::ostream &operator<<(std::ostream &os, const ZDOCK &obj) {
  char s[256];
  if (obj.ismzdock_) {
    // M-ZDOCK header 1, 2
    snprintf(s, sizeof(s), "%d\t%.1f\t%d\n%.6f\t%.6f\t%.6f\n", obj.boxsize_,
             obj.spacing_, obj.symmetry_, obj.receptor_.rotation[0],
             obj.receptor_.rotation[1], obj.receptor_.rotation[2]);
  } else if (obj.isfixed_) {
    // ZDOCK (fixed) header 1, 2
    snprintf(s, sizeof(s), "%d\t%.1f\n%.6f\t%.6f\t%.6f\n", obj.boxsize_,
             obj.spacing_, obj.ligand_.rotation[0], obj.ligand_.rotation[1],
             obj.ligand_.rotation[2]);
  } else {
    // ZDOCK (new format / not fixed) header 1, 2
    snprintf(s, sizeof(s), "%d\t%.1f\t%d\n%.6f\t%.6f\t%.6f\n%.6f\t%.6f\t%.6f\n",
             obj.boxsize_, obj.spacing_, obj.isswitched_,
             obj.receptor_.rotation[0], obj.receptor_.rotation[1],
             obj.receptor_.rotation[2], obj.ligand_.rotation[0],
             obj.ligand_.rotation[1], obj.ligand_.rotation[2]);
  }
  os << s;
  // receptor / ligand filename and translation
  os << (obj.isswitched_ ? obj.ligand_ : obj.receptor_);
  if (!obj.ismzdock_) {
    // ligand / receptor filename and translation
    os << '\n' << (obj.isswitched_ ? obj.receptor_ : obj.ligand_);
  }
  // predictions
  for (const auto &x : obj.predictions_) {
    os << '\n' << x;
  }
  return os;
}

} // namespace zdock
