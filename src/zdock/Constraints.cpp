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

#include "Constraints.hpp"
#include "Utils.hpp"
#include <fstream>
#include <regex>

namespace zdock {

Constraints::Constraints(const std::string &filename)
    : filename_(Utils::realpath(filename)) {
  std::ifstream infile(filename_);
  std::string line;
  int linenum = 1;
  if (infile.is_open()) {
    Constraint c;
    try {
      while (infile >> c) {
        cons_.push_back(c);
        linenum++;
      }
    } catch (const ConstraintException &e) {
      throw ConstraintException(
          "Error reading constraints (line: " + std::to_string(linenum) + ")");
    }
  } else {
    throw ConstraintException("Error opening file '" + filename_ + "'");
  }
}

const std::vector<Constraint> &Constraints::constraints() const {
  return cons_;
}

std::istream &operator>>(std::istream &s, Constraint &c) {
  const std::regex r("^([0-9]+)\\s+([A-Z0-9]+)\\s+([A-Z0-9]+)\\s+(.)\\s+([0-9]+"
                     ")\\s+([0-9]+)\\s+([A-Z0-9]+)\\s+([A-Z0-9]+)\\s+(.)\\s+(["
                     "0-9]+)\\s+([0-9.]+)(\\s+(MIN|MAX))?$");
  std::smatch cm;
  std::string line;
  std::string ctype;
  if (std::getline(s, line)) {
    if (std::regex_match(line, cm, r)) {
      try {
        c.recCoord.serialNum = std::stoi(cm[1]);
        c.recCoord.atomName = cm[2];
        c.recCoord.resName = cm[3];
        c.recCoord.chain = cm[4].str().c_str()[0];
        c.recCoord.resNum = std::stoi(cm[5]);
        c.ligCoord.serialNum = std::stoi(cm[6]);
        c.ligCoord.atomName = cm[7];
        c.ligCoord.resName = cm[8];
        c.ligCoord.chain = cm[9].str().c_str()[0];
        c.ligCoord.resNum = std::stoi(cm[10]);
        c.distance = std::stod(cm[11]);
        c.constraintType =
            ("MIN" == std::string(cm[13]) ? Constraint::MIN : Constraint::MAX);
      } catch (const std::invalid_argument& e) {
        throw ConstraintException("Constraint format error; invalid argument");
      }
    } else {
      throw ConstraintException("Constraint format error");
    }
  }
  return s;
}

} // namespace zdock
