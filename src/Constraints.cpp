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
      } catch (const std::invalid_argument e) {
        throw ConstraintException("Constraint format error; invalid argument");
      }
    } else {
      throw ConstraintException("Constraint format error");
    }
  }
  return s;
}

} // namespace zdock
