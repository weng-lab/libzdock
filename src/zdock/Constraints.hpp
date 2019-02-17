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
#include "PDB.hpp"
#include <sstream>
#include <string>
#include <vector>

namespace zdock {

class Constraint {
public:
  Coord recCoord;
  Coord ligCoord;
  double distance;
  enum { MIN, MAX } constraintType;
  Constraint() : recCoord(), ligCoord(), distance(0.0), constraintType(MAX) {}
};

class Constraints {
private:
  const std::string filename_;
  std::vector<Constraint> cons_;

public:
  Constraints() {}
  Constraints(const std::string &filename);
  const std::vector<Constraint> &constraints() const;
};

std::istream &operator>>(std::istream &s, Constraint &c);

inline std::ostream &operator<<(std::ostream &s, const Constraint &c) {
  std::ostringstream os;
  os << c.recCoord << '\t';
  os << c.ligCoord << '\t';
  os << c.distance << '\t';
  os << (Constraint::MAX == c.constraintType ? "MAX" : "MIN");
  s << os.str();
  return s;
}

} // namespace zdock
