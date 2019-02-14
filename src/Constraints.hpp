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

inline std::istream &operator>>(std::istream &s, Constraint &c) {
  std::string line;
  std::string ctype;
  if (std::getline(s, line)) {
    std::istringstream ss(line);
    if (!(ss >> c.recCoord.serialNum >> c.recCoord.atomName >>
          c.recCoord.resName >> c.recCoord.chain >> c.recCoord.resNum >>
          c.ligCoord.serialNum >> c.ligCoord.atomName >> c.ligCoord.resName >>
          c.ligCoord.chain >> c.ligCoord.resNum >> c.distance)) {
      throw ConstraintException("Error reading constraint.");
    }
    if (ss >> ctype && "MIN" == ctype) {
      c.constraintType = Constraint::MIN;
    } else {
      c.constraintType = Constraint::MAX;
    }
  }
  return s;
}

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
