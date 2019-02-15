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
