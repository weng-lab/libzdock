#include "Constraints.hpp"
#include "Utils.hpp"
#include <fstream>

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
}
