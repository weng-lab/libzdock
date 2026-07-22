// Copyright (c) 2019-2026 Arjan van der Velde, Weng Lab
// SPDX-License-Identifier: BSD-2-Clause

#pragma once

#include "Exception.hpp"
#include "PDB.hpp"
#include <sstream>
#include <string>
#include <vector>

namespace zdock {

  /**
   * @brief Distance constraint
   */
class Constraint {
public:
  //! Receptor (first structure) RecordCoord
  RecordCoord recCoord;
  //! Ligand (second structure) RecordCoord
  RecordCoord ligCoord;
  //! Constraint distance
  double distance;
  //! Possible constraint types (MIN for minimum and MAX for maximum distance)
  enum ConstraintType {
    MIN, //!< minimum distance constraint
    MAX  //!< maximum distance constraint
  } constraintType;
  //! Constructor
  Constraint() : recCoord(), ligCoord(), distance(0.0), constraintType(MAX) {}
};

/**
 * @brief Constraints file parser
 *
 * Constraint files are line based. Each line contains a distance constraint
 * (either minimum distance, or maximum distance) between two atoms in two
 * structures. For ZDOCK these represent the "receptor" and "ligand" structures
 * and for M-ZDOCK they refer to two atoms in the same structure.
 *
 * The format whitespace separated and looks as follows:
 *
 * ```
 * 13  OE2 GLU A   5    101  OD1 ASP b  12 7.3 MIN
 * 13  OE2 GLU A   5    101  OD1 ASP b  12 7.5
 * ```
 *
 * Column 1-5 represent:
 *
 *   - ATOM/HETATM serial (integer)
 *   - ATOM name
 *   - Residue name
 *   - Chain identifier (exactly one character)
 *   - Residue sequence number
 *
 *   Column 6-10 represent the second atom
 *   Column 11 contains the distance (double)
 *   Column 12 is optionally "MIN" or "MAX". If none is specified, "MAX" is assumed.
 *
 */
class Constraints {
private:
  //! constraints file file name
  const std::string filename_;
  //! vector of constraints
  std::vector<Constraint> cons_;

public:
  /**
   * @brief Constructor
   */
  Constraints() {}
  /**
   * @brief Constructor
   *
   * @param filename file to read from
   */
  Constraints(const std::string &filename);
  /**
   * @brief Get constraints vector
   *
   * @return vector of constraints
   */
  const std::vector<Constraint> &constraints() const;
};

//! parser of input stream for constraint (see Constraints)
std::istream &operator>>(std::istream &s, Constraint &c);

//! output stream representation of constraint
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
