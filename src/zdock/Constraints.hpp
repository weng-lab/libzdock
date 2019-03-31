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
  //! Possible constraint types (MIN for mimimum and MAX for maximum distance)
  enum {
    MIN, /**< minimum distance constraint */
    MAX  /**< maximum distance constraint */
  } constraintType;
  //! Constructor
  Constraint() : recCoord(), ligCoord(), distance(0.0), constraintType(MAX) {}
};

/**
 * @brief Constraints parser
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

//! input stream for constraint
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
