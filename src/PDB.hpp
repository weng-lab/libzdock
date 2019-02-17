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

#include "pdb++.h"
#include <Eigen/Dense>
#include <iostream>
#include <mutex>
#include <vector>

namespace zdock {

class Coord {
public:
  int serialNum;
  std::string atomName;
  std::string resName;
  char chain;
  int resNum;
  Coord() : serialNum(0), atomName(""), resName(""), chain('\0'), resNum(0) {}
};

class PDB {
private:
  typedef Eigen::Matrix<double, 3, Eigen::Dynamic> Matrix;
  typedef Eigen::Transform<double, 3, Eigen::Affine> Transform;

  std::vector<size_t> atoms_;
  std::vector<libpdb::PDB> records_;
  Matrix m_;
  void read_(const std::string &fn, const int model = MODEL_FIRST,
             std::function<bool(const libpdb::PDB &)> filter =
                 [](const libpdb::PDB &) { return true; });
  std::mutex insertmtx_;

public:
  static const int MODEL_ALL = -1;
  static const int MODEL_FIRST = 0;
  PDB();
  PDB(const PDB &p);
  PDB(const std::string &fn, const int model = MODEL_FIRST,
      std::function<bool(const libpdb::PDB &)> filter =
          [](const libpdb::PDB &) { return true; });
  PDB &operator=(const PDB &p);
  const Matrix &matrix() const;
  const Matrix &transform(const Transform &t);
  const Matrix &setMatrix(const Matrix &m);
  const std::vector<libpdb::PDB> &records();
  const Eigen::Vector3d centroid() const;
  const libpdb::PDB &operator[](const int serial) const;
  const libpdb::PDB &operator[](const Coord &coord) const;
  void append(const libpdb::PDB &);
};

inline std::ostream &operator<<(std::ostream &s, const Coord &c) {
  std::ostringstream os;
  os << c.serialNum << '\t';
  os << c.atomName << '\t';
  os << c.resName << '\t';
  os << c.chain << '\t';
  os << c.resNum;
  s << os.str();
  return s;
}

} // namespace zdock
