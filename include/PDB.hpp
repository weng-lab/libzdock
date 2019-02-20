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
#include <memory>
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

class Model;

class PDB {
private:
  void read_(const std::string &fn);

public:
  typedef Eigen::Matrix<double, 3, Eigen::Dynamic> Matrix;
  typedef Eigen::Transform<double, 3, Eigen::Affine> Transform;
  typedef std::shared_ptr<libpdb::PDB> Record;
  typedef std::shared_ptr<zdock::Model> Model;

protected:
  std::vector<Model> models_;   // zero or more models
  std::vector<Record> records_; // all records
  std::vector<Record> atoms_;   // just atoms
  Matrix matrix_;               // eigen matrix w/ atom coords
  // atomic inserts...
  std::mutex lock_;
  // atom filter
  const std::function<bool(const libpdb::PDB &)> filter_;

public:
  PDB();
  PDB(const PDB &p);
  PDB(const std::string &filename,
      std::function<bool(const libpdb::PDB &)> filter =
          [](const libpdb::PDB &) { return true; });
  PDB &operator=(const PDB &p);
  const Matrix &matrix() const;
  const Matrix &setMatrix(const Matrix &m);
  const std::vector<Model> &models() const;
  size_t nmodels() const;
  const std::vector<Record> &records() const;
  const std::vector<Record> &atoms() const;
  const Record &operator[](const int serial) const;
  const Record &operator[](const Coord &coord) const;
  void append(const libpdb::PDB &, const int model = 0);
  void append(const Record &, const int model = 0);
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

class Model : public PDB {
private:
  const std::vector<Model> &models() const = delete;
  const std::vector<Record> &records() const = delete;
  int modelNum_;

public:
  int modelNum() const { return modelNum_; }
  friend void PDB::append(const Record &r, const int model);
};

} // namespace zdock
