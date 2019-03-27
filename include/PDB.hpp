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

#include "pdb++.h"
#include <Eigen/Dense>
#include <iostream>
#include <memory>
#include <mutex>
#include <vector>

namespace zdock {

class RecordCoord {
public:
  int serialNum;
  std::string atomName;
  std::string resName;
  char chain;
  int resNum;
  RecordCoord() : serialNum(0), atomName(""), resName(""), chain('\0'), resNum(0) {}
};

class Model;

class PDB {
private:
  void read_(const std::string &fn);

public:
  typedef Eigen::Matrix<double, 3, Eigen::Dynamic> Matrix;
  typedef Eigen::Transform<double, 3, Eigen::Affine> Transform;
  typedef Eigen::Vector3d Coord;
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
  const Record &operator[](const RecordCoord &coord) const;
  void append(const libpdb::PDB &, const int model = 0);
  void append(const Record &, const int model = 0);
  Coord centroid() const;
};

inline std::ostream &operator<<(std::ostream &s, const RecordCoord &c) {
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
