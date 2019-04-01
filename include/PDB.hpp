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

//! Exact 'pointer' to a record in a PDB file
class RecordCoord {
public:
  //! ATOM/HETATM serial number
  int serialNum;
  //! Atom name
  std::string atomName;
  //! Residue name
  std::string resName;
  //! Chain ID
  char chain;
  //! Residue number
  int resNum;
  /**
   * @brief Constructor
   */
  RecordCoord() : serialNum(0), atomName(""), resName(""), chain('\0'), resNum(0) {}
};

class Model;

/**
 * @brief Collection of PDB records, representing a PDB file
 */
class PDB {
private:
  /**
   * @brief Read from file
   *
   * @param fn file name
   */
  void read_(const std::string &fn);

public:
  //! PDB coordinate matrix type
  typedef Eigen::Matrix<double, 3, Eigen::Dynamic> Matrix;
  //! Shorthand for eigen Transformation type
  typedef Eigen::Transform<double, 3, Eigen::Affine> Transform;
  //! Coordinate (x, y, z)
  typedef Eigen::Vector3d Coord;
  //! Shared Pointer type for PDB record
  typedef std::shared_ptr<libpdb::PDB> Record;
  //! Shared Pointer type for Model
  typedef std::shared_ptr<zdock::Model> Model;

protected:
  std::vector<Model> models_;   //!< zero or more models
  std::vector<Record> records_; //!< all records
  std::vector<Record> atoms_;   //!< just atoms
  Matrix matrix_;               //!< eigen matrix w/ atom coords
  // atomic inserts...
  std::mutex lock_; //!< lock for atomic updates
  //! Atom filter function
  const std::function<bool(const libpdb::PDB &)> filter_;

public:
  PDB(); //!< Constructor
  PDB(const PDB &p); //!< Copy constructor
  /**
   * @brief Constructor
   * @param filename PDB file name to read from
   * @param filter filter function for ATOM/HETATM records
   */
  PDB(const std::string &filename,
      std::function<bool(const libpdb::PDB &)> filter =
          [](const libpdb::PDB &) { return true; });
  /**
   * @brief Assignement operator
   * @param p other PDB object
   * @return reference to *this, updated from p
   */
  PDB &operator=(const PDB &p);
  //! get coordinate matrix
  const Matrix &matrix() const;
  //! set coordinate matrix
  const Matrix &setMatrix(const Matrix &m);
  //! get models
  const std::vector<Model> &models() const;
  //! get number of models
  size_t nmodels() const;
  //! get all records (see Record)
  const std::vector<Record> &records() const;
  //! get atom records (see Record)
  const std::vector<Record> &atoms() const;
  //! get ATOM/HETATM by serial
  const Record &operator[](const int serial) const;
  //! get ATOM/HETATM by RecordCoord coordinate
  const Record &operator[](const RecordCoord &coord) const;
  //! append record, from actual object, optionally to model by number
  void append(const libpdb::PDB &, const int model = 0);
  //! append record, from shared pointer, optionally to model by number
  void append(const Record &, const int model = 0);
  //! get centroid (i.e. mean x, y, z) of strcuture
  Coord centroid() const;
};

// output stream representation of RecordCoord
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


/**
 * @brief Model, a sub-PDB structure
 */
class Model : public PDB {
private:
  //! a Model cannot itself have more models
  const std::vector<Model> &models() const = delete;
  //! a Model contains only atom records
  const std::vector<Record> &records() const = delete;
  //! model number of this model
  int modelNum_;

public:
  int modelNum() const { return modelNum_; }
  friend void PDB::append(const Record &r, const int model);
};

} // namespace zdock
