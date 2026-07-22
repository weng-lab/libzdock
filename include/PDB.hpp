// Copyright (c) 2019-2026 Arjan van der Velde, Weng Lab
// SPDX-License-Identifier: BSD-2-Clause

#pragma once

#include "pdb++.h"
#include <Eigen/Dense>
#include <iostream>
#include <memory>
#include <mutex>
#include <unordered_map>
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
  typedef std::unordered_map<const libpdb::PDB *, Record> RecordMap;

  void copyFrom_(const PDB &p, RecordMap &records);

  std::vector<Model> models_;   //!< zero or more models
  std::vector<Record> records_; //!< all records
  std::vector<Record> atoms_;   //!< just atoms
  Matrix matrix_;               //!< eigen matrix w/ atom coords
  // atomic inserts...
  std::mutex lock_; //!< lock for atomic updates
  //! Atom filter function
  std::function<bool(const libpdb::PDB &)> filter_;

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
   * @brief Assignment operator
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
  //! get centroid (i.e. mean x, y, z) of structure
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
  friend class PDB;
  friend void PDB::append(const Record &r, const int model);
};

} // namespace zdock
