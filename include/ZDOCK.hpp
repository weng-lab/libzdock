// Copyright (c) 2019-2026 Arjan van der Velde, Weng Lab
// SPDX-License-Identifier: BSD-2-Clause

#pragma once

#include <string>
#include <vector>

namespace zdock {

/**
 * @brief ZDOCK structure records (i.e. recepter, ligand, structure)
 */
class Structure {
public:
  /**
   * @brief Constructor
   */
  Structure()
      : filename(""), translation{0.0, 0.0, 0.0}, rotation{0.0, 0.0, 0.0} {}
  /**
   * @brief Structure PDB file name
   */
  std::string filename;
  /**
   * @brief Structure initial translation (actual coordinate units)
   */
  double translation[3]; // actual
  /**
   * @brief Initial rotation
   */
  double rotation[3];
  /**
   * @brief stream representation of Structure
   *
   * @param os output stream
   * @param obj Structure object
   *
   * @return output stream
   */
  friend std::ostream &operator<<(std::ostream &os, const Structure &obj);
};

/**
 * @brief (M-)ZDOCK prediction record
 */
class Prediction {
public:
  /**
   * @brief Constructor
   */
  Prediction()
      : rotation{0.0, 0.0, 0.0}, translation{0, 0, 0}, score(0.0),
        ismzdock(false) {}
  /**
   * @brief Rotation
   */
  double rotation[3];
  /**
   * @brief Translation (in grid units)
   */
  int translation[3]; // grid cells
  /**
   * @brief score
   */
  double score;
  /**
   * @brief Flag to indicate this is an M-ZDOCK prediction
   */
  bool ismzdock;
  /**
   * @brief Stream representation of Prediction
   *
   * @param os output stream
   * @param obj Prediction object
   *
   * @return output stream
   */
  friend std::ostream &operator<<(std::ostream &os, const Prediction &obj);
};

/**
 * @brief (M-)ZDOCK output file parser
 */
class ZDOCK {
private:
  // structures and predictions
  /**
   * @brief Receptor structure (or structure in case of M-ZDOCK)
   */
  Structure receptor_;
  /**
   * @brief Ligand structure (in case of M-ZDOCK)
   */
  Structure ligand_;
  /**
   * @brief Vector of predictions
   */
  std::vector<Prediction> predictions_;

  // metadata (i.e. from header)
  //! box size
  int boxsize_; // box size
  //! grid spacing
  double spacing_; // grid size
  //! indicates switched receptor and ligand
  bool isswitched_; // receptor and ligand are switched
  //! indicates file is M-ZDOCK
  bool ismzdock_; // mzdock format
  //! flag whether receptor was pre-rotated or fixed
  bool isfixed_; // -F was used (receptor fixed)
  //! ZDOCK output format 'version'
  int version_; // format version; (old/mzdock:0; new: 1)
  //! symmetry (M-ZDOCK only)
  int symmetry_; // m-zdock symmetry
  //! (M-)ZDOCK file name
  const std::string filename_;

  // private methods
  /**
   * @brief read from (M-)ZDOCK output file
   */
  void read_();

public:
  /**
   * @brief Constructor
   *
   * @param fn (M-)ZDOCK output file name
   */
  ZDOCK(const std::string &fn);

  // simple getters
  /**
   * @brief Get box size
   *
   * @return box size
   */
  int boxsize() const { return boxsize_; }
  /**
   * @brief Get spacing
   *
   * @return grid spacing
   */
  double spacing() const { return spacing_; }
  /**
   * @brief Is switched?
   *
   * @return true if receptor and ligand entries are switched
   */
  bool isswitched() const { return isswitched_; }
  /**
   * @brief Is M-ZDOCK?
   *
   * @return true if file is M-ZDOCK
   */
  bool ismzdock() const { return ismzdock_; }
  /**
   * @brief is ZDOCK?
   *
   * @return true if file is ZDOCK
   */
  bool iszdock() const { return !ismzdock_; }
  /**
   * @brief Is fixed?
   *
   * @return true if receptor fixed
   */
  bool isfixed() const { return isfixed_; }
  /**
   * @brief Get version
   *
   * @return 0 for old-style and M-ZDOCK, 1 for new style
   */
  int version() const { return version_; }
  /**
   * @brief Get symmetry (M-ZDOCK only)
   *
   * @return symmetry
   */
  int symmetry() const;
  /**
   * @brief Get file name
   *
   * @return (M-)ZDOCK file name
   */
  const std::string &filename() const { return filename_; }

  // structures
  /**
   * @brief Get receptor (only for ZDOCK)
   *
   * @return receptor structure
   */
  Structure &receptor();
  /**
   * @brief Get receptor (only for ZDOCK)
   *
   * @return receptor structure
   */
  const Structure &receptor() const;
  /**
   * @brief Get ligand (only for ZDOCK)
   *
   * @return ligand structure
   */
  Structure &ligand();
  /**
   * @brief Get ligand (only for ZDOCK)
   *
   * @return ligand structure
   */
  const Structure &ligand() const;
  /**
   * @brief Get structure (only for M-ZDOCK)
   *
   * @return structure
   */
  Structure &structure();
  /**
   * @brief Get structure (only for M-ZDOCK)
   *
   * @return structure
   */
  const Structure &structure() const;

  // predictions
  /**
   * @brief Get number of predictions
   *
   * @return number of predictions
   */
  size_t npredictions() const { return predictions_.size(); }
  /**
   * @brief Get predictions
   *
   * @return vector of predictions
   */
  std::vector<Prediction> &predictions() { return predictions_; }
  /**
   * @brief Get predictions
   *
   * @return vector of predictions
   */
  const std::vector<Prediction> &predictions() const { return predictions_; }

  // stream out
  friend std::ostream &operator<<(std::ostream &os, const ZDOCK &obj);
};

} // namespace zdock
