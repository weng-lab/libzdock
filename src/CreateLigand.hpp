// Copyright (c) 2019-2026 Arjan van der Velde, Weng Lab
// SPDX-License-Identifier: BSD-2-Clause

#pragma once

#include "Exception.hpp"
#include "ZDOCK.hpp"

namespace zdock {

/**
 * @brief Create transformed ligand, or whole structure, for a ZDOCK prediction.
 */
class CreateLigand {
private:
  /**
   * @brief ZDOCK output file name
   */
  const std::string zdockfn_;
  /**
   * @brief Ligand PDB file name
   */
  const std::string ligandfn_;
  /**
   * @brief Receptor PDB file name
   */
  const std::string receptorfn_;
  /**
   * @brief Prediction number in ZDOCK file (1-based)
   */
  const size_t n_;
  /**
   * @brief Toggle generation of full complex (requires #receptorfn_)
   */
  const bool complex_;
  /**
   * @brief Toggle return of all input PDB records vs only ATOM/HETATM records
   */
  const bool allrecords_;

public:
  /**
   * @brief Constructor
   *
   * @param zdockoutput ZDOCK output file name
   * @param ligand Ligand PDB file name
   * @param receptor Receptor PDB file name
   * @param n Prediction number in ZDOCK file (1-based)
   * @param cmplx Toggle generation of full complex
   * @param allrecords Toggle return of all input PDB records
   */
  CreateLigand(const std::string &zdockoutput, const std::string &ligand,
               const std::string &receptor, const size_t n, const bool cmplx,
               const bool allrecords);
  /**
   * @brief Actually perform generation of ligand/complex
   */
  void doCreate();
};

/**
 * @brief General exception for CreateLigand
 */
class CreateLigandException : public Exception {
private:
  const std::string what_;

public:
  CreateLigandException(const std::string &msg) : Exception(msg) {}
};

} // namespace zdock
