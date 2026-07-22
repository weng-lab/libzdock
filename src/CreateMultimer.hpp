// Copyright (c) 2019-2026 Arjan van der Velde, Weng Lab
// SPDX-License-Identifier: BSD-2-Clause

#pragma once

#include "Exception.hpp"
#include "ZDOCK.hpp"

namespace zdock {

/**
 * @brief Create multimer, or component thereof, for a M-ZDOCK prediction
 */
class CreateMultimer {
private:
  /**
   * @brief M-ZDOCK output file name
   */
  const std::string zdockfn_;
  /**
   * @brief Receptor PDB file name
   */
  const std::string structurefn_;
  /**
   * @brief Prediction number in M-ZDOCK output file (1-based)
   */
  const size_t n_;
  /**
   * @brief Component to generate (0-based)
   */
  const int mer_;
  /**
   * @brief Toggle whether to return all PDB records rather than just ATOM/HETATM
   */
  const bool allrecords_;

public:
  /**
   * @brief Constructor
   *
   * @param zdockoutput M-ZDOCK output file name
   * @param structure Structure PDB file name
   * @param n Prediction number in M-ZDOCK output file (1-based)
   * @param mer Component number if single component required
   * @param allrecords Toggle whether to return all PDB records rather than just ATOM/HETATM
   */
  CreateMultimer(const std::string &zdockoutput, const std::string &structure,
                 const size_t n, const int mer, const bool allrecords);
  /**
   * @brief Actually perform multimer creation
   */
  void doCreate();
};

/**
 * @brief General exception in CreateMultimer
 */
class CreateMultimerException : public Exception {
private:
  const std::string what_;

public:
  CreateMultimerException(const std::string &msg) : Exception(msg) {}
};

} // namespace zdock
