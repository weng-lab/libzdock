// Copyright (c) 2019-2026 Arjan van der Velde, Weng Lab
// SPDX-License-Identifier: BSD-2-Clause

#pragma once

#include "Exception.hpp"
#include "ZDOCK.hpp"

#include <cstdint>

namespace zdock {

/**
 * @brief Split a ZDOCK or M-ZDOCK file into smaller pieces
 */
class Split {
private:
  const std::string zdockfn_; //!< zdock file name
  const int chunksize_;       //!< number of lines per chunk
  const std::string prefix_;  //!< prefix for output file names
  /**
   * @brief generate suffix file number
   * @param i file number
   * @return a suffix
   */
  std::string suffix(const uint64_t i) const;
public:
  /**
   * @brief Constructor
   *
   * @param zdockfn ZDOCK file name
   * @param chunksize Number of lines per output file
   * @param prefix prefix for output files
   */
  Split(const std::string &zdockfn, const int chunksize,
        const std::string &prefix);
  //! actually perform the split
  void split();
};

/**
 * @brief General exception in Split
 */
class SplitException : public Exception {
public:
  SplitException(const std::string &msg) : Exception(msg) {}
};

} // namespace zdock
