// Copyright (c) 2019-2026 Arjan van der Velde, Weng Lab
// SPDX-License-Identifier: BSD-2-Clause

#pragma once

#include "Exception.hpp"
#include "ZDOCK.hpp"
#include <string>
#include <vector>

namespace zdock {

/**
 * @brief Reconstitute a ZDOCK output file from multiple smaller files
 */
class UnSplit {
private:
  //! list of file names
  const std::vector<std::string> files_;

public:
  /**
   * @brief Constructor
   *
   * @param files list of file names
   */
  UnSplit(const std::vector<std::string> &files);
  //! actually perform the concatenation
  void unsplit();
};

} // namespace zdock
