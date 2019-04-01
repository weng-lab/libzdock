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
   * @brief Constructur
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
