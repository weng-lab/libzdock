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
