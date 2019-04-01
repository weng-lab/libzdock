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
