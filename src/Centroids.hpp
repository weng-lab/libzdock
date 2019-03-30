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
#include "PDB.hpp"
#include <string>

namespace zdock {

/**
 * @brief calculate center of mass for top-N ZDOCK predictions
 */
class Centroids {
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
   * @brief Chain ID to use for output
   */
  const std::string chain_;
  /**
   * @brief top-N centroids are produced
   */
  size_t n_;
  /**
   * @brief Template ATOM for centroids
   */
  const libpdb::PDB::Atom templateAtom_ = {
      0,    // serialNum
      "N",  // name
      '\0', // altLoc
      {
          "HOH",       // name,
          'Z',         // chainId,
          0,           // seqNum
          '\0'         // insertCode
      },               // residue
      {0.0, 0.0, 0.0}, // xyz
      1.0,             // occupancy
      0.0,             // tempFactor
      "N",             // element
      "\0",            // charge
      false,           // iszdatom
      0,               // type
      0,               // surface
      0.0,             // rad
      "\0",            // segid
      0.0              // chg
  };

public:
  /**
   * @brief Constructur
   *
   * @param zdockoutput ZDOCK output file name
   * @param ligand Ligand PDB file name
   * @param n Top-N centroids are produced
   * @param chain Chain ID to use for output
   */
  Centroids(const std::string &zdockoutput, const std::string &ligand,
            const size_t n, const std::string &chain);
  /**
   * @brief Actually perform centroids generation
   */
  void doCentroids();
};

/**
 * @brief General exception during centroid generation
 */
class CentroidsException : public Exception {
public:
  CentroidsException(const std::string &msg) : Exception(msg) {}
};

} // namespace zdock
