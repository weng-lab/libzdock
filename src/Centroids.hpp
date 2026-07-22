// Copyright (c) 2019-2026 Arjan van der Velde, Weng Lab
// SPDX-License-Identifier: BSD-2-Clause

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
   * @brief Constructor
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
