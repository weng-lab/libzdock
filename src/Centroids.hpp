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
  const zdock::pdb::AtomRecord templateAtom_ = [] {
    zdock::pdb::AtomRecord atom;
    atom.serialNum = 0;
    atom.name = "N";
    atom.altLoc = '\0';
    atom.residue.name = "HOH";
    atom.residue.chainId = 'Z';
    atom.residue.seqNum = 0;
    atom.residue.insertCode = '\0';
    atom.xyz = {0.0, 0.0, 0.0};
    atom.occupancy = 1.0;
    atom.tempFactor = 0.0;
    atom.element = "N";
    atom.charge = "";
    return atom;
  }();

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
