/**
 * Copyright 2019 Arjan van der Velde, vandervelde.ag [at] gmail
 * Licensed under the Apache License, Version 2.0 (the "License");
 * you may not use this file except in compliance with the License.
 * You may obtain a copy of the License at
 * http://www.apache.org/licenses/LICENSE-2.0
 *
 * Unless required by applicable law or agreed to in writing, software
 * distributed under the License is distributed on an "AS IS" BASIS,
 * WITHOUT WARRANTIES OR CONDITIONS OF ANY KIND, either express or implied.
 * See the License for the specific language governing permissions and
 * limitations under the License.
 */

#include "Exception.hpp"
#include "PDB.hpp"
#include <string>

namespace zdock {

class Centroids {
private:
  const std::string zdockfn_;
  const std::string ligandfn_;
  const size_t n_;
  const libpdb::PDB::Atom templateAtom_ = {
      0,    // serialNum
      "N",  // name
      '\0', // altLoc
      {
          "ALA",       // name,
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
  Centroids(const std::string &zdockoutput, const std::string &ligand,
            const size_t n);
  void doCentroids();
};

class CentroidsException : public Exception {
public:
  CentroidsException(const std::string &msg) : Exception(msg) {}
};

} // namespace zdock
