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

#include <string>
#include <vector>

namespace zdock {

class Structure {
public:
  Structure()
      : filename(""), translation{0.0, 0.0, 0.0}, rotation{0.0, 0.0, 0.0} {}
  std::string filename;
  double translation[3]; // actual
  double rotation[3];
  friend std::ostream &operator<<(std::ostream &os, const Structure &obj);
};

class Prediction {
public:
  Prediction()
      : rotation{0.0, 0.0, 0.0}, translation{0, 0, 0}, score(0.0),
        ismzdock(false) {}
  double rotation[3];
  int translation[3]; // grid cells
  double score;
  bool ismzdock;
  friend std::ostream &operator<<(std::ostream &os, const Prediction &obj);
};

class ZDOCK {
private:
  // structures and predictions
  Structure receptor_;
  Structure ligand_;
  std::vector<Prediction> predictions_;

  // metadata (i.e. from header)
  int boxsize_;     // box size
  double spacing_;  // grid size
  bool isswitched_; // receptor and ligand are switched
  bool ismzdock_;   // mzdock format
  bool isfixed_;    // -F was used (receptor fixed)
  int version_;     // format version; (old/mzdock:0; new: 1)
  int symmetry_;    // m-zdock symmetry
  const std::string filename_;

  // private methods
  void read_();

public:
  ZDOCK(const std::string &fn);

  // simple getters
  int boxsize() const { return boxsize_; }
  double spacing() const { return spacing_; }
  bool isswitched() const { return isswitched_; }
  bool ismzdock() const { return ismzdock_; }
  bool iszdock() const { return !ismzdock_; }
  bool isfixed() const { return isfixed_; }
  int version() const { return version_; }
  int symmetry() const;
  const std::string &filename() const { return filename_; }

  // structures
  Structure &receptor() { return receptor_; }
  const Structure &receptor() const { return receptor_; }
  Structure &ligand();
  const Structure &ligand() const;

  // predictions
  size_t npredictions() const { return predictions_.size(); }
  std::vector<Prediction> &predictions() { return predictions_; }
  const std::vector<Prediction> &predictions() const { return predictions_; }

  // stream out
  friend std::ostream &operator<<(std::ostream &os, const ZDOCK &obj);
};

} // namespace zdock
