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
  bool isfixed() const { return isfixed_; }
  int version() const { return version_; }
  int symmetry() const;
  const std::string &filename() const { return filename_; }

  // structures
  Structure &receptor() { return receptor_; }
  Structure &ligand();

  // predictions
  size_t npredictions() const { return predictions_.size(); }
  std::vector<Prediction> &predictions() { return predictions_; }

  // stream out
  friend std::ostream &operator<<(std::ostream &os, const ZDOCK &obj);
};

} // namespace zdock