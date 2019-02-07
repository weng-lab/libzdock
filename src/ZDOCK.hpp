#pragma once

#include <string>
#include <vector>

namespace zlab {

class structure {
public:
  structure()
      : filename(""), translation{0.0, 0.0, 0.0}, rotation{0.0, 0.0, 0.0} {}
  std::string filename;
  double translation[3]; // actual
  double rotation[3];
  friend std::ostream &operator<<(std::ostream &os, const structure &obj);
};

struct prediction {
public:
  prediction()
      : rotation{0.0, 0.0, 0.0}, translation{0, 0, 0}, score(0.0),
        ismzdock(false) {}
  double rotation[3];
  int translation[3]; // grid cells
  double score;
  bool ismzdock;
  friend std::ostream &operator<<(std::ostream &os, const prediction &obj);
};

class ZDOCK {
private:
  // structures and predictions
  structure receptor_;
  structure ligand_;
  std::vector<prediction> predictions_;

  // metadata (i.e. from header)
  int boxsize_;     // box size
  double spacing_;  // grid size
  bool isswitched_; // receptor and ligand are switched
  bool ismzdock_;   // mzdock format
  bool isfixed_;    // -F was used (receptor fixed)
  int version_;     // format version; (old/mzdock:0; new: 1)
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
  const std::string &filename() const { return filename_; }

  // structures
  const structure &receptor() const { return receptor_; }
  const structure &ligand() const;

  // predictions
  size_t npredictions() const { return predictions_.size(); }
  const std::vector<prediction> &predictions() const { return predictions_; }
};

} // namespace zlab
