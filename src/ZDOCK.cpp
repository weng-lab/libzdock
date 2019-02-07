#include "ZDOCK.hpp"
#include "Exception.hpp"
#include <fstream>
#include <iostream>
#include <sstream>

namespace zlab {

ZDOCK::ZDOCK(const std::string &fn)
    : boxsize_(0), spacing_(0.0), isswitched_(false), ismzdock_(false),
      isfixed_(false), version_(0), filename_(fn) {
  read_();
}

void ZDOCK::read_() {
  std::vector<std::string> header;
  std::ifstream infile(filename_);
  std::string line;
  if (infile.is_open()) {
    predictions_.clear();
    while (std::getline(infile, line)) {
      prediction p;
      p.ismzdock = false;
      if (7 == std::sscanf(line.c_str(), "%lf\t%lf\t%lf\t%d\t%d\t%d\t%lf",
                           &p.rotation[0], &p.rotation[1], &p.rotation[2],
                           &p.translation[0], &p.translation[1],
                           &p.translation[2], &p.score)) {
        if (ismzdock_) {
          // M-ZDOCK established but 7-column prediction encountered
          throw ZDOCKInvalidFormat(filename_, "Invalid M-ZDOCK prediction");
        }
        predictions_.push_back(p);
      } else if (5 == std::sscanf(line.c_str(), "%lf\t%lf\t%d\t%d\t%lf",
                                  &p.rotation[0], &p.rotation[1],
                                  &p.translation[0], &p.translation[1],
                                  &p.score)) {
        if (!ismzdock_ && predictions_.size() > 0) {
          // ZDOCK established but 5-column prediction encountered
          throw ZDOCKInvalidFormat(filename_, "Invalid ZDOCK prediction");
        }
        ismzdock_ = true;
        p.ismzdock = true;
        predictions_.push_back(p);
      } else {
        header.push_back(line);
      }
    }
  }

  // figure out header
  if (!ismzdock_) {
    if (5 == header.size()) {
      int sw;
      version_ = 1; // new style has five header rows
      if (3 != std::sscanf(header[0].c_str(), "%d\t%lf\t%d", &boxsize_,
                           &spacing_, &sw)) {
        throw ZDOCKInvalidFormat(filename_, "ZDOCK header error");
      }
      isswitched_ = !!sw;
      isfixed_ = false;
    } else if (4 == header.size()) {
      version_ = 0; // old style has 4 header rows
      isfixed_ = true;
      isswitched_ = false;
      if (2 !=
          std::sscanf(header[0].c_str(), "%d\t%lf", &boxsize_, &spacing_)) {
        throw ZDOCKInvalidFormat(filename_, "ZDOCK header error");
      }
    } else {
      // ZDOCK _must_ have 4 or 5 header rows
      throw ZDOCKInvalidFormat(filename_, "ZDOCK header error");
    }
  } else {
    if (3 != header.size()) {
      // M-ZDOCK _must_ have 3 header rows
      throw ZDOCKInvalidFormat(filename_, "M-ZDOCK header error");
    }
  }

  // receptor
  int count =
      std::sscanf(header[1].c_str(), "%lf\t%lf\t%lf", &receptor_.rotation[0],
                  &receptor_.rotation[1], &receptor_.rotation[2]);
  if (3 != count) {
    throw ZDOCKInvalidFormat(filename_,
                             "Unable to obtain receptor initial rotation");
  }
  try {
    std::stringstream ss(header[(isswitched_ ? 4 : 3) - (!version_)]);
    ss >> receptor_.filename;
    ss >> receptor_.translation[0];
    ss >> receptor_.translation[1];
    ss >> receptor_.translation[2];
  } catch (const std::exception &e) {
    throw ZDOCKInvalidFormat(filename_,
                             "Unable to obtain receptor initial translation");
  }

  // ligand
  if (!ismzdock_) {
    if (0 != version_) {
      if (3 != std::sscanf(header[2].c_str(), "%lf\t%lf\t%lf",
                           &ligand_.rotation[0], &ligand_.rotation[1],
                           &ligand_.rotation[2])) {
        throw ZDOCKInvalidFormat(filename_,
                                 "Unable to obtain receptor initial rotation");
      }
    }
    try {
      std::stringstream ss(header[(isswitched_ ? 3 : 4) - (!version_)]);
      ss >> ligand_.filename;
      ss >> ligand_.translation[0];
      ss >> ligand_.translation[1];
      ss >> ligand_.translation[2];
    } catch (const std::exception &e) {
      throw ZDOCKInvalidFormat(filename_,
                               "Unable to obtain receptor initial translation");
    }
  }
}

// get ligand (unsupported for M-ZDOCK)
const structure &ZDOCK::ligand() const {
  if (ismzdock_) {
    throw ZDOCKUnsupported("ligand() not supported for M-ZDOCK output");
  }
  return ligand_;
}

// text representation of structure (in (m-)zdock format)
std::ostream &operator<<(std::ostream &os, const structure &obj) {
  char s[256];
  snprintf(s, sizeof(s), "%s\t%.3f\t%.3f\t%.3f", obj.filename.c_str(),
           obj.translation[0], obj.translation[1], obj.translation[2]);
  os << s;
  return os;
}

// text representation of prediction (actual (m-)zdock.out format)
std::ostream &operator<<(std::ostream &os, const prediction &obj) {
  char s[256];
  if (obj.ismzdock) {
    snprintf(s, sizeof(s), "%.3f\t%.3f\t%d\t%d\t%.2f", obj.rotation[0],
             obj.rotation[1], obj.translation[0], obj.translation[1],
             obj.score);
  } else {
    snprintf(s, sizeof(s), "%.3f\t%.3f\t%.3f\t%d\t%d\t%d\t%.3f",
             obj.rotation[0], obj.rotation[1], obj.rotation[2],
             obj.translation[0], obj.translation[1], obj.translation[2],
             obj.score);
  }
  os << s;
  return os;
}

} // namespace zlab

