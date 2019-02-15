#include "Constraints.hpp"
#include "PDB.hpp"
#include "TransformLigand.hpp"
#include "ZDOCK.hpp"

namespace zdock {

class Pruning {
private:
  typedef Eigen::Transform<double, 3, Eigen::Affine> Transform;
  typedef Eigen::Matrix<double, 3, Eigen::Dynamic> Matrix;

  ZDOCK zdock_;                        // zdock output
  const TransformLigand txl_;          // ligand tranfomation class
  std::string recfn_, ligfn_;          // receptor and ligand filenames

public:
  Pruning(const std::string &zdockoutput,
          const std::string &receptorpdb = "", // or grab from zdock.out
          const std::string &ligandpdb = ""    // or grab from zdock.out
  );

  // perform pruning
  void prune(const double cutoff);

  // ligand pdb record to stdout
  void makeComplex(const size_t n);

  // filter predictions based on constraints
  void filterConstraints(const std::string &fn);
};

} // namespace zdock
