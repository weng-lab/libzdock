#include "Pruning.hpp"

int main(int argc, char **argv) {
  try {
    if (argc > 2) {
      const std::string zfile = std::string(argv[1]);
      const double cutoff = std::stod(argv[2]);
      zdock::Pruning p(zfile); //, receptor, ligand);
      // p.makeComplex(0);
      p.prune(cutoff);
    } else {
      std::cerr << "Usage: " << argv[0] << " <zdock-output-file> <cut-off>" << std::endl;
      exit(1);
    }
  } catch (const std::exception &e) {
    std::cerr << "Exception: " << e.what() << std::endl;
    exit(1);
  }
}

