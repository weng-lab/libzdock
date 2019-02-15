#include "Pruning.hpp"
#include "Constraints.hpp"
#include <iostream>

int main(int argc, char **argv) {
  try {
    if (argc > 3) {
      const std::string zfile = std::string(argv[1]);
      const double cutoff = std::stod(argv[2]);
      const std::string cfile = std::string(argv[3]);
      zdock::Pruning p(zfile);
//      p.prune(cutoff);
      p.filterConstraints(cfile);
 //     std::cout << cutoff << std::endl;
    } else {
      std::cerr << "Usage: " << argv[0] << " <zdock-output-file> <cut-off> <constrains-file>" << std::endl;
      exit(1);
    }
  } catch (const std::exception &e) {
    std::cerr << "Exception: " << e.what() << std::endl;
    exit(1);
  }
}

