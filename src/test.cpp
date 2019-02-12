#include "Pruning.hpp"

int main() {
  const std::string zfile = "572b42d9-8406-41df-8200-fabd347a9548/zdock.out.ranked";
  const std::string receptor = "572b42d9-8406-41df-8200-fabd347a9548/receptor.pdb";
  const std::string ligand = "572b42d9-8406-41df-8200-fabd347a9548/ligand.pdb";

  try {
    zdock::Pruning p(zfile); //, receptor, ligand);
    //p.makeComplex(0);
    p.prune(6.0);
  } catch (const std::exception &e) {
    std::cerr << "Exception: " << e.what() << std::endl;
    exit(1);
  }
}

