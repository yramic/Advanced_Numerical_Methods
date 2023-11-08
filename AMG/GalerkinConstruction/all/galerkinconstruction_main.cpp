/**
 * @ file galerkinconstruction_main.cpp
 * @ brief NPDE homework GalerkinConstruction MAIN FILE
 * @ author
 * @ date
 * @ copyright Developed at SAM, ETH Zurich
 */

#include "galerkinconstruction.h"

int main(int /*argc*/, char** /*argv*/) {
  std::cout << "Running code for ADVNCSE HW GalerkinConstruction" << std::endl;

  GalerkinConstruction::tabulateRuntimes(std::vector<double>({8,16,32,64,128,256,512,1024,2048,4096,9192}));
  
  return 0;
}
