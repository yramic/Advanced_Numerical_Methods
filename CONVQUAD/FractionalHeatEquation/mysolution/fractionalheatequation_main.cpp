/**
 * @ file fractionalheatequation_main.cpp
 * @ brief NPDE homework FractionalHeatEquation MAIN FILE
 * @ author Jörg Nick
 * @ date October 2023
 * @ copyright Developed at SAM, ETH Zurich
 */

#include "fractionalheatequation.h"

int main(int /*argc*/, char** /*argv*/) {
  FractionalHeatEquation::SqrtsMplusA Amat(2, std::complex<double>(1,1));

  std::cout << "Running code for ADVNCSE HW FractionalHeatEquation" << std::endl;
       //     << Amat.p_matrix_;
  return 0;
}
