/**
 * @ file fractionalheatequation_main.cpp
 * @ brief NPDE homework FractionalHeatEquation MAIN FILE
 * @ author Jörg Nick
 * @ date October 2023
 * @ copyright Developed at SAM, ETH Zurich
 */

#include "fractionalheatequation.h"

int main(int /*argc*/, char** /*argv*/) {
  Eigen::VectorXd w = FractionalHeatEquation::cqWeights(5, 0.1);
  std::cout << "Running code for ADVNCSE HW FractionalHeatEquation" << std::endl
            << w;
  return 0;
}
