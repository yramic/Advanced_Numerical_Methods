/**
 * @ file mfmg_main.cpp
 * @ brief NPDE homework MFMG MAIN FILE
 * @ author
 * @ date October 2023
 * @ copyright Developed at SAM, ETH Zurich
 */

#include "mfmg.h"

int main(int /*argc*/, char** /*argv*/) {
  std::cout << "Running code for ADVNCSE HW MFMG" << std::endl;
  MFMG::tabulateMGConvergenceRate();
  return 0;
}
