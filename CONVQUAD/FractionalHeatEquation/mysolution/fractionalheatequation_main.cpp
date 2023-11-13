/**
 * @ file fractionalheatequation_main.cpp
 * @ brief NPDE homework FractionalHeatEquation MAIN FILE
 * @ author Jörg Nick, Bob Schreiner
 * @ date October 2023
 * @ copyright Developed at SAM, ETH Zurich
 */

#include <chrono>
#include <fstream>

#include "fractionalheatequation.h"

/* SAM_LISTING_BEGIN_0 */
int main(int /*argc*/, char** /*argv*/) {
// **********************************************************************
// Your Solution here
// **********************************************************************/
=======
  Eigen::VectorXd mu_MOT =
      FractionalHeatEquation::evlMOT(f, n, T, std::pow(2, L) - 1);
  Eigen::VectorXd mu_Toep =
      FractionalHeatEquation::evlTriangToeplitz(f, n, T, L);
  Eigen::VectorXd mu_ASAO = FractionalHeatEquation::evlASAOCQ(f, n, T, L);
  std::cout << "mu_MOT: \n" << mu_MOT << std::endl << std::endl;
  std::cout << "mu_Toep: \n" << mu_Toep << std::endl << std::endl;
  std::cout << "mu_ASAO: \n" << mu_ASAO << std::endl << std::endl;
/* SAM_LISTING_BEGIN_0 */
// **************************
// Code for runtime mesaurement     
// **************************
/* SAM_LISTING_END_0 */
  
>>>>>>> 9f23489a85abae1a4183f950200f027aadd5db8b
  return 0;
}
/* SAM_LISTING_END_0 */