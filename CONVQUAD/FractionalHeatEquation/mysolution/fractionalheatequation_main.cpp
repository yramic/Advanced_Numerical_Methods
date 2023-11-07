/**
 * @ file fractionalheatequation_main.cpp
 * @ brief NPDE homework FractionalHeatEquation MAIN FILE
 * @ author Jörg Nick, Bob Schreiner
 * @ date October 2023
 * @ copyright Developed at SAM, ETH Zurich
 */

#include "fractionalheatequation.h"

int main(int /*argc*/, char** /*argv*/) {
  FractionalHeatEquation::SqrtsMplusA Amat(2, std::complex<double>(1, 1));

  std::cout << "Running code for ADVNCSE HW FractionalHeatEquation"
            << std::endl;
  double T = 1.0;
  int L = 3;
  double tau;
  int n = 4;
  std::function<double(double, Eigen::Vector2d)> f = [](double t , Eigen::Vector2d){return 1.0;};
  Eigen::VectorXd mu_MOT = FractionalHeatEquation::evlMOT(f, n , T , std::pow(2,L));
  Eigen::VectorXd mu_Toep = FractionalHeatEquation::evlTriangToeplitz(f, n , T ,L);
  Eigen::VectorXd mu_ASAO = FractionalHeatEquation::evlASAOCQ(f, n , T ,L);
  std::cout<< "mu_MOT: \n" << mu_MOT <<std::endl;
  std::cout<< "mu_Toep: \n" << mu_Toep <<std::endl;
  std::cout<< "mu_ASAO: \n" << mu_ASAO <<std::endl;

  //     << Amat.p_matrix_;
  return 0;
}
