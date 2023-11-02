#include <Eigen/Dense>
#include <Eigen/Sparse>
#include <Eigen/SparseLU>
#include <cmath>
#include <iomanip>
#include <iostream>

#include "absorbingboundarycondition.h"
using namespace Eigen;
using namespace std;

void test_cqweights() {
  auto F = [](std::complex<double> s) { return std::pow(s, -1); };
  auto delta = [](std::complex<double> z) {
    return 1.0 / 2.0 * z * z - 2.0 * z + 3.0 / 2.0;
  };
  double tau = 0.1;
  size_t N = 10;
  Eigen::VectorXd w = AbsorbingBoundaryCondition::cqweights_by_dft(F, delta, tau, N);
  // std::cout << tau*w <<std::endl;
}

int main() {
  /* SAM_LISTING_BEGIN_3 */
  // TODO: Tabulate the H1-error of the Galerkin discretization + convolution
  // quadrature
  /* SAM_LISTING_END_3 */
}
