#include <Eigen/Dense>
#include <cmath>
#include <fstream>
#include <iostream>
#include <istream>

#include "genLaguerreRule.hpp"

// Quadrature rule struct
struct QuadRule {
  std::size_t dim;    // dimension of space
  std::size_t n;      // number of nodes/weights
  Eigen::MatrixXd x;  // quadrature nodes (columns of a matrix with dim rows)
  Eigen::VectorXd w;  // vector of quadrature weights
};

static const int KIND = 5;

//------------------------------------------------------------------------------
/* SAM_LISTING_BEGIN_0 */
bool testGaussLaguerre(int N) {
  // get Gauss-Laguerre quadrature points and weights
  double *w = new double[N];
  double *x = new double[N];
  cgqf(N, KIND, 0, 0, 0, 1, x, w);

  // TODO: Implement your code
  return true;
}
/* SAM_LISTING_END_0 */

//------------------------------------------------------------------------------
/* SAM_LISTING_BEGIN_1 */
Eigen::VectorXd testGaussLaguerreConvergence(int N = 20) {
  Eigen::VectorXd error(N);
  // TODO: Implement your code
  return error;
}
/* SAM_LISTING_END_1 */

//------------------------------------------------------------------------------
/* SAM_LISTING_BEGIN_2 */
QuadRule getLogWeightQR(double a, int n) {
  // get Gauss-Laguerre quadrature points and weights
  double *w = new double[n];
  double *x = new double[n];
  cgqf(n, KIND, 0, 0, 0, 1, x, w);

  // create new Quadrature rule
  QuadRule logWeightQR;
  // TODO: Implement your code
  return logWeightQR;
}
/* SAM_LISTING_END_2 */

//------------------------------------------------------------------------------
/* SAM_LISTING_BEGIN_3 */
Eigen::VectorXd testLogWeightQRConvergence(int n) {
  Eigen::VectorXd error(n);
  // TODO: Implement your code
  return error;
}
/* SAM_LISTING_END_3 */

//------------------------------------------------------------------------------
/* SAM_LISTING_BEGIN_4 */
/* SAM_LISTING_END_4 */

//------------------------------------------------------------------------------
int main() {
  //-------------------------------------------------------------
  /* SAM_LISTING_BEGIN_5 */
  std::cout << "TESTING GAUSS-LAGUERRE QUADRATURE" << std::endl;
  for (int j = 1; j < 30; j++) {
    bool pass = testGaussLaguerre(j);
    if (pass) {
      std::cout << "Test passed for N = " << j << std::endl;
    } else {
      std::cout << "Test for N = " << j << " failed " << std::endl;
    }
  }
  /* SAM_LISTING_END_5 */
  std::cout << std::endl;

  //-------------------------------------------------------------
  /* SAM_LISTING_BEGIN_6 */
  std::cout << "TESTING QUADRATURES' CONVERGENCE" << std::endl << std::endl;
  Eigen::VectorXd errorGL = testGaussLaguerreConvergence(30);
  Eigen::VectorXd errorLW = testLogWeightQRConvergence(30);
  Eigen::VectorXi N = Eigen::VectorXi::LinSpaced(30, 1, 30);
  // output for plot
  std::ofstream out_errorGL("GLQR_errors.txt");
  out_errorGL << std::setprecision(18) << errorGL;
  out_errorGL.close();
  std::ofstream out_errorLW("WLQR_errors.txt");
  out_errorLW << std::setprecision(18) << errorLW;
  out_errorLW.close();
  std::ofstream out_N("QR_N.txt");
  out_N << N;
  out_N.close();
  /* SAM_LISTING_END_6 */

  //-------------------------------------------------------------
  /* SAM_LISTING_BEGIN_7 */
  std::cout << "TESTING LOG-WEIGHT QUADRATURE" << std::endl;
  // TODO: Implement your code

  return 0;
}
