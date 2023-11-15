#include <Eigen/Dense>
#include <cmath>
#include <ctime>
#include <iomanip>
#include <iostream>
#include <unsupported/Eigen/FFT>

#include "utilities.h"

#ifndef ABELINTEGRALEQUATION_H_
#define ABELINTEGRALEQUATION_H_

using namespace Eigen;
using namespace std;

namespace AbelIntegralEquation {
/* @brief Compute Gaussian quadrature nodes and weights for n nodes over
 * interval [a,b] \param[in] a,b Interval [a,b] endpoints \param[in] n Number of
 * quadrature points \param[out] xq,wq Pair of Gauss-Legendre quadrature points
 * and weights
 */
std::pair<Eigen::RowVectorXd, Eigen::RowVectorXd> gauleg(double a, double b,
                                                         int n,
                                                         double eps = 1.e-13) {
  const double PI = 4. * std::atan(1.);  // PI
  int i, j, m;
  double xmid, xlen, p1, p2, p3, dp1, z, z1, wqi;
  Eigen::RowVectorXd xq(n), wq(n);

  m = (n + 1) / 2;
  xmid = 0.5 * (a + b);
  xlen = 0.5 * (b - a);

  // get roots
  for (i = 0; i < m; i++) {
    // i-th root guess
    z = std::cos(PI * (i + 1 - 0.25) / (n + 0.5));

    // get i-th root
    do {
      p1 = 1.;
      p2 = 0.;
      for (j = 1; j <= n; j++) {
        p3 = p2;
        p2 = p1;
        p1 = ((2. * j - 1.) * z * p2 - (j - 1.) * p3) / j;
      }
      dp1 = n * (z * p1 - p2) / (z * z - 1.0);
      z1 = z;
      z = z1 - p1 / dp1;
    } while (std::abs(z - z1) > eps);

    // set nodes
    xq(i) = xmid - xlen * z;
    xq(n - 1 - i) = xmid + xlen * z;

    // set weights
    wqi = 2. * xlen / ((1. - z * z) * dp1 * dp1);
    wq(i) = wqi;
    wq(n - 1 - i) = wqi;
  }

  return std::make_pair(xq, wq);
}

/* @brief Find the unknown function u in the Abel integral equation
 * using Galerkin discretization with a polynomial basis.
 * \param y Template function for the right-hand side
 * \param p Maximum degree of the polynomial basis and
 * order of the quadrature rule to compute the righ-hand side
 * \param tau Meshwidth of the grid where to compute the values of u
 * \\return Values of u on a grid in [0,1] with meshwidth tau
 */
/* SAM_LISTING_BEGIN_0 */
template <typename FUNC>
VectorXd poly_spec_abel(const FUNC& y, size_t p, double tau) {
  MatrixXd A = MatrixXd::Zero(p + 1, p + 1);
  VectorXd b = VectorXd::Zero(p + 1);

  // generate Gauss-Legendre points and weights
  Eigen::RowVectorXd gauss_pts_p, gauss_wht_p;
  std::tie(gauss_pts_p, gauss_wht_p) = gauleg(0., 1., p);

  // set-up the Galerkin matrix and rhs vector

// **********************************************************************
// PROBLEM 3-2f:

  // We take the solution of 3-2e and now we need to implement it here.
  // Note that for the gamma function std::tgamma can be used!
  for (int i{0}; i <= p; ++i) {
    for (int j{0}; j <= p; ++j) {
      // Set up for the Galerkin Matrix based on subproblem 3-2e:
      A(i,j) = std::sqrt(M_PI) * std::tgamma(j + 1) / ((i + j + 3./2.) * std::tgamma(j + 3./2.));
    }
    // Next: set up of the rhs vector of the LSE (Linear System of Equations) based on 
    // Gauss-Legendre quadrature
    for (int k{0}; k < p; ++k) {
      const double t_k = gauss_pts_p[k]; 
      const double w_k = gauss_wht_p[k];
      // Now based on the formula 3.2.9, also note that b here is equivalent to phi!
      b(i) += w_k * std::pow(t_k, i) * y(t_k);
    }
  }

// **********************************************************************/

  // linear system solve using QR decomposition
  VectorXd x = A.colPivHouseholderQr().solve(b);

  size_t N = round(1. / tau);
  VectorXd grid = VectorXd::LinSpaced(N + 1, 0., 1.);
  VectorXd u = VectorXd::Zero(N + 1);

  // generate solution at grid points
  for (int i = 0; i <= N; i++) {
    for (int j = 0; j <= p; j++) {
      u(i) += x(j) * pow(grid(i), j);
    }
  }

  return u;
}
/* SAM_LISTING_END_0 */
/* @brief Find the unknown function u in the Abel integral equation
 * using convolution quadrature (implicit Euler)
 * \param y Template function for the right-hand side
 * \param N Number of discretization steps
 * \\return Values of u from convolution quadrature
 */
/* SAM_LISTING_BEGIN_2 */
template <typename FUNC>
Eigen::VectorXd cq_ieul_abel(const FUNC& y, size_t N) {
  Eigen::VectorXd u(N + 1);
// **********************************************************************
// PROBLEM 3-2j:
// Convolution Quadrature Implicit Euler
// Calculation weights of Convolution Quadrature based on 3-2h!
  Eigen::VectorXd w(N + 1);
  w(0) = 1; // Initial Value (Fromula 3-2h)
  // Note: achieved through Taylor Series Expansion!
  for (int l{1}; l < N + 1; ++l) {
    w(l) = w(l - 1) * (static_cast<double>(l) - 0.5) / static_cast<double>(l);
  }
  // Note that the sqrt(tau * pi) where tau = 1/N is missing:
  w *= M_PI / N;
  // Now: To solve the CQ Grid Points have to be defined first!
  Eigen::VectorXd grid = Eigen::VectorXd::LinSpaced(N+1, 0., 1.);
  // Set up the RHS Vector of the LSE (Linear System of Equations)
  Eigen::VectorXd y_N(N + 1);
  for (int i{0}; i < N + 1; ++i) {
    y_N(i) = y(grid(i));
  }
  // Set up the coefficient Matrix:
  Eigen::MatrixXd T = toeplitz_triangular(w);
  // Solve the Linear System Equation with Eigen's build in triangular solver!
  u = T.triangularView<Lower>().solve(y_N);
// **********************************************************************/
  return u;
}
/* SAM_LISTING_END_2 */

/* @brief Find the unknown function u in the Abel integral equation
 * using convolution quadrature (BDF-2)
 * \param y Template function for the right-hand side
 * \param N Number of discretization steps
 * \\return Values of u from convolution quadrature
 */
/* SAM_LISTING_BEGIN_3 */
template <typename FUNC>
Eigen::VectorXd cq_bdf2_abel(const FUNC& y, size_t N) {
  Eigen::VectorXd u(N + 1);
// **********************************************************************
// PROBLEM 3-2k:
// Note: Calculations rely on the subproblem 3-2i !
  Eigen::VectorXd w_1(N + 1); // First Factor
  Eigen::VectorXd w_2(N + 1); // Second Factor
  // Same approach for the computation of the weights as in example 3-2j
  // Taylor Series Expansion
  w_1(0) = w_2(0) = 1;
  for (int i{1}; i < N + 1; ++i) {
    w_1(i) = w_1(i - 1) * (static_cast<double>(i) - 0.5) / static_cast<double>(i);
    w_2(i) = w_1(i) / std::pow(3, i);
  }
  // Full Expansion by Cauchy Product:
  // The weights can be computed by a discrete convolution of the coefficients
  // of the two Taylor Series Expansions.
  Eigen::VectorXd w = myconv(w_1, w_2).head(N+1).real();
  // Note that the sqrt(2/3) and sqrt(pi/N) has to be included
  w *= std::sqrt(2./3.) * std::sqrt(M_PI/N);

  // Again a Grid is required to solve the LSE (Linear System of Equations)
  Eigen::VectorXd grid = Eigen::VectorXd::LinSpaced(N + 1, 0., 1.);
  // Also again setting up the rhs vector of the LSE
  Eigen::VectorXd y_N(N + 1);
  for (int i{0}; i < N + 1; ++i) {
    y_N(i) = y(grid(i));
  }
  // Set up the coefficient matrix
  Eigen::MatrixXd T = toeplitz_triangular(w);
  // Also here the LSE gets solved with Eigen's built in triangular eliminator solver:
  u = T.triangularView<Lower>().solve(y_N);
// **********************************************************************/
  return u;
}
/* SAM_LISTING_END_3 */

}  // namespace AbelIntegralEquation
#endif
