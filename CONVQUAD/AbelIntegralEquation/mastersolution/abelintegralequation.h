#include <Eigen/Dense>
#include <cmath>
#include <ctime>
#include <iomanip>
#include <iostream>
#include <unsupported/Eigen/FFT>

#include "utilities.h"

#ifndef ABELINTEGRALEQUATION_H_
#define ABELINTEGRALEQUATION_H_

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
Eigen::VectorXd poly_spec_abel(const FUNC& y, std::size_t p, double tau) {
  Eigen::MatrixXd A = MatrixXd::Zero(p + 1, p + 1);
  Eigen::VectorXd b = Eigen::VectorXd::Zero(p + 1);

  // generate Gauss-Legendre quadrature points and weights
  const auto [gauss_pts_p, gauss_wht_p] = gauleg(0., 1., p);

  // set up the Galerkin matrix and rhs vector

  // Precompute factors related to gamma function
  // Watch out for integer division
  const Eigen::VectorXd gamma_factor =
      Eigen::VectorXd::NullaryExpr(p + 1, [](Eigen::Index j) {
        return std::tgamma(j + 1) / std::tgamma(j + 3. / 2.);
      });
  const double sqrt_pi = std::sqrt(M_PI);

  for (int i = 0; i <= p; i++) {
    // Set up Galerkin matrix based on \prbcref{subprb:sg2}
    for (int j = 0; j <= p; j++) {
      A(i, j) = sqrt_pi * gamma_factor[j] / (1.5 + i + j);
    }

    // Set up the rhs vector of the LSE based on Gauss-Legendre quadrature
    for (int k = 0; k < p; k++) {
      const double tk = gauss_pts_p(k);
      const double wk = gauss_wht_p(k);
      // Based on \prbeqref{eq:aieAent}
      b(i) += wk * std::pow(tk, i) * y(tk);
    }
  }

  // linear system solve using QR decomposition
  const Eigen::ColPivHouseholderQR<Eigen::MatrixXd> solver(A);
  const Eigen::VectorXd x = solver.solve(b);

  // generate points on the grid
  const std::size_t N = std::round(1. / tau);
  const Eigen::VectorXd grid = Eigen::VectorXd::LinSpaced(N + 1, 0., 1.);
  Eigen::VectorXd u = Eigen::VectorXd::Zero(N + 1);

  // generate solution at grid points
  for (int j = 0; j <= p; j++) {
    u += x(j) * grid.array().pow(j).matrix();
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
  Eigen::VectorXd w(N + 1);
  w(0) = 1.;
  // Calculate weights of convolution quadrature based on \prbcref{subprb:cq1}
  for (int l = 1; l < N + 1; ++l) {
    w(l) = w(l - 1) * (l - 0.5) / l;  // denominator is factorial
  }
  w *= std::sqrt(M_PI / N);

  // Solve the convolution quadrature:
  // Generate points on the grid
  Eigen::VectorXd grid = Eigen::VectorXd::LinSpaced(N + 1, 0., 1.);
  // Set up the rhs vector of the LSE
  Eigen::VectorXd y_N(N + 1);
  for (int i = 0; i < N + 1; ++i) {
    y_N(i) = y(grid(i));
  }
  // Set up the coefficient matrix
  Eigen::MatrixXd T = toeplitz_triangular(w);
  // Solve the lse with Eigen's build-in triangular elimination solver
  Eigen::VectorXd u = T.triangularView<Lower>().solve(y_N);
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
  Eigen::VectorXd w1(N + 1);
  w1(0) = 1.;
  // Calculate weights of convolution quadrature based on \prbcref{subprb:cq2}
  // Taylor series expansion of the first factor
  for (int l = 1; l < N + 1; ++l) {
    w1(l) = w1(l - 1) * (l - 0.5) / l;  // denominator is factorial
  }

  // Taylor series expansion of the second factor
  Eigen::VectorXd w2 = w1;
  for (int l = 1; l < N + 1; ++l) {
    w2(l) /= pow(3, l);
  }

  // Full expansion by Cauchy product
  Eigen::VectorXd w = myconv(w1, w2).head(N + 1).real();
  w *= std::sqrt(M_PI / N) * std::sqrt(2. / 3.);

  // Solve the convolution quadrature:
  // Generate points on the grid
  Eigen::VectorXd grid = Eigen::VectorXd::LinSpaced(N + 1, 0., 1.);
  // Set up the rhs vector of the LSE
  Eigen::VectorXd y_N(N + 1);
  for (int i = 0; i < N + 1; ++i) {
    y_N(i) = y(grid(i));
  }

  // Set up the coefficient matrix
  Eigen::MatrixXd T = toeplitz_triangular(w);
  // Solve the lse with Eigen's build-in triangular elimination solver
  Eigen::VectorXd u = T.triangularView<Lower>().solve(y_N);

  return u;
}
/* SAM_LISTING_END_3 */

}  // namespace AbelIntegralEquation
#endif
