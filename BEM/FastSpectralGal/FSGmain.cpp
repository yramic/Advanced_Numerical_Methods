#include <Eigen/Dense>
#include <cassert>
#include <cmath>
#include <fstream>
#include <iomanip>
#include <iostream>
#include <istream>

using namespace Eigen;

#if SOLUTION
#include <boost/math/special_functions/bessel.hpp>  // For Bessel Function
#include <boost/math/special_functions/chebyshev.hpp>

#include "chebyshev_gauss_quadrature.hpp"
#include "gauleg.hpp"
#endif  // SOLUTION

/* @brief Compute matrix M using analytic expression.
 *
 * \param[in] N Discretization parameter indicating number of basis functions.
 */
/* SAM_LISTING_BEGIN_0 */
MatrixXd computeM(unsigned int N) {
  Eigen::MatrixXd M(N, N);
  M.setZero();

#if SOLUTION
  // Matrix assembly
  for (unsigned int i = 0; i < N; ++i) M(i, i) = M_PI / 4. / (i + 1);

#else   // TEMPLATE
  // TODO: Compute M
#endif  // TEMPLATE
  return M;
}
/* SAM_LISTING_END_0 */

/* @brief Compute right hand side using Chebyshev Gauss Quadrature
 *
 * \tparam FUNC Type for the Dirichlet B.C. supporting : g(const double&)
 * \param[in] g Dirichlet boundary condition; Signature : double(const double&)
 * \param[in] N Discretization parameter indicating number of basis functions.
 */
/* SAM_LISTING_BEGIN_1 */
template <typename FUNC>
VectorXd computeG(const FUNC &g, int N) {
  // Initialize right hand side vector
  VectorXd RHS(N);
  RHS.setZero();
#if SOLUTION
  // Get Chebyshev Gauss Quadrature points and weight
  unsigned int order = 2 * N;
  Eigen::RowVectorXd points;
  double weight;
  std::tie(points, weight) = ChebyshevGaussQuad(order);
  // Fill vector entries
  for (unsigned int i = 0; i < N; ++i) {
    // Evaluate integral by using quadrature
    for (unsigned int qp = 0; qp < order; ++qp) {
      double x = points(qp);
      RHS(i) += weight * g(x) * boost::math::chebyshev_t(i + 1, x);
    }  // end iteration over quadrature points
  }

#else   // TEMPLATE
  // TODO: Compute RHS
#endif  // TEMPLATE
  return RHS;
}
/* SAM_LISTING_END_1 */

/* @brief Build and solve boundary integral equation V rho = g
 *
 * \tparam FUNC Type for the Dirichlet B.C. supporting : g(const double&)
 * \param[in] g Dirichlet boundary condition; Signature : double(const double&)
 * \param[in] N Discretization parameter indicating number of basis functions.
 */
/* SAM_LISTING_BEGIN_2 */
template <typename FUNC>
VectorXd solveBIE(const FUNC &g, int N) {
#if SOLUTION
  // Assemble Galerkin Matrix M
  Eigen::MatrixXd M = computeM(N);
  // Build RHS
  VectorXd RHS = computeG(g, N);
  // Use direct solver
  MatrixXd LHS = M.eval();
  VectorXd sol = LHS.lu().solve(RHS);
  return sol;
#else   // TEMPLATE
  // TODO: Build BIE system and solve it
#endif  // TEMPLATE
}

#if SOLUTION
/* SAM_LISTING_END_2 */

/* @brief Reconstruct function UN from its coefficients and evaluate it at t.
 *
 * \param[in] coeffs coefficients of UN
 * \param[in] t evaluation point [-1,1]
 */
/* SAM_LISTING_BEGIN_3a */
double reconstructRho(const VectorXd &coeffs, double t) {
  assert(t >= -1 && t <= 1);  // Asserting evaluation is within the domain
  int N = coeffs.rows();
  double rho_N = 0.;
  for (unsigned int i = 0; i < N; ++i)
    // Coefficients start from $T_1(x)$
    rho_N += coeffs(i) * boost::math::chebyshev_t(i + 1, t);
  return rho_N / std::sqrt(1 - t * t);
}
#endif

/* SAM_LISTING_END_3a */

/* @brief Compute L2 norm of UN from its coefficients using Gauss Legendre
 *        Quadrature rule
 *
 * \param[in] coeffs coefficients of UN
 */
/* SAM_LISTING_BEGIN_3b */
double L2norm(const VectorXd &coeffs) {
  double norm = 0.;
#if SOLUTION
  int N = coeffs.rows();
  // Get quadrature points and weight for Gauss Legendre Quadrature
  unsigned int order = 2 * N;  // Quadrature order
  Eigen::RowVectorXd weights, points;
  std::tie(points, weights) = gauleg(-1, 1, order);
  // Iterating over quadrature points
  for (int qp = 0; qp < order; qp++) {
    auto z = points(qp);
    // evaluating the function
    double rho = reconstructRho(coeffs, z);
    norm += weights(qp) * rho * rho;
  }
#else   // TEMPLATE
  // TODO: implement the computation of the L2norm
#endif  // TEMPLATE
  return std::sqrt(norm);
}
/* SAM_LISTING_END_3b */

/* SAM_LISTING_BEGIN_4 */
int main() {
  std::cout << "Test for source term g1(x) = sin(2*Pi*x)" << std::endl;
  std::cout << "N" << std::setw(15) << "L2error" << std::endl;
  std::cout << "############################" << std::endl;

  // std::function object for source term g1 using a lambda expression
  std::function<double(const double &)> g1 = [](const double &x) {
    return sin(2 * M_PI * x);
  };
#if SOLUTION
  // Analytically evaluating the coefficients for the truncated infinite series
  unsigned int trunc_num = 100;  // Number of terms in truncated series
  Eigen::VectorXd exact_coeffs(trunc_num);
  // Evaluating the exact coefficients
  for (unsigned int i = 0; i < trunc_num; ++i)
    // Coefficients starting from $T_1(x)$
    exact_coeffs(i) = 4 * (i + 1) *
                      boost::math::cyl_bessel_j((i + 1), 2 * M_PI) *
                      sin((i + 1) * M_PI / 2.);

  // Varying the discretization parameter
  for (unsigned int N = 2; N < 50; ++N) {
    // Solving the BIE to get the coefficients
    Eigen::MatrixXd M = computeM(N);
    Eigen::VectorXd coeffs = solveBIE(g1, N);
    // Extending the evaluated coefficients to truncted series' length by zeros
    Eigen::VectorXd extended_coeffs(trunc_num);
    Eigen::VectorXd zeros(trunc_num - N);
    zeros.setZero();
    extended_coeffs << coeffs, zeros;  // Coefficients extended with zeros
    // Error coefficients
    Eigen::VectorXd error_coeffs = extended_coeffs - exact_coeffs;
    // Evaluating error
    double l2error = L2norm(error_coeffs);

    std::cout << N << std::setw(15) << l2error << std::endl;
  }
#else   // TEMPLATE
  // TODO: implement the convergence test
#endif  // TEMPLATE
  return 0;
}
/* SAM_LISTING_END_4 */
