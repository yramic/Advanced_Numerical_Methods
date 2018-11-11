#include <Eigen/Dense>
#include <cassert>
#include <cmath>
#include <fstream>
#include <iostream>
#include <istream>

using namespace Eigen;


/* @brief Compute matrix M using analytic expression.
 *
 * \param[in] N Discretization parameter indicating number of basis functions.
 */
/* SAM_LISTING_BEGIN_0 */
MatrixXd computeM(unsigned int N) {
  Eigen::MatrixXd M(N, N);
  M.setZero();

  // TODO: Compute M
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
  // TODO: Compute RHS
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
  // TODO: Build BIE system and solve it
}

 
/* SAM_LISTING_END_3a */

/* @brief Compute L2 norm of UN from its coefficients using Gauss Legendre
 *        Quadrature rule
 *
 * \param[in] coeffs coefficients of UN
 */
/* SAM_LISTING_BEGIN_3b */
double L2norm(const VectorXd &coeffs) {
  // TODO: implement the computation of the L2norm
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
  // TODO: implement the convergence test
  return 0;
}
/* SAM_LISTING_END_4 */
