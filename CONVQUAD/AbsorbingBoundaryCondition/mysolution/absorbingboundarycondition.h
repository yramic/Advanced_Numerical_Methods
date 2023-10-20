#include <Eigen/Dense>
#include <Eigen/Sparse>
#include <Eigen/SparseLU>
#include <cmath>
#include <iomanip>
#include <iostream>
#include <unsupported/Eigen/FFT>
using namespace Eigen;
using namespace std;

/* @brief Compute the convolution quadrature weights for Laplace transform F
 * \param F     Template function for the Laplace transform
 * \param delta Template function determining the multistep method
 * \param tau   Time step size  
 * \param N     Number of time steps
 * \\return convolution quadrature weights
 */
/* SAM_LISTING_BEGIN_0 */
template <typename FFUNC, typename DFUNC>
VectorXd cqweights_by_dft(const FFUNC& F, const DFUNC& delta, double tau,
                          size_t N) {
  Eigen::VectorXcd w = Eigen::VectorXd::Zero(N + 1);
  // Setting $r = EPS^{\frac{1}{2N+2}}$, as discussed in Rem. 3.4.3.20
  double r = std::pow(10, -16.0 / (2 * N + 2));
  // TODO: Compute the convolution weights for Laplace transform F

  return w.real();
}
/* SAM_LISTING_END_0 */

/* @brief Find the unknown function u at final time t = 1 in the evolution problem
 * using Galerkin discretization and convolution quadrature (BDF-2)
 * \param g Template function for the right-hand side
 * \param M Number of discretization intervals in time
 * \param N Number of discretization intervals in space
 * \param p Order of quadrature rule
 * \\return Values of u at final time t = 1
 */
/* SAM_LISTING_BEGIN_2 */
template <typename FUNC>
VectorXd solve_IBVP(const FUNC& g, size_t M, size_t N, int p) {
  // TODO: Find the unknown function u at final time t = 1 in the evolution problem
}
/* SAM_LISTING_END_2 */
