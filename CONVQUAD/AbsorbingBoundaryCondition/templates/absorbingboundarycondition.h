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



<<<<<<< HEAD:CONVQUAD/AbsorbingBoundaryCondition/all/absorbingboundarycondition.h
=======
/* @brief Compute the BDF-2 convolution weights for Laplace transform F
 * \param F Template function for the Laplace transform
 * \param n Number of convolution weights, plus 1
 * \param p Order of quadrature rule
 * \param r Radius of circumference of integration (default = 1.0E-7)
 * \\return BDF-2 convolution weights
 */
/* SAM_LISTING_BEGIN_0 */
template <typename FFUNC>
VectorXd conv_wght_bdf2(const FFUNC& F, size_t n, int p, double r = 1.0E-7) {
  // TODO: Compute the convolution weights for Laplace transform F
}
/* SAM_LISTING_END_0 */

>>>>>>> origin:CONVQUAD/AbsorbingBoundaryCondition/AbsorbingBoundaryCondition.cpp

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

<<<<<<< HEAD:CONVQUAD/AbsorbingBoundaryCondition/all/absorbingboundarycondition.h
=======
int main() {
  /* SAM_LISTING_BEGIN_3 */
  // TODO: Tabulate the H1-error of the Galerkin discretization + convolution quadrature
  /* SAM_LISTING_END_3 */
}
>>>>>>> origin:CONVQUAD/AbsorbingBoundaryCondition/AbsorbingBoundaryCondition.cpp
