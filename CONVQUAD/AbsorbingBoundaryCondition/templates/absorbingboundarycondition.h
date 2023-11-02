#ifndef ABC_H_
#define ABC_H_


#ifndef ABC_H_
#define ABC_H_


#include <Eigen/Dense>
#include <Eigen/Sparse>
#include <Eigen/SparseLU>
#include <cmath>
#include <fstream>
#include <iomanip>
#include <iostream>
#include <unsupported/Eigen/FFT>

using namespace Eigen;
using namespace std;

namespace AbsorbingBoundaryCondition {
/* @brief Compute the convolution quadrature weights for Laplace transform F
 * \param F     Template function for the Laplace transform
 * \param delta Template function determining the multistep method
 * \param tau   Time step size
 * \param M     Number of time steps
 * \\return convolution quadrature weights
 */
/* SAM_LISTING_BEGIN_0 */
template <typename FFUNC, typename DFUNC>
VectorXd cqweights_by_dft(const FFUNC& F, const DFUNC& delta, double tau,
                          size_t M) {
  Eigen::VectorXcd w = Eigen::VectorXd::Zero(M + 1);
  // **********************************************************************
  // Your Solution here
  // **********************************************************************/
  return w.real();
}
/* SAM_LISTING_END_0 */

/* @brief Build the sparse symmetric tri-diagonal matrix
 * \param N Number of discretization intervals in space
 * \\return SparseMatrix A
 */
/* SAM_LISTING_BEGIN_1 */
SparseMatrix<double> compute_matA(size_t N) {
  double h = 1. / N;
  VectorXd grid = VectorXd::LinSpaced(N + 1, 0., 1.);

  SparseMatrix<double> A(N + 1, N + 1);
  A.reserve(3 * N + 1);  // 3(N+1) - 2
// **********************************************************************
// Your Solution here
// **********************************************************************/
  return A;
}
/* SAM_LISTING_END_1 */

/* @brief Find the unknown function u at final time t = 1 in the evolution problem
 * using Galerkin discretization and convolution quadrature (BDF-2)
 * \param g Template function for the right-hand side
 * \param M Number of discretization intervals in time
 * \param N Number of discretization intervals in space
 * \param T Final Time
 * \\return Values of u at final time t = T
 */
/* SAM_LISTING_BEGIN_2 */
template <typename FUNC>
VectorXd solve_IBVP(const FUNC& g, size_t M, size_t N, double T) {
  MatrixXcd u = MatrixXcd::Zero(N + 1, M + 1);
  auto F = [](complex<double> s) { return log(s) / (s * s + 1.); };
  // matrix A
  SparseMatrix<double> A(N + 1, N + 1);
  A = compute_matA(N);
  // convolution weights
  auto delta = [](std::complex<double> z) {
    return 1.0 / 2.0 * z * z - 2.0 * z + 3.0 / 2.0;
  };
  double tau = T / M;
  Eigen::VectorXd w = cqweights_by_dft(F, delta, tau, M);
// **********************************************************************
// Your Solution here
// **********************************************************************/
  return u.col(M).real();
}
/* SAM_LISTING_END_2 */
}  // namespace AbsorbingBoundaryCondition
#endif  // Macro ABC_H_
