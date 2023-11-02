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
#if SOLUTION
  // Setting $r = EPS^{\frac{1}{2M+2}}$, as discussed in Rem. 3.4.3.20
  double r = std::pow(10, -16.0 / (2 * M + 2));
  // Initialize vector for the evaluations in the Laplace domain
  Eigen::VectorXcd f = Eigen::VectorXd::Zero(M + 1);
  // Imaginary unit stored in imag
  std::complex<double> imag = std::complex<double>(0, 1);
  // Variable for the frequencies along the contour
  std::complex<double> s_k;
  for (int k = 0; k < M + 1; k++) {
    // Integrate on circumference centered in (0,0) with radius r
    s_k = delta(r * std::exp(2 * M_PI * imag * ((double)k) / (double)(M + 1))) /
          tau;
    f[k] = F(s_k);
  }
  // Integrate with trapezoidal rule at M+1 points
  Eigen::FFT<double> fft;
  w = fft.fwd(f) / (M + 1);
  for (int k = 0; k < M + 1; k++) {
    // Rescale by the radius of the circle, which arise from the $z^l$ -factor
    // in the integrand
    w[k] = w[k] / std::pow(r, k);
  }
#else
  // **********************************************************************
  // Your Solution here
  // **********************************************************************/
#endif
  return w.real();
}
/* SAM_LISTING_END_0 */

/* @brief Build the sparse symmetric tri-diagonal matrix
 * @param N Number of discretization intervals in space
 * @return SparseMatrix A
 */
/* SAM_LISTING_BEGIN_1 */
SparseMatrix<double> compute_matA(size_t N) {
  double h = 1. / N;
  VectorXd grid = VectorXd::LinSpaced(N + 1, 0., 1.);

  SparseMatrix<double> A(N + 1, N + 1);
  A.reserve(3 * N + 1);  // 3(N+1) - 2
#if SOLUTION
  // Inserting endpoints.
  A.insert(0, 0) =
      1. / h + 1. / (h * h * pow(M_PI, 3)) *
                   ((pow(M_PI * h, 2) - 2.) + 2. * cos(M_PI * h));  // A(0,0)
  A.insert(N, N) = 1. / h + 1. / (h * h * pow(M_PI, 3)) *
                                ((pow(M_PI * h, 2) - 2.) -
                                 2. * cos(M_PI * (1 - h)));  // A(N,N)

  for (int i = 1; i <= N; ++i) {
    if (i < N) {
      // Inserting diagonal entries
      A.insert(i, i) =
          2. / h + 2. / (h * h * pow(M_PI, 3)) *
                       (2 * M_PI * h * sin(M_PI * grid(i)) +
                        cos(M_PI * grid(i + 1)) - cos(M_PI * grid(i - 1)));
    }
    // Inserting sub- and super-diagonal entries
    // A(i-1,i) = A(i,i-1)
    A.insert(i - 1, i) = A.insert(i, i - 1) =
        -1. / h +
        1. / (h * h * pow(M_PI, 3)) *
            (2. * (cos(M_PI * grid(i - 1)) - cos(M_PI * grid(i))) -
             M_PI * h * (sin(M_PI * grid(i - 1)) + sin(M_PI * grid(i))));
  }
#else
// **********************************************************************
// Your Solution here
// **********************************************************************/
#endif
  return A;
}
/* SAM_LISTING_END_1 */

/* @brief Find the unknown function u at final time t = 1 in the evolution
 * problem using Galerkin discretization and convolution quadrature (BDF-2)
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
#if SOLUTION
  // Aw <- A + lowToeplitz(w)*B; B(N,N) = 1, else B(i,j) = 0
  SparseMatrix<complex<double> > Aw = A.cast<complex<double> >();
  Aw.coeffRef(N, N) += w(0);
  SparseLU<SparseMatrix<complex<double> > > solver;
  solver.compute(Aw);
  // For visualization
  ofstream f_ref;
  if (M == 4096 && N == 4096) f_ref.open(CURRENT_BINARY_DIR "/u_ref.txt");
  // run solver from t=0 -> t=T
  for (int i = 1; i <= M; ++i) {
    // rhs\_cq: from the convolution weights $l=0,1,\ldots,n-1$
    complex<double> rhs_cq = 0.;
    for (int l = 0; l < i; ++l) {
      rhs_cq += w(i - l) * u(N, l);
    }

    // function $\phi$
    VectorXcd phi = VectorXcd::Zero(N + 1);
    phi(0) = -(complex<double>)g(i * tau);

    VectorXcd rhs = phi;
    rhs(N) -= rhs_cq;              // rhs
    u.col(i) = solver.solve(rhs);  // solution at $t = t_n$
  }
  // For visualization
  if (f_ref.is_open()) {
    f_ref << u.real();
    f_ref.close();
  }
#else
// **********************************************************************
// Your Solution here
// **********************************************************************/
#endif
  return u.col(M).real();
}
/* SAM_LISTING_END_2 */

}  // namespace AbsorbingBoundaryCondition

#endif
