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
  // Initialize vector for the evaluations in the Laplace domain
  Eigen::VectorXcd f = Eigen::VectorXd::Zero(N + 1);
  // Imaginary unit stored in imag
  std::complex<double> imag = std::complex<double>(0, 1);
  // Variable for the frequencies along the contour
  std::complex<double> s_k;
  for (int k = 0; k < N + 1; k++) {
    // Integrate on circumference centered in (0,0) with radius r
    s_k = delta(r * std::exp(2 * M_PI * imag * ((double)k) / (double)(N + 1))) /
          tau;
    f[k] = F(s_k);
  }
  // Integrate with trapezoidal rule at N+1 points
  Eigen::FFT<double> fft;
  w = fft.fwd(f) / (N + 1);
  for (int k = 0; k < N + 1; k++) {
    // Rescale by the radius of the circle, which arise from the z^l -factor in the integrand
    w[k] = w[k] / std::pow(r, k);
  }

  return w.real();
}
/* SAM_LISTING_END_0 */

/* @brief Compute the IE convolution weights for Laplace transform F
 * \param F Template function for the Laplace transform
 * \param n Number of convolution weights, plus 1
 * \param p Order of quadrature rule
 * \param r Radius of circumference of integration (default = 1.0E-7)
 * \\return IE convolution weights
 */
template <typename FFUNC>
VectorXd conv_wght_ieu(const FFUNC& F, size_t n, int p, double r = 1.0E-7) {
  double tau = 1. / n;
  VectorXcd w = VectorXd::Zero(n + 1);

  // integrate on circumference centered in (0,0) with radius r

  // define integrand after transformation: $s = r \exp(2 \pi \varphi)$, where parameter is $\varphi$
  auto integrand = [](const FFUNC& F, double tau, complex<double> s, int l) {
    return s / pow(1. - tau * s, l + 1) * F(s);
  };

  for (int i = 0; i < p; ++i) {
    complex<double> s =
        r * complex<double>(cos(2. * M_PI * i / p), sin(2. * M_PI * i / p));

    // integrate with trapezoidal rule
    double coeff = 2.0;
    if (i == 0 || i == p - 1) coeff = 1.0;

    for (int l = 0; l < n + 1; ++l) {
      w(l) += coeff * integrand(F, tau, s, l);
    }
  }
  w *= -0.5 * tau / p;  // $\tau$ times quadrature weight $1./(2*p)$

  return w.real();
}

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

  A.insert(0, 0) =
      1. / h + 1. / (h * h * pow(M_PI, 3)) *
                   ((pow(M_PI * h, 2) - 2.) + 2. * cos(M_PI * h));  //A(0,0)
  A.insert(N, N) = 1. / h + 1. / (h * h * pow(M_PI, 3)) *
                                ((pow(M_PI * h, 2) - 2.) -
                                 2. * cos(M_PI * (1 - h)));  //A(N,N)

  for (int i = 1; i <= N; ++i) {
    if (i < N) {  // A(i,i)
      A.insert(i, i) =
          2. / h + 2. / (h * h * pow(M_PI, 3)) *
                       (2 * M_PI * h * sin(M_PI * grid(i)) +
                        cos(M_PI * grid(i + 1)) - cos(M_PI * grid(i - 1)));
    }
    // A(i-1,i) = A(i,i-1)
    A.insert(i - 1, i) = A.insert(i, i - 1) =
        -1. / h +
        1. / (h * h * pow(M_PI, 3)) *
            (2. * (cos(M_PI * grid(i - 1)) - cos(M_PI * grid(i))) -
             M_PI * h * (sin(M_PI * grid(i - 1)) + sin(M_PI * grid(i))));
  }

  return A;
}
/* SAM_LISTING_END_1 */

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
  //auto F = [](complex<double> s) { return s / (std::exp(-s) + 1.); };
  auto F = [](complex<double> s) { return log(s) / (s * s + 1.); };

  // matrix A
  SparseMatrix<double> A(N + 1, N + 1);
  A = compute_matA(N);

  // convolution weights
  auto delta = [](std::complex<double> z) {
    return 1.0 / 2.0 * z * z - 2.0 * z + 3.0 / 2.0;
  };
  // Assume the final time to be $t=1$
  double tau = 1.0 / M;
  Eigen::VectorXd w = cqweights_by_dft(F, delta, tau, M);
  std::cout << (wv - w).squaredNorm() << std::endl;
  // Aw <- A + lowToeplitz(w)*B; B(N,N) = 1, else B(i,j) = 0
  SparseMatrix<complex<double> > Aw = A.cast<complex<double> >();
  Aw.coeffRef(N, N) += w(0);
  SparseLU<SparseMatrix<complex<double> > > solver;
  solver.compute(Aw);
  // run solver from t=0 -> t=1
  //double tau = 1. / M;
  MatrixXcd u = MatrixXcd::Zero(N + 1, M + 1);
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
    rhs(N) -= rhs_cq;              //rhs
    u.col(i) = solver.solve(rhs);  // solution at $t = t_n$
  }

  return u.col(M).real();
}
/* SAM_LISTING_END_2 */
