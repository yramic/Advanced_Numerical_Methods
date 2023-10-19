#include <Eigen/Dense>
#include <Eigen/Sparse>
#include <Eigen/SparseLU>
#include <cmath>
#include <iomanip>
#include <iostream>

using namespace Eigen;
using namespace std;

#if SOLUTION
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
#endif

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
#if SOLUTION
  double tau = 1. / n;
  VectorXcd w = VectorXd::Zero(n + 1);

  // integrate on circumference centered in (0,0) with radius r

  // define integrand after transformation: $s = r \exp(2 \pi \varphi)$, where parameter is $\varphi$
  auto integrand = [](const FFUNC& F, double tau, complex<double> s, int l) {
    return s / sqrt(1. + 2. * tau * s) *
           (1. / pow(2. + sqrt(1. + 2. * tau * s), l + 1) -
            1. / pow(2. - sqrt(1. + 2. * tau * s), l + 1)) *
           F(s);
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
#else   // TEMPLATE
  // TODO: Compute the convolution weights for Laplace transform F
#endif  // TEMPLATE
}
/* SAM_LISTING_END_0 */

#if SOLUTION
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
#endif

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
#if SOLUTION
  auto F = [](complex<double> s) { return log(s) / (s * s + 1.); };

  // matrix A
  SparseMatrix<double> A(N + 1, N + 1);
  A = compute_matA(N);

  // convolution weights
  VectorXd w = conv_wght_bdf2(F, M, p);

  // Aw <- A + lowToeplitz(w)*B; B(N,N) = 1, else B(i,j) = 0
  SparseMatrix<complex<double> > Aw = A.cast<complex<double> >();
  Aw.coeffRef(N, N) += w(0);
  SparseLU<SparseMatrix<complex<double> > > solver;
  solver.compute(Aw);

  // run solver from t=0 -> t=1
  double tau = 1. / M;
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
#else   // TEMPLATE
  // TODO: Find the unknown function u at final time t = 1 in the evolution problem
#endif  // TEMPLATE
}
/* SAM_LISTING_END_2 */

int main() {
  /* SAM_LISTING_BEGIN_3 */
#if SOLUTION
  auto g = [](double t) { return sin(M_PI * t); };

  // compute reference solution
  int M_ref = 4096;
  int N_ref = 4096;
  double h_ref = 1. / N_ref;
  VectorXd u_ref = solve_IBVP(g, M_ref, N_ref, 20);

  // compute H1 norm of reference solution
  double norm_u_ref = 0.;
  for (int i = 1; i <= N_ref; ++i) {
    norm_u_ref += pow((u_ref(i) - u_ref(i - 1)), 2);
  }
  norm_u_ref = sqrt(norm_u_ref / h_ref);

  cout << "\nConvergence wrt spatial discretisation" << endl;
  cout << "N\tRelative H1-error" << endl;
  for (int N = 16; N <= N_ref / 4; N *= 2) {
    double h = 1. / N;
    VectorXd u_tmp = solve_IBVP(g, M_ref, N, 20);
    double error = 0.;
    double ratio = N_ref / N;
    for (int i = 1; i <= N_ref; ++i) {
      int j = ceil(i / ratio);
      error += pow(
          (u_ref(i) - u_ref(i - 1)) / h_ref - (u_tmp(j) - u_tmp(j - 1)) / h, 2);
    }
    error = sqrt(error * h_ref) / norm_u_ref;
    cout << N << "\t" << scientific << setprecision(10) << error << endl;
  }

  cout << "\nConvergence wrt time discretisation" << endl;
  cout << "M\tRelative H1-error" << endl;
  for (int M = 16; M <= M_ref / 4; M *= 2) {
    VectorXd u_tmp = solve_IBVP(g, M, N_ref, 20);
    double error = 0.;
    for (int i = 1; i <= N_ref; ++i) {
      error += pow((u_ref(i) - u_ref(i - 1)) - (u_tmp(i) - u_tmp(i - 1)), 2);
    }
    error = sqrt(error / h_ref) / norm_u_ref;
    cout << M << "\t" << scientific << setprecision(10) << error << endl;
  }
#else   // TEMPLATE
  // TODO: Tabulate the H1-error of the Galerkin discretization + convolution quadrature
#endif  // TEMPLATE
        /* SAM_LISTING_END_3 */
}
