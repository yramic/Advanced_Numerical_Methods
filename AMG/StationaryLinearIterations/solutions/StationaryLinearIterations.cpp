#include <Eigen/Dense>
#include <Eigen/Sparse>
#include <cmath>
#include <ctime>
#include <iomanip>
#include <iostream>
#include <unsupported/Eigen/FFT>

using namespace Eigen;
using namespace std;

/* @brief Poisson matrix for a uniform triangulation on a unit square domain
 * \param N matrix size
 * \\return sparse Poisson matrix
 */
/* SAM_LISTING_BEGIN_0 */
SparseMatrix<double> poissonMatrix(const int n) {
  int N = (n - 1) * (n - 1);

  // define vector of triplets and reserve memory
  typedef Triplet<double> triplet;
  vector<triplet> entries;
  int nnz = 0;
  if (n == 2) {
    nnz = 1;
  } else if (n > 2) {
    nnz = 5 * n * n - 14 * n + 9;
  }
  entries.reserve(nnz);

  // set the vector of triplets
  for (int block_id = 0; block_id < (n - 1); block_id++) {
    int start_id = block_id * (n - 1);
    int end_id = (block_id + 1) * (n - 1) - 1;

    // tri-diagonal matrix T
    for (int i = start_id; i <= end_id; i++) {
      entries.push_back(triplet(i, i, 4));
      if (i > start_id) {  // T(i-1,i) = T(i,i-1)
        entries.push_back(triplet(i, i - 1, -1));
        entries.push_back(triplet(i - 1, i, -1));
      }
    }

    // identity blocks
    if (block_id > 0) {
      for (int i = start_id; i <= end_id; i++) {
        entries.push_back(triplet(i, i - (n - 1), -1));
        entries.push_back(triplet(i - (n - 1), i, -1));
      }
    }
  }

  // create the sparse matrix
  SparseMatrix<double> A(N, N);
  A.setFromTriplets(entries.begin(), entries.end());

  return A;
}
/* SAM_LISTING_END_0 */

/* @brief The Gauss-Seidel iterative method
 * to solve sparse linear system of equations
 * \param A Sparse matrix, $\IR^{N \times N}$
 * \param phi right-hand-side vector, $\IR^{N}$
 * \param mu initial guess and output solution, $\IR^{N \times N}$
 * \param TOL iteration error tolerance for termination criteria
 */
/* SAM_LISTING_BEGIN_1 */
void gaussSeidel(const SparseMatrix<double> &A, const VectorXd &phi,
                 VectorXd &mu, double TOL = 1.0E-06) {
  VectorXd delta(A.rows());
  do {
    delta = A.triangularView<Lower>().solve(phi - A * mu);
    mu += delta;
  } while (delta.norm() > TOL * mu.norm());
}
/* SAM_LISTING_END_1 */

/* SAM_LISTING_BEGIN_2 */
double comp_lmax_gaussSeidel(const SparseMatrix<double> &X,
                             double TOL = 1.0E-03) {
  int N = X.rows();

  VectorXd v = MatrixXd::Random(N, 1);
  double lambda_new = 0;
  double lambda_old = 1;

  int itr = 0;
  do {
    lambda_old = lambda_new;
    v /= v.norm();
    v -= X.triangularView<Lower>().solve(X *
                                         v);  // E = I - M*A; M = inv(tril(A))
    lambda_new = v.norm();
  } while (fabs(lambda_new - lambda_old) > TOL * lambda_new);

  return lambda_new;
}
/* SAM_LISTING_END_2 */

/* @brief Determines the asymptotic convergence rate of the Gauss-Seidel
 * iterative method using power iteration,
 * for the matrix $\mathbf{X} = \mathbf{A} + c \mathbf{I}_{N}$.
 * Here $\mathbf{A}$ is the Poisson matrix on a unit square domain.
 * \param N matrix size
 * \param c linear combination coefficient
 * \\return asymptotic convergence rate
 */
/* SAM_LISTING_BEGIN_3 */
double gaussSeidelRate(const int n, double c, double TOL = 1.0E-03) {
  int N = (n - 1) * (n - 1);

  SparseMatrix<double> X = poissonMatrix(n);
  SparseMatrix<double> I(N, N);
  I.setIdentity();
  X += c * I;
  double lambda_max = comp_lmax_gaussSeidel(X, TOL);

  return lambda_max;
}
/* SAM_LISTING_END_3 */

double comp_lmax_jacobi(const SparseMatrix<double> &X, double w,
                        double TOL = 1.0E-03) {
  int N = X.rows();

  VectorXd v = MatrixXd::Random(N, 1);
  double lambda_new = 0;
  double lambda_old = 1;

  VectorXd D = X.diagonal();
  int itr = 0;
  do {
    lambda_old = lambda_new;
    v /= v.norm();
    VectorXd X_mult_v = X * v;
    for (int i = 0; i < N; i++) {  // E = I - M*A; M = w*inv(diag(A))
      v(i) -= w * X_mult_v(i) / D(i);
    }
    lambda_new = v.norm();
  } while (fabs(lambda_new - lambda_old) > TOL * lambda_new);

  return lambda_new;
}

/* @brief Determines the asymptotic convergence rate of the Jacobi
 * iterative method using power iteration,
 * for the Poisson matrix $\mathbf{A}$ on a unit square domain.
 * \param N matrix size
 * \param w relaxation parameter
 * \\return asymptotic convergence rate
 */
double jacobiRate(const int n, double w, double TOL = 1.0E-03) {
  int N = (n - 1) * (n - 1);
  SparseMatrix<double> A = poissonMatrix(n);
  double lambda_max = comp_lmax_jacobi(A, w, TOL);
  return lambda_max;
}

int main() {
  /* SAM_LISTING_BEGIN_4 */
  {
    cout << "\n\nAccuracy test of the Gauss-Seidel method" << endl;
    for (int l = 2; l <= 3; l++) {
      int n = pow(2, l);
      int N = (n - 1) * (n - 1);
      double TOL = 1.0E-08;

      SparseMatrix<double> A = poissonMatrix(n);

      VectorXd mu_exact = MatrixXd::Random(N, 1);
      VectorXd phi = A * mu_exact;

      VectorXd mu = MatrixXd::Identity(N, 1);
      gaussSeidel(A, phi, mu, TOL);
      cout << n << "\t" << (mu_exact - mu).norm() / mu_exact.norm() << endl;
    }
  }
  /* SAM_LISTING_END_4 */

  /* SAM_LISTING_BEGIN_5 */
  {
    cout
        << "\n\nCompute asymptotic convergence rates of the Gauss-Seidel method"
        << endl;

    VectorXd c(3);
    c << 0, 1, 10;
    double TOL = 1.0E-03;

    for (int k = 0; k < c.size(); k++) {
      cout << "\nFor c = " << c(k) << ":" << endl;
      for (int l = 2; l <= 10; l++) {
        int n = pow(2, l);
        double lambda_max = gaussSeidelRate(n, c(k), TOL);
        cout << n << "\t" << lambda_max << endl;
      }
    }
  }
  /* SAM_LISTING_END_5 */

  {
    cout << "\n\nCompute asymptotic convergence rates of the Jacobi method"
         << endl;

    VectorXd w(10);
    for (int i = 0; i < 10; i++) w(i) = (i + 1) * 0.1;
    double TOL = 1.0E-06;

    for (int k = 0; k < w.size(); k++) {
      cout << "\nFor w = " << w(k) << ":" << endl;
      for (int l = 2; l <= 6; l++) {
        int n = pow(2, l);
        double lambda_max = jacobiRate(n, w(k), TOL);
        cout << n << "\t" << lambda_max << endl;
      }
    }
  }
}
