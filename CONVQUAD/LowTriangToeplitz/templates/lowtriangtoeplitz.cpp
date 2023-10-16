#include <Eigen/Dense>
#include <chrono>
#include <cmath>
#include <iostream>
#include <unsupported/Eigen/FFT>

namespace LowTriangToeplitz {

/* @brief Generate a Toeplitz matrix
 * \param c Vector of entries of first column of the Toeplitz matrix
 * \param r Vector of entries of first row of the Toeplitz matrix
 * \\return Toeplitz matrix
 */
Eigen::MatrixXcd toeplitz(const Eigen::VectorXcd& c,
                          const Eigen::VectorXcd& r) {
  if (c(0) != r(0)) {
    std::cerr << "First entries of c and r are different!" << std::endl
              << "We assign the first entry of c to the diagonal" << std::endl;
  }

  // Initialization
  const std::size_t m = c.size();
  const std::size_t n = r.size();
  Eigen::MatrixXcd T(m, n);

  for (int i = 0; i < n; ++i) {
    T.col(i).tail(m - i) = c.head(m - i);
  }
  for (int i = 0; i < m; ++i) {
    T.row(i).tail(n - i - 1) = r.segment(1, n - i - 1);
  }  // Do not reassign the diagonal!

  return T;
}

/* @brief Multiply a Circulant matrix with a vector, using FFT
 * \param u Generating vector for the Circulant matrix
 * \param x Vector
 * \\return Circulant(u)*x
 */
Eigen::VectorXcd pconvfft(const Eigen::VectorXcd& u,
                          const Eigen::VectorXcd& x) {
  Eigen::FFT<double> fft;
  const Eigen::VectorXcd tmp = (fft.fwd(u)).cwiseProduct(fft.fwd(x));
  return fft.inv(tmp);
}

/* @brief Multiply two lower triangular Toeplitz matrices
 * \param f Vector of entries of first lower triangular Toeplitz matrix
 * \param g Vector of entries of second lower triangular Toeplitz matrix
 * \\return Vector of entries of output lower triangular Toeplitz matrix
 */
/* SAM_LISTING_BEGIN_0 */
Eigen::VectorXcd ltpMult(const Eigen::VectorXcd& f, const Eigen::VectorXcd& g) {
  assert(f.size() == g.size() && "f and g vectors must have the same length!");
  const std::size_t n = f.size();
  Eigen::VectorXcd res(n);
  // **********************************************************************
  // Your Solution here
  // **********************************************************************
  return res;
}
/* SAM_LISTING_END_0 */

/* SAM_LISTING_BEGIN_1 */
std::tuple<double, double, double> runtimes_ltpMult(unsigned int N) {
  // Runtime of matrix-matrix, matrix-vector and vector-vector multiplication
  // in seconds
  double s_dense, s_mv, s_ltp;

  // **********************************************************************
  // Code to be supplemented
  // **********************************************************************
  return {s_dense, s_mv, s_ltp};
}
/* SAM_LISTING_END_1 */

/* @brief Multiply a Toeplitz matrix with a vector, uses pconvfft
 * \param c Vector of entries of first column of the Toeplitz matrix
 * \param r Vector of entries of first row of the Toeplitz matrix
 * \param x Vector
 * \\return toeplitz(c,r)*x
 */
/* SAM_LISTING_BEGIN_2 */
Eigen::VectorXcd toepMatVecMult(const Eigen::VectorXcd& c,
                                const Eigen::VectorXcd& r,
                                const Eigen::VectorXcd& x) {
  assert(c.size() == x.size() && r.size() == x.size() &&
         "c, r, x have different lengths!");

  const std::size_t n = c.size();
  Eigen::VectorXcd y(2 * n);
  // **********************************************************************
  // Code to be supplemented
  // **********************************************************************
  return y.head(n);
}
/* SAM_LISTING_END_2 */

/* @brief Solve a linear problem involving a lower triangular Toeplitz matrix
 * \param f Vector of entries of lower triangular Toeplitz matrix
 * \param y Right-hand side of linear problem
 * \\return Solution of linear problem
 */
/* SAM_LISTING_BEGIN_3 */
Eigen::VectorXcd ltpSolve(const Eigen::VectorXcd& f,
                          const Eigen::VectorXcd& y) {
  assert(f.size() == y.size() && "f and y vectors must have the same length!");
  assert(abs(f(0)) > 1e-10 &&
         "Lower triangular Toeplitz matrix must be invertible!");
  assert(log2(f.size()) == floor(log2(f.size())) &&
         "Size of f must be a power of 2!");

  const std::size_t n = f.size();
  Eigen::VectorXcd u(n);
  // **********************************************************************
  // Code to be supplemented
  // **********************************************************************
  return u;
}
/* SAM_LISTING_END_3 */

/* SAM_LISTING_BEGIN_4 */
std::pair<double, double> runtimes_ltpSolve(unsigned int N) {
  // Runtime of Eigen's triangular solver and ltpSolve() in seconds
  double s_tria, s_ltp;

  // **********************************************************************
  // Code to be supplemented
  // **********************************************************************
  return {s_tria, s_ltp};
}
/* SAM_LISTING_END_4 */

// check accuracy ltpMult
void test_accuracy_ltpMult() {
  const std::size_t n = 4;
  Eigen::VectorXcd c1(n), c2(n), r1(n), r2(n), y(n);
  c1 << 1, 2, 3, 4;
  r1 << 1, 0, 0, 0;
  c2 << 5, 6, 7, 8;
  r2 << 5, 0, 0, 0;
  const Eigen::MatrixXcd T1 = toeplitz(c1, r1);
  const Eigen::MatrixXcd T2 = toeplitz(c2, r2);

  std::cout << "\nCheck that ltpMult is correct" << std::endl;
  const Eigen::VectorXcd c1c2 = ltpMult(c1, c2);
  const Eigen::MatrixXcd T1T2 = T1 * T2;
  std::cout << "Error = " << (c1c2 - T1T2.col(0)).norm() << std::endl;
}

// check accuracy ltpSolve
void test_accuracy_ltpSolve() {
  const std::size_t n = 4;
  Eigen::VectorXcd c(n), r(n), y(n);
  c << 1, 2, 3, 4;
  r.setZero();
  r(0) = c(0);
  y << 9, 10, 11, 12;
  Eigen::MatrixXcd T = toeplitz(c, r);

  std::cout << "\nCheck that ltpSolve is correct" << std::endl;
  const Eigen::VectorXcd u_rec = ltpSolve(c, y);
  const Eigen::VectorXcd u_sol = T.triangularView<Eigen::Lower>().solve(y);
  std::cout << "Error = " << (u_rec - u_sol).norm() << std::endl;
}

}  // namespace LowTriangToeplitz

// End of file
