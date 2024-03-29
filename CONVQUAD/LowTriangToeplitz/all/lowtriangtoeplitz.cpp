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
  // Circulant matrix $(C)_{l,j}=p_{l-j}$, Fourier matrix $F_n$
  // $C = F^{-1}_n diag(F_nu)F_n$   ref{eq:circmv}
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
#if SOLUTION
  // The sequence of the product matrix can be obtained through discrete convolution
  // Discrete convolution can be treated as multiplication of circulant matrix, which
  // can be represented by Fourier matrix ref{par:circul}

  // Zero padding
  Eigen::VectorXcd f_long = Eigen::VectorXcd::Zero(2 * n);
  Eigen::VectorXcd g_long = Eigen::VectorXcd::Zero(2 * n);
  f_long.head(n) = f;
  g_long.head(n) = g;
  // Periodic discrete convolution using FFT
  res = pconvfft(f_long, g_long).head(n);
#else
  // **********************************************************************
  // Your Solution here
  // **********************************************************************
#endif
  return res;
}
/* SAM_LISTING_END_0 */

/* SAM_LISTING_BEGIN_1 */
std::tuple<double, double, double> runtimes_ltpMult(unsigned int N) {
  // Runtime of matrix-matrix, matrix-vector and vector-vector multiplication
  // in seconds
  double s_dense, s_mv, s_ltp;

#if SOLUTION
  // Measure runtime several times
  const int num_repetitions = 6;

  // Sequence of Toeplitz matrices
  Eigen::VectorXcd c(N), r(N), v(N);
  // Initialization of the sequence of Toeplitz matrices
  c = Eigen::VectorXcd::Random(N);
  v = Eigen::VectorXcd::Constant(N, 1.0);

  // Generate dense representation of Toeplitz matrices
  r.setZero();
  r(0) = c(0);
  const Eigen::MatrixXcd T = toeplitz(c, r);
  r(0) = v(0);
  const Eigen::MatrixXcd V = toeplitz(v, r);

  // Runtime when using Eigen's built-in multiplication of dense matrices
  s_dense = std::numeric_limits<double>::max();
  Eigen::MatrixXcd T_mult_V;
  for (int k = 0; k < num_repetitions; k++) {
    // Use C++ chrono library to measure runtimes
    auto t1 = std::chrono::high_resolution_clock::now();
    T_mult_V = T * V;
    auto t2 = std::chrono::high_resolution_clock::now();

    // Getting number of seconds as a double
    std::chrono::duration<double> ms_double = (t2 - t1);

    // Taking the minimal measured time as the result
    s_dense = std::min(s_dense, ms_double.count());
  }

  // Runtime when multiplying one of the Toeplitz matrix
  // with the vector defining the second
  s_mv = std::numeric_limits<double>::max();
  Eigen::VectorXcd T_mult_v;
  for (int k = 0; k < num_repetitions; k++) {
    // Use C++ chrono library to measure runtimes
    auto t1 = std::chrono::high_resolution_clock::now();
    T_mult_v = T * v;
    auto t2 = std::chrono::high_resolution_clock::now();

    // Getting number of seconds as a double
    std::chrono::duration<double> ms_double = (t2 - t1);

    // Taking the minimal measured time as the result
    s_mv = std::min(s_mv, ms_double.count());
  }

  // Runtime when using ltpMult() from \prbcref{subprb:tp3}
  s_ltp = std::numeric_limits<double>::max();
  Eigen::VectorXcd c_conv_v;
  for (int k = 0; k < num_repetitions; k++) {
    // Use C++ chrono library to measure runtimes
    auto t1 = std::chrono::high_resolution_clock::now();
    c_conv_v = ltpMult(c, v);
    auto t2 = std::chrono::high_resolution_clock::now();

    // Getting number of seconds as a double
    std::chrono::duration<double> ms_double = (t2 - t1);

    // Taking the minimal measured time as the result
    s_ltp = std::min(s_ltp, ms_double.count());
  }
#else
  // **********************************************************************
  // Code to be supplemented
  // **********************************************************************
#endif
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
#if SOLUTION
  // The multiplication of a Toeplitz matrix with a vector can be reduced
  // to the multiplication of circulant matrix with a vector ref{emp:tpmv},
  // which can be computed by FFT as in ref{par:circul}
  Eigen::VectorXcd cr_tmp(2 * n), x_tmp(2 * n);

  // Assemble the sequence of Toeplitz matrix
  cr_tmp.head(n) = c;
  cr_tmp.tail(n) = Eigen::VectorXcd::Zero(n);
  // Need to reverse the row sequence and add it to the end
  cr_tmp.tail(n - 1) = r.tail(n - 1).reverse();

  x_tmp.head(n) = x;
  // Zero padding
  x_tmp.tail(n) = Eigen::VectorXcd::Zero(n);

  // Periodic discrete convolution
  y = pconvfft(cr_tmp, x_tmp);
#else
  // **********************************************************************
  // Code to be supplemented
  // **********************************************************************
#endif
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
#if SOLUTION
  // When it reduces to scalar, solve directly
  if (n == 1) {
    return y.cwiseQuotient(f);
  }

  // Solve by recursion:
  // Solve for the first half of $\Vu$
  const Eigen::VectorXcd u_head = ltpSolve(f.head(n / 2), y.head(n / 2));
  // Update the right hand side by subtracting the product of the lower left part
  // of Toeplitz matrix and the first half of $\Vu$
  const Eigen::VectorXcd t =
      y.tail(n / 2) -
      toepMatVecMult(f.tail(n / 2), f.segment(1, n / 2).reverse(), u_head);
  // Solve for the second half of $\Vu$
  // (upper left part of the matrix is identical to the lower right part)
  const Eigen::VectorXcd u_tail = ltpSolve(f.head(n / 2), t);
  // Assemble results
  u << u_head, u_tail;
#else
  // **********************************************************************
  // Code to be supplemented
  // **********************************************************************
#endif
  return u;
}
/* SAM_LISTING_END_3 */

/* SAM_LISTING_BEGIN_4 */
std::pair<double, double> runtimes_ltpSolve(unsigned int N) {
  // Runtime of Eigen's triangular solver and ltpSolve() in seconds
  double s_tria, s_ltp;

#if SOLUTION
  // Measure runtime several times
  int num_repititions = 6;

  // Initialization of Toeplitz matrix and vector
  Eigen::VectorXcd c(N), r(N), v(N);
  for (int i = 0; i < N; i++) {
    c(i) = i + 1;
  }
  v = Eigen::VectorXcd::Constant(N, 1.0);
  r.setZero();
  r(0) = c(0);

  // Generate dense representation of Toeplitz and rhs vector of LSE
  Eigen::MatrixXcd T = toeplitz(c, r);
  Eigen::VectorXcd T_mult_v = ltpMult(c, v);

  // Runtime when using Eigen's built-in triangular solver
  Eigen::VectorXcd u_sol;
  s_tria = std::numeric_limits<double>::max();
  for (int k = 0; k < num_repititions; k++) {
    // Use C++ chrono library to measure runtimes
    auto t1 = std::chrono::high_resolution_clock::now();
    u_sol = T.triangularView<Eigen::Lower>().solve(T_mult_v);
    auto t2 = std::chrono::high_resolution_clock::now();

    // Getting number of seconds as a double
    std::chrono::duration<double> ms_double = (t2 - t1);

    // Taking the minimal measured time as the result
    s_tria = std::min(s_tria, ms_double.count());
  }

  // Runtime when using ltpSolve() from \prbcref{subprb:tp5}
  Eigen::VectorXcd u_rec;
  s_ltp = std::numeric_limits<double>::max();
  for (int k = 0; k < num_repititions; k++) {
    // Use C++ chrono library to measure runtimes
    auto t1 = std::chrono::high_resolution_clock::now();
    u_rec = ltpSolve(c, T_mult_v);
    auto t2 = std::chrono::high_resolution_clock::now();

    // Getting number of seconds as a double
    std::chrono::duration<double> ms_double = (t2 - t1);

    // Taking the minimal measured time as the result
    s_ltp = std::min(s_ltp, ms_double.count());
  }

#else
  // **********************************************************************
  // Code to be supplemented
  // **********************************************************************
#endif
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
