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
  // Eigen::FFT<double> ... Object to perform FFT operations is created
  Eigen::FFT<double> fft; 
  // .fwd(x) ... Perform the forward FFT on the signal input x
  // forward FFT: transforms time-domain signals into frequency-
  // domain signals
  // .cwiseProduct ... elementwise product of two input signals 
  Eigen::VectorXcd tmp = (fft.fwd(u)).cwiseProduct(fft.fwd(x));
  // .inv ... transformation from the frequency domain to the time domain
  return fft.inv(tmp);
}

/* @brief Multiply a Toeplitz matrix with a vector, uses pconvfft
 * \param c Vector of entries of first column of the Toeplitz matrix
 * \param r Vector of entries of first row of the Toeplitz matrix
 * \param x Vector
 * \\return toeplitz(c,r)*x
 */
Eigen::VectorXcd toepMatVecMult(const Eigen::VectorXcd& c,
                                const Eigen::VectorXcd& r,
                                const Eigen::VectorXcd& x) {
  assert(c.size() == x.size() && r.size() == x.size() &&
         "c, r, x have different lengths!");

  std::size_t n = c.size();
  // The multiplication of a Toeplitz matrix with a vector can be reduced
  // We want a multiplication of a circulant matrix with a vector
  // This is because we can use then FFT!
  Eigen::VectorXcd cr_tmp(2 * n), x_tmp(2 * n);

  // Assemble the sequence of Toeplitz Matrix:
  // First elements have the entries of the columnvector c
  cr_tmp.head(n) = c; 
  // last elements set to 0 from N to N^2-1
  cr_tmp.tail(n) = Eigen::VectorXcd::Zero(n);
  // Now we take the row vector r and plug them in, in reversed order 
  cr_tmp.tail(n - 1) = r.tail(n - 1).reverse();

  // Set the values to x from 0 to N-1
  x_tmp.head(n) = x;
  // Zero Padding for values from N to N^2-1
  x_tmp.tail(n) = Eigen::VectorXcd::Zero(n);

  Eigen::VectorXcd y = pconvfft(cr_tmp, x_tmp);

  return y.head(n);
}

Eigen::VectorXcd ltpMultold(const Eigen::VectorXcd& f,
                            const Eigen::VectorXcd& g) {
  assert(f.size() == g.size() && "f and g vectors must have the same length!");

  std::size_t n = f.size();
  return toepMatVecMult(f, Eigen::VectorXcd::Zero(n), g);
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
  // Your Solution here PROBLEM 3-1C:

  // Note: VectorXcd explained
  // - X ... declares the dynamic size
  // - c ... matrices of complex numbers
  // - d ... real and imaginary parts are represented by doubles

  // The sequence of the product matrix can be obtained through discrete convolution
  // Discrete convolution can be treated as multiplication of circulant matrix, which
  // can be represented by Fourier matrix refpar:circul

  // First step according to Remark 4.1.4.15 is the need for a zero padding!
  Eigen::VectorXcd f_long = Eigen::VectorXcd::Zero(2*n); // Increased f!
  Eigen::VectorXcd g_long = Eigen::VectorXcd::Zero(2*n); // Increased g!

  f_long.head(n) = f; // Assign first n elements!
  g_long.head(n) = g;

  // Now the Periodic Discrete Convolution can be computed
  // For further information look at the function pconvfft()
  res = pconvfft(f_long, g_long).head(n);
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
  // Code to be supplemented PROBLEM 3-1D:
  
  // As done previously we will do the computation several times to have a
  // better approximation of the runtime duration!
  int n_run{6}; // 6 runs for each Operation!

  // First, we want to compute LOWER triangular Toeplitz matrices:
  // Hence, the row vector consists out of zeros and the column vector out
  // of random numbers:
  Eigen::VectorXcd r_vec(N), c1_vec(N), c2_vec(N);
  r_vec.setZero();
  c1_vec = Eigen::VectorXcd::Random(N);
  c2_vec = Eigen::VectorXcd::Constant(N, 1.0);

  // For the Teoplitz Matrix since it is set up with the first row represented
  // by r and the first column represented by c, the following condition:
  // r(0) == c(0) needs to be fullfilled!
  r_vec(0) = c1_vec(0);
  Eigen::MatrixXcd T_1 = toeplitz(c1_vec, r_vec);
  r_vec(0) = c2_vec(0);
  Eigen::MatrixXcd T_2 = toeplitz(c2_vec, r_vec);

  // 1) Matrix x Matrix:

  // Initially, we don't have a measured value, so we need to set s_dense 
  // to a very large number as a placeholder!
  s_dense = std::numeric_limits<double>::max(); // max possible value

  Eigen::MatrixXcd Test_MxM;
  for (int i{0}; i < n_run; ++i) {
    auto t_0 = std::chrono::high_resolution_clock::now();
    Test_MxM = T_1 * T_2;
    auto t_1 = std::chrono::high_resolution_clock::now();
    std::chrono::duration<double> ms_dur = (t_1 - t_0);
    // Minimal measured Time:
    // Note: .count() is necessary since two doubles need to be compared!
    s_dense = (ms_dur.count() < s_dense) ? ms_dur.count() : s_dense;
  }

  // 2) Matrix x Vector

  // Same approach here!
  s_mv = std::numeric_limits<double>::max();

  Eigen::VectorXcd Test_MxV;
  for (int i{0}; i < n_run; ++i) {
    auto t_0 = std::chrono::high_resolution_clock::now();
    Test_MxV = T_1 * c2_vec;
    auto t_1 = std::chrono::high_resolution_clock::now();
    std::chrono::duration<double> ms_dur = (t_1 - t_0);
    s_mv = (ms_dur.count() < s_mv) ? ms_dur.count() : s_mv;
  }

  // 3) Usage of ltpMult (Vector x Vector)

  // Same approach here again!
  s_ltp = std::numeric_limits<double>::max();

  Eigen::VectorXcd Test_VxV;
  for (int i{0}; i < n_run; ++i) {
    auto t_0 = std::chrono::high_resolution_clock::now();
    Test_VxV = ltpMult(c1_vec, c2_vec);
    auto t_1 = std::chrono::high_resolution_clock::now();
    std::chrono::duration<double> ms_dur = (t_1 - t_0);
    s_ltp = (ms_dur.count() < s_ltp) ? ms_dur.count() : s_ltp;
  }
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

  std::size_t n = f.size();
  // When the problem reduces to a scalar, solve directly!
  if (n == 1) {
    return y.cwiseQuotient(f);
  }

  Eigen::VectorXcd u_head = ltpSolve(f.head(n / 2), y.head(n / 2));
  Eigen::VectorXcd t =
      y.tail(n / 2) -
      toepMatVecMult(f.tail(n / 2), f.segment(1, n / 2).reverse(), u_head);
  Eigen::VectorXcd u_tail = ltpSolve(f.head(n / 2), t);
  Eigen::VectorXcd u(n);
  return u;
}
/* SAM_LISTING_END_3 */

/* SAM_LISTING_BEGIN_4 */
std::pair<double, double> runtimes_ltpSolve(unsigned int N) {
  // Runtime of Eigen's triangular solver and ltpSolve() in seconds
  double s_tria, s_ltp;

  // **********************************************************************
  // Problem 3-1f:

  // runtimes:
  const int runs {6};

  // Generate Toeplitz Matrix and Vector:
  Eigen::VectorXcd c(N), r(N), v(N);
  // Note: c... column, r... row, v... vector
  for (unsigned int i{0}; i < N; ++i) {
    if (i == 0) {
      c(i) = 1;
    }
    c(i) = 1/i;
  }
  r.setZero();
  r(0) = c(0); // Condition for the Toeplitz Matrix!
  v = Eigen::VectorXcd::Constant(N, 1.0);

  // Build the actual Toeplitz Matrix:
  Eigen::MatrixXcd T = toeplitz(c, r);
  // Build the right hand side vector:
  Eigen::VectorXcd v_rhs = ltpMult(c, v);

  // Vector for the solutions:
  Eigen::VectorXcd sol_direct(N);
  Eigen::VectorXcd sol_ltp(N);

  // Initialize s_tria and s_ltp:
  s_tria = std::numeric_limits<double>::max();
  s_ltp = std::numeric_limits<double>::max();

  // 1) Runtime estimation for direct lower triangular method:
  for (unsigned int i{0}; i < runs; ++i) {
    auto t0 = std::chrono::high_resolution_clock::now();
    sol_direct = T.triangularView<Eigen::Lower>().solve(v_rhs);
    auto t1 = std::chrono::high_resolution_clock::now();
    std::chrono::duration<double> dur = t1 - t0;
    s_tria = (dur.count() < s_tria) ? dur.count() : s_tria; 
  }

  // 2) Runtime estimation for ltp method:
  for (unsigned int i {0}; i < runs; ++i) {
    auto t0 = std::chrono::high_resolution_clock::now();
    sol_ltp = ltpSolve(c, v_rhs);
    auto t1 = std::chrono::high_resolution_clock::now();

    std::chrono::duration<double> dur = t1 - t0;
    s_ltp = (dur.count() < s_ltp) ? dur.count() : s_ltp;
  }


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
