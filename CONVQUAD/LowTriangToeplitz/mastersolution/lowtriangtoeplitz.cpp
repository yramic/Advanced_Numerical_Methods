#include <Eigen/Dense>
#include <cmath>
#include <ctime>
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
  const auto m = c.size();
  const auto n = r.size();
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

  const auto n = c.size();
  Eigen::VectorXcd cr_tmp(2 * n);
  Eigen::VectorXcd x_tmp(2 * n);

  cr_tmp.head(n) = c;
  cr_tmp[n] = 0;
  cr_tmp.tail(n - 1) = r.tail(n - 1).reverse();

  x_tmp.head(n) = x;
  x_tmp.tail(n) = Eigen::VectorXcd::Zero(n);

  Eigen::VectorXcd y = pconvfft(cr_tmp, x_tmp);

  return y.head(n);
}

Eigen::VectorXcd ltpMultold(const Eigen::VectorXcd& f,
                            const Eigen::VectorXcd& g) {
  assert(f.size() == g.size() && "f and g vectors must have the same length!");

  const auto n = f.size();
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
  const auto n = f.size();
  Eigen::VectorXcd res(n);
  Eigen::VectorXcd f_long = Eigen::VectorXcd::Zero(2 * n);
  Eigen::VectorXcd g_long = Eigen::VectorXcd::Zero(2 * n);
  f_long.head(n) = f;
  g_long.head(n) = g;
  res = pconvfft(f_long, g_long).head(n);
  return res;
}
/* SAM_LISTING_END_0 */

/* @brief Solve a linear problem involving a lower triangular Toeplitz matrix
 * \param f Vector of entries of lower triangular Toeplitz matrix
 * \param y Right-hand side of linear problem
 * \\return Solution of linear problem
 */
/* SAM_LISTING_BEGIN_1 */
Eigen::VectorXcd ltpSolve(const Eigen::VectorXcd& f,
                          const Eigen::VectorXcd& y) {
  assert(f.size() == y.size() && "f and y vectors must have the same length!");
  assert(abs(f(0)) > 1e-10 &&
         "Lower triangular Toeplitz matrix must be invertible!");
  assert(log2(f.size()) == floor(log2(f.size())) &&
         "Size of f must be a power of 2!");

  const auto n = f.size();
  if (n == 1) {
    return y.cwiseQuotient(f);
  }
  Eigen::VectorXcd u(n);

  // implements the algorithm presented in \lref{par:dcealg}
  const Eigen::VectorXcd u_head = ltpSolve(f.head(n / 2), y.head(n / 2)); // step 1
  const Eigen::VectorXcd t = y.tail(n / 2) -
      toepMatVecMult(f.tail(n / 2), f.segment(1, n / 2).reverse(), u_head); // step 2
  const Eigen::VectorXcd u_tail = ltpSolve(f.head(n / 2), t); // step 3
  u << u_head, u_tail;
  return u;
}
/* SAM_LISTING_END_1 */

// check accuracy ltpMult
void test_accuracy_ltpMult() {
  std::size_t n = 4;
  Eigen::VectorXcd c1(n), c2(n), r1(n), r2(n), y(n);
  c1 << 1, 2, 3, 4;
  r1 << 1, 0, 0, 0;
  c2 << 5, 6, 7, 8;
  r2 << 5, 0, 0, 0;
  Eigen::MatrixXcd T1 = toeplitz(c1, r1);
  Eigen::MatrixXcd T2 = toeplitz(c2, r2);

  std::cout << "\nCheck that ltpMult is correct" << std::endl;
  Eigen::VectorXcd c1c2 = ltpMult(c1, c2);
  Eigen::MatrixXcd T1T2 = T1 * T2;
  std::cout << "Error = " << (c1c2 - T1T2.col(0)).norm() << std::endl;
}

// measure time ltpMult
void time_measure_ltpMult() {
  std::size_t nl = 12;
  std::size_t n_start = 4;
  std::size_t n_end = n_start * pow(2, nl - 1);
  int num_repititions = 6;

  Eigen::VectorXd error(nl), et_slow(nl), et_fast(nl);
  std::clock_t start_time, end_time;
  double et_sum;

  std::cout << "\nMatrix size, start: " << n_start << std::endl;
  std::cout << "Matrix size, end: " << n_end << std::endl;
  std::cout << "Number of matrices: " << nl << "\n" << std::endl;

  std::size_t n = n_start;
  for (int l = 0; l < nl; l++) {
    Eigen::VectorXcd c(n), r(n), v(n);
    c = Eigen::VectorXcd::Random(n);
    v = Eigen::VectorXcd::Constant(n, 1.0);
    r.setZero();
    r(0) = c(0);

    Eigen::MatrixXcd T = toeplitz(c, r);

    et_sum = 0;
    Eigen::VectorXcd T_mult_v;
    for (int k = 0; k < num_repititions; k++) {
      start_time = clock();
      T_mult_v = T * v;
      end_time = clock();
      if (k > 0) et_sum += static_cast<double>(end_time - start_time) / CLOCKS_PER_SEC;
    }
    et_slow(l) = et_sum / (num_repititions - 1);

    et_sum = 0;
    Eigen::VectorXcd c_conv_v;
    for (int k = 0; k < num_repititions; k++) {
      start_time = clock();
      c_conv_v = ltpMult(c, v);
      end_time = clock();
      if (k > 0) et_sum += static_cast<double>(end_time - start_time) / CLOCKS_PER_SEC;
    }
    et_fast(l) = et_sum / (num_repititions - 1);

    error(l) = (c_conv_v - T_mult_v.col(0)).norm();
    std::cout << l << "\t" << n << "\t" << error(l) << "\t" << et_slow(l)
              << "\t" << et_fast(l) << std::endl;

    n *= 2;
  }
}

// check accuracy ltpSolve
void test_accuracy_ltpSolve() {
  std::size_t n = 4;
  Eigen::VectorXcd c(n), r(n), y(n);
  c << 1, 2, 3, 4;
  r.setZero();
  r(0) = c(0);
  y << 9, 10, 11, 12;
  Eigen::MatrixXcd T = toeplitz(c, r);

  std::cout << "\nCheck that ltpSolve is correct" << std::endl;
  Eigen::VectorXcd u_rec = ltpSolve(c, y);
  Eigen::VectorXcd u_sol = T.triangularView<Eigen::Lower>().solve(y);
  std::cout << "Error = " << (u_rec - u_sol).norm() << std::endl;
}

// measure time ltpSolve
void time_measure_ltpSolve() {
  std::size_t nl = 13;
  std::size_t n_start = 4;
  std::size_t n_end = n_start * std::pow(2, nl - 1);
  int num_repititions = 6;

  Eigen::VectorXd error(nl), et_slow(nl), et_fast(nl);
  std::clock_t start_time, end_time;
  double et_sum;

  std::cout << "\nMatrix size, start: " << n_start << std::endl;
  std::cout << "Matrix size, end: " << n_end << std::endl;
  std::cout << "Number of matrices: " << nl << "\n" << std::endl;

  std::size_t n = n_start;
  for (int l = 0; l < nl; l++) {
    Eigen::VectorXcd c(n), r(n), v(n);
    for (int i = 0; i < n; i++) {
      c(i) = i + 1;
    }
    v = Eigen::VectorXcd::Constant(n, 1.0);
    r.setZero();
    r(0) = c(0);

    Eigen::MatrixXcd T = toeplitz(c, r);
    Eigen::VectorXcd T_mult_v = ltpMult(c, v);

    et_sum = 0;
    Eigen::VectorXcd u_sol;
    for (int k = 0; k < num_repititions; k++) {
      start_time = clock();
      u_sol = T.triangularView<Eigen::Lower>().solve(T_mult_v);
      end_time = clock();
      if (k > 0) et_sum += double(end_time - start_time) / CLOCKS_PER_SEC;
    }
    et_slow(l) = et_sum / (num_repititions - 1);

    et_sum = 0;
    Eigen::VectorXcd u_rec;
    for (int k = 0; k < num_repititions; k++) {
      start_time = clock();
      u_rec = ltpSolve(c, T_mult_v);
      end_time = clock();
      if (k > 0) et_sum += double(end_time - start_time) / CLOCKS_PER_SEC;
    }
    et_fast(l) = et_sum / (num_repititions - 1);

    error(l) = (u_sol - u_rec).norm() / (T_mult_v).norm();
    std::cout << l << "\t" << n << "\t" << error(l) << "\t" << et_slow(l)
              << "\t" << et_fast(l) << std::endl;

    n *= 2;
  }
}

}  // namespace LowTriangToeplitz

// End of file
