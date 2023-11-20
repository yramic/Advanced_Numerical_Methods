/**
 * @file fractionalheatequation.h
 * @brief NPDE homework FractionalHeatEquation code
 * @author Jörg Nick, Bob Schreiner
 * @date October 2023
 * @copyright Developed at SAM, ETH Zurich
 */

#ifndef FRACTIONALHEATEQUATION_H_
#define FRACTIONALHEATEQUATION_H_

#include <Eigen/Dense>
#include <Eigen/Sparse>
#include <cassert>
#include <complex>
#include <iomanip>
#include <iostream>
#include <memory>
#include <stdexcept>
#include <unsupported/Eigen/FFT>
#include <unsupported/Eigen/KroneckerProduct>
#define assertm(exp, msg) assert(((void)msg, exp))

namespace FractionalHeatEquation {

/** @brief Encoding sparse matrix $\sqrt(s)*M + A$  */
/* SAM_LISTING_BEGIN_1 */
class SqrtsMplusA {
 public:
  // Constructor initializes sparse matrix
  SqrtsMplusA(unsigned int n, std::complex<double> s);
  SqrtsMplusA(const SqrtsMplusA &) = delete;
  SqrtsMplusA &operator=(const SqrtsMplusA &) = delete;
  SqrtsMplusA(SqrtsMplusA &&) = default;
  SqrtsMplusA &operator=(SqrtsMplusA &&) = default;
  virtual ~SqrtsMplusA() = default;

  // Matrix x vector
  template <typename SCALAR>
  [[nodiscard]] Eigen::VectorXcd eval(
      const Eigen::Matrix<SCALAR, Eigen::Dynamic, 1> &v) const;
  // Solve sparse linear system, LU decomposition on demand!
  template <typename SCALAR>
  [[nodiscard]] Eigen::VectorXcd solve(
      const Eigen::Matrix<SCALAR, Eigen::Dynamic, 1> &y);

 private:
  std::unique_ptr<Eigen::SparseMatrix<std::complex<double>>> p_matrix_;
  std::unique_ptr<Eigen::SparseLU<Eigen::SparseMatrix<std::complex<double>>>>
      p_LUdec_{nullptr};

 public:
  static unsigned int solve_cnt;
  static unsigned int ludec_cnt;
};
/* SAM_LISTING_END_1 */

template <typename SCALAR>
Eigen::VectorXcd SqrtsMplusA::eval(
    const Eigen::Matrix<SCALAR, Eigen::Dynamic, 1> &v) const {
  if (!p_matrix_) {
    throw std::runtime_error("eval: matrix not initialized!");
  }
  // Matrix-vector product
  return (*p_matrix_) * v.template cast<std::complex<double>>();
}

template <typename SCALAR>
Eigen::VectorXcd SqrtsMplusA::solve(
    const Eigen::Matrix<SCALAR, Eigen::Dynamic, 1> &v) {
  if (!p_LUdec_) {
    // Compute LU decomposition on demand
    p_LUdec_ =
        std::move(std::make_unique<
                  Eigen::SparseLU<Eigen::SparseMatrix<std::complex<double>>>>(
            *p_matrix_));
    if ((*p_LUdec_).info() != Eigen::Success) {
      throw std::runtime_error("LU factorization failed!");
    }
    ludec_cnt++;
  }
  solve_cnt++;
  return (*p_LUdec_).solve(v.template cast<std::complex<double>>());
}

/** @brief Compute CQ weights for IE-CQ and F(s) = sqrt(s) */
Eigen::VectorXd cqWeights(unsigned int M, double tau);

/** @brief Generate uniform grid of square with n elements along the edge.*/
std::vector<Eigen::Vector2d> generateGrid(unsigned n);

/** @brief Solve fully discrete evolution by MOT */
/* SAM_LISTING_BEGIN_2 */
template <typename SOURCEFN,
          typename RECORDER = std::function<void(const Eigen::VectorXd &)>>
Eigen::VectorXd evlMOT(
    SOURCEFN &&f, unsigned int n, double T, unsigned int M,
    RECORDER rec = [](const Eigen::VectorXd &mu_n) {}) {
  const unsigned int N = n * n;  // Number of FE d.o.f.s
  // Vector storing all the states; big memory consumption
  std::vector<Eigen::VectorXd> mu_vecs{M + 1, Eigen::VectorXd(N)};
  double tau = T * 1.0 / M;  // timestep size
  double h = 1.0 / (n + 1);  // meshwidth
  // See \prbcref{sp:1}
  Eigen::VectorXd w = cqWeights(M, tau);
// **********************************************************************
// PROBLEM 3-4b:

// The Task is to implement the MOT (Marching On in Time) Algorithm!

// First I need to initialize a matrix to invert at every timestep gridpoints
// and the rhs!
  SqrtsMplusA w0MplusA(n, std::complex<double>(std::pow(w[0], 2), 0.0));
  // Generate the Grid Points:
  std::vector<Eigen::Vector2d> gridpoints = generateGrid(n);
  // Generate RHS!
  Eigen::VectorXd rhs(N);
  // Note: idx_t is a time dependent indices
  for (int idx_t{0}; idx_t < (M + 1); ++idx_t) {
    // Next we need to evaluate the rhs and use a loop over the spacial component:
    for (int idx_x{0}; idx_x < N; ++idx_x) {
      rhs[idx_x] = f(idx_t * tau, gridpoints[idx_x]);
    }
    // Add memory tail from convolution:
    for (int eval_t{0}; eval_t < idx_t; ++eval_t) {
      rhs += -w[idx_t- eval_t] * mu_vecs[eval_t];
    }
    rhs *= h * h; // * Meshwidth squared!
    // Now the forward elimination algorithm can start:
    mu_vecs[idx_t] = w0MplusA.solve(rhs).real(); // Check if this is correct!
    // Print the result:
    std::cout << mu_vecs[idx_t] << std::endl;
    rec(mu_vecs[idx_t]);
  }
// **********************************************************************/
  return mu_vecs.back();
}
/* SAM_LISTING_END_2 */

/** @brief Real square Toeplitz matrix operator */
/* SAM_LISTING_BEGIN_4 */
class ToeplitzOp {
 public:
  explicit ToeplitzOp(const Eigen::VectorXd &);
  [[nodiscard]] Eigen::VectorXd eval(const Eigen::VectorXd &);

 private:
  Eigen::VectorXcd u_;
  Eigen::VectorXcd tmp_;
  Eigen::FFT<double> fft_;

 public:
  static unsigned int init_cnt;
  static unsigned int eval_cnt;
};
/* SAM_LISTING_END_4 */

/** @brief Solve fully discrete evolution by recursive algorithm for triangular
   Toeplitz system */
/* SAM_LISTING_BEGIN_3 */
template <typename SOURCEFN,
          typename RECORDER = std::function<void(const Eigen::VectorXd &)>>
Eigen::VectorXd evlTriangToeplitz(
    SOURCEFN &&f, unsigned int n, double T, unsigned int L,
    RECORDER rec = [](const Eigen::Vector2d &mu_n) {}) {
  const unsigned int N = n * n;
  const unsigned int M = std::pow(2, L) - 1;
  Eigen::MatrixXd mu_vecs(N, M + 1);
  const double tau = T * 1.0 / M;
  const double h = 1.0 / (n + 1);
  Eigen::VectorXd cq_weights = cqWeights(M, tau);
  // Initialise matrix to invert at every timestep, gridpoints and rhs
  SqrtsMplusA A(n, std::complex<double>(std::pow(cq_weights[0], 2), 0.0));
  //Generate rhs vector
  std::vector<Eigen::Vector2d> gridpoints = generateGrid(n);
  Eigen::MatrixXd rhs(N, (M + 1));
  for (int time_ind = 0; time_ind < M + 1; time_ind++) {
    // Evaluate rhs
    for (int space_ind = 0; space_ind < N; space_ind++) {
      rhs(space_ind, time_ind) =
          h * h * f(time_ind * tau, gridpoints[space_ind]);
    }
  }
  // ************************************************************
  // TO BE SUPPLEMENTED
  // Use recursive lambda function, see
  // https://gitlab.math.ethz.ch/NumCSE/NumCSE/-/blob/master/CppTutorial/lambdarecurse.cpp

  // PROBLEM 3-4d:
  std::function<Eigen::MatrixXd(unsigned int, Eigen::VectorXd, Eigen::MatrixXd, 
                Eigen::MatrixXd)> rec_solver = [&](unsigned int L, 
                                                  Eigen::VectorXd weights, 
                                                  Eigen::MatrixXd mu,
                                                  Eigen::MatrixXd phi) {
    // This lambda function takes the input by reference [&]
    // Note that weights are a vector according to the shape of cq_weights above
    // mu and phi are matrices since they are also defined above

    // Note that in 3-1 (not mandatory) we worked with vectors to compute convolutions
    // and when we worked with Toeplitz matrices! This could be done here also instead
    // of taking mu and phi in as matrices!

    // First we want to check if the dimension of mu and phi allign:
    assert(mu.size() == phi.size());


    // Solve directly if L = 0. Note that the number of timesteps is defined as
    // M = 2^L - 1. Hence, L = 0 would lead to M = 0 and no recursion is required!
    if (L == 0) {
      // .template cast<std::complex<double>>... deals with template types and ensures a
      // correct template conversion.
      // In particular cast<std::complex<double>> converts elements of the vector phi to
      // complex numbers (doubles)!
      const Eigen::VectorXcd y_direct = phi.template cast<std::complex<double>>();
      // Now return the direct solution if M and L = 0!
      return mu = (A.solve(y_direct)).real();
    } else {
      // else: Recursion required!
      const unsigned int n = mu.cols(), k = n/2; // 2*k - 1 entries that have to be multiplied!

      // Now we solve the UPPER LEFT part of the matrix recursively:
      mu.leftCols(k) = rec_solver(L-1, weights.head(k), mu.leftCols(k), phi.leftCols(k));

      // LOWER LEFT part can be solved with fft as in Problem 3-1 (LowTriangToeplitz)!
      ToeplitzOp T(h * h * weights.tail(n-1));
      for (unsigned l{0}; l < N; ++l) {
        phi.row(l).tail(k) = phi.row(l).tail(k) - T.eval(mu.row(l).head(k)).transpose();
      }

      // LOWER RIGHT part needs to be solved recursively again!
      mu.rightCols(k) = rec_solver(L-1, weights.head(k), mu.rightCols(k), phi.rightCols(k));
    }
    return mu;
  };
  // Start recursion:
  mu_vecs = rec_solver(L, cq_weights, mu_vecs, rhs);
  // ************************************************************
  return mu_vecs.col(M);
}
/* SAM_LISTING_END_3 */

/** @brief Solve fully discrete evolution by all-steps-in-one forward CQ */
/* SAM_LISTING_BEGIN_5 */
template <typename SOURCEFN,
          typename RECORDER = std::function<void(const Eigen::VectorXd &)>>
Eigen::VectorXd evlASAOCQ(
    SOURCEFN &&f, unsigned int n, double T, unsigned int L,
    RECORDER rec = [](const Eigen::VectorXd &mu_n) {}) {
  const unsigned int N = n * n;
  const unsigned int M = std::pow(2, L) - 1;
  double tau = T * 1.0 / M;
  double h = 1.0 / (n + 1.0);
  auto delta = [](std::complex<double> z) {
    return 1.0 - z;
    //return 1.0 / 2.0 * z * z - 2.0 * z + 3.0 / 2.0;
  };
  // Initialize the numerical solution. This implementation is, for the sake 
  // of clarity, not memory-efficient.
  Eigen::MatrixXd mu_vecs(N, M + 1);
  // **********************************************************************
  // PROBLEM 3-4e:
  // Convolution Quadrature Algorithm has to be implemented!

  // Note: for the FFT it's important to note that we work with arrays and 
  // zero padding has to be included, similar to 3-1!

  // Initialize a set of arrays for the whole rhs (all timepoints):
  Eigen::MatrixXd phi(N, M+1);
  // // Initialize an array for the rhs (single timepoint):
  // Eigen::VectorXd phi_part(N); // NOT NECESSARY!
  // Generate Grid Points:
  const std::vector<Eigen::Vector2d> gridpoints {generateGrid(N)};

  // Loop over Time idx_t
  for (int idx_t{0}; idx_t < (M + 1); ++ idx_t) {
    // Evaluate RHS!
    // Looper of Space idx_x:
    for (int idx_x{0}; idx_x < N; ++idx_x) {
      phi(idx_x, idx_t) = h * h * f(idx_t * tau, gridpoints[idx_x]);
    }
  }
  Eigen::FFT<double> fft;
  // Radius of Integral Contour
  double r {std::pow(1e-12, 1./(2.*M + 2.))};

  // Now we want to transform the rhs from the time domain to the frequency domain
  // Initialize Matrix for the frequency domain:
  Eigen::MatrixXcd phi_hat(N, M+1);
  for (unsigned int idx_t{0}; idx_t < (M+1); ++idx_t) {
    phi.col(idx_t) *= std::pow(r,idx_t);
  }
  for (unsigned int idx_x{0}; idx_x < N; ++idx_x) {
    phi_hat.row(idx_x) = fft.fwd(phi.row(idx_x));
  }

  // Next: Initialization of the numerical solution inside the frequency domain!
  Eigen::MatrixXcd mu_hat(N, M+1);

  for (unsigned int idx_freq{0}; idx_freq < (M+1); ++ idx_freq) {
    // Complex variable containing the imaginary unit and the discrete frequencies
    // Applying the time-harmonic operator G(s) where G(s) = (\sqrt(s)M + A)^{-1}
    std::complex<double> s = delta(r*std::exp(std::complex<double>(0, -2*M_PI*idx_freq/(M + 1))))/tau;
    SqrtsMplusA sMplusA(n, s);
    // Initialize Helper:
    Eigen::VectorXcd phi_hat_part(N);
    Eigen::VectorXcd mu_hat_part(N);
    phi_hat_part = phi_hat.col(idx_freq);
    mu_hat_part = sMplusA.solve(phi_hat_part);
    // Store the solution!
    mu_hat.col(idx_freq) = mu_hat_part;
  }
  // Now we need to transform the numerical solution from the frequency domain
  // into the time domain!
  for (unsigned int idx_x{0}; idx_x < N; ++idx_x) {
    mu_vecs.row(idx_x) = (fft.inv(mu_hat.row(idx_x))).real();
  }
  // Rescaling the numerical solution:
  for (unsigned int idx_t{0}; idx_t < (M+1); ++idx_t) {
    mu_vecs.col(idx_t) *= std::pow(r, -idx_t);
  }
  // **********************************************************************/
  return mu_vecs.col(M);
}
/* SAM_LISTING_END_5 */

}  // namespace FractionalHeatEquation

#endif
