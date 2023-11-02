/**
 * @file fractionalheatequation.h
 * @brief NPDE homework FractionalHeatEquation code
 * @author Jörg Nick
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
#define assertm(exp, msg) assert(((void)msg, exp))

namespace FractionalHeatEquation {

/** @brief Encoding sparse matrix \sqrt(s)*M + A  */
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
    RECORDER rec = [](const Eigen::Vector2d &mu_n) {}) {
  const unsigned int N = n * n;
  std::vector<Eigen::VectorXd> mu_vecs{M + 1, Eigen::VectorXd(N)};
  double tau = T * 1.0 / M;
  double h = 1.0 / (n + 1);
  Eigen::VectorXd w = cqWeights(M, tau);
#if SOLUTION
  // Initialise matrix to invert at every timestep, gridpoints and rhs
  SqrtsMplusA w0MplusA(n, std::complex<double>(std::pow(w[0], 2), 0.0));
  std::vector<Eigen::Vector2d> gridpoints = generateGrid(n);
  Eigen::VectorXd rhs(N);
  for (int time_ind = 0; time_ind < M + 1; time_ind++) {
    // Evaluate rhs
    for (int space_ind = 0; space_ind < N; space_ind++) {
      rhs[space_ind] = f(time_ind * tau, gridpoints[space_ind]);
    }
    // Add memory tail from convolution
    for (int l = 0; l < time_ind; l++) {
      rhs += -w[time_ind - l] * mu_vecs[l];
    }
    rhs *= h * h;
    // Solve timestep
    mu_vecs[time_ind] = w0MplusA.solve(rhs);
    rec(mu_vecs[time_ind]);
  }
#else
// **********************************************************************
// Your Solution here
// **********************************************************************/
#endif
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

/** @brief Solve fully discrete evolution my recursive algorithm for triangular
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
  Eigen::VectorXd cq_weights(M + 1);
  // ************************************************************
  // TO BE SUPPLEMENTED
  // Use recursive lambda function, see
  // https://gitlab.math.ethz.ch/NumCSE/NumCSE/-/blob/master/CppTutorial/lambdarecurse.cpp
  // ************************************************************
  return mu_vecs.col(M);
}
/* SAM_LISTING_END_3 */

/** @brief Solve fully discrete evolution my all-steps-in-one forward CQ */
/* SAM_LISTING_BEGIN_X */
template <typename SOURCEFN,
          typename RECORDER = std::function<void(const Eigen::VectorXd &)>>
Eigen::VectorXd evlASAOCQ(
    SOURCEFN &&f, unsigned int n, double T, unsigned int L,
    RECORDER rec = [](const Eigen::VectorXd &mu_n) {}) {
  const unsigned int N = n * n;
  const unsigned int M = std::pow(2, L) - 1;
  double tau = T * 1.0 / M;
  auto delta = [](std::complex<double> z) {
    return 1.0 / 2.0 * z * z - 2.0 * z + 3.0 / 2.0;
  };
  // Initialize the numerical solution. This implementation is, for the sake of clarity, not memory-efficient.
  Eigen::MatrixXd mu_vecs(N, M + 1);
  // Initialise array for the whole right hand side (all timepoints)
  Eigen::MatrixXd phi(N, M + 1);
  // Initialise array for the right hand side at a single timepoint
  Eigen::MatrixXd phi_slice(N);
  // Set radius of integral contour
  double r = std::pow(10, -16.0 / (2 * M + 2));
  // Set gridpoints
  std::vector<Eigen::Vector2d> gridpoints = generateGrid(n);
  for (int time_ind = 0; time_ind < M + 1; time_ind++) {
    // Evaluate rhs
    for (int space_ind = 0; space_ind < N; space_ind++) {
      phi_slice[space_ind] = f(time_ind * tau, gridpoints[space_ind]);
    }
    phi[time_ind] = std::pow(r, time_ind) * phi_slice;
  }
  // Transform the right-hand side from the time domain into the frequency domain
  Eigen::MatrixXcd phi_hat(N, M + 1);
  Eigen::FFT<double> fft;
  for (int space_ind = 0; space_ind < N; space_ind++) {
    Eigen::VectorXcd in = phi.row(space_ind);
    Eigen::VectorXcd out(N);
    out = fft.fwd(in);
    phi_hat.row(space_ind) = out;
  }
  // Initializing the frequency domain numerical solution
  Eigen::MatrixXcd mu_hat(N, M + 1);
  // Complex variable containing the imaginary unit and the discrete frequencies
  std::complex<double> s_l;
  std::complex<double> imag(0, 1);
  for (int freq_ind = 0; freq_ind < M + 1; freq_ind++) {
    s_l = delta(r * std::exp(2 * M_PI * imag * ((double)freq_ind) /
                             (double)(M + 1))) /
          tau;
    // Applying the time-harmonic operator $G(s_l)^{-1}= (\sqrt{s_l}M+A)^{-1}$
    SqrtsMplusA slMplusA(n, s_l);
    Eigen::VectorXcd phi_hat_slice(N);
    Eigen::VectorXcd mu_hat_slice(N);
    phi_hat_slice = phi_hat.col(freq_ind);
    mu_hat_slice = slMplusA.solve(phi_hat_slice);
    mu_hat.col(freq_ind) = mu_hat_slice;
  }
  // Transform the numerical solution from the frequency domain to the time domain
  for (int space_ind = 0; space_ind < N; space_ind++) {
    Eigen::VectorXcd in = mu_hat.row(space_ind);
    Eigen::VectorXcd out(N);
    out = fft.inv(in);
    mu_vecs.row(space_ind) = out;
  }
  // Rescaling of the numerical solution
  for (int time_ind = 0; time_ind < M + 1; time_ind++) {
    mu_vecs.col(time_ind) *= std::pow(r, -time_ind);
  }
  // ************************************************************
  // TO BE SUPPLEMENTED
  // ************************************************************
  return mu_vecs.col(M);
}
/* SAM_LISTING_END_X */

}  // namespace FractionalHeatEquation

#endif
