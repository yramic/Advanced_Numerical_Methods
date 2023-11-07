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
  const unsigned int N = n * n;
  std::vector<Eigen::VectorXd> mu_vecs{M + 1, Eigen::VectorXd(N)};
  double tau = T * 1.0 / M;
  double h = 1.0 / (n + 1);
  Eigen::VectorXd w = cqWeights(M, tau);
// **********************************************************************
// Your Solution here
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
  double h = 1.0 / (n + 1);
  const unsigned int M = std::pow(2, L) - 1;
  double tau = T * 1.0 / M;
  auto delta = [](std::complex<double> z) {
    return 1.0 - z;
    //return 1.0 / 2.0 * z * z - 2.0 * z + 3.0 / 2.0;
  };
  // Initialize the numerical solution. This implementation is, for the sake of clarity, not memory-efficient.
  Eigen::MatrixXd mu_vecs(N, M + 1);
  // **********************************************************************
  // Your Solution here
  // **********************************************************************/
  return mu_vecs.col(M);
}
/* SAM_LISTING_END_5 */

}  // namespace FractionalHeatEquation

#endif
