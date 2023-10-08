/**
 * @file fractionalheatequation.h
 * @brief NPDE homework FractionalHeatEquation code
 * @author Jörg Nick
 * @date October 2023
 * @copyright Developed at SAM, ETH Zurich
 */

#ifndef FHE_H_
#define FHE_H_

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

/** @brief Solve fully discrete evolution my MOT */
/* SAM_LISTING_BEGIN_2 */
template <typename SOURCEFN,
          typename RECORDER = std::function<void(const Eigen::VectorXd &)>>
Eigen::VectorXd evlMOT(
    SOURCEFN &&f, unsigned int n, double T, unsigned int M,
    RECORDER rec = [](const Eigen::Vector2d &mu_n) {}) {
  const unsigned int N = n * n;
  std::vector<Eigen::VectorXd> mu_vecs{M + 1, Eigen::VectorXd(N)};
  // ************************************************************
  // TO BE SUPPLEMENTED
  // ************************************************************
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
  // https://gitlab.math.ethz.ch/NumCSE/NumCSE/-/blob/master/CppTutorial/lambdarecurse.cpp?ref_type=heads
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
    RECORDER rec = [](const Eigen::Vector2d &mu_n) {}) {
  const unsigned int N = n * n;
  const unsigned int M = std::pow(2, L) - 1;
  Eigen::MatrixXd mu_vecs(N, M + 1);
  // ************************************************************
  // TO BE SUPPLEMENTED
  // ************************************************************
  return mu_vecs.col(M);
}
/* SAM_LISTING_END_X */

}  // namespace FractionalHeatEquation

#endif
