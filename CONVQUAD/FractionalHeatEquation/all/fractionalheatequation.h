/**
 * @file fractionalheatequation.h
 * @brief NPDE homework FractionalHeatEquation code
 * @author Jörg Nick
 * @date October 2023
 * @copyright Developed at SAM, ETH Zurich
 */

#include <Eigen/Dense>
#include <Eigen/Sparse>
#include <cassert>
#include <complex>
#include <exception>
#include <iomanip>
#include <iostream>
#include <memory>
#include <stdexcept>

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
  Eigen::VectorXcd eval(
      const Eigen::Matrix<SCALAR, Eigen::Dynamic, 1> &v) const;
  // Solve sparse linear system, LU decomposition on demand!
  template <typename SCALAR>
  Eigen::VectorXcd solve(const Eigen::Matrix<SCALAR, Eigen::Dynamic, 1> &y);

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
  // ************************************************************
  return mu_vecs.col(M);
}
/* SAM_LISTING_END_3 */

}  // namespace FractionalHeatEquation
