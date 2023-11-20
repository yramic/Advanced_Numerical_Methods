/**
 * @file fractionalheatequation.cpp
 * @brief NPDE homework FractionalHeatEquation code
 * @author Jörg Nick
 * @date October 2023
 * @copyright Developed at SAM, ETH Zurich
 */

#include "fractionalheatequation.h"

#include <cassert>
#include <memory>

/* SAM_LISTING_BEGIN_9 */
namespace FractionalHeatEquation {
std::vector<Eigen::Vector2d> generateGrid(unsigned n) {
  std::vector<Eigen::Vector2d> gridpoints{n * n, Eigen::Vector2d()};
  double h = 1.0 / (n + 1);
  int rown = 0;
  for (int yind = 1; yind < n + 1; yind++) {
    for (int xind = 1; xind < n + 1; xind++, rown++) {
      Eigen::Vector2d tmp;
      tmp << xind * h, yind * h;
      gridpoints[rown] = tmp;
    }
  }
  return gridpoints;
}
/* SAM_LISTING_END_9 */

unsigned int SqrtsMplusA::solve_cnt{0};
unsigned int SqrtsMplusA::ludec_cnt{0};

SqrtsMplusA::SqrtsMplusA(unsigned int n, std::complex<double> s)
    : p_matrix_(std::make_unique<Eigen::SparseMatrix<std::complex<double>>>(
          n * n, n * n)) {
  (*p_matrix_).reserve(Eigen::RowVectorXi::Constant(n * n, 5));
  const double h = 1.0 / (1 + n);
  unsigned int rown = 0;
  for (int i = 0; i < n; ++i) {
    for (int j = 0; j < n; ++j, ++rown) {
      (*p_matrix_).insert(rown, rown) = 4.0 + h * h * std::sqrt(s);
      if (j > 0) (*p_matrix_).insert(rown, rown - 1) = -1.0;
      if (j < n - 1) (*p_matrix_).insert(rown, rown + 1) = -1.0;
      if (i > 0) (*p_matrix_).insert(rown, rown - n) = -1.0;
      if (i < n - 1) (*p_matrix_).insert(rown, rown + n) = -1.0;
    }
  }
}

/* SAM_LISTING_BEGIN_1 */
Eigen::VectorXd cqWeights(unsigned int M, double tau) {
  Eigen::VectorXd w(M + 1);
  // *************************************************************************
  // PROBLEM 3-4a:
  // Note that this part was already implemented!
  w(0) = std::pow(tau, -0.5);
  // Calculate weights of convolution quadrature based on \prbcref{eq:smuw}
  for (int l = 1; l < M + 1; ++l) {
    w(l) = w(l - 1) * (-1) * (0.5 - (l - 1)) / l;  // denominator is factorial
  }
  // *************************************************************************
  return w;
}
/* SAM_LISTING_END_1 */

unsigned int ToeplitzOp::eval_cnt{0};
unsigned int ToeplitzOp::init_cnt{0};

/* SAM_LISTING_BEGIN_3 */
ToeplitzOp::ToeplitzOp(const Eigen::VectorXd &v)
    : u_(Eigen::VectorXcd::Zero(v.size() + 1)),
      tmp_(Eigen::VectorXcd::Zero(v.size() + 1)) {
  const unsigned int n = v.size() + 1;
  assertm(n > 1, "Vector v missing");
  assertm((n % 2) == 0,
          "Only odd-length vectors can define square Toeplitz matrices");
  // Initialize vector defining circulant extension of Toeplitz matrix
  tmp_.head(n / 2) =
      (v.segment((n - 2) / 2, n / 2)).template cast<std::complex<double>>();
  tmp_[n / 2] = std::complex<double>(0.0, 0.0);
  tmp_.tail((n - 2) / 2) =
      (v.head((n - 2) / 2)).template cast<std::complex<double>>();
  // Store DFT of circulant-defining vector
  u_ = fft_.fwd(tmp_);
  init_cnt++;
}
/* SAM_LISTING_END_3 */

/* SAM_LISTING_BEGIN_2 */
Eigen::VectorXd ToeplitzOp::eval(const Eigen::VectorXd &x) {
  const unsigned int d = x.size();
  const unsigned int m = u_.size();
  assertm(d == m / 2, "Vector length mismatch");
  eval_cnt++;
  tmp_.head(d) = x.template cast<std::complex<double>>();
  tmp_.tail(d) = Eigen::VectorXcd::Zero(d);
  return (fft_.inv(u_.cwiseProduct(fft_.fwd(tmp_)).eval()).real()).head(d);
}
/* SAM_LISTING_END_2 */

}  // namespace FractionalHeatEquation
