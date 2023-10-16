/**
 * @file lowrankmerge.cpp
 * @brief NPDE homework LowRankMerge code
 * @author
 * @date
 * @copyright Developed at SAM, ETH Zurich
 */

#include "lowrankmerge.h"

#include <Eigen/QR>
#include <Eigen/SVD>

namespace LowRankMerge {

/* SAM_LISTING_BEGIN_0 */
std::pair<Eigen::MatrixXd, Eigen::MatrixXd> low_rank_merge(
    const Eigen::MatrixXd &A1, const Eigen::MatrixXd &B1,
    const Eigen::MatrixXd &A2, const Eigen::MatrixXd &B2) {
  assert(A1.cols() == B1.cols() && A2.cols() == B2.cols() &&
         A1.cols() == A2.cols() && "All no.s of cols should be equal to q");

#if SOLUTION
  size_t m = B1.rows();
  size_t n = B1.cols();
  // Find low-rank factors of B1 and B2 using QR decomposition as in (2.4.2.25)
  // \cref{par:mergetrunc}
  Eigen::HouseholderQR<Eigen::MatrixXd> QR1 = B1.householderQr();
  Eigen::HouseholderQR<Eigen::MatrixXd> QR2 = B2.householderQr();

  // Build the thin matrix R
  Eigen::MatrixXd R1 = Eigen::MatrixXd::Identity(std::min(m, n), m) *
                       QR1.matrixQR().triangularView<Eigen::Upper>();
  Eigen::MatrixXd R2 = Eigen::MatrixXd::Identity(std::min(m, n), m) *
                       QR2.matrixQR().triangularView<Eigen::Upper>();

  // About QR decomposition with Eigen:
  // If $\VB_1: m \times n$, then $\VQ_1: m \times m$ and $\VR_1: m \times n$.
  // If $m > n$, however, the extra columns of $\VQ_1$ and extra rows of $\VR_1$
  // are not needed. Matlab returns this "economy-size" format calling
  // "qr(A,0)", which does not compute these extra entries. With the code above,
  // Eigen is smart enough not to compute the discarded vectors.

  // Build $\hat Z$ from $A_i$ and $R_i$
  Eigen::MatrixXd Z(A1.rows(), R1.rows() + R2.rows());
  Z << A1 * R1.transpose(), A2 * R2.transpose();

  // Compute SVD of $\hat Z$ as in (2.4.2.25) \cref{par:mergetrunc}
  Eigen::JacobiSVD<Eigen::MatrixXd> SVD(
      Z, Eigen::ComputeThinU | Eigen::ComputeThinV);
  Eigen::VectorXd s = SVD.singularValues();

  // Only consider first q singular values
  Eigen::MatrixXd S;
  S.setZero(A1.cols(), A1.cols());
  S.diagonal() = s.head(A1.cols());

  // Only consider first q columns of U and V
  Eigen::MatrixXd U = SVD.matrixU().leftCols(A1.cols());
  Eigen::MatrixXd V = SVD.matrixV().leftCols(A1.cols());

  // Split V to be compatible with Q1 and Q2 in compressed format
  Eigen::MatrixXd V1 =
      Eigen::MatrixXd::Identity(m, std::min(m, n)) * V.topRows(std::min(m, n));
  Eigen::MatrixXd V2 = Eigen::MatrixXd::Identity(m, std::min(m, n)) *
                       V.bottomRows(std::min(m, n));

  // About SVD decomposition with Eigen:
  // With Eigen::JacobiSVD you can ask for thin $\VU$ or $\VV$ to be computed.
  // In case of a rectangular $m \times n$ matrix,
  // with $j$ the smaller value among $m$ and $n$,
  // there can only be at most $j$ singular values.
  // The remaining columns of $\VU$ and $\VV$ do not correspond
  // to actual singular vectors and are not computed in thin format.

  // Compute $\tilde A$ as in \eqref{eq:lrfac1}
  Eigen::MatrixXd Atilde = U * S;

  // Compute $\tilde B$ while avoiding recovering Q as a dense matrix
  Eigen::MatrixXd Btilde1 = QR1.householderQ() * V1;
  Eigen::MatrixXd Btilde2 = QR2.householderQ() * V2;
  Eigen::MatrixXd Btilde(Btilde1.rows() + Btilde2.rows(), Btilde1.cols());
  Btilde << Btilde1, Btilde2;

  // Return factors of $\tilde Z$
  return {Atilde, Btilde};

#else
  // TODO: Compute {Atilde,Btilde} as in \eqref{eq:lrfac1}/\eqref{eq:lrfac2}

  // Dummy solution
  return {Eigen::MatrixXd::Zero(3, 3), Eigen::MatrixXd::Zero(3, 3)};
#endif
}
/* SAM_LISTING_END_0 */

/* SAM_LISTING_BEGIN_1 */
std::pair<double, double> test_low_rank_merge(size_t n) {
#if SOLUTION
  double nf = static_cast<double>(n);  // convert data type
  // Build $K_1$, $K_2$ and $Z$ defined in \prbcref{subprb:1}
  Eigen::MatrixXd X1(n, n), X2(n, n);
  for (int i = 0.; i < n; ++i) {
    for (int j = 0.; j < n; ++j) {
      X1(i, j) = std::sin((i - j) / nf);
      X2(i, j) = std::cos((i - j - 0.5) / nf);
    }
  }
  Eigen::MatrixXd Z(n, 2 * n);
  Z << X1, X2;

  // Build A1, B1, A2, B2 using \prbcref{subprb:1}
  Eigen::MatrixXd A1(n, 2), B1(n, 2), A2(n, 2), B2(n, 2);
  for (int i = 0.; i < n; ++i) {
    A1(i, 0) = std::sin(i / nf);
    A1(i, 1) = std::cos(i / nf);
    B1(i, 0) = std::cos(i / nf);
    B1(i, 1) = -std::sin(i / nf);
    A2(i, 0) = std::cos(i / nf);
    A2(i, 1) = std::sin(i / nf);
    B2(i, 0) = std::cos((i + 0.5) / nf);
    B2(i, 1) = std::sin((i + 0.5) / nf);
  }

  // Call the function implemented in \prbcref<prb:lrm:subprb:2>
  std::pair<Eigen::MatrixXd, Eigen::MatrixXd> AB =
      low_rank_merge(A1, B1, A2, B2);

  // Compute low-rank approximation error
  Eigen::MatrixXd Ztilde = AB.first * AB.second.transpose();
  Eigen::MatrixXd diff = Z - Ztilde;

  double err_Frob = diff.norm() / n;            // scaled Frobenius norm
  double err_max = diff.cwiseAbs().maxCoeff();  // maximum norm

  // Return scaled Frobunius norm and maximum norm of approximation error
  return {err_Frob, err_max};

#else
  // TODO: Compute {err_Frob,err_max}, approximation error in
  // scaled Frobunius norm and maximum norm

  // Dummy solution
  return {0, 0};
#endif
}
/* SAM_LISTING_END_1 */

/* SAM_LISTING_BEGIN_2 */
std::pair<Eigen::MatrixXd, Eigen::MatrixXd> adap_rank_merge(
    const Eigen::MatrixXd &A1, const Eigen::MatrixXd &B1,
    const Eigen::MatrixXd &A2, const Eigen::MatrixXd &B2, double rtol,
    double atol) {
  assert(A1.cols() == B1.cols() && A2.cols() == B2.cols() &&
         A1.cols() == A2.cols() && "All no.s of cols should be equal to q");

#if SOLUTION
  size_t m = B1.rows();
  size_t n = B1.cols();
  // Find low-rank factors of B1 and B2 using QR decomposition as in (2.4.2.25)
  // \cref{par:mergetrunc}
  Eigen::HouseholderQR<Eigen::MatrixXd> QR1 = B1.householderQr();
  Eigen::HouseholderQR<Eigen::MatrixXd> QR2 = B2.householderQr();

  // Build the thin matrix R
  Eigen::MatrixXd R1 = Eigen::MatrixXd::Identity(std::min(m, n), m) *
                       QR1.matrixQR().triangularView<Eigen::Upper>();
  Eigen::MatrixXd R2 = Eigen::MatrixXd::Identity(std::min(m, n), m) *
                       QR2.matrixQR().triangularView<Eigen::Upper>();

  // Construct $\hat Z$ from $A_i$ and $R_i$
  Eigen::MatrixXd Z(A1.rows(), R1.rows() + R2.rows());
  Z << A1 * R1.transpose(), A2 * R2.transpose();

  // Compute SVD of $\hat Z$
  Eigen::JacobiSVD<Eigen::MatrixXd> SVD(
      Z, Eigen::ComputeThinU | Eigen::ComputeThinV);
  Eigen::VectorXd s = SVD.singularValues();

  // Find singular values larger than atol and rtol as in \eqref{eq:adaptrunc}
  unsigned p = s.size();
  for (unsigned q = 1; q < s.size(); ++q) {
    // $q \in \{1 , ..., p-1\}$ : $\sigma_{q} \le \text{rtol} \cdot \sigma_0$
    if ((s(q) <= s(0) * rtol) or (s(q) <= atol)) {
      p = q;
      break;
    }
  }

  // Only consider the largest p singular values
  Eigen::MatrixXd S;
  S.setZero(p, p);
  S.diagonal() = s.head(p);

  // Only consider the first p columns of U and V
  Eigen::MatrixXd U = SVD.matrixU().leftCols(p);
  Eigen::MatrixXd V = SVD.matrixV().leftCols(p);

  // Split V to be compatible with Q1 and Q2 in compressed format
  Eigen::MatrixXd V1 =
      Eigen::MatrixXd::Identity(m, std::min(m, n)) * V.topRows(std::min(m, n));
  Eigen::MatrixXd V2 = Eigen::MatrixXd::Identity(m, std::min(m, n)) *
                       V.bottomRows(std::min(m, n));

  // Compute $\tilde A$ as in \eqref{eq:lrfac1}
  Eigen::MatrixXd Atilde = U * S;

  // Compute $\tilde B$ while avoiding recovering Q as a dense matrix
  Eigen::MatrixXd Btilde1 = QR1.householderQ() * V1;
  Eigen::MatrixXd Btilde2 = QR2.householderQ() * V2;
  Eigen::MatrixXd Btilde(Btilde1.rows() + Btilde2.rows(), Btilde1.cols());
  Btilde << Btilde1, Btilde2;

  // Return factors of $\tilde Z$
  return {Atilde, Btilde};

#else
  // TODO: Compute {Atilde,Btilde} as in \eqref{eq:lrfac1}/\eqref{eq:lrfac2}, given \eqref{eq:adaptrunc}

  // Dummy solution, to be replaced
  return {Eigen::MatrixXd::Zero(3, 3), Eigen::MatrixXd::Zero(3, 3)};
#endif
}
/* SAM_LISTING_END_2 */

/* SAM_LISTING_BEGIN_3 */
std::pair<double, size_t> test_adap_rank_merge(size_t n, double rtol) {
#if SOLUTION
  // Make  sure that trigonometric functions receive float arguments
  // "Hidden integer arithmetic" is a trap of C/C++
  const double nf = static_cast<double>(n);
  // Build $K_1$, $K_2$ and $Z$ defined in \prbcref{prb:lrm:subprb:1}
  Eigen::MatrixXd X1(n, n), X2(n, n);
  for (int i = 0.; i < n; ++i) {
    for (int j = 0.; j < n; ++j) {
      X1(i, j) = std::sin((i - j) / nf);
      X2(i, j) = std::cos((i - j - 0.5) / nf);
    }
  }
  Eigen::MatrixXd Z(n, 2 * n);
  Z << X1, X2;

  // Build A1, B1, A2, B2 using \prbcref{subprb:1}
  Eigen::MatrixXd A1(n, 2), B1(n, 2), A2(n, 2), B2(n, 2);
  for (int i = 0.; i < n; ++i) {
    A1(i, 0) = std::sin(i / nf);
    A1(i, 1) = std::cos(i / nf);
    B1(i, 0) = std::cos(i / nf);
    B1(i, 1) = -std::sin(i / nf);
    A2(i, 0) = std::cos(i / nf);
    A2(i, 1) = std::sin(i / nf);
    B2(i, 0) = std::cos((i + 0.5) / nf);
    B2(i, 1) = std::sin((i + 0.5) / nf);
  }

  // Call the function from \prbcref{subprb:4}
  std::pair<Eigen::MatrixXd, Eigen::MatrixXd> AB =
      adap_rank_merge(A1, B1, A2, B2, rtol, __DBL_MIN__);

  // Compute low-rank approximation error
  Eigen::MatrixXd Ztilde = AB.first * AB.second.transpose();
  Eigen::MatrixXd diff = Z - Ztilde;

  double err_Frob = diff.norm() / n;  // scaled Frobenius norm

  // Return scaled Frobenius norm of approximation error and
  // the rank required to achieve the relative tolerance rtol
  return {err_Frob, AB.first.cols()};

#else
  // TODO: Compute {err_Frob,p}, with p := no. of singular values larger than
  // tolerance

  // Dummy solution, to be replaced
  return {0, 0};
#endif
}
/* SAM_LISTING_END_3 */

}  // namespace LowRankMerge
