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
std::pair<Eigen::MatrixXd,Eigen::MatrixXd> low_rank_merge(
        const Eigen::MatrixXd &A1, const Eigen::MatrixXd &B1,
        const Eigen::MatrixXd &A2, const Eigen::MatrixXd &B2)
{
    assert(A1.cols() == B1.cols() &&
           A2.cols() == B2.cols() &&
           A1.cols() == A2.cols() &&
           "All no.s of cols should be equal to q");

    size_t m = B1.rows();
    size_t n = B1.cols();
    Eigen::HouseholderQR<Eigen::MatrixXd> QR1 = B1.householderQr();
    Eigen::MatrixXd Q1 = QR1.householderQ() * Eigen::MatrixXd::Identity(m, std::min(m, n));
    Eigen::MatrixXd R1 = Eigen::MatrixXd::Identity(std::min(m, n), m) * QR1.matrixQR().triangularView<Eigen::Upper>();
    // About QR  decomposition with Eigen:
    // If $\VB_1: m \times n$, then $\VQ_1: m \times m$ and $\VR_1: m \times n$.
    // If $m > n$, however, the extra columns of $\VQ_1$ and extra rows of $\VR_1$ are not needed.
    // Matlab returns this "economy-size" format calling "qr(A,0)",
    // which does not compute these extra entries.
    // With the code above, Eigen is smart enough not to compute the discarded vectors.

    m = B2.rows();
    n = B2.cols();
    Eigen::HouseholderQR<Eigen::MatrixXd> QR2 = B2.householderQr();
    Eigen::MatrixXd Q2 = QR2.householderQ() * Eigen::MatrixXd::Identity(m, std::min(m, n));
    Eigen::MatrixXd R2 = Eigen::MatrixXd::Identity(std::min(m, n), m) * QR2.matrixQR().triangularView<Eigen::Upper>();

    Eigen::MatrixXd Z(A1.rows(), R1.rows()+R2.rows());
    Z << A1 * R1.transpose(), A2 * R2.transpose();
    Eigen::JacobiSVD<Eigen::MatrixXd> SVD(Z, Eigen::ComputeThinU | Eigen::ComputeThinV);
    Eigen::VectorXd s = SVD.singularValues();
    Eigen::MatrixXd S; S.setZero(A1.cols(), A1.cols());    // only consider first q singular values
    S.diagonal() = s.head(A1.cols());
    Eigen::MatrixXd U = SVD.matrixU().leftCols(A1.cols()); // only consider first q columns
    Eigen::MatrixXd V = SVD.matrixV();
    // About SVD decomposition with Eigen:
    // With Eigen::JacobiSVD you can ask for thin $\VU$ or $\VV$ to be computed.
    // In case of a rectangular $m \times n$ matrix,
    // with $j$ the smaller value among $m$ and $n$,
    // there can only be at most $j$ singular values.
    // The remaining columns of $\VU$ and $\VV$ do not correspond
    // to actual singular vectors and are not computed in thin format.

    Eigen::MatrixXd Q = Eigen::MatrixXd::Zero(Q1.rows()+Q2.rows(), Q1.cols()+Q2.cols());
    Q.block(0,0,
            Q1.rows(),Q1.cols()) = Q1;
    Q.block(Q1.rows(),Q1.cols(),
            Q2.rows(),Q2.cols()) = Q2;

    Eigen::MatrixXd Atilde =  U * S; // (2.4.38a)
    Eigen::MatrixXd Btilde = (Q * V).leftCols(A1.cols()); // (2.4.38b)

    return {Atilde,Btilde};
}
/* SAM_LISTING_END_0 */

/* SAM_LISTING_BEGIN_1 */
std::pair<double,double> test_low_rank_merge(size_t n)
{

    // Build X1, X2 defined in sub-problem (2.4.a)
    Eigen::MatrixXd X1(n,n), X2(n,n);
    for(double i=0.; i<n; ++i) {
        for(double j=0.; j<n; ++j) {
            X1(i,j) = std::sin((i-j)/n);
            X2(i,j) = std::cos((i-j-0.5)/n);
        }
    }

    // Build A1, B1, A2, B2 using HINT 1 (2-6.c)
    Eigen::MatrixXd A1(n,2), B1(n,2), A2(n,2), B2(n,2);
    for(double i=0.; i<n; ++i) {
        A1(i,0) = std::sin(i/n);
        A1(i,1) = std::cos(i/n);
        B1(i,0) = std::cos(i/n);
        B1(i,1) = -std::sin(i/n);
        A2(i,0) = std::cos(i/n);
        A2(i,1) = std::sin(i/n);
        B2(i,0) = std::cos((i+0.5)/n);
        B2(i,1) = std::sin((i+0.5)/n);
    }

    Eigen::MatrixXd Z(n,2*n);
    Z << X1, X2;

    std::pair<Eigen::MatrixXd,Eigen::MatrixXd> AB = low_rank_merge(A1, B1, A2, B2); // Call the function implemented in sub-problem (2.4.b)
    Eigen::MatrixXd Ztilde = AB.first * AB.second.transpose();

    Eigen::MatrixXd diff = Z - Ztilde;
    double err_Frob = diff.norm()/n; // $n^{-1} \| Z - \widetilde{Z} \|_F$
    double err_max  = diff.cwiseAbs().maxCoeff();

    return {err_Frob,err_max};
}
/* SAM_LISTING_END_1 */

/* SAM_LISTING_BEGIN_2 */
std::pair<Eigen::MatrixXd,Eigen::MatrixXd> adap_rank_merge(
        const Eigen::MatrixXd &A1, const Eigen::MatrixXd &B1,
        const Eigen::MatrixXd &A2, const Eigen::MatrixXd &B2, double rtol)
{
    assert(A1.cols() == B1.cols() &&
           A2.cols() == B2.cols() &&
           A1.cols() == A2.cols() &&
           "All no.s of cols should be equal to q");

    size_t m = B1.rows();
    size_t n = B1.cols();
    Eigen::HouseholderQR<Eigen::MatrixXd> QR1 = B1.householderQr();
    Eigen::MatrixXd Q1 = QR1.householderQ() * Eigen::MatrixXd::Identity(m, std::min(m, n));
    Eigen::MatrixXd R1 = Eigen::MatrixXd::Identity(std::min(m, n), m) * QR1.matrixQR().triangularView<Eigen::Upper>();

    m = B2.rows();
    n = B2.cols();
    Eigen::HouseholderQR<Eigen::MatrixXd> QR2 = B2.householderQr();
    Eigen::MatrixXd Q2 = QR2.householderQ() * Eigen::MatrixXd::Identity(m, std::min(m, n));
    Eigen::MatrixXd R2 = Eigen::MatrixXd::Identity(std::min(m, n), m) * QR2.matrixQR().triangularView<Eigen::Upper>();

    Eigen::MatrixXd Z(A1.rows(), R1.rows()+R2.rows());
    Z << A1 * R1.transpose(), A2 * R2.transpose();
    Eigen::JacobiSVD<Eigen::MatrixXd> SVD(Z, Eigen::ComputeThinU | Eigen::ComputeThinV);
    Eigen::VectorXd s = SVD.singularValues();

    unsigned p = s.size();
    for(unsigned q=1; q<s.size(); ++q) {
    // $q \in \{0 , ..., p-1\}$ : $\sigma_{q} \le \text{rtol} \cdot \sigma_0$
        if(s(q) <= s(0) * rtol) { 
            p = q;
            break;
        }
    }

    Eigen::MatrixXd S; S.setZero(p, p);
    S.diagonal() = s.head(p);
    Eigen::MatrixXd U = SVD.matrixU().leftCols(p);
    Eigen::MatrixXd V = SVD.matrixV().leftCols(p);

    Eigen::MatrixXd Q = Eigen::MatrixXd::Zero(Q1.rows()+Q2.rows(), Q1.cols()+Q2.cols());
    Q.block(0,0,
            Q1.rows(),Q1.cols()) = Q1;
    Q.block(Q1.rows(),Q1.cols(),
            Q2.rows(),Q2.cols()) = Q2;

    Eigen::MatrixXd Atilde = U * S;
    Eigen::MatrixXd Btilde = Q * V;

    return {Atilde,Btilde};

}
/* SAM_LISTING_END_2 */

/* SAM_LISTING_BEGIN_3 */
std::pair<double,size_t> test_adap_rank_merge(size_t n, double rtol) {

    Eigen::MatrixXd X1(n,n), X2(n,n);
    for(double i=0.; i<n; ++i) {
        for(double j=0.; j<n; ++j) {
            X1(i,j) = std::sin((i-j)/n);
            X2(i,j) = std::cos((i-j-0.5)/n);
        }
    }
    Eigen::MatrixXd Z(n,2*n);
    Z << X1, X2;

    // Build A1, B1, A2, B2 using HINT 1 (2-6.c)
    Eigen::MatrixXd A1(n,2), B1(n,2), A2(n,2), B2(n,2);
    for(double i=0.; i<n; ++i) {
        A1(i,0) = std::sin(i/n);
        A1(i,1) = std::cos(i/n);
        B1(i,0) = std::cos(i/n);
        B1(i,1) = -std::sin(i/n);
        A2(i,0) = std::cos(i/n);
        A2(i,1) = std::sin(i/n);
        B2(i,0) = std::cos((i+0.5)/n);
        B2(i,1) = std::sin((i+0.5)/n);
    }

    std::pair<Eigen::MatrixXd,Eigen::MatrixXd> AB = adap_rank_merge(A1, B1, A2, B2, rtol); // Call the function from sub-problem (2.4.d)
    Eigen::MatrixXd Ztilde = AB.first * AB.second.transpose();

    Eigen::MatrixXd diff = Z - Ztilde;
    double err_Frob = diff.norm()/n;

    return {err_Frob,AB.first.cols()};

}
/* SAM_LISTING_END_3 */

}  // namespace LowRankMerge
