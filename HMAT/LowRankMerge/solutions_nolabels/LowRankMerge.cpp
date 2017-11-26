//// 
//// Copyright (C) 2017 SAM (D-MATH) @ ETH Zurich
//// Author(s): curzuato < > 
//// Contributors:  dcasati 
//// This file is part of the AdvNumCSE repository.
////
#include <Eigen/Dense>
#include <Eigen/QR>
#include <Eigen/SVD>
#include <cmath>


/* @brief Compute the rank-q best approximation of [A1*B1 A2*B2]
 */
std::pair<Eigen::MatrixXd,Eigen::MatrixXd> low_rank_merge(
        const Eigen::MatrixXd &A1, const Eigen::MatrixXd &B1,
        const Eigen::MatrixXd &A2, const Eigen::MatrixXd &B2) {

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

    Eigen::MatrixXd Zhat(A1.rows(), R2.rows());
    Zhat << A1 * R1.transpose(), A2 * R2.transpose();
    Eigen::JacobiSVD<Eigen::MatrixXd> SVD(Zhat, Eigen::ComputeThinU | Eigen::ComputeThinV);
    Eigen::VectorXd s = SVD.singularValues();
    Eigen::MatrixXd S; S.setZero(s.size(), s.size());
    S.diagonal() = s;
    Eigen::MatrixXd U = SVD.matrixU();
    Eigen::MatrixXd V = SVD.matrixV();

    Eigen::MatrixXd Q = Eigen::MatrixXd::Zero(Q1.rows()+Q2.rows(), Q1.cols()+Q2.cols());
    Q.block(1,1,
            Q1.rows(),Q1.cols()) = Q1;
    Q.block(Q1.rows(),Q1.cols(),
            Q2.rows(),Q2.cols()) = Q2;

    Eigen::MatrixXd Atilde = U * S;
    Eigen::MatrixXd Btilde = Q * V;

    return {Atilde,Btilde};

}

// About QR  decomposition with Eigen:
// If $\VB_1: m \times n$, then $\VQ_1: m \times m$ and $\VR_1: m \times n$.
// If $m > n$, however, the extra columns of $\VQ_1$ and extra rows of $\VR_1$ are not needed.
// Matlab returns this ''economy-size´´ format calling ''qr(A,0)´´,
// which does not compute these extra entries.
// With the code above, Eigen is smart enough to not compute the discarded vectors.

// About SVD decomposition with Eigen:
// With Eigen::JacobiSVD you can ask for thin $\VU$ or $\VV$ to be computed.
// In case of a rectangular $m \times n$ matrix,
// with $j$ the smaller value among $m$ and $n$,
// there can only be at most $j$ singular values.
// The remaining columns of $\VU$ and $\VV$ do not correspond
// to actual singular vectors and are not computed in thin format.


/* @brief Compute errors between $\VZ$ and $\tilde{\VZ}$ in scaled Frobenius and max norms.
 * \param[in] n Number of rows/columns of matrices $VX_1$ and $VX_2$ defined as in Subproblem (2.1.a)
 */
std::pair<double,double> test_low_rank_merge(size_t n) {

    Eigen::MatrixXd X1(n,n), X2(n,n);
    for(unsigned i=0; i<n; ++i) {
        for(unsigned j=0; j<n; ++j) {
            X1(i,j) = std::sin((i-j)/n);
            X2(i,j) = std::cos((i-j-1.)/n);
        }
    }
    Eigen::MatrixXd Z(n,2*n);
    Z << X1, X2;



    std::pair<Eigen::MatrixXd,Eigen::MatrixXd> AB = low_rank_merge();
    Eigen::MatrixXd Ztilde = AB.first * AB.second.transpose();

    Eigen::MatrixXd diff = Z - Ztilde;
    double err_Frob = diff.norm()/n;
    double err_max  = diff.cwiseAbs().maxCoeff();

    return {err_Frob,err_max};

}


int main() {

}
