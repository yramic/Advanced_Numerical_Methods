#include <Eigen/Dense>
#include <Eigen/QR>
#include <Eigen/SVD>
#include <cmath>
#include <iomanip>
#include <iostream>


/* @brief Compute the rank-q best approximation of [A1*B1 A2*B2]
 */
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
    // Matlab returns this ''economy-size´´ format calling ''qr(A,0)´´,
    // which does not compute these extra entries.
    // With the code above, Eigen is smart enough to not compute the discarded vectors.

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
    Q.block(1,1,
            Q1.rows(),Q1.cols()) = Q1;
    Q.block(Q1.rows(),Q1.cols(),
            Q2.rows(),Q2.cols()) = Q2;

    Eigen::MatrixXd Atilde =  U * S; // (2.4.38a)
    Eigen::MatrixXd Btilde = (Q * V).leftCols(A1.cols()); // (2.4.38b)

    return {Atilde,Btilde};
}
/* SAM_LISTING_END_0 */


/* @brief Compute errors between $\VZ$ and $\tilde{\VZ}$ in scaled Frobenius and max norms.
 * \param[in] n Number of rows/columns of matrices $VX_1$ and $VX_2$ defined as in Subproblem (2.1.a)
 */
/* SAM_LISTING_BEGIN_1 */
std::pair<double,double> test_low_rank_merge(size_t n)
{
    Eigen::MatrixXd X1(n,n), X2(n,n);
    for(double i=0.; i<n; ++i) {
        for(double j=0.; j<n; ++j) {
            X1(i,j) = std::sin((i-j)/n);
            X2(i,j) = std::cos((i-j-0.5)/n);
        }
    }
    Eigen::MatrixXd Z(n,2*n);
    Z << X1, X2;

    Eigen::JacobiSVD<Eigen::MatrixXd> SVD1(X1, Eigen::ComputeThinU | Eigen::ComputeThinV);
    Eigen::VectorXd s1 = SVD1.singularValues();
    Eigen::MatrixXd S1; S1.setZero(SVD1.rank(), SVD1.rank());
    S1.diagonal() = s1.head(SVD1.rank());
    Eigen::MatrixXd A1 = SVD1.matrixU().leftCols(SVD1.rank()) * S1;
    Eigen::MatrixXd B1 = SVD1.matrixV().leftCols(SVD1.rank());
    Eigen::JacobiSVD<Eigen::MatrixXd> SVD2(X2, Eigen::ComputeThinU | Eigen::ComputeThinV);
    Eigen::VectorXd s2 = SVD2.singularValues();
    Eigen::MatrixXd S2; S2.setZero(SVD2.rank(), SVD2.rank());
    S2.diagonal() = s2.head(SVD2.rank());
    Eigen::MatrixXd A2 = SVD2.matrixU().leftCols(SVD2.rank()) * S2;
    Eigen::MatrixXd B2 = SVD2.matrixV().leftCols(SVD2.rank());

    std::pair<Eigen::MatrixXd,Eigen::MatrixXd> AB = low_rank_merge(A1, B1, A2, B2);
    Eigen::MatrixXd Ztilde = AB.first * AB.second.transpose();

    Eigen::MatrixXd diff = Z - Ztilde;
    double err_Frob = diff.norm()/n;
    double err_max  = diff.cwiseAbs().maxCoeff();

    return {err_Frob,err_max};
}
/* SAM_LISTING_END_1 */


/* @brief Compute the rank-q best approximation of [A1*B1 A2*B2]
 * by setting all singular values <= rtol * s_1 to 0.
 */
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
    Q.block(1,1,
            Q1.rows(),Q1.cols()) = Q1;
    Q.block(Q1.rows(),Q1.cols(),
            Q2.rows(),Q2.cols()) = Q2;

    Eigen::MatrixXd Atilde = U * S;
    Eigen::MatrixXd Btilde = Q * V;

    return {Atilde,Btilde};

}
/* SAM_LISTING_END_2 */


/* @brief Compute the error between $\VZ$ and $\tilde{\VZ}$ in scaled Frobenius norm
 * and the number of singular values different from 0, given rtol.
 * \param[in] n Number of rows/columns of matrices $VX_1$ and $VX_2$ defined as in Subproblem (2.1.a)
 * \param[in] rtol Relative tolerance such that all singular values <= s_1 * rtol are set to 0.
 */
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

    Eigen::JacobiSVD<Eigen::MatrixXd> SVD1(X1, Eigen::ComputeThinU | Eigen::ComputeThinV);
    Eigen::VectorXd s1 = SVD1.singularValues();
    Eigen::MatrixXd S1; S1.setZero(SVD1.rank(), SVD1.rank());
    S1.diagonal() = s1.head(SVD1.rank());
    Eigen::MatrixXd A1 = SVD1.matrixU().leftCols(SVD1.rank()) * S1;
    Eigen::MatrixXd B1 = SVD1.matrixV().leftCols(SVD1.rank());
    Eigen::JacobiSVD<Eigen::MatrixXd> SVD2(X2, Eigen::ComputeThinU | Eigen::ComputeThinV);
    Eigen::VectorXd s2 = SVD2.singularValues();
    Eigen::MatrixXd S2; S2.setZero(SVD2.rank(), SVD2.rank());
    S2.diagonal() = s2.head(SVD2.rank());
    Eigen::MatrixXd A2 = SVD2.matrixU().leftCols(SVD2.rank()) * S2;
    Eigen::MatrixXd B2 = SVD2.matrixV().leftCols(SVD2.rank());

    std::pair<Eigen::MatrixXd,Eigen::MatrixXd> AB = adap_rank_merge(A1, B1, A2, B2, rtol);
    Eigen::MatrixXd Ztilde = AB.first * AB.second.transpose();

    Eigen::MatrixXd diff = Z - Ztilde;
    double err_Frob = diff.norm()/n;

    return {err_Frob,AB.first.cols()};

}
/* SAM_LISTING_END_3 */


int main() {

    std::cout << "Problem 2.3.c" << std::endl;
    for(int p=3; p<=12; ++p) {
        unsigned n = std::pow(2,p);
        std::pair<double,double> errs = test_low_rank_merge(n);
        std::cout << "n = " << n << std::setw(15)
                  << "Frobenius = "
                  << std::scientific << std::setprecision(3)
                  << errs.first  << std::setw(10)
                  << "Max = "
                  << std::scientific << std::setprecision(3)
                  << errs.second << std::endl;
    }

    std::cout << "Problem 2.3.e"       << std::endl;
    std::cout << "Fixed rtol = 0.0001" << std::endl;
    double rtol = 0.0001;
    for(int p=3; p<=12; ++p) {
        unsigned n = std::pow(2,p);
        std::pair<double,size_t> errs = test_adap_rank_merge(n, rtol);
        std::cout << "n = " << n << std::setw(15)
                  << "Frobenius = "
                  << std::scientific << std::setprecision(3)
                  << errs.first  << std::setw(10)
                  << "Rank = "
                  << errs.second << std::endl;
    }
    std::cout << "Problem 2.3.e" << std::endl;
    std::cout << "Fixed n = 500" << std::endl;
    unsigned n = 500;
    for(int p=1; p<=8; ++p) {
        double rtol = std::pow(10,-p);
        std::pair<double,size_t> errs = test_adap_rank_merge(n, rtol);
        std::cout << "rtol = " << rtol << std::setw(15)
                  << "Frobenius = "
                  << std::scientific << std::setprecision(3)
                  << errs.first  << std::setw(10)
                  << "Rank = "
                  << errs.second << std::endl;
    }
}
