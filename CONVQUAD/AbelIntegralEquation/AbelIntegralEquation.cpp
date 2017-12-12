#include <Eigen/Dense>
#include <cmath>
//#include <iomanip>
//#include <iostream>
extern "C" {
#include "../../../BEM/CppHilbert/Library/source/gaussQuadrature.h"
}


/* @brief Find the unknown function u in the Abel integral equation
 * using Galerkin discretization with a polynomial basis.
 * \param y Template function for the right-hand side
 * \param p Maximum degree of the polynomial basis and
 * order of the quadrature rule to compute the righ-hand side
 * \param tau Meshwidth of the grid where to compute the values of u
 * \\return Values of u on a grid in [0,1] with meshwidth tau
 */
/* SAM_LISTING_BEGIN_0 */
template<typename FUNC>
Eigen::VectorXd poly_spec_abel(const FUNC& y, size_t p, double tau)
{
#if SOLUTION
    Eigen::MatrixXd A = Eigen::MatrixXd::Zero(p,p);
    Eigen::VectorXd b = Eigen::VectorXd::Zero(p);

    const double* gauss_pts_pp = getGaussPoints(p*p);
    const double* gauss_wht_pp = getGaussWeights(p*p);

    const double* gauss_pts_p = getGaussPoints(p);
    const double* gauss_wht_p = getGaussWeights(p);

    for(int i=0; i<p; ++i) {

        for(int j=0; j<p; ++j) {

            for(int k=0; k<p*p; ++k) {

                double tk = 0.5 * (gauss_pts_pp[k] + 1.);
                double wk = 0.5 *  gauss_wht_pp[k];

                for(int l=0; l<p*p; ++l) {

                    double xl = 0.5 * std::sqrt(tk) * (gauss_pts_pp[k] + 1.);
                    double wl = 0.5 * std::sqrt(tk) *  gauss_wht_pp[k];

                    A(i,j) += wk*wl * 2. * std::pow(tk,i) * std::pow(tk - xl*xl,j);
                }
            }
        }

        for(int k=0; k<p; ++k) {

            double tk = 0.5 * (gauss_pts_p[k] + 1.);
            double wk = 0.5 *  gauss_wht_p[k];

            b(i) += wk * std::pow(tk,i) * y(tk);
        }
    }

    Eigen::VectorXd x = A.colPivHouseholderQr().solve(b);

    size_t N = std::round(1./tau) + 1;
    Eigen::VectorXd grid = Eigen::LinSpaced(N, 0., 1.);
    Eigen::VectorXd u    = Eigen::VectorXd::Zero(N);

    for(int i=0; i<N; ++i) {
        for(int j=0; j<p; ++j) {
            u(i) += x(j) * std::pow(grid(i),j);
        }
    }

    return grid;
    #else // TEMPLATE
    // TODO: Find the unknown function u in the Abel integral equation
    #endif // TEMPLATE
}
/* SAM_LISTING_END_0 */


/* @brief Compute errors between $\VZ$ and $\tilde{\VZ}$ in scaled Frobenius and max norms.
 * \param[in] n Number of rows/columns of matrices $VX_1$ and $VX_2$ defined as in Subproblem (2.1.a)
 */
/* SAM_LISTING_BEGIN_1 */
std::pair<double,double> test_low_rank_merge(size_t n)
{
#if SOLUTION
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
    #else // TEMPLATE
    // TODO: Compute {err_Frob,err_max}
    #endif // TEMPLATE
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

#if SOLUTION
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

    #else // TEMPLATE
    // TODO: Compute {Atilde,Btilde} as in (2.4.38a)/(2.4.38b), given (2.4.32)
    #endif // TEMPLATE
}
/* SAM_LISTING_END_2 */


/* @brief Compute the error between $\VZ$ and $\tilde{\VZ}$ in scaled Frobenius norm
 * and the number of singular values different from 0, given rtol.
 * \param[in] n Number of rows/columns of matrices $VX_1$ and $VX_2$ defined as in Subproblem (2.1.a)
 * \param[in] rtol Relative tolerance such that all singular values <= s_1 * rtol are set to 0.
 */
/* SAM_LISTING_BEGIN_3 */
std::pair<double,size_t> test_adap_rank_merge(size_t n, double rtol) {
#if SOLUTION

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

    #else // TEMPLATE
    // TODO: Compute {err_Frob,p}, with p := no. of singular values != 0
    #endif // TEMPLATE
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
    std::cout << "Problem 2.3.c" << std::endl;
    std::cout << "Fixed n = 500" << std::endl;
    unsigned n = 500;
    for(int p=1; p<=8; ++p) {
        double rtol = std::pow(10,-p);
        std::pair<double,size_t> errs = test_adap_rank_merge(n, rtol);
        std::cout << "n = " << n << std::setw(15)
                  << "Frobenius = "
                  << std::scientific << std::setprecision(3)
                  << errs.first  << std::setw(10)
                  << "Rank = "
                  << errs.second << std::endl;
    }
}
