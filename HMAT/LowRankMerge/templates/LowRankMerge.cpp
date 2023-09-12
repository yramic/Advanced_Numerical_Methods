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

    // TODO: Compute {Atilde,Btilde} as in (2.4.38a)/(2.4.38b)
}
/* SAM_LISTING_END_0 */


/* @brief Compute errors between $\VZ$ and $\tilde{\VZ}$ in scaled Frobenius and max norms.
 * \param[in] n Number of rows/columns of matrices $VX_1$ and $VX_2$ defined as in Subproblem (2.1.a)
 */
/* SAM_LISTING_BEGIN_1 */
std::pair<double,double> test_low_rank_merge(size_t n)
{
    // TODO: Compute {err_Frob,err_max}
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

    // TODO: Compute {Atilde,Btilde} as in (2.4.38a)/(2.4.38b), given (2.4.32)
}
/* SAM_LISTING_END_2 */


/* @brief Compute the error between $\VZ$ and $\tilde{\VZ}$ in scaled Frobenius norm
 * and the number of singular values different from 0, given rtol.
 * \param[in] n Number of rows/columns of matrices $VX_1$ and $VX_2$ defined as in Subproblem (2.1.a)
 * \param[in] rtol Relative tolerance such that all singular values <= s_1 * rtol are set to 0.
 */
/* SAM_LISTING_BEGIN_3 */
std::pair<double,size_t> test_adap_rank_merge(size_t n, double rtol) {
    // TODO: Compute {err_Frob,p}, with p := no. of singular values != 0
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
