#include <Eigen/Dense>
#include <Eigen/QR>
#include <Eigen/SVD>
#include <cmath>


/* @brief Compute the rank-q best approximation of [A1*B1 A2*B2]
 */
/* SAM_LISTING_BEGIN_0 */
std::pair<Eigen::MatrixXd,Eigen::MatrixXd> low_rank_merge(
        const Eigen::MatrixXd &A1, const Eigen::MatrixXd &B1,
        const Eigen::MatrixXd &A2, const Eigen::MatrixXd &B2) {
    // TODO: Compute {Atilde,Btilde} as in (2.4.38a)/(2.4.38b)
}
/* SAM_LISTING_END_0 */

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
/* SAM_LISTING_BEGIN_1 */
std::pair<double,double> test_low_rank_merge(size_t n) {
    // TODO: Compute {err_Frob,err_max}
}
/* SAM_LISTING_END_1 */


/* @brief Compute the rank-q best approximation of [A1*B1 A2*B2]
 * by setting all singular values <= rtol * s_1 to 0.
 */
/* SAM_LISTING_BEGIN_2 */
std::pair<Eigen::MatrixXd,Eigen::MatrixXd> adap_rank_merge(
        const Eigen::MatrixXd &A1, const Eigen::MatrixXd &B1,
        const Eigen::MatrixXd &A2, const Eigen::MatrixXd &B2, double rtol) {
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

    for(unsigned p=3; p<=12; ++p) {
        unsigned n = std::pow(2,p);
        std::pair<double,double> errs = test_low_rank_merge(n);
        std::cout << "n = " << n << std::setw(10)
                  << "Frobenius = "
                  << std::scientific << std::setprecision(3)
                  << errs.first  << std::setw(10)
                  << "Max = "
                  << std::scientific << std::setprecision(3)
                  << errs.second << std::endl;
    }

    double rtol = 0.0001;
    for(unsigned p=3; p<=12; ++p) {
        unsigned n = std::pow(2,p);
        std::pair<double,size_t> errs = test_adap_rank_merge(n, rtol);
        std::cout << "n = " << n << std::setw(10)
                  << "Frobenius = "
                  << std::scientific << std::setprecision(3)
                  << errs.first  << std::setw(10)
                  << "Rank = "
                  << errs.second << std::endl;
    }

    unsigned n = 500;
    for(unsigned p=1; p<=8; ++p) {
        double rtol = std::pow(10,-p);
        std::pair<double,size_t> errs = test_adap_rank_merge(n, rtol);
        std::cout << "n = " << n << std::setw(10)
                  << "Frobenius = "
                  << std::scientific << std::setprecision(3)
                  << errs.first  << std::setw(10)
                  << "Rank = "
                  << errs.second << std::endl;
    }
}
