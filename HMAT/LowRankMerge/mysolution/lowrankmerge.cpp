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

  // TODO: Compute {Atilde,Btilde} as in \eqref{eq:lrfac1}/\eqref{eq:lrfac2}

  // Dummy solution
  return {Eigen::MatrixXd::Zero(3,3), Eigen::MatrixXd::Zero(3,3)};
}
/* SAM_LISTING_END_0 */

/* SAM_LISTING_BEGIN_1 */
std::pair<double, double> test_low_rank_merge(size_t n) {

  // TODO: Compute {err_Frob,err_max}, approximation error in 
  // scaled Frobunius norm and maximum norm

  // Dummy solution
  return {0, 0};
}
/* SAM_LISTING_END_1 */

/* SAM_LISTING_BEGIN_2 */
std::pair<Eigen::MatrixXd, Eigen::MatrixXd> adap_rank_merge(
    const Eigen::MatrixXd &A1, const Eigen::MatrixXd &B1,
    const Eigen::MatrixXd &A2, const Eigen::MatrixXd &B2, double rtol,
    double atol) {
  assert(A1.cols() == B1.cols() && A2.cols() == B2.cols() &&
         A1.cols() == A2.cols() && "All no.s of cols should be equal to q");

  // TODO: Compute {Atilde,Btilde} as in \eqref{eq:lrfac1}/\eqref{eq:lrfac2}, given \eqref{eq:adaptrunc}
  
  // Dummy solution, to be replaced
  return {Eigen::MatrixXd::Zero(3,3), Eigen::MatrixXd::Zero(3,3)};
}
/* SAM_LISTING_END_2 */

/* SAM_LISTING_BEGIN_3 */
std::pair<double, size_t> test_adap_rank_merge(size_t n, double rtol) {

  // TODO: Compute {err_Frob,p}, with p := no. of singular values larger than tolerance

  // Dummy solution, to be replaced
  return {0, 0};
}
/* SAM_LISTING_END_3 */

}  // namespace LowRankMerge
