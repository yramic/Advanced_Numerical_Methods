/**
 * @file lowrankmerge.h
 * @brief NPDE homework LowRankMerge code
 * @author
 * @date
 * @copyright Developed at SAM, ETH Zurich
 */

#ifndef LOWRANKMERGE_H_
#define LOWRANKMERGE_H_

#include <Eigen/Dense>

namespace LowRankMerge {

/** @brief Compute the rank-q best approximation of Z = [X1 X2]
 *
 * @param Ai, Bi the factorized form of Xi: $X_i = A_i B_i^T$
 */
std::pair<Eigen::MatrixXd,Eigen::MatrixXd> low_rank_merge(
        const Eigen::MatrixXd &A1, const Eigen::MatrixXd &B1,
        const Eigen::MatrixXd &A2, const Eigen::MatrixXd &B2);


/** @brief Compute errors between $\VZ$ and $\tilde{\VZ}$ in scaled Frobenius and max norms.
 *
 * @param n Number of rows/columns of matrices $VK_1$ and $VK_2$ defined as in Subproblem (2-6.a)
 */
std::pair<double,double> test_low_rank_merge(size_t n);

/** @brief Compute the rank-q best approximation of Z = [X1 X2] 
 *  by setting all singular values <= rtol * s_1 (largest singular value) to 0.
 *
 * @param Ai, Bi the factorized form of Xi: $X_i = A_i B_i^T$
 * @param rtol Relative tolerance such that all singular values <= s_1 * rtol are set to 0.
 * @param atol Absolute tolerance such that all singular values <= atol are set to 0.
 */
std::pair<Eigen::MatrixXd,Eigen::MatrixXd> adap_rank_merge(
        const Eigen::MatrixXd &A1, const Eigen::MatrixXd &B1,
        const Eigen::MatrixXd &A2, const Eigen::MatrixXd &B2, double rtol, double atol);

/** @brief Compute the error between $\VZ$ and $\tilde{\VZ}$ in scaled Frobenius norm
 *  and the number of singular values different from 0, given rtol.
 * 
 * @param n Number of rows/columns of matrices $VK_1$ and $VK_2$ defined as in Subproblem (2-6.a)
 * @param rtol Relative tolerance such that all singular values <= s_1 * rtol are set to 0.
 */
std::pair<double,size_t> test_adap_rank_merge(size_t n, double rtol);

}  // namespace LowRankMerge

#endif
