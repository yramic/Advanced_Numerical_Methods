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
    // TODO: Compute {Atilde,Btilde} as in (2.4.38a)/(2.4.38b)
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
    // TODO: Compute {err_Frob,err_max}
}


int main() {

}
