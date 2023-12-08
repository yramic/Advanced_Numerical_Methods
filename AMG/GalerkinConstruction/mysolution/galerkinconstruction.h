/**
 * @file galerkinconstruction.h
 * @brief NPDE homework GalerkinConstruction code
 * @author R. Hiptmair
 * @date Novermber 2023
 * @copyright Developed at SAM, ETH Zurich
 */

#ifndef GC_H_
#define GC_H_

#include <Eigen/Core>
#include <Eigen/SparseCore>
#include <cassert>
#include <iomanip>
#include <iostream>
#include <vector>
#include <chrono>

#define assertm(exp, msg) assert(((void)msg, exp))

namespace GalerkinConstruction {

/** @brief Initialization of Poisson matrix as sparse matrix */
Eigen::SparseMatrix<double> poissonMatrix(unsigned int n);

/** @brief Initialization of prolongation matrix P as sparse matrix
 * @param M *number of cells* in each direction of the fine mesh
 *        (must be even)
 */
Eigen::SparseMatrix<double, Eigen::RowMajor> prolongationMatrix(
    unsigned int M, bool bilinear = true);

/** @brief Computation of coarse-grid matrix using Eigen's built-in sparse
 *         matrix arithmetic
 * @param A Galerkin matrix
 * @param P prolongation matrix
 *
 * Implementation based on the multiplication of Eigen sparse matrices
 */
Eigen::SparseMatrix<double> buildAH_eigen(
    const Eigen::SparseMatrix<double> &A,
    const Eigen::SparseMatrix<double, Eigen::RowMajor> &P);

/** @brief Computation of coarse-grid matrix using Eigen's built-in sparse
 *         matrix arithmetic
 * @param A Galerkin matrix
 * @param P prolongation matrix
 *
 * Loop and triplet-based implementation
 */
Eigen::SparseMatrix<double> buildAH(
    const Eigen::SparseMatrix<double> &A,
    const Eigen::SparseMatrix<double, Eigen::RowMajor> &P);

/* SAM_LISTING_BEGIN_X */
template <typename SEQ>
void tabulateRuntimes(SEQ &&M_vals) {
// **********************************************************************
// To be supplemented
// **********************************************************************
}
/* SAM_LISTING_END_X */

}  // namespace GalerkinConstruction

#endif
