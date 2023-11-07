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
 *
 */
Eigen::SparseMatrix<double> buildAH(const Eigen::SparseMatrix<double> &A,
                                    const Eigen::SparseMatrix<double> &P);

}  // namespace GalerkinConstruction

#endif
