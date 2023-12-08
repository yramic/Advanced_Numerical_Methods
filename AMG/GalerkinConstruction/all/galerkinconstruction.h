/**
 * @file galerkinconstruction.h
 * @brief NPDE homework GalerkinConstruction code
 * @author R. Hiptmair
 * @date Novermber 2023
 * @copyright Developed at SAM, ETH Zurich
 */

#ifndef GC_H_
#define GC_H_

#include <Eigen/Dense>
#include <Eigen/Sparse>
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
#if SOLUTION
  // Take the minimum over several runs as runtime
  int n_runs = 5;
  std::cout << "M\tEigen\t\tMine" << std::endl;

  for (auto M : M_vals) {
    // Obtain the fine-grid matrix
    const Eigen::SparseMatrix<double> Ah = poissonMatrix(M);
    // Generate the prolongation matrix
    const Eigen::SparseMatrix<double, Eigen::RowMajor> P =
        prolongationMatrix(M, true);
    // Measure runtime of Eigen's multiplication of sparse matrices
    double ms_eigen = std::numeric_limits<double>::max();
    for (int r = 0; r < n_runs; r++) {
      auto t1_eigen = std::chrono::high_resolution_clock::now();
      const Eigen::SparseMatrix<double> AH_eigen = buildAH_eigen(Ah, P);
      auto t2_eigen = std::chrono::high_resolution_clock::now();
      std::chrono::duration<double, std::milli> ms_double = t2_eigen - t1_eigen;
      ms_eigen = std::min(ms_eigen, ms_double.count());
    }

    // Measure runtime of \prbcref{sp:2a}
    double ms_mine = std::numeric_limits<double>::max();
    for (int r = 0; r < n_runs; r++) {
      auto t1_mine = std::chrono::high_resolution_clock::now();
      const Eigen::SparseMatrix<double> AH_mine = buildAH(Ah, P);
      auto t2_mine = std::chrono::high_resolution_clock::now();
      std::chrono::duration<double, std::milli> ms_double = t2_mine - t1_mine;
      ms_mine = std::min(ms_mine, ms_double.count());
    }

    std::cout << M << "\t" << ms_eigen << "\t\t" << ms_mine << std::endl;
  }
#else
// **********************************************************************
// To be supplemented
// **********************************************************************
#endif
}
/* SAM_LISTING_END_X */

}  // namespace GalerkinConstruction

#endif
