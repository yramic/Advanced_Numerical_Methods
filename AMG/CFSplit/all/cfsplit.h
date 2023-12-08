/**
 * @file cfsplit.h
 * @brief NPDE homework CFSplit code
 * @author R. Hiptmair
 * @date November 2023
 * @copyright Developed at SAM, ETH Zurich
 */

#ifndef CFSplit_H_
#define CFSplit_H_

#include <stdio.h>

#include <Eigen/Dense>
#include <Eigen/Sparse>
#include <functional>
#include <iomanip>
#include <iostream>
#include <utility>

#include "mmio.h"

// Routines from Matrix Market C library
// Implementation in mmio.c

extern "C" {
int mm_read_banner(FILE *f, MM_typecode *matcode);
int mm_read_mtx_crd_size(FILE *f, int *M, int *N, int *nz);
char *mm_typecode_to_str(MM_typecode matcode);
}

namespace GalerkinConstruction {
/** @brief Computation of coarse-grid matrix using Eigen's built-in sparse
 *         matrix arithmetic
 * @param A Galerkin matrix
 * @param P prolongation matrix
 *
 * Loop and triplet-based implementation.
 * From ADVNCSE homework project "GalerkinConsrtruction"
 */
Eigen::SparseMatrix<double> buildAH(
    const Eigen::SparseMatrix<double> &A,
    const Eigen::SparseMatrix<double, Eigen::RowMajor> &P);
}  // namespace GalerkinConstruction

namespace CFSplit {
// Reading a Eigen sparse matrix from file in MM format
Eigen::SparseMatrix<double> readMMatrixFromFile(
    std::string filename = "galmat.mm");

/** @brief Class handling AMG iterative solution of sparse linear systems of
 * equations
 *
 * The constructor is in charge of the setup phase
 * - Coarse-fine splitting
 * - Construction of prolongation matrices
 */
/* SAM_LISTING_BEGIN_1 */
class RugeStuebenAMG {
 public:
  using nodeidx = unsigned int;  // Number of a degree of freedom = node
  enum NodeFlag : int { UNDECIDED = 0, FINE = 1, COARSE = -1 };
  // Standard constructors
  RugeStuebenAMG() = delete;
  RugeStuebenAMG(const RugeStuebenAMG &) = delete;
  RugeStuebenAMG(RugeStuebenAMG &&) = default;
  RugeStuebenAMG &operator=(const RugeStuebenAMG &) = delete;
  RugeStuebenAMG &operator=(RugeStuebenAMG &&) = default;
  // Constructor performing the setup phase for the AMG method
  template <typename RECORDER =
                std::function<void(const Eigen::SparseMatrix<double> &,
                                   const Eigen::SparseMatrix<double> &,
                                   const std::vector<NodeFlag> &)>>
  explicit RugeStuebenAMG(
      const Eigen::SparseMatrix<double> &A, unsigned int min_mat_size = 10,
      double tauC = 0.25, double sigma = 1.5,
      RECORDER &&rec = [](const Eigen::SparseMatrix<double> & /*AH*/,
                          const Eigen::SparseMatrix<double> & /*P*/,
                          const std::vector<NodeFlag> & /*flags*/) -> void {});
  virtual ~RugeStuebenAMG() = default;

  // Iteration operator for AMG
  void iterate(const Eigen::VectorXd &phi, Eigen::VectorXd &mu) const;

 private:
  // Determine sets of strongly connected nodes
  [[nodiscard]] std::pair<std::vector<std::vector<nodeidx>>,
                          std::vector<std::vector<nodeidx>>>
  getStrongConn(const Eigen::SparseMatrix<double> &A, double tauC) const;
  // Set smoothable node to F-nodes
  void setSmoothable(const Eigen::SparseMatrix<double> &A, double sigma,
                     std::vector<NodeFlag> &nodeflags) const;
  // Compute connectivity measure for all nodes
  [[nodiscard]] std::vector<unsigned int> connectivityMeasure(
      const std::vector<NodeFlag> &nodeflags,
      const std::vector<std::vector<nodeidx>> &Sj_star) const;
  // Set flags for nodes defining C/F splitting
  [[nodiscard]] std::vector<NodeFlag> CFsplit(
      const Eigen::SparseMatrix<double> &A,
      const std::vector<std::vector<nodeidx>> &Sj_star) const;
  // Compute the AMG prolongation matrix
  [[nodiscard]] Eigen::SparseMatrix<double, Eigen::RowMajor> computeP(
      const Eigen::SparseMatrix<double> &A,
      const std::vector<std::vector<nodeidx>> &Si,
      const std::vector<NodeFlag> &nodeflags) const;

  // Finest level: level of linear system of equations
  unsigned int L_;
  // Array (length L) of system matrices on all levels
  std::vector<Eigen::SparseMatrix<double>> Amats_;
  // Array (length L-1) of prolongation matrices
  std::vector<Eigen::SparseMatrix<double>> Pmats_;
  // LU decomposition of matrix on coarsest level
  Eigen::SparseLU<Eigen::SparseMatrix<double>, Eigen::COLAMDOrdering<int>>
      CoarseLU_;
};
/* SAM_LISTING_END_1 */

/* SAM_LISTING_BEGIN_C */
template <typename RECORDER>
RugeStuebenAMG::RugeStuebenAMG(const Eigen::SparseMatrix<double> &A,
                               unsigned int min_mat_size, double tauC,
                               double sigma, RECORDER &&rec) {
#if SOLUTION

#else
  // ************************************************************
  // Code to be supplemented
  // ************************************************************
#endif
}
/* SAM_LISTING_END_C */

/** @brief Measuring the convergence rate */
double cvgRateAMG(const Eigen::SparseMatrix<double> &A, double tol = 1.0E-4);

}  // namespace CFSplit

#endif
