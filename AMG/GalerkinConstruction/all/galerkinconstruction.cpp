/**
 * @file galerkinconstruction.cpp
 * @brief NPDE homework GalerkinConstruction code
 * @author
 * @date
 * @copyright Developed at SAM, ETH Zurich
 */

#include "galerkinconstruction.h"

#include <Eigen/src/Core/util/Constants.h>

namespace GalerkinConstruction {
/* SAM_LISTING_BEGIN_1 */
using triplet = Eigen::Triplet<double>;
Eigen::SparseMatrix<double> poissonMatrix(unsigned int n) {
  const int N = (n - 1) * (n - 1);  // Matrix size
  // define vector of triplets and reserve memory
  std::vector<triplet> entries;  // For temporary COO format
  // set the vector of triplets
  for (int block_id = 0; block_id < (n - 1); block_id++) {
    const int start_id = block_id * (n - 1);
    const int end_id = (block_id + 1) * (n - 1) - 1;
    // tri-diagonal matrix T
    for (int i = start_id; i <= end_id; i++) {
      entries.emplace_back(i, i, 4);  // Diagonal entry
      if (i > start_id) {             // T(i-1,i) = T(i,i-1)
        entries.emplace_back(i, i - 1, -1.0);
        entries.emplace_back(i - 1, i, -1.0);
      }
    }
    // identity blocks
    if (block_id > 0) {
      for (int i = start_id; i <= end_id; i++) {
        entries.emplace_back(i, i - (n - 1), -1.0);
        entries.emplace_back(i - (n - 1), i, -1.0);
      }
    }
  }
  // create the sparse matrix in CCS format
  Eigen::SparseMatrix<double> A(N, N);
  A.setFromTriplets(entries.begin(), entries.end());
  return A;
}
/* SAM_LISTING_END_1 */

Eigen::SparseMatrix<double, Eigen::RowMajor> prolongationMatrix(
    unsigned int M, bool bilinear) {
  assertm(((M > 3) and (M % 2 == 0)), "prolongationMatrix: M must be even!");
  const unsigned int N = (M - 1) * (M - 1);  // No of inter nodes of fine grid
  const unsigned int m =
      M / 2;  // Number of cells in each direction of coarse grid
  const unsigned int n =
      (m - 1) * (m - 1);  // No of interior nodes of fine grid
  // define vector of triplets and reserve memory
  std::vector<triplet> P_trp;  // For temporary COO format
  // Conversion of node positions to indices: lexikographic ordering
  auto idxH = [&m](unsigned int I, unsigned int J) -> unsigned int {
    assertm(((I > 0) and (I < m) and (J > 0) and (J < m)),
            "idxH: index out of range");
    return (I - 1) + (J - 1) * (m - 1);
  };
  auto idxh = [&M](unsigned int i, unsigned int j) -> unsigned int {
    assertm(((i > 0) and (i < M) and (j > 0) and (j < M)),
            "idxh: index out of range");
    return (i - 1) + (j - 1) * (M - 1);
  };
  // Traverse all nodes of the coarse mesh
  for (unsigned int I = 1; I < m; ++I) {
    for (unsigned int J = 1; J < m; ++J) {
      // Set values according to prolongation stencil
      const unsigned int i = 2 * I;
      const unsigned int j = 2 * J;
      const unsigned int idx_H = idxH(I, J);
      P_trp.emplace_back(idxh(i, j), idx_H, 1.0);
      P_trp.emplace_back(idxh(i + 1, j), idx_H, 0.5);
      P_trp.emplace_back(idxh(i - 1, j), idx_H, 0.5);
      P_trp.emplace_back(idxh(i, j + 1), idx_H, 0.5);
      P_trp.emplace_back(idxh(i, j - 1), idx_H, 0.5);
      if (bilinear) {
      P_trp.emplace_back(idxh(i + 1, j + 1), idx_H, 0.25);
      P_trp.emplace_back(idxh(i + 1, j - 1), idx_H, 0.25);
      P_trp.emplace_back(idxh(i - 1, j + 1), idx_H, 0.25);
      P_trp.emplace_back(idxh(i - 1, j - 1), idx_H, 0.25);
      }
      else {
      P_trp.emplace_back(idxh(i + 1, j + 1), idx_H, 0.5);
      P_trp.emplace_back(idxh(i - 1, j - 1), idx_H, 0.5);
      }
    }
  }
  // create the sparse matrix in CCS format
  Eigen::SparseMatrix<double> P(N, n);
  P.setFromTriplets(P_trp.begin(), P_trp.end());
  return P;
}

Eigen::SparseMatrix<double> buildAH(const Eigen::SparseMatrix<double> &A,
                                    const Eigen::SparseMatrix<double> &P) {
  Eigen::SparseMatrix<double> AH;
  AH = P.transpose() * A * P;
  return AH;
}

}  // namespace GalerkinConstruction
