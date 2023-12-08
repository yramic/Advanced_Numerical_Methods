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

/* SAM_LISTING_BEGIN_2 */
Eigen::SparseMatrix<double, Eigen::RowMajor> prolongationMatrix(unsigned int M,
                                                                bool bilinear) {
  assertm(((M > 3) and (M % 2 == 0)), "prolongationMatrix: M must be even!");
  const unsigned int N = (M - 1) * (M - 1);  // No of inter nodes of fine grid
  const unsigned int m =
      M / 2;  // Number of cells in each direction of coarse grid
  const unsigned int n =
      (m - 1) * (m - 1);  // No of interior nodes of fine grid
  // Sparse matrix in CCS format
  Eigen::SparseMatrix<double, Eigen::RowMajor> P(N, n);
  // **********************************************************************
  // To be supplemented
  // **********************************************************************
  return P;
}
/* SAM_LISTING_END_2 */

/* SAM_LISTING_BEGIN_3 */
Eigen::SparseMatrix<double> buildAH_eigen(
    const Eigen::SparseMatrix<double> &A,
    const Eigen::SparseMatrix<double, Eigen::RowMajor> &P) {
  Eigen::SparseMatrix<double> AH;
  AH = P.transpose() * A * P;
  return AH;
}
/* SAM_LISTING_END_3 */

/* SAM_LISTING_BEGIN_4 */
Eigen::SparseMatrix<double> buildAH(
    const Eigen::SparseMatrix<double> &A,
    const Eigen::SparseMatrix<double, Eigen::RowMajor> &P) {
  Eigen::SparseMatrix<double> AH(P.cols(), P.cols());
// **********************************************************************
// To be supplemented
// **********************************************************************
  return AH;
}
/* SAM_LISTING_END_4 */

}  // namespace GalerkinConstruction
