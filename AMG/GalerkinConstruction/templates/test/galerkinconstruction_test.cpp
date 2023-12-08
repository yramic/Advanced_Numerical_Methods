/**
 * @file galerkinconstruction_test.cc
 * @brief NPDE homework GalerkinConstruction test code
 * @author
 * @date
 * @copyright Developed at SAM, ETH Zurich
 */

#include "../galerkinconstruction.h"

#include <gtest/gtest.h>

namespace GalerkinConstruction::test {
TEST(GalConst, test_PTAP) {
  /* SAM_LISTING_BEGIN_1 */
  const unsigned int N = 8;
  const unsigned int n = N / 2;

  const Eigen::SparseMatrix<double> Ah = poissonMatrix(N);
  const Eigen::SparseMatrix<double> AH = poissonMatrix(n);
  const auto P_L = prolongationMatrix(N, false);
  const auto P_B = prolongationMatrix(N, true);
  const Eigen::SparseMatrix<double> AH_L = buildAH(Ah, P_L);
  const Eigen::SparseMatrix<double> AH_B = buildAH(Ah, P_B);

  const double diff_L = (AH_L - AH).norm();
  const double diff_B = (AH_B - AH).norm();
  std::cout << "(AH_L - AH).norm() = " << diff_L
            << ",\n(AH_B - AH).norm() = " << diff_B << std::endl;
  /* SAM_LISTING_END_1 */

  ASSERT_NEAR(diff_L, 0.0, 1.0E-6);
}

}  // namespace GalerkinConstruction::test
