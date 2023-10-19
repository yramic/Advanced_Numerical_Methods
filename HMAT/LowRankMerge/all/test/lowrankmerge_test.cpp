/**
* @file lowrankmerge_test.cpp
* @brief NPDE homework LowRankMerge test code
* @author
* @date
* @copyright Developed at SAM, ETH Zurich
*/

#include "../lowrankmerge.h"

#include <gtest/gtest.h>

namespace LowRankMerge::test {
TEST(LowRankMerge, low_rank_merge) {
  // Simple test case for low_rank_merge(): If A1 = A2, then Atilde = A1 = A2, Btilde = [B1 B2]
  Eigen::MatrixXd A1(3, 3), B1(3, 3), A2(3, 3), B2(3, 3), X1(3, 3), X2(3, 3);
  A1 << 1, 0, 0, 0, 0, 1, 0, 1, 0;
  B1 << 1, 0, 0, 2, 3, 0, 4, 6, 5;
  X1 << 1, 2, 4, 0, 0, 5, 0, 3, 6;
  A2 = A1;
  B2 = B1;
  X2 = X1;

  std::pair<Eigen::MatrixXd, Eigen::MatrixXd> AB =
      LowRankMerge::low_rank_merge(A1, B1, A2, B2);
  Eigen::MatrixXd Ztilde = AB.first * AB.second.transpose();
  ASSERT_EQ(Ztilde.rows(), 3);
  ASSERT_EQ(Ztilde.cols(), 6);
  Eigen::MatrixXd Z(3, 6);
  Z << X1, X2;
  ASSERT_NEAR((Z - Ztilde).norm(), 0.0, 1.0E-4);
}
}  // namespace LowRankMerge::test
