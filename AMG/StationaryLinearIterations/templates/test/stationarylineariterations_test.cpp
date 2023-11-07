/**
 * @file stationarylineariterations_test.cc
 * @brief ADVNCSE homework StationaryLinearIterations test code
 * @author Bob Schreiner
 * @date October 2023
 * @copyright Developed at SAM, ETH Zurich
 */

#include "../stationarylineariterations.h"

#include <gtest/gtest.h>

namespace StationaryLinearIterations::test {
TEST(StationaryLinearIterations, test_poissonMatrix) {
  const int n = 20;
  const Eigen::VectorXd c = Eigen::VectorXd::Constant((n - 1) * (n - 1), 1);
  const Eigen::SparseMatrix<double> M = poissonMatrix(n);
  Eigen::VectorXd res = (M * c).segment(n - 1, (n - 1) * (n - 3));
  res[0] -= 1;
  for (int i = 1; i < n - 3; ++i) {
    res[i * (n - 1)] -= 1;
    res[i * (n - 1) - 1] -= 1;
  }
  res[(n - 1) * (n - 3) - 1] -= 1;
  ASSERT_NEAR((res).norm(), 0, 1e-5);
}
}  // namespace StationaryLinearIterations::test
