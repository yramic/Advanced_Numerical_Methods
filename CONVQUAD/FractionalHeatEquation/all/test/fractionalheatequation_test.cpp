/**
 * @file fractionalheatequation_test.cc
 * @brief NPDE homework FractionalHeatEquation test code
 * @author Dr. Jörg Nick
 * @date October 2023
 * @copyright Developed at SAM, ETH Zurich
 */

#include "../fractionalheatequation.h"

#include <gtest/gtest.h>

namespace FractionalHeatEquation::test {

TEST(FractionalHeatEquation, SqrtsMPlusA) {
  const int n = 3;
  FractionalHeatEquation::SqrtsMplusA mat(n, std::complex<double>(0, 1.0));
  Eigen::VectorXd v(n * n);
  v.setZero();
  Eigen::MatrixXcd mat_d(n * n, n * n);
  for (int i = 0; i < n * n; ++i) {
    if (i > 0) {
      v[i - 1] = 0.0;
    }
    v[i] = 1.0;
    mat_d.col(i) = mat.eval(v);
  }
  std::cout << "n = " << n << ", matrix = " << std::endl << mat_d << std::endl;

  for (int i = 0; i < 10; ++i) {
    v = Eigen::VectorXd::Random(n * n);
    EXPECT_NEAR((mat.eval(mat.solve(v)) - v).norm(), 0.0, 1.0E-10);
  }
  std::cout << FractionalHeatEquation::SqrtsMplusA::solve_cnt << " solves, "
            << FractionalHeatEquation::SqrtsMplusA::ludec_cnt
            << " lu decompositions" << std::endl;
}

}  // namespace FractionalHeatEquation::test
