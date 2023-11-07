/**
 * @file fractionalheatequation_test.cc
 * @brief NPDE homework FractionalHeatEquation test code
 * @author Jörg Nick, Bob Schreiner
 * @date October 2023
 * @copyright Developed at SAM, ETH Zurich
 */

#include "../fractionalheatequation.h"

#include <gtest/gtest.h>

#include "../../../AbsorbingBoundaryCondition/mastersolution/absorbingboundarycondition.h"

namespace FractionalHeatEquation::test {

TEST(FractionalHeatEquation, generateGrid) {
  std::vector<Eigen::Vector2d> gridpoints = generateGrid(2);
  ASSERT_NEAR(gridpoints[0][0], 1.0 / (3), 1E-12);
}

TEST(FractionalHeatEquation, cqWeights) {
  size_t M = 5;
  double tau = 0.1;
  Eigen::VectorXd w_ex = cqWeights(M, tau);
  auto F = [](std::complex<double> s) { return std::pow(s, 0.5); };
  auto delta = [](std::complex<double> z) {
    return 1.0 - z;
    //return 1.0 / 2.0 * z * z - 2.0 * z + 3.0 / 2.0;
  };
  Eigen::VectorXd w_tr = AbsorbingBoundaryCondition::cqweights_by_dft(F, delta, tau, M);
  //AbsorbingBoundaryCondition::cqweights_by_dft(F, delta, tau, M);
  std::cout << w_tr - w_ex << std::endl;
  ASSERT_NEAR((w_tr - w_ex).norm(), 0, 1E-8);
}

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

TEST(FHE, Toeplitzop) {
  Eigen::VectorXd v = (Eigen::VectorXd(7) << 1, 2, 3, 4, 5, 6, 7).finished();
  ToeplitzOp T(v);
  Eigen::MatrixXd T_mat(4, 4);
  T_mat << 4, 3, 2, 1, 5, 4, 3, 2, 6, 5, 4, 3, 7, 6, 5, 4;
  Eigen::VectorXd x = (Eigen::VectorXd(4) << 0.5, 2.5, 3.5, 7.5).finished();
  std::cout << "Teoplitz matrix = " << std::endl << T_mat << std::endl;
  std::cout << "T_mat * x = " << (T_mat * x).transpose() << std::endl;
  std::cout << "T.eval(x) = " << T.eval(x).transpose() << std::endl;
  EXPECT_NEAR((T_mat * x - T.eval(x)).norm(), 0.0, 10E-10);
}

}  // namespace FractionalHeatEquation::test
