/**
 * @file tensorproductchebintp_test.cc
 * @brief ADVNCSE homework TensorProductChebIntp test code
 * @author Bob Schreiner
 * @date August 2023
 * @copyright Developed at SAM, ETH Zurich
 */

#include "../lowtriangtoeplitz.h"

#include <gtest/gtest.h>

namespace LowTriangToeplitz::test {

Eigen::VectorXcd ltpMultold(const Eigen::VectorXcd& f,
                            const Eigen::VectorXcd& g) {
  assert(f.size() == g.size() && "f and g vectors must have the same length!");

  const std::size_t n = f.size();
  return toepMatVecMult(f, Eigen::VectorXcd::Zero(n), g);
}

TEST(LowTriangToeplitz, ltpMult) {
  const Eigen::VectorXcd u = Eigen::VectorXcd::Random(10);
  const Eigen::VectorXcd v = Eigen::VectorXcd::Random(10);
  ASSERT_NEAR((ltpMultold(u, v) - ltpMult(u, v)).real().norm(), 0.0, 1e-8);
}

TEST(LowTriangToeplitz, ltpSolve) {
  const Eigen::VectorXcd u = Eigen::VectorXcd::Random(16);
  const Eigen::VectorXcd v = Eigen::VectorXcd::Random(16);
  Eigen::VectorXcd zero(16);
  zero[0] = u[0];
  zero.tail(15) = Eigen::VectorXcd::Zero(15);
  const Eigen::MatrixXcd K = toeplitz(u, zero);
  const Eigen::VectorXcd res = K.triangularView<Eigen::Lower>().solve(v);
  const Eigen::VectorXcd stud_res = ltpSolve(u, v);

  ASSERT_NEAR((res - stud_res).norm(), 0.0, 1e-8);
}

}  // namespace LowTriangToeplitz::test
