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

// Test ltpMult()
TEST(LowTriangToeplitz, test_ltpMult) {
  // Test size
  int size = 16;
  // Random initialization
  Eigen::VectorXcd f = Eigen::VectorXcd::Random(size);
  Eigen::VectorXcd g = Eigen::VectorXcd::Random(size);
  // Result vector
  Eigen::VectorXcd rv = ltpMult(f, g);
  // Build Toeplitz matrices from sequences
  Eigen::MatrixXcd Mf, Mg, gt, rm;
  Mf.setZero(size, size);
  Mg.setZero(size, size);
  rm.setZero(size, size);
  for (int i = 0; i < size; i++) {
    Mf.diagonal(-i) = f(i) * Eigen::VectorXcd::Ones(size - i);
    Mg.diagonal(-i) = g(i) * Eigen::VectorXcd::Ones(size - i);
    rm.diagonal(-i) = rv(i) * Eigen::VectorXcd::Ones(size - i);
  }
  // Ground truth
  gt = Mf * Mg;

  ASSERT_NEAR((gt - rm).real().norm(), 0.0, 1e-8);
}

// Test ltpSolve()
TEST(LowTriangToeplitz, test_ltpSolve) {
  // Test size
  int size = 16;
  // Random initialization
  Eigen::VectorXcd f = Eigen::VectorXcd::Random(size);
  Eigen::VectorXcd y = Eigen::VectorXcd::Random(size);
  // Result vector
  Eigen::VectorXcd x = ltpSolve(f, y);
  // Build Toeplitz matrices from sequences
  Eigen::MatrixXcd Mf;
  Mf.setZero(size, size);
  for (int i = 0; i < size; i++) {
    Mf.diagonal(-i) = f(i) * Eigen::VectorXcd::Ones(size - i);
  }
  // Ground truth
  Eigen::VectorXcd gt = Mf.colPivHouseholderQr().solve(y);

  ASSERT_NEAR((gt - x).real().norm(), 0.0, 1e-8);
}

}  // namespace LowTriangToeplitz::test
