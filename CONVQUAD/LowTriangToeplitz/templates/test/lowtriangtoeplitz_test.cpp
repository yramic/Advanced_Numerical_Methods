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

TEST(LowTriangToeplitz, toeplitz) {
    Eigen::VectorXcd u(3);
    u << 1,2,3;


     ASSERT_NEAR((ltpMultold(u,u)-ltpMult(u,u)).real().norm(),0.0, 1e-8); }

}  // namespace LowTriangToeplitz::test
