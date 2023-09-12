/**
* @file gravitationalforces_test.cpp
* @brief NPDE homework GravitationalForces test code
* @author
* @date
* @copyright Developed at SAM, ETH Zurich
*/


#include <gtest/gtest.h>
#include "../gravitationalforces.h"

namespace GravitationalForces::test {
  TEST(GravitationalForces,FourStars) {
    // Simple test case for the computation of exact forces
    std::vector<Eigen::Vector2d> pos(4, Eigen::Vector2d::Zero());
    const std::vector<double> mass(4, 1.0);
    pos[0] << 0.0, 0.0;
    pos[1] << 0.0, 1.0;
    pos[2] << 1.0, 0.0;
    pos[3] << 1.0, 1.0;
    std::vector<Eigen::Vector2d> forces =
        GravitationalForces::computeForces_direct(pos, mass);
    ASSERT_EQ(forces.size(),4);
    Eigen::Vector2d f0(0.107712, 0.107712);
    Eigen::Vector2d f1(0.107712, -0.107712);
    Eigen::Vector2d f2(-0.107712,  0.107712);
    Eigen::Vector2d f3(-0.107712, -0.107712);
    ASSERT_NEAR((forces[0]-f0).norm(),0.0,1.0E-4);
    ASSERT_NEAR((forces[1]-f1).norm(),0.0,1.0E-4);
    ASSERT_NEAR((forces[2]-f2).norm(),0.0,1.0E-4);
    ASSERT_NEAR((forces[3]-f3).norm(),0.0,1.0E-4);
  }
} // namespace GravitationalForces::test
