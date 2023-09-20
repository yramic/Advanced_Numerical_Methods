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

TEST(GravitationalForces, MinDist) {
  // tests the minimum distance condition
  const int n = 10;
  const double mindist = 1. / (2. * std::sqrt(n));
  const auto pos = initStarPositions(n, mindist);

  for (auto p1 = pos.begin(); p1 != pos.end(); ++p1) {
    for (auto p2 = p1 + 1; p2 != pos.end(); ++p2) {
      ASSERT_GE((*p1 - *p2).norm(), mindist);
    }
  }
}

/////////////////////////
// assumes up then right ordering of son nodes though other ordering are possible
/////////////////////////
TEST(GravitationalForces, QuadTree) {
  // Simple test for the QuadTree structure
  const int n = 4;
  std::vector<Eigen::Vector2d> pos(4, Eigen::Vector2d::Zero());
  const std::vector<double> mass(4, 1.0);
  pos[0] << 0.25, 0.25;
  pos[1] << 0.25, 0.75;
  pos[2] << 0.3, 0.8;
  pos[3] << 0.75, 0.6;
  const GravitationalForces::StarQuadTreeClustering clustering(pos, mass);

  ASSERT_EQ(clustering.root_->star_idx_.size(), 4);
  const auto& sons = clustering.root_->sons_;
  ASSERT_EQ(sons[0]->star_idx_.size(), 1);
  ASSERT_EQ(sons[1]->star_idx_.size(), 2);
  ASSERT_EQ(sons[2], nullptr);
  ASSERT_EQ(sons[3]->star_idx_.size(), 1);

  ASSERT_EQ(clustering.no_leaves_, n);
}

TEST(GravitationalForces, isAdmissible) {
  // tests whether the isAdmissible function works properly
  const int n = 1;
  const auto pos = GravitationalForces::initStarPositions(n);
  const std::vector<double> mass(n, 1.0);
  const GravitationalForces::StarQuadTreeClustering clustering(pos, mass);

  const Eigen::Vector2d center{-0.5, 0.5};
  ASSERT_TRUE(clustering.isAdmissible(*clustering.root_, center, 0.2));
  ASSERT_FALSE(clustering.isAdmissible(*clustering.root_, center, 0.6));
}

TEST(GravitationalForces, forceOnStar) {
  // Simple test case for the computation of approximate forces
  std::vector<Eigen::Vector2d> pos(4, Eigen::Vector2d::Zero());
  const std::vector<double> mass(4, 1.0);
  pos[0] << 0.0, 0.0;
  pos[1] << 0.0, 1.0;
  pos[2] << 1.0, 0.0;
  pos[3] << 1.0, 1.0;
  const GravitationalForces::StarQuadTreeClustering clustering(pos, mass);
  const double eta = 1.;

  Eigen::MatrixXd f(2, 4);
  f << 0.107712, 0.107712, -0.107712, -0.107712,
      0.107712, -0.107712, 0.107712, -0.107712;
  for (int i = 0; i < 4; ++i) {
    ASSERT_NEAR((clustering.forceOnStar(i, eta)-f.col(i)).norm(),0.0,1.0E-4);
  }
}


} // namespace GravitationalForces::test
