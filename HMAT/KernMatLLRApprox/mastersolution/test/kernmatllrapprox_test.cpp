/**
 * @file kernmatllrapprox_test.cc
 * @brief NPDE homework KernMatLLRApprox test code
 * @author
 * @date
 * @copyright Developed at SAM, ETH Zurich
 */

#include "../kernmatllrapprox.h"

#include <gtest/gtest.h>

namespace KernMatLLRApprox::test {

// Creating of 2D points for testing
// Copied from HMAT/CLUSTERING/clustertree_main.cpp
std::vector<HMAT::Point<2>> initPoints(unsigned int npts) {
  std::vector<HMAT::Point<2>> pts;
  for (int i = 0; i < npts; ++i) {
    for (int j = 0; j < npts; ++j) {
      HMAT::Point<2> p;
      p.idx = i * npts + j;
      const double pos_x = static_cast<double>(i) / (npts - 1);
      const double pos_y = static_cast<double>(j) / (npts - 1);
      p.x = Eigen::Vector2d(pos_x * pos_x, pos_y * pos_y);
      pts.push_back(p);
    }
  }
  return pts;
}

TEST(KernMatLLRApprox, check_clustertree) {
  std::cout << "TEST: check_clustertree" << std::endl;
  // Building a cluster tree for testing
  HMAT::ClusterTree<HMAT::CtNode<2>> T;
  T.init(initPoints(5));
  std::cout << "TEST: check_clustertree: tree built" << std::endl;
  ASSERT_TRUE(KernMatLLRApprox::checkClusterTree(T));
}

TEST(KernMatLLRApprox, check_matrixpartition) {
  std::cout << "TEST: check_matrixpartition" << std::endl;
  // Build a partitioned matrix, copied from
  // HMAT/CLUSTERING/matrixpartition_main.cpp
  const unsigned int npts = 16;
  const double eta = 2;

  // Build cluster tree
  std::vector<HMAT::Point<1>> pts;
  for (int n = 0; n < npts; n++) {
    HMAT::Point<1> p;
    p.idx = n;
    p.x[0] = static_cast<double>(n) / (npts - 1);
    pts.push_back(p);
  }
  auto T = std::make_shared<HMAT::ClusterTree<HMAT::CtNode<1>>>();
  T->init(pts);
  // Build block cluster tree
  HMAT::BlockPartition<HMAT::CtNode<1>, HMAT::IndexBlock<HMAT::CtNode<1>>,
                       HMAT::IndexBlock<HMAT::CtNode<1>>>
      bP(T, T);
  bP.init(eta);  // Admissibility parameter eta = 2.00
  std::cout << "TEST: check_matrixpartition: Call test" << std::endl;
  ASSERT_TRUE(KernMatLLRApprox::checkMatrixPartition(bP));
}

}  // namespace KernMatLLRApprox::test
