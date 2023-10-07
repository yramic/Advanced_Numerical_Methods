/**
 * @file kernmatllrapprox_test.cc
 * @brief NPDE homework KernMatLLRApprox test code
 * @author
 * @date
 * @copyright Developed at SAM, ETH Zurich
 */

#include "../kernmatllrapprox.h"

#include <gtest/gtest.h>

#include <memory>

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

template <typename TRANSFORM = std::function<double(double)>>
std::shared_ptr<KernMatLLRApprox::LLRClusterTree<HMAT::CtNode<1>>>
make1DClusterTree(
    unsigned int q, unsigned int npts,
    TRANSFORM trf = [](double x) -> double { return x; }) {
  // Build cluster tree
  std::vector<HMAT::Point<1>> pts;
  for (int n = 0; n < npts; n++) {
    HMAT::Point<1> p;
    p.idx = n;
    p.x[0] = trf(static_cast<double>(n) / (npts - 1));
    pts.push_back(p);
  }
  auto T =
      std::make_shared<KernMatLLRApprox::LLRClusterTree<HMAT::CtNode<1>>>(q);
  T->init(pts);
  return T;
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
  HMAT::BlockPartition<HMAT::ClusterTree<HMAT::CtNode<1>>> bP(T, T);
  bP.init(eta);  // Admissibility parameter eta = 2.00
  std::cout << "TEST: check_matrixpartition: Call test" << std::endl;
  ASSERT_TRUE(KernMatLLRApprox::checkMatrixPartition(bP));
}

TEST(KernMatLLRApprox, check_V) {
  // Test whether row sums of interpolation matrices are equal to 1
  auto T = make1DClusterTree(5, 16);
  std::function<bool(const HMAT::CtNode<1> *)> rec_check =
      [&](const HMAT::CtNode<1> *node) -> bool {
    if (node) {
      const Eigen::MatrixXd &V{T->Vs[node->nodeNumber]};
      const Eigen::VectorXd row_sums = V.rowwise().sum();
      double dev = (row_sums - Eigen::VectorXd::Constant(V.rows(), 1.0)).norm();
      if (dev > 1.0E-10) {
        std::cout << "V with row sum different from 1" << std::endl
                  << V << std::endl
                  << "row sums = " << row_sums.transpose() << std::endl;
        ;
        return false;
      }
      return (rec_check(node->sons[0]) and rec_check(node->sons[1]));
    }
    return true;
  };
  ASSERT_TRUE(rec_check(T->root));
}

TEST(KernMatLLRApprox, validateLLR) {
  ASSERT_TRUE(KernMatLLRApprox::validateLLR(5));
}

}  // namespace KernMatLLRApprox::test
