/**
 * @ file kernmatllrapprox_main.cpp
 * @ brief NPDE homework XXX MAIN FILE
 * @ author R. Hiptmair
 * @ date September 2023
 * @ copyright Developed at SAM, ETH Zurich
 */

#include "kernmatllrapprox.h"

int main(int /*argc*/, char** /*argv*/) {
  std::cout
      << "ADVNCSE HW KernMatLLRApprox: Compressing 1D kernel collocation matrix"
      << " with bi-directional Chebychev interpolation " << std::endl;
  // Build cluster tree
  /* SAM_LISTING_BEGIN_1 */
  const int q = 5;      // Degree of Chebychev interpolation
  const int npts = 16;  // Number of collocation points
  std::vector<HMAT::Point<1>> pts;
  for (int n = 0; n < npts; n++) {
    HMAT::Point<1> p;
    p.idx = n;
    p.x[0] = static_cast<double>(n) / (npts - 1);
    pts.push_back(p);
  }
  // Allocate cluster tree object (the same for both directions)
  auto T_row = std::make_shared<
      KernMatLLRApprox::LLRClusterTree<KernMatLLRApprox::InterpNode<1>>>(q);
  T_row->init(pts);
  auto T_col = std::make_shared<
      KernMatLLRApprox::LLRClusterTree<KernMatLLRApprox::InterpNode<1>>>(q);
  T_col->init(pts);
  // Kernel for 1D collocation
  struct Kernel {
    Kernel() = default;
    double operator()(double x, double y) const {
      return ((x != y) ? -std::log(std::abs(x - y)) : 0.0);
    }
  } G;
  KernMatLLRApprox::BiDirChebPartMat1D<Kernel> Mt(T_row, T_col, G, q, 2.0);
  /* SAM_LISTING_END_1 */

  // // Output row cluster tree (the column tree is the same)
  // std::cout << "CLUSTER TREE" << std::endl << *(Mt.rowT->root) << std::endl;
  // // Output block partition
  // std::cout << Mt;

  // std::cout << "Computing sparsity measure" << std::endl;
  // unsigned int spm = KernMatLLRApprox::computeSparsityMeasure(Mt,
  // &std::cout); std::cout << "Sparsity measure = " << spm << std::endl;

  std::cout << "Matrix x vector" << std::endl;
  Eigen::VectorXd x = Eigen::VectorXd::Constant(16, 1.0);
  Eigen::VectorXd y = KernMatLLRApprox::mvLLRPartMat(Mt, x);
  std::cout << y << std::endl;
  // KernMatLLRApprox::tabulateConvergenceLLR(
  //     {10, 20, 40, 60, 80, 160}, {3, 4, 5, 6});

  // KernMatLLRApprox::runtimeMatVec(
  //     {1024, 2048, 4096, 8192, 16384, 32768, 65536, 131072, 262144});
  return 0;
}
