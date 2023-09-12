/**
 * @ file kernmatllrapprox_main.cpp
 * @ brief NPDE homework XXX MAIN FILE
 * @ author R. Hiptmair
 * @ date September 2023
 * @ copyright Developed at SAM, ETH Zurich
 */

#include "kernmatllrapprox.h"

int main(int /*argc*/, char** /*argv*/) {
  std::cout << "ADVNCSE HW KernMatLLRApprox: Compressing 1D kernel collocation matrix"
            << " with bi-directional Chebychev interpolation " << std::endl;
  // Build cluster tree
  /* SAM_LISTING_BEGIN_1 */
  const int q = 5;      // Degree of Chebychev interpolation
  const int npts = 64;  // Number of collocation points
  std::vector<HMAT::Point<1>> pts;
  for (int n = 0; n < npts; n++) {
    HMAT::Point<1> p;
    p.idx = n;
    p.x[0] = static_cast<double>(n) / (npts - 1);
    pts.push_back(p);
  }
  // Allocate cluster tree object (the same for both directions)
  auto T = std::make_shared<KernMatLLRApprox::LLRClusterTree<KernMatLLRApprox::InterpNode<1>>>(q);
  T->init(pts);

  struct LogKernel {
    LogKernel() = default;
    double operator()(double x, double y) {
      return ((x != y) ? -std::log(std::abs(x - y)) : 0.0);
    }
  } G;
  KernMatLLRApprox::BiDirChebPartMat1D<LogKernel> Mt(T, T, G, q);
  /* SAM_LISTING_END_1 */

  KernMatLLRApprox::tabulateConvergenceLLR({10, 20, 40, 60, 80, 160, 320, 640, 1280},
                               {3, 4, 5, 6, 7});
  KernMatLLRApprox::runtimeMatVec(
      {1024, 2048, 4096, 8192, 16384, 32768, 65536, 131072, 262144});
  return 0;
}
