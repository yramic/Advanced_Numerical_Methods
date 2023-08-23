/***********************************************************************
 *                                                                     *
 * Code for Course "Advanced Numerical Methods for CSE"                *
 * (Prof. Dr. R. Hiptmair)                                             *
 * Author: Ralf Hiptmair                                               *
 * Date: August 2023                                                   *
 * (C) Seminar for Applied Mathematics, ETH Zurich                     *
 * This code can be freely used for non-commercial purposes as long    *
 * as this header is left intact.                                      *
 ***********************************************************************/

#include "locallowrank.h"

/** Main program */
int main(int /*argc*/, char** /*argv*/) {
  std::cout << "Compressing 1D kernel collocation matrix"
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
  auto T = std::make_shared<HMAT::LLRClusterTree<HMAT::InterpNode<1>>>(q);
  T->init(pts);

  struct LogKernel {
    LogKernel() = default;
    double operator()(double x, double y) {
      return ((x != y) ? -std::log(std::abs(x - y)) : 0.0);
    }
  } G;
  HMAT::BiDirChebPartMat1D<LogKernel> Mt(T, T, G, q);
  /* SAM_LISTING_END_1 */
  exit(0);
}
