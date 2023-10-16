/***********************************************************************
 *                                                                     *
 * Code for Course "Advanced Numerical Methods for CSE"                *
 * (Prof. Dr. R. Hiptmair)                                             *
 * Author: Ralf Hiptmair                                               *
 * Date: Augsut 2023
 * (C) Seminar for Applied Mathematics, ETH Zurich                     *
 * This code can be freely used for non-commercial purposes as long    *
 * as this header is left intact.                                      *
 ***********************************************************************/

#include <memory>

#include "matrixpartition.h"

// Auxiliary function for building a matrix partitioning
template <typename TRANSFORM = std::function<double(double)>>
HMAT::BlockPartition<HMAT::ClusterTree<HMAT::CtNode<1>>> matrixPartition1d(
    unsigned int npts, double eta = 2.0,
    TRANSFORM trf = [](double x) -> double { return x; }) {
  // Build cluster tree
  std::vector<HMAT::Point<1>> pts;
  for (int n = 0; n < npts; n++) {
    HMAT::Point<1> p;
    p.idx = n;
    p.x[0] = trf(static_cast<double>(n) / (npts - 1));
    pts.push_back(p);
  }
  auto T = std::make_shared<HMAT::ClusterTree<HMAT::CtNode<1>>>(pts);

  HMAT::BlockPartition<HMAT::ClusterTree<HMAT::CtNode<1>>> bP(T, T);
  bP.init(eta);  // Admissibility parameter eta = 0.5
  return bP;
}

/** Main program */
int main(int /*argc*/, char** /*argv*/) {
  std::cout << "Construction of block cluster tree" << std::endl;
  // Construct 1D cluster tree based on equidistant points
  auto matpart = matrixPartition1d(16);
  // Output cluster tree
  std::cout << matpart << "Part I: Normal termination" << std::endl;
  (void)HMAT::printGeometricPartition(matpart);
  exit(0);
}
