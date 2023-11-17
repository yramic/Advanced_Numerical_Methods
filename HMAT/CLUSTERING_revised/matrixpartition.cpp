/***********************************************************************
 *                                                                     *
 * Code for Course "Advanced Numerical Methods for CSE"                *
 * (Prof. Dr. R. Hiptmair)                                             *
 * Author: R.H.                                                        *
 * Date: Nov 18, 2017                                                  *
 * (C) Seminar for Applied Mathematics, ETH Zurich                     *
 * This code can be freely used for non-commercial purposes as long    *
 * as this header is left intact.                                      *
 ***********************************************************************/

#include "matrixpartition.h"

namespace HMAT {
// Auxiliary function for 1D case:
// output of geometric partition of unit square
void printGeometricPartition(
    const BlockPartition<ClusterTree<CtNode<1>>> &partmat) {
  std::cout << "NF = [";
  for (const IndexBlock<CtNode<1>> &b : partmat.nearField) {
    const BBox<1> bbx{partmat.rowT->getBBox(&(b.nx))};
    const BBox<1> bby{partmat.colT->getBBox(&(b.ny))};
    std::cout << bbx.minc << " , " << bbx.maxc << " ," << bby.minc << " , "
              << bby.maxc << "; ";
  }
  std::cout << "]" << std::endl;
  std::cout << "FF = [";
  for (const IndexBlock<CtNode<1>> &b : partmat.farField) {
    const BBox<1> bbx{partmat.rowT->getBBox(&(b.nx))};
    const BBox<1> bby{partmat.colT->getBBox(&(b.ny))};
    std::cout << bbx.minc << " , " << bbx.maxc << " , " << bby.minc << " , "
              << bby.maxc << "; ";
  }
  std::cout << "]" << std::endl;
}
}  // namespace HMAT
