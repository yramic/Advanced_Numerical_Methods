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

#include "clustertree.h"

namespace HMAT {
// distance of 1D intervals \cob{$\cintv{a,b}$} and \cob{$\cintv{c,d}$}
double dist(double a, double b, double c, double d) {
  if (b < a) {
    std::swap(a, b);
  }
  if (d < c) {
    std::swap(c, d);
  }
  if (c < a) {
    std::swap(a, c);
    std::swap(b, d);
  }
  return (c < b) ? 0.0 : c - b;
}

// Non-recursive output operator for 1D cluster, for debugging
void outCTNode(const CtNode<1> &ctnode, std::ostream &o) {
  o << "CTNode: ";
  for (const HMAT::Point<1> pt : ctnode.pts) {
    o << "(" << pt.idx << ", " << pt.x[0] << ") ";
  }
  o << std::endl;
}

}  // namespace HMAT
