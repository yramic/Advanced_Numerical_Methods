/***********************************************************************
 *                                                                     *
 * Code for Course "Advanced Numerical Methods for CSE"                *
 * (Prof. Dr. R. Hiptmair)                                             *
 * Author: Ralf Hiptmair * Date: August 2023 (C) Seminar for Applied
 *Mathematics, ETH Zurich                     * This code can be freely used for
 *non-commercial purposes as long    * as this header is left intact. *
 ***********************************************************************/

#include <cstddef>

#include "clustertree.h"

// Auxiliary function for outputting a cluster tree
template <int DIM>
void printClusterTree(const HMAT::ClusterTree<HMAT::CtNode<DIM>> &T,
                      const std::unique_ptr<HMAT::CtNode<DIM>> &node_ptr,
                      int level = 0) {
  if (node_ptr != nullptr) {
    std::cout << "Node level " << level << ": " << node_ptr->noIdx()
              << " points, " << T.getBBox(node_ptr.get()) << ", idx = [";
    const int n_sons = ((node_ptr->sons[0] != nullptr) ? 1 : 0) +
                       ((node_ptr->sons[1] != nullptr) ? 1 : 0);
    const std::vector<std::size_t> node_idxs{node_ptr->I};
    for (std::size_t idx : node_idxs) {
      std::cout << idx << ' ';
    }
    std::cout << "]: " << n_sons << " sons" << std::endl;
    if (node_ptr->sons[0] != nullptr) {
      printClusterTree(T, node_ptr->sons[0], level + 1);
    }
    if (node_ptr->sons[1] != nullptr) {
      printClusterTree(T, node_ptr->sons[1], level + 1);
    }
  }
}

/** Main program */
int main(int /*argc*/, char ** /*argv*/) {
  std::cout << "Construction of cluster tree" << std::endl;
  {
    std::cout << "\t >> 1D:" << std::endl;
    const int d = 1;
    const std::size_t npts = 16;
    std::vector<HMAT::Point<d>> pts;
    for (int n = 0; n < npts; n++) {
      HMAT::Point<d> p;
      p.idx = n;
      p.x[0] = static_cast<double>(n) / (npts - 1);
      pts.push_back(p);
    }
    HMAT::ClusterTree<HMAT::CtNode<d>> T;
    T.init(pts);
    std::cout << "1D Cluster tree" << T << std::endl;
    printClusterTree(T, T.root);
  }
  {
    std::cout << "\t >> 2D:" << std::endl;
    std::vector<HMAT::Point<2>> pts;
    const std::size_t npts_1d = 8;
    for (int i = 0; i < npts_1d; ++i) {
      for (int j = 0; j < npts_1d; ++j) {
        HMAT::Point<2> p;
        p.idx = i * npts_1d + j;
        const double pos_x = static_cast<double>(i) / (npts_1d - 1);
        const double pos_y = static_cast<double>(j) / (npts_1d - 1);
        p.x = Eigen::Vector2d(pos_x * pos_x, pos_y * pos_y);
        pts.push_back(p);
      }
    }
    HMAT::ClusterTree<HMAT::CtNode<2>> T;
    T.init(pts);
    std::cout << "2D Cluster tree" << T << std::endl;
    printClusterTree(T, T.root);
  }
  exit(0);
}
