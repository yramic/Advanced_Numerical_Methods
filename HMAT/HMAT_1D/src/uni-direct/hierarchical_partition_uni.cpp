/***********************************************************************
 *                                                                     *
 * Code for Course "Advanced Numerical Methods for CSE"                *
 * (Prof. Dr. R. Hiptmair)                                             *
 * Author: Daniele Casati                                              *
 * Date: 11/2017                                                       *
 * (C) Seminar for Applied Mathematics, ETH Zurich                     *
 * This code can be freely used for non-commercial purposes as long    *
 * as this header is left intact.                                      *
 ***********************************************************************/
#include <Eigen/Dense>
#include <algorithm>
#include <iostream>
#include <vector>

#include "../../include/block_nearf.hpp"
#include "../../include/cheby.hpp"
#include "../../include/ctree.hpp"
#include "../../include/hierarchical_partition.hpp"
#include "../../include/is_admissible.hpp"
#include "../../include/point.hpp"
#include "../../include/uni-direct/block_cluster_Y.hpp"
#include "../../include/uni-direct/node_Y.hpp"

// add pairs of pointers to Node of the Cluster Tree in the Near and Far Field
// Vectors of pairs
template <>
void HierarchicalPartitioning<BlockCluster_Y, Node_Y>::setNearFar_recursion(
    Node* xnode, Node* ynode, double eta) {
  // if *xnode or *ynode is a leaf, we add the pair (*xnode,*ynode) to the near
  // field vector
  if ((*xnode).getLChild() == NULL || (*ynode).getRChild() == NULL) {
    NearField_.push_back(new BlockNearF(xnode, ynode));
  } else {
    double xl = (*xnode).getPoints().front().getX();
    double xr = (*xnode).getPoints().back().getX();
    double yl = (*ynode).getPoints().front().getX();
    double yr = (*ynode).getPoints().back().getX();
    // if the cluster corresponding to *xnode and *ynode is admissible, we add
    // the pair (*xnode,*ynode) to the far field vector
    if (is_admissible(xl, xr, yl, yr, eta)) {
      // the line above checks the admissibility condition (eta)
      FarField_.push_back(new BlockCluster_Y(xnode, ynode));
      FarField_xnds_.push_back(xnode);
      FarField_ynds_.push_back(ynode);
    } else {  // else we consider the children of *xnode and *ynode and check
              // whether their clusters are admissible
      setNearFar_recursion((*xnode).getLChild(), (*ynode).getLChild(), eta);
      setNearFar_recursion((*xnode).getRChild(), (*ynode).getLChild(), eta);
      setNearFar_recursion((*xnode).getLChild(), (*ynode).getRChild(), eta);
      setNearFar_recursion((*xnode).getRChild(), (*ynode).getRChild(), eta);
    }
  }
}

// compute the Far and Near Field pairs
template <>
void HierarchicalPartitioning<BlockCluster_Y, Node_Y>::setNearFar() {
  setNearFar_recursion(Tx_.getRoot(), Ty_.getRoot(), eta_);
  auto checkpointers = [](Node* x, Node* y) -> bool { return x < y; };
  std::sort(FarField_xnds_.begin(), FarField_xnds_.end(), checkpointers);
  FarField_xnds_.erase(
      std::unique(FarField_xnds_.begin(), FarField_xnds_.end()),
      FarField_xnds_.end());
  std::sort(FarField_ynds_.begin(), FarField_ynds_.end(), checkpointers);
  FarField_ynds_.erase(
      std::unique(FarField_ynds_.begin(), FarField_ynds_.end()),
      FarField_ynds_.end());
}

// return the bounding box corresponding to index i of the far-field vector
template <>
std::pair<std::pair<double, double>, std::pair<double, double> >
HierarchicalPartitioning<BlockCluster_Y, Node_Y>::getBB(int i) {
  Node* xnode = FarField_[i]->getXNode();
  std::vector<Point> xpts = xnode->getPoints();
  double xmin = xpts[0].getX();
  double xmax = xpts[0].getX();
  for (int j = 1; j < xpts.size(); j++) {
    if (xpts[j].getX() < xmin) {
      xmin = xpts[j].getX();
    }
    if (xpts[j].getX() > xmax) {
      xmax = xpts[j].getX();
    }
  }
  Node* ynode = FarField_[i]->getYNode();
  std::vector<Point> ypts = ynode->getPoints();
  double ymin = ypts[0].getX();
  double ymax = ypts[0].getX();
  for (int j = 1; j < ypts.size(); j++) {
    if (ypts[j].getX() < ymin) {
      ymin = ypts[j].getX();
    }
    if (ypts[j].getX() > ymax) {
      ymax = ypts[j].getX();
    }
  }
  return {{xmin, ymin}, {xmax, ymax}};
}
