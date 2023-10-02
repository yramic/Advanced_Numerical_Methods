/***********************************************************************
 *                                                                     *
 * Code for Course "Advanced Numerical Methods for CSE"                *
 * (Prof. Dr. R. Hiptmair)                                             *
 * Author: Ioannis Magkanaris                                          *
 * Date: 11/2017                                                       *
 * (C) Seminar for Applied Mathematics, ETH Zurich                     *
 * This code can be freely used for non-commercial purposes as long    *
 * as this header is left intact.                                      *
 ***********************************************************************/
#include <iostream>

#include "include/ctree.hpp"
#include "include/hierarchical_partition.hpp"
#include "include/node.hpp"
#include "include/point.hpp"

void setNodeIds(Node* x, int& id) {
  if (x == NULL) return;
  x->setNodeId(id);
  id++;
  // std::cout << id << std::endl;
  setNodeIds(x->getTl_Child(), id);
  setNodeIds(x->getTr_Child(), id);
  setNodeIds(x->getBl_Child(), id);
  setNodeIds(x->getBr_Child(), id);
}

int main() {
  unsigned n = 16;  // number of points
  std::vector<int> x = {24, 22, 73, 63, 14, 17, 39, 99,
                        83, 41, 40, 4,  83, 30, 65, 23},
                   y = {24, 55, 87, 91, 30, 1,  9,  28,
                        67, 85, 10, 84, 57, 72, 86, 56},
                   v = {66, 63, 27, 46, 83, 58, 46, 8,
                        95, 57, 2,  79, 34, 21, 64, 95};
  std::vector<Point> PPoints;  // initalizing Polygon Points properties
  PPoints.reserve(n);
  for (int i = 0; i < n; i++) {
    Point p;
    p.setId(i);
    p.setX(x[i]);
    p.setY(y[i]);
    p.setV(v[i]);
    PPoints.push_back(p);
  }
  unsigned d = 2;
  double eta = 2;
  HierarchicalPartitioning HP(PPoints, eta, d);
  // cTree test_tree(PPoints, d);           // creating a tree
  int id = 0;
  setNodeIds(HP.getTx().getRoot(), id);  // setting ID in each node of the tree
  Eigen::MatrixXd cmatrix = Eigen::MatrixXd::Zero(
      id + 1, id + 1);  // matrix that rows and columns represent the IDs of the
                        // nodes that exist on the tree
  HP.setNearFar(cmatrix);  // set the near and far field of each node and add 1
                           // for every compression that is done in this process
                           // in the referring cell of the matrix
  std::cout << cmatrix << std::endl;
}
