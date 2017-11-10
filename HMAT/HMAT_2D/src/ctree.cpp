/*!
  * \author Daniele Casati, Ioannis Magkanaris
  * \date 11/2017
  * \mainpage Low Rank Approximation for BEM
  */
#include "../include/ctree.hpp"
#include "../include/point.hpp"
#include "../include/is_admissible.hpp"
#include "../include/node.hpp"
#include <Eigen/Dense>
#include <vector>
#include <iostream>

// actual constructor
cTree::cTree(const std::vector<Point> &PPoints, unsigned deg):
    root_(NULL), PPointsTree_(PPoints)
{
    unsigned n = PPoints.size();
    if(n > 0) { // build the tree
        root_ = new Node(PPoints,deg); // root is a node, leaves are added
    }
}
