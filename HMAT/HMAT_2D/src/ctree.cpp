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
#include "../include/ctree.hpp"
#include "../include/is_admissible.hpp"
#include "../include/node.hpp"
#include "../include/point.hpp"
#include <Eigen/Dense>
#include <vector>
#include <iostream>

// actual constructor
cTree::cTree(const std::vector<Point>& PPoints, unsigned deg):
    root_(NULL), PPointsTree_(PPoints)
{
    unsigned n = PPoints.size();
    if(n > 0) { // build the tree
        root_ = new Node(PPoints,deg); // root is a node, leaves are added
    }
}
