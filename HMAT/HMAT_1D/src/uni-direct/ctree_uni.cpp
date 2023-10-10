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
#include "../../include/ctree.hpp"
#include "../../include/uni-direct/node_Y.hpp"
#include <Eigen/Dense>
#include <vector>

// Actual constructor
template<>
cTree<Node_Y>::cTree(const std::vector<Point>& GPoints, unsigned deg):
    root_(NULL)
{
    unsigned n = GPoints.size();
    if(n > 1) { // build the tree
        int node_id = 0;
        root_ = new Node_Y(GPoints, node_id, deg); // root is a node, leaves are added
    }
}
