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
#include "../../../include/uni-direct/node_Y.hpp"
#include "../../../include/cheby.hpp"
#include "../../../include/point.hpp"
#include <Eigen/Dense>
#include <iostream>

// actual constructor: adds a tree below the node if left_index != right_index
Node_Y::Node_Y(std::vector<Point> points, int& id, unsigned deg)
{
    l_child_ = NULL; r_child_ = NULL;
    node_points_ = points; nodeId_ = id; deg_ = deg;
    CVc_node_ = Eigen::VectorXd::Zero(deg+1);
    setSons(id);
}

// destructor
Node_Y::~Node_Y()
{
    if((l_child_ == NULL) && (r_child_ == NULL))
        std::cout << "leaves destroyed" << std::endl;
    if(l_child_ != NULL) delete l_child_;
    if(r_child_ != NULL) delete r_child_;
}

// build tree recursively
void Node_Y::setSons(int& id)
{
    if(node_points_.size() > 1) {
        std::vector<Point>::iterator it;                    // we assume that the points are in ascending order
        it=node_points_.begin()+(node_points_.size()+1)/2;  // set iterator in the middle of the vector of points in this node
        std::vector<Point> lc,rc;                           // left and right child´s vectors
        lc.assign(node_points_.begin(), it);                // division of node´s points into it´s childs
        rc.assign(it, node_points_.end());
        id++;                                               // increase the identification number of the nodes of the tree
        l_child_ = new Node_Y(lc,id,deg_);                  // recursivly do the same for the left child
        id++;                                               // increase the identification number of the nodes of the tree
        r_child_ = new Node_Y(rc,id,deg_);                  // recursivly do the same for the right child
        double xmin = node_points_.front().getX();          // find min coordinate of the axis
        double xmax = node_points_.back().getX();           // find max coordinate of the axis
        Cheby cb(xmin, xmax, deg_);                         // Checyshev interpolation
        tk_ = cb.getNodes();                                // Chebyshev interpolation nodes
        wk_ = cb.getWghts();                                // Weights of Lagrange polynomial
    }
}

// compute fake V-matrix of ynode with uni-directional interpolation
//(just the identity)
/* SAM_LISTING_BEGIN_0 */
unsigned Node_Y::setV()
{
    int n = node_points_.size();
    V_node_ = Eigen::MatrixXd::Identity(n, n);
    return 1; // return no. of 'operations' performed
}
/* SAM_LISTING_END_0 */
