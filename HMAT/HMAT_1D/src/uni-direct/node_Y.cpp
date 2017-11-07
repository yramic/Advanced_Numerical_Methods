/***********************************************************************
 *                                                                     *
 * Code for Course "Advanced Numerical Methods for CSE"                *
 * (Prof. Dr. R. Hiptmair)                                             * 
 * Author:                                                             *
 * Date:                                                               *
 * (C) Seminar for Applied Mathematics, ETH Zurich                     *
 * This code can be freely used for non-commercial purposes as long    *
 * as this header is left intact.                                      *
 ***********************************************************************/
#include <Eigen/Dense>
#include "../../include/uni-direct/node_Y.hpp"
#include "../../include/cheby.hpp"
#include <iostream>

// destructor
Node_Y::~Node_Y()
{
    if((l_child_ == NULL) && (r_child_ == NULL))
        std::cout << "leaves destroyed" << std::endl;
    if(l_child_ != NULL) delete l_child_;
    if(r_child_ != NULL) delete r_child_;
}

// compute fake V-matrix of ynode with uni-directional interpolation
//(just the identity)
void Node_Y::setV()
{
    int n = node_points_.size();
    V_node_ = Eigen::MatrixXd::Identity(n, n);
}
