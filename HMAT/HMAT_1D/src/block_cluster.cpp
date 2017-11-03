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
#include "../include/block_cluster.hpp"
#include "../include/cheby.hpp"

// Constructor
BlockCluster::BlockCluster(Eigen::VectorXd tkx, Eigen::VectorXd tky, unsigned deg, Kernel G):
    tkx_(tkx), tky_(tky), deg_(deg), G_(G), X_(Eigen::MatrixXd::Zero(deg+1,deg+1))
{ setMatrix(); }

// compute matrix $X_{\sigma,\mu}$
void BlockCluster::setMatrix()
{
    // Compute collocation matrix for Chebychev nodes
    for(unsigned i=0; i<=deg_; ++i)
        for(unsigned j=0; j<=deg_; ++j)
            X_(i,j) = G_(tkx_[i],tky_[j]);
}

