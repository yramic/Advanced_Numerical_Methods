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
#include "../../../include/uni-direct/block_cluster_Y.hpp"
#include <Eigen/Dense>

// compute matrix $C_{\sigma,\mu}$ for uni-directional interpolation
/* SAM_LISTING_BEGIN_0 */
unsigned BlockCluster_Y::setMatrix(Kernel* G)
{
    Eigen::VectorXd tkx = pair_.first->getTK();
    C_.resize(tkx.size(),
              pair_.second->getPoints().size());
    // Compute collocation matrix
    // for Chebychev nodes along x-dimension
    // and for grid points along y-dimension
    for(unsigned i=0; i<tkx.size(); ++i)
        for(unsigned j=0; j<pair_.second->getPoints().size(); ++j)
            C_(i,j) = (*G)(tkx(i),
                           pair_.second->getPoints()[j].getX());
    return C_.rows()*C_.cols(); // return no. of 'operations' performed
}
/* SAM_LISTING_END_0 */
