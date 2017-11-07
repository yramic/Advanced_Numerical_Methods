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
#ifndef BLOCK_CLUSTER_UNI_HPP
#define BLOCK_CLUSTER_UNI_HPP

#include <Eigen/Dense>
#include "../block_cluster.hpp"

/*!
 * \brief Block cluster class to compute matrix \f$X_{\sigma,\mu}\f$ for uni-directional interpolation
 */
class BlockClusterUni : public BlockCluster
{
    using BlockCluster::BlockCluster; // C++11 inheritance of constructors

public:
    /*!
     * \brief compute matrix $C_{\sigma,\mu}$ for uni-directional interpolation
     */
    void setMatrix();
};
#endif // BLOCK_CLUSTER_UNI_HPP
