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
#ifndef BLOCK_CLUSTER_Y_HPP
#define BLOCK_CLUSTER_Y_HPP

#include <Eigen/Dense>

#include "../block_cluster.hpp"

/*!
 * \brief Block cluster class to compute matrix \f$X_{\sigma,\mu}\f$ for uni-directional interpolation
 */
class BlockCluster_Y : public BlockCluster {
  using BlockCluster::BlockCluster;  // C++11 inheritance of constructors

 public:
  /*!
     * \brief compute matrix $C_{\sigma,\mu}$ for uni-directional interpolation
     * \param G any kind of derived kernel from base class Kernel
     * \return no. of 'operations' performed
     */
  unsigned setMatrix(Kernel* G);
};
#endif  // BLOCK_CLUSTER_Y_HPP
