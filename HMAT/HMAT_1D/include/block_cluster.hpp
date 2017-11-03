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
#ifndef BLOCK_CLUSTER_HPP
#define BLOCK_CLUSTER_HPP

#include <Eigen/Dense>
#include "kernel.hpp"

/*!
 * \brief Block cluster class to compute matrix \f$X_{\sigma,\mu}\f$
 */
class BlockCluster
{
public:

    /*!
     * \brief Constructor
     * \param tkx Chebyshev nodes for the x axis
     * \param tky Chebyshev nodes for the y axis
     * \param deg degree of interpolation
     * \param G Kernel Function
     */
    BlockCluster(Eigen::VectorXd tkx, Eigen::VectorXd tky, unsigned deg, Kernel G);

    /*!
     * \brief return matrix \f$X_{\sigma,\mu}\f$, where \f$\sigma\f$ and \f$\mu\f$ denote the clusters
     */
    Eigen::MatrixXd getMatrix() const {
        return X_;
    }

    /*!
     * \brief compute matrix \f$X_{\sigma,\mu}\f$
     */
    void setMatrix();

private:
    Eigen::VectorXd tkx_;   //!< Chebychev nodes for x axis
    Eigen::VectorXd tky_;   //!< Chebyshev nodes for y axis
    unsigned deg_; //!< degree of interpolating polynomial
    Kernel G_; //!< kernel
    Eigen::MatrixXd X_; //!< matrix \f$X_{\sigma,\mu}\f$, where \f$\sigma\f$ and \f$\mu\f$ denote the clusters
};
#endif // BLOCK_CLUSTER_HPP
