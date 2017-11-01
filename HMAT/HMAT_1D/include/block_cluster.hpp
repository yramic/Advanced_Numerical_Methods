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
//#define ver1
#define ver2
#include <Eigen/Dense>
#include "kernel.hpp"

#ifdef ver2
/**
* \brief Block cluster class to compute matrix $X_{\sigma,\mu}$
*/
class BlockCluster
{
public:

    /**
    * \brief Constructor
    */
    BlockCluster(Eigen::VectorXd tkx, Eigen::VectorXd tky, unsigned deg, Kernel G);

    /**
    * \brief Getter
    */
    // return matrix $X_{\sigma,\mu}$, where $\sigma$ and $\mu$ denote the clusters
    Eigen::MatrixXd getMatrix() const {
        return X_;
    }

    /**
    * \brief Setter
    */
    // compute matrix $X_{\sigma,\mu}$
    void setMatrix();

private:
    Eigen::VectorXd tkx_;   // Chebychev nodes for x axis
    Eigen::VectorXd tky_;   // Chebyshev nodes for y axis
    unsigned deg_; // degree of interpolating polynomial
    Kernel G_; // kernel
    Eigen::MatrixXd X_; // matrix $X_{\sigma,\mu}$, where $\sigma$ and $\mu$ denote the clusters
};
#endif
#ifdef ver1
    /**
    * \brief Block cluster class to compute matrix $X_{\sigma,\mu}$
    */
    class BlockCluster
    {
    public:

        /**
        * \brief Constructor
        */
        BlockCluster(double xl, double xr, double yl, double yr, unsigned deg, Kernel G);

        /**
        * \brief Getter
        */
        // return matrix $X_{\sigma,\mu}$, where $\sigma$ and $\mu$ denote the clusters
        Eigen::MatrixXd getMatrix() const {
            return X_;
        }

        /**
        * \brief Setter
        */
        // compute matrix $X_{\sigma,\mu}$
        void setMatrix();

    private:

        double xl_;	// left  boundary of bounding box of *xcluster
        double xr_;	// right boundary of bounding box of *xcluster
        double yl_;	// left  boundary of bounding box of *ycluster
        double yr_;	// right boundary of bounding box of *ycluster
        unsigned deg_; // degree of interpolating polynomial
        Kernel G_; // kernel
        Eigen::MatrixXd X_; // matrix $X_{\sigma,\mu}$, where $\sigma$ and $\mu$ denote the clusters
    };
#endif
#endif // BLOCK_CLUSTER_HPP
