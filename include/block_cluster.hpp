#ifndef BLOCK_CLUSTER_HPP
#define BLOCK_CLUSTER_HPP

#include <Eigen/Dense>
#include "kernel.hpp"


/**
* \brief Block cluster class to compute matrix $X_{\sigma,\mu}$
*/
class BlockCluster
{
public:

    /**
    * \brief Constructor
    */
    BlockCluster(double xl, double xr, double yl, double yr, unsigned deg_, Kernel G);

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

#endif // BLOCK_CLUSTER_HPP
