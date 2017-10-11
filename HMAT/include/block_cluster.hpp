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
    BlockCluster(double x1l, double x1r, double y1l, double y1r, double x2l, double x2r, double y2l, double y2r, unsigned deg, Kernel G);
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
    double x1l_; // left  x coordinate of first bounding box
    double x1r_; // right x coordinate of first bounding box
    double y1l_; // left  y coordinate of first bounding box
    double y1r_; // right y coordinate of first bounding box
    double x2l_; // left  x coordinate of second bounding box
    double x2r_; // right x coordinate of second bounding box
    double y2l_; // left  y coordinate of second bounding box
    double y2r_; // right y coordinate of second bounding box
    unsigned deg_; // degree of interpolating polynomial
    Kernel G_; // kernel
    Eigen::MatrixXd X_; // matrix $X_{\sigma,\mu}$, where $\sigma$ and $\mu$ denote the clusters
};

#endif // BLOCK_CLUSTER_HPP
