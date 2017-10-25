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
    * \brief Constructor for 4D
    */
    //BlockCluster(double x1l, double x1r, double y1l, double y1r, double x2l, double x2r, double y2l, double y2r, unsigned deg, Kernel4D* G);
    //BlockCluster(double x1l, double x1r, double y1l, double y1r, double x2l, double x2r, double y2l, double y2r, unsigned deg, PolynomialKernel* G);
    //BlockCluster(double x1l, double x1r, double y1l, double y1r, double x2l, double x2r, double y2l, double y2r, unsigned deg, ConstantKernel* G);
    BlockCluster(double x1l, double x1r, double y1l, double y1r, double x2l, double x2r, double y2l, double y2r, unsigned deg, Kernel* G);
    /*!
    * \brief Constructor for 2D
    */
    BlockCluster(double xl, double xr, double yl, double yr, unsigned deg, Kernel2D G);


    /*!
    * \brief return matrix \f$X_{\sigma,\mu}\f$, where \f$\sigma\f$ and \f$\mu\f$ denote the clusters
    */
    // return matrix $X_{\sigma,\mu}$, where $\sigma$ and $\mu$ denote the clusters
    Eigen::MatrixXd getMatrix() const {
        return X_;
    }

    /*!
    * \brief compute matrix \f$X_{\sigma,\mu}\f$ for 2D
    */
    // compute matrix $X_{\sigma,\mu}$
    void setMatrix2D();
    /*!
    * \brief compute matrix \f$X_{\sigma,\mu}\f$ for 4D
    */
    void setMatrix4D(Kernel* G);
private:

    double xl_;     //!< left  boundary of bounding box of *xcluster
    double xr_;     //!< right boundary of bounding box of *xcluster
    double yl_;     //!< left  boundary of bounding box of *ycluster
    double yr_;     //!< right boundary of bounding box of *ycluster
    double x1l_;    //!< left  x coordinate of first bounding box
    double x1r_;    //!< right x coordinate of first bounding box
    double y1l_;    //!< left  y coordinate of first bounding box
    double y1r_;    //!< right y coordinate of first bounding box
    double x2l_;    //!< left  x coordinate of second bounding box
    double x2r_;    //!< right x coordinate of second bounding box
    double y2l_;    //!< left  y coordinate of second bounding box
    double y2r_;    //!< right y coordinate of second bounding box
    unsigned deg_;  //!< degree of interpolating polynomial
    Kernel2D G2_;           //!< 2D kernel
    Kernel4D G4_;           //!< 4D Galerkin kernel
    PolynomialKernel GP_;   //!< POlynomial Kernel
    ConstantKernel GC_;     //!< Constant Kernel
    Eigen::MatrixXd X_;     //!< matrix \f$X_{\sigma,\mu}\f$, where \f$\sigma\f$ and \f$\mu\f$ denote the clusters
};

#endif // BLOCK_CLUSTER_HPP
