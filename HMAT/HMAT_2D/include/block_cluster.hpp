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
    * \brief Constructor for 2D
    * \param x1l x coordinate of left edge of the first bounind box
    * \param x1r x coordinate of right edge of the first bounind box
    * \param y1l y coordinate of bottom edge of the first bounind box
    * \param y1r y coordinate of top edge of the first bounind box
    * \param x2l x coordinate of left edge of the second bounind box
    * \param x2r x coordinate of right edge of the second bounind box
    * \param y2l y coordinate of bottom edge of the second bounind box
    * \param y2r y coordinate of top edge of the second bounind box
    * \param deg degree of interpolation
    * \param G any kind of derived kernel from base class Kernel
    */ 
    BlockCluster(double x1l, double x1r, double y1l, double y1r, double x2l, double x2r, double y2l, double y2r, unsigned deg, Kernel* G);

    /*!
    * \brief return matrix \f$X_{\sigma,\mu}\f$, where \f$\sigma\f$ and \f$\mu\f$ denote the clusters
    */
    Eigen::MatrixXd getMatrix() const {
        return X_;
    }

    /*!
    * \brief compute matrix \f$X_{\sigma,\mu}\f$ for 2D
    * \param G any kind of derived kernel from base class Kernel
    */
    void setMatrix2D(Kernel* G);
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
    KernelGalerkin G4_;           //!< 2D Galerkin kernel
    PolynomialKernel GP_;   //!< Polynomial Kernel
    ConstantKernel GC_;     //!< Constant Kernel
    Eigen::MatrixXd X_;     //!< Matrix \f$X_{\sigma,\mu}\f$, where \f$\sigma\f$ and \f$\mu\f$ denote the clusters
};

#endif // BLOCK_CLUSTER_HPP
