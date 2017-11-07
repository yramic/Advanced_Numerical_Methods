#ifndef BLOCK_CLUSTER_HPP
#define BLOCK_CLUSTER_HPP

#include <Eigen/Dense>
#include "kernel.hpp"
#include "node.hpp"

/*!
* \brief Block cluster class to compute matrix \f$X_{\sigma,\mu}\f$
*/
class BlockCluster
{
public:
    /*!
    * \brief Constructor for 2D
    * \param tk1x Chebyshev nodes for the x axis of the first bounding box
    * \param tk1y Chebyshev nodes for the y axis of the first bounding box
    * \param tk2x Chebyshev nodes for the x axis of the second bounding box
    * \param tk2y Chebyshev nodes for the y axis of the second bounding box
    * \param deg degree of interpolation
    * \param G any kind of derived kernel from base class Kernel
    */
    //BlockCluster(Eigen::VectorXd tk1x, Eigen::VectorXd tk1y, Eigen::VectorXd tk2x, Eigen::VectorXd tk2y, unsigned deg, Kernel* G);

    /*!
     * \brief Constructor
     * \param ndx Node* x
     * \param ndy Node* y
     */
    BlockCluster(Node* ndx, Node* ndy);

    /*!
     * \brief Constructor
     * \param ndx Node* x
     * \param ndy Node* y
     * \param G Kernel Function
     */
    BlockCluster(Node* ndx, Node* ndy, Kernel* G);

    /*!
    * \brief return matrix \f$C_{\sigma,\mu}\f$, where \f$\sigma\f$ and \f$\mu\f$ denote the clusters
    */
    Eigen::MatrixXd getMatrix() const {
        return C_;
    }
    /*!
    * \brief compute matrix \f$C_{\sigma,\mu}\f$ for 2D
    * \param G any kind of derived kernel from base class Kernel
    */
    void setMatrix2D(Kernel* G);

    /*!
     * \brief compute CVc vector and store it in xnode
     * \param c vector for multiplication with the matrix
     */
    void setCVc(Kernel* G);

    /*!
     * \brief return pointer to xnode
     */
    Node* getXNode() {
        return pair_.first;
    }

    /*!
     * \brief return pointer to ynode
     */
    Node* getYNode() {
        return pair_.second;
    }
private:
    /*Eigen::VectorXd tk1x_;  //!< Chebychev nodes for x axis
    Eigen::VectorXd tk1y_;  //!< Chebyshev nodes for y axis
    Eigen::VectorXd tk2x_;  //!< Chebychev nodes for x axis
    Eigen::VectorXd tk2y_;  //!< Chebyshev nodes for y axis*/
    //unsigned deg_;          //!< degree of interpolating polynomial
    std::pair<Node*,Node*> pair_;
    Eigen::MatrixXd C_;     //!< Matrix \f$C_{\sigma,\mu}\f$, where \f$\sigma\f$ and \f$\mu\f$ denote the clusters
};

#endif // BLOCK_CLUSTER_HPP
