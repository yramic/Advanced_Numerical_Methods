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

#include "kernel.hpp"
#include "node.hpp"
#include <Eigen/Dense>
#include <fstream>

/*!
 * \brief Block cluster class to compute matrix \f$X_{\sigma,\mu}\f$
 */
class BlockCluster
{
public:

    /*!
     * \brief Constructor
     * \param ndx XNode
     * \param ndy YNode
     * \param G Kernel Function
     */
    BlockCluster(Node* xnode, Node* ynode);

    /*!
     * \brief Constructor
     * \param ndx XNode
     * \param ndy YNode
     * \param G Kernel Function
     */
    BlockCluster(Node* xnode, Node* ynode, Kernel G);

    /*!
     * \brief return matrix \f$C_{\sigma,\mu}\f$, where \f$\sigma\f$ and \f$\mu\f$ denote the clusters
     */
    Eigen::MatrixXd getMatrix() const {
        return C_;
    }

    /*!
     * \brief compute matrix \f$C_{\sigma,\mu}\f$
     * \param G any kind of derived kernel from base class Kernel
     * \return no. of 'operations' performed
     */
    unsigned setMatrix(Kernel* G);

    /*!
     * \brief compute CVc vector and store it in xnode
     * \return no. of 'operations' performed
     */
    unsigned setCVc();

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

    /*!
     * \brief return matrix \f$V_{\sigma}C_{\sigma,\mu}V_{\mu}^\top\f$,
     * where \f$\sigma\f$ and \f$\mu\f$ denote the clusters
     */
    Eigen::MatrixXd getVCV() const;

protected:
    std::pair<Node*,Node*> pair_;
    Eigen::MatrixXd C_; //!< matrix \f$C_{\sigma,\mu}\f$, where \f$\sigma\f$ and \f$\mu\f$ denote the clusters
};
#endif // BLOCK_CLUSTER_HPP
