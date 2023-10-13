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
#ifndef BLOCK_CLUSTER_HPP
#define BLOCK_CLUSTER_HPP

#include <Eigen/Dense>

#include "kernel.hpp"
#include "node.hpp"

/*!
* \brief Block cluster class to compute matrix \f$X_{\sigma,\mu}\f$
*/
class BlockCluster {
 public:
  /*!
     * \brief Constructor
     * \param ndx Node* x
     * \param ndy Node* y
     */
  BlockCluster(Node* ndx, Node* ndy);

  /*!
    * \brief return matrix \f$C_{\sigma,\mu}\f$, where \f$\sigma\f$ and \f$\mu\f$ denote the clusters
    */
  Eigen::MatrixXd getMatrix() const { return C_; }

  /*!
    * \brief compute matrix \f$C_{\sigma,\mu}\f$ for 2D
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
  Node* getXNode() { return pair_.first; }

  /*!
     * \brief return pointer to ynode
     */
  Node* getYNode() { return pair_.second; }

  /*!
     * \brief return matrix \f$V_{\sigma}C_{\sigma,\mu}V_{\mu}^\top\f$,
     * where \f$\sigma\f$ and \f$\mu\f$ denote the clusters
     */
  Eigen::MatrixXd getVCV() const;

 private:
  std::pair<Node*, Node*> pair_;
  Eigen::MatrixXd
      C_;  //!< Matrix \f$C_{\sigma,\mu}\f$, where \f$\sigma\f$ and \f$\mu\f$ denote the clusters
};

#endif  // BLOCK_CLUSTER_HPP
