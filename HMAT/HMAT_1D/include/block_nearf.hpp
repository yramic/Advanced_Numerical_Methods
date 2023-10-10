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
#ifndef BLOCK_NEARF_HPP
#define BLOCK_NEARF_HPP

#include <Eigen/Dense>

#include "kernel.hpp"
#include "node.hpp"

/*!
 * \brief Helper class to compute near field block matrices
 */
class BlockNearF {
 public:
  /*!
     * \brief Constructor
     * \param ndx XNode
     * \param ndy YNode
     * \param G Kernel Function
     */
  BlockNearF(Node* xnode, Node* ynode);

  /*!
     * \brief return near field block matrix
     */
  Eigen::MatrixXd getMatrix() const { return C_; }

  /*!
     * \brief compute near field block matrix
     * \return no. of 'operations' performed
     */
  unsigned setMatrix(Kernel* G);

  /*!
     * \brief return pointer to xnode
     */
  Node* getXNode() { return pair_.first; }

  /*!
     * \brief return pointer to ynode
     */
  Node* getYNode() { return pair_.second; }

 private:
  std::pair<Node*, Node*> pair_;
  Eigen::MatrixXd C_;  //!< near field block matrix
};
#endif  // BLOCK_NEARF_HPP
