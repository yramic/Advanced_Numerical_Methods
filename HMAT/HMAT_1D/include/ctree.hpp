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
#ifndef CTREE_HPP
#define CTREE_HPP

#include <Eigen/Dense>
#include <vector>

#include "point.hpp"

// forward declarations to avoid cross-referencing
class Node;
class NodeY;

/**
* \brief Cluster tree class
*/
template <typename NODE_Y = Node>
class cTree {
 public:
  /*!
    * \brief Default Constructor
    */
  cTree() : root_(NULL) {}
  /*!
    * \brief Actual  constructor
    * \param x Vector of oints of the grid
    * \param deg Degree of interpolation
    */
  cTree(const std::vector<Point>& x, unsigned deg);
  /*!
    * \brief copy the root of a tree
    * \param T Tree to copy the root from
    */
  cTree(const cTree& T) : root_(T.getRoot()) {}
  /*!
    * \brief return a pointer to node-root of "cTree"
    */
  Node* getRoot() const { return root_; }

 private:
  Node* root_;  //!< pointer to node-root of ``cTree''
  friend class Node;
  friend class NodeY;
};
#endif  // CTREE_HPP
