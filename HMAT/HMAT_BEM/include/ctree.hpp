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

#include "segment.hpp"

// forward declaration to avoid cross-referencing
class Node;

/*!
* \brief Cluster tree class
*/
class cTree {
 public:
  /*!
    * \brief Default Constructor
    */
  cTree() : root_(NULL) {}
  /*!
    * \brief Actual Constructor for 2D
    * \param segments Vector of segments
    * \param deg Degree of interpolation
    */
  cTree(const std::vector<Segment>& segments, unsigned deg);
  /*!
    * \brief Copy Constructor for 2D
    * \param T The cTree to copy
    */
  cTree(const cTree& T) : root_(T.getRoot()) {}
  /*!
    * \brief Return a pointer to node-root of "cTree"
    */
  Node* getRoot() const { return root_; }

 private:
  Node* root_;  //!< pointer to node-root of "cTree"
  const std::vector<Segment>
      segments_;  //!< vector that has all the Polygon Points
  friend class Node;
};

#endif  // CTREE_HPP
