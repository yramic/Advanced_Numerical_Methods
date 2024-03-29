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
#include "../include/block_nearf.hpp"

#include <Eigen/Dense>

#include "../include/kernel.hpp"
#include "../include/node.hpp"

// Constructor
BlockNearF::BlockNearF(Node* xnode, Node* ynode)
    : pair_(std::make_pair(xnode, ynode)) {}

// Constructor
BlockNearF::BlockNearF(Node* xnode, Node* ynode, Kernel* G)
    : pair_(std::make_pair(xnode, ynode)) {}

// compute near field block matrix
unsigned BlockNearF::setMatrix(Kernel* G) {
  C_.resize(pair_.first->getSegments().size(),
            pair_.second->getSegments().size());
  // Compute collocation matrix for near field points
  for (unsigned i = 0; i < pair_.first->getSegments().size(); ++i)
    for (unsigned j = 0; j < pair_.second->getSegments().size(); ++j)
      C_(i, j) = (*G)(pair_.first->getSegments()[i].getA(),
                      pair_.first->getSegments()[i].getB(),
                      pair_.second->getSegments()[j].getA(),
                      pair_.second->getSegments()[j].getB());
  return C_.rows() * C_.cols();  // return no. of 'operations' performed
}
