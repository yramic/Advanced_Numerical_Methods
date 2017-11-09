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
#include "../include/block_nearf.hpp"
#include "../include/kernel.hpp"
#include "../include/node.hpp"
#include <Eigen/Dense>

// Constructor
BlockNearF::BlockNearF(Node* xnode, Node* ynode):
    pair_(std::make_pair(xnode,ynode))
{ }

// Constructor
BlockNearF::BlockNearF(Node* xnode, Node* ynode, Kernel* G):
    pair_(std::make_pair(xnode,ynode))
{ }

// compute near field block matrix
unsigned BlockNearF::setMatrix(Kernel* G)
{
    C_.resize(pair_.first->getPPoints().size(),
              pair_.second->getPPoints().size());
    // Compute collocation matrix for near field points
    for(unsigned i=0; i<pair_.first->getPPoints().size(); ++i)
        for(unsigned j=0; j<pair_.second->getPPoints().size(); ++j)
            C_(i,j) = (*G)(pair_.first->getPPoints()[i].getX(),pair_.first->getPPoints()[i].getY(),
                           pair_.second->getPPoints()[j].getX(),pair_.second->getPPoints()[j].getY());
    return C_.rows()*C_.cols(); // return no. of 'operations' performed
}

