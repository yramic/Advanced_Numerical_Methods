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
#include <Eigen/Dense>
#include "../include/block_nearf.hpp"
#include "../include/kernel.hpp"
#include "../include/node.hpp"

// Constructor
BlockNearF::BlockNearF(Node* xnode, Node* ynode):
    pair_(std::make_pair(xnode,ynode))
{ }

// Constructor
BlockNearF::BlockNearF(Node* xnode, Node* ynode, Kernel G):
    pair_(std::make_pair(xnode,ynode)), G_(G)
{
    setMatrix();
}

// compute near field block matrix
void BlockNearF::setMatrix()
{
    C_.resize(pair_.first->getPoints().size(),
              pair_.second->getPoints().size());
    // Compute collocation matrix for near field points
    for(unsigned i=0; i<pair_.first->getPoints().size(); ++i)
        for(unsigned j=0; j<pair_.second->getPoints().size(); ++j)
            C_(i,j) = G_(pair_.first->getPoints()[i].getX(),
                         pair_.second->getPoints()[j].getX());
}

// set kernel
void BlockNearF::setKernel(Kernel G) {
    G_ = G;
}
