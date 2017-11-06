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
#include "../include/kernel.hpp"
#include "../include/node.hpp"

// Constructor
BlockCluster::BlockCluster(Node* ndx, Node* ndy):
    pair_(std::make_pair(ndx,ndy)), deg_(deg)
{ }

// Constructor
BlockCluster::BlockCluster(Node* ndx, Node* ndy, Kernel G):
    pair_(std::make_pair(ndx,ndy)), deg_(deg), G_(G)
{
    setMatrix();
}

// compute matrix $C_{\sigma,\mu}$
void BlockCluster::setMatrix()
{
    Eigen::VectorXd tkx = pair_.first->getTK();
    Eigen::VectorXd tky = pair_.second->getTK();
    C_.resize(tkx.size(),tky.size());
    // Compute collocation matrix for Chebychev nodes
    for(unsigned i=0; i<=tkx.size(); ++i)
        for(unsigned j=0; j<=tky.size(); ++j)
            C_(i,j) = G_(tkx[i],tky[j]);
}

// set kernel
void BlockCluster::setKernel(Kernel G) {
    G_ = G;
}

// compute V matrix
void BlockCluster::setV() {
    pair_.first->setV(); // xnode
}

// compute Vc vector
void BlockCluster::setVc(const Eigen::VectorXd& c) {
    pair_.second->setVc(c); // ynode
}

// compute CVc vector and store it into xnode
void BlockCluster::setCVc() {
    setMatrix(); // here because needed for each pair of nodes,
                 // cannot be moved to pre-processing
    Eigen::VectorXd  Vc = pair_.second->getVc_Node();
    Eigen::VectorXd CVc = C_ * Vc;
    pair_.first->setCVc(CVc); // xnode
}
