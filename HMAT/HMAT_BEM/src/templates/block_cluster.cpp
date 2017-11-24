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
#include "../../include/block_cluster.hpp"
#include "../../include/cheby.hpp"
#include <iostream>

// Constructor
BlockCluster::BlockCluster(Node* ndx, Node* ndy):
    pair_(std::make_pair(ndx,ndy))
{ }

// compute matrix $C_{\sigma,\mu}$
unsigned BlockCluster::setMatrix(Kernel* G)
{
    // TODO
    return C_.rows()*C_.cols(); // return no. of 'operations' performed
}

// compute CVc vector and store it into xnode
unsigned BlockCluster::setCVc()
{
    Eigen::VectorXd  Vc = pair_.second->getVc_Node();
    Eigen::VectorXd CVc = C_ * Vc;
    unsigned nops = C_.rows()*C_.cols();
    nops += pair_.first->setCVc(CVc); // xnode
    return nops; // return no. of 'operations' performed
}

// return matrix $V_{\sigma}C_{\sigma,\mu}V_{\mu}^\top$
Eigen::MatrixXd BlockCluster::getVCV() const
{
    Eigen::MatrixXd Vsigma = pair_.first->getV_Node();
    Eigen::MatrixXd Vmu    = pair_.second->getV_Node();
    Eigen::MatrixXd VCV = Vsigma * C_ * Vmu.transpose();
    Eigen::MatrixXd VCVord(VCV.rows(), VCV.cols());
    for(unsigned i=0; i<pair_.first->getSegments().size(); ++i)
        for(unsigned j=0; j<pair_.second->getSegments().size(); ++j)
            VCVord(pair_.first->getSegments()[i].getId(),
                   pair_.second->getSegments()[j].getId()) = VCV(i,j);
    return VCVord;
}
