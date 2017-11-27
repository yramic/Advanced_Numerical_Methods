/***********************************************************************
 *                                                                     *
 * Code for Course "Advanced Numerical Methods for CSE"                *
 * (Prof. Dr. R. Hiptmair)                                             *
 * Author: Ioannis Magkanaris                                          *
 * Date: 11/2017                                                       *
 * (C) Seminar for Applied Mathematics, ETH Zurich                     *
 * This code can be freely used for non-commercial purposes as long    *
 * as this header is left intact.                                      *
 ***********************************************************************/
#include "../../include/block_cluster.hpp"
#include "../../include/cheby.hpp"

// Constructor
BlockCluster::BlockCluster(Node* ndx, Node* ndy):
    pair_(std::make_pair(ndx,ndy))
{ }

// compute matrix $C_{\sigma,\mu}$
unsigned BlockCluster::setMatrix(Kernel* G)
{
    /* SAM_LISTING_BEGIN_0 */
    Eigen::VectorXd tk1x = pair_.first->getTkx();
    Eigen::VectorXd tk1y = pair_.first->getTky();
    Eigen::VectorXd tk2x = pair_.second->getTkx();
    Eigen::VectorXd tk2y = pair_.second->getTky();
    C_.resize(tk1x.size()*tk1x.size(),tk1y.size()*tk1y.size());
    for(int i=0; i<tk1x.size(); i++){
        for(int j=0; j<tk1y.size(); j++){
            for(int k=0; k<tk2x.size(); k++){
                for(int l=0; l<tk2y.size(); l++){
                    C_(i*(tk1y.size())+j,
                       k*(tk2y.size())+l) = (*G)(tk1x[i],tk1y[j],tk2x[k],tk2y[l]);
                }
            }
        }
    }
    /* SAM_LISTING_END_0 */
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
    return Vsigma * C_ * Vmu.transpose();
}
