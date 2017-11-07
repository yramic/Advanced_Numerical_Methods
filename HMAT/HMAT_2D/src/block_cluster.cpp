#include "../include/block_cluster.hpp"
#include "../include/cheby.hpp"
#include <iostream>

// constructor for BEM(2D)
/*BlockCluster::BlockCluster(Eigen::VectorXd tk1x, Eigen::VectorXd tk1y, Eigen::VectorXd tk2x, Eigen::VectorXd tk2y, unsigned deg, Kernel* G):
    tk1x_(tk1x), tk1y_(tk1y), tk2x_(tk2x), tk2y_(tk2y), deg_(deg), X_(Eigen::MatrixXd::Zero((deg+1)*(deg+1),(deg+1)*(deg+1)))
{
    setMatrix2D(G);
}*/

// Constructor
BlockCluster::BlockCluster(Node* ndx, Node* ndy):
    pair_(std::make_pair(ndx,ndy))
{ }

// Constructor
BlockCluster::BlockCluster(Node* ndx, Node* ndy, Kernel* G):
    pair_(std::make_pair(ndx,ndy))
{ }

// compute matrix $C_{\sigma,\mu}$
void BlockCluster::setMatrix2D(Kernel* G)
{
    Eigen::VectorXd tk1x = pair_.first->getTkx();
    Eigen::VectorXd tk1y = pair_.first->getTky();
    Eigen::VectorXd tk2x = pair_.second->getTkx();
    Eigen::VectorXd tk2y = pair_.second->getTky();
    C_.resize(tk1x.size()*tk1x.size(),tk1y.size()*tk1y.size());
    for(int i=0; i<tk1x.size(); i++){
        for(int j=0; j<tk1y.size(); j++){
            for(int k=0; k<tk2x.size(); k++){
                for(int l=0; l<tk2y.size(); l++){
                    C_(i*(tk1y.size())+j,k*(tk2y.size())+l) = (*G)(tk1x[i],tk1y[j],tk2x[k],tk2y[l]);
                }
            }
        }
    }
}

// compute CVc vector and store it into xnode
void BlockCluster::setCVc(Kernel* G) {
    setMatrix2D(G); // here because needed for each pair of nodes,
                  // cannot be moved to pre-processing
    Eigen::VectorXd  Vc = pair_.second->getVc_Node();
    Eigen::VectorXd CVc = C_ * Vc;
    pair_.first->setCVc(CVc); // xnode
}
