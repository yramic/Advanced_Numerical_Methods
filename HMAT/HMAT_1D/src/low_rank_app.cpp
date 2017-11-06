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
#include "../include/low_rank_app.hpp"
#include "../include/block_cluster.hpp"
#include "../include/ctree.hpp"
#include "../include/kernel.hpp"
#include "../include/node.hpp"
#include "../include/point.hpp"
#include <iostream>

// constructor
LowRankApp::LowRankApp(Kernel kernel, const std::vector<Point> &GPoints, double eta, unsigned deg):
    kernel_(kernel), GPoints_(GPoints), HP_(GPoints,eta,deg), deg_(deg)
{ }

// approximate matrix-vector multiplication
Eigen::VectorXd LowRankApp::mvProd(const Eigen::VectorXd& c)
{
    // compute Far and Near Field relationships between Nodes of the Cluster Tree
    HP_.setNearFar();
    // compute far field contribution
    Eigen::VectorXd f_approx = Eigen::VectorXd::Zero(c.size());
    ff_contribution(HP_.getFF(), c, f_approx);
    // compute near-field contribution
    nf_contribution(HP_.getNF(), c, f_approx);

    return f_approx;
}

// pre-processing: initialize matrix V and vector Vc for all far field nodes
// do these steps once for each node, not every time the node appears in a pair
void LowRankApp::preProcess(std::vector<BlockCluster> ff_v, const Eigen::VectorXd& c)
{
    for(auto& pair : ff_v){ // iterate for all the pairs of far field nodes
        pair.setV();
        pair.setVc(c);
    }
}

// block-processing: compute vector CVc for all far field pairs and store it into xnode
// all vectors CVc of an xnode can already be summed together
void LowRankApp::blockProcess(std::vector<BlockCluster> ff_v)
{
    for(auto& pair : ff_v){ // iterate for all the pairs of far field nodes
        pair.setKernel(kernel_);
        pair.setCVc();
    }
}

// post-processing: compute vector Vx*CVc for all far field xnodes and add it to vector f in the right place
void LowRankApp::postProcess(std::vector<BlockCluster> ff_v, Eigen::VectorXd& f)
{
    for(auto& pair : ff_v){ // iterate for all the pairs of far field nodes
        Node* xnode = pair.getXNode();
        Eigen::VectorXd CVc = xnode->getCVc_Node();
        Eigen::MatrixXd  Vx = xnode->getV_Node();
        Eigen::VectorXd f_seg = Vx * CVc;
        for(int i=0; i<xnode->getPoints().size(); i++){
            f[xnode->getPoints()[i].getId()] += f_seg[i]; // add contribution of far field to ``f''
        }
    }
}

// compute far field contribution
void LowRankApp::ff_contribution(std::vector<BlockCluster> ff_v, const Eigen::VectorXd& c, Eigen::VectorXd& f)
{
    preProcess(ff_v, c);
    blockProcess(ff_v);
    postProcess(ff_v, f);
}

// compute near-field contribution
void LowRankApp::nf_contribution(std::vector<std::pair<Node*,Node*> > nf_v, const Eigen::VectorXd& c, Eigen::VectorXd& f)
{
    int n = nf_v.size();
    for(int i = 0; i<n; i++){
        Node* xnode = nf_v[i].first;
        Node* ynode = nf_v[i].second;
        for(int j=0; j<xnode->getPoints().size(); j++){
            for(int k=0; k<ynode->getPoints().size(); k++){
                f(xnode->getPoints()[j].getId()) += kernel_(xnode->getPoints()[j].getX(), ynode->getPoints()[k].getX()) * c(ynode->getPoints()[k].getId()); // add near field contribution to ``f''
            }
        }
    }
}
