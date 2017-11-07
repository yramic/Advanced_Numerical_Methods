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
#include "../../include/low_rank_app.hpp"
#include "../../include/uni-direct/block_cluster_Y.hpp"
#include "../../include/block_nearf.hpp"
#include "../../include/ctree.hpp"
#include "../../include/kernel.hpp"
#include "../../include/uni-direct/node_Y.hpp"
#include "../../include/point.hpp"
#include <iostream>

// constructor
template<>
LowRankApp<BlockCluster_Y,Node_Y>::LowRankApp(Kernel kernel, const std::vector<Point> &GPoints, double eta, unsigned deg):
    kernel_(kernel), GPoints_(GPoints), HP_(GPoints,eta,deg), deg_(deg)
{ }

// pre-processing: initialize matrix V and vector Vc for all far field nodes
// do these steps once for each node, not every time the node appears in a pair
template<>
void LowRankApp<BlockCluster_Y,Node_Y>::preProcess(std::vector<Node*> ff_v_x, std::vector<Node*> ff_v_y, const Eigen::VectorXd& c)
{
    for(auto& xnode : ff_v_x){ // iterate for all the far field xnodes
        xnode->setV();
    }
    for(auto& ynode : ff_v_y){ // iterate for all the far field ynodes
        Node_Y* ynode_ = static_cast<Node_Y*>(ynode);
        ynode_->setV();
        ynode_->setVc(c);
    }
}

// block-processing: compute vector CVc for all far field pairs and store it into xnode
// all vectors CVc of an xnode can already be summed together
template<>
void LowRankApp<BlockCluster_Y,Node_Y>::blockProcess(std::vector<BlockCluster*> ff_v)
{
    for(auto& pair : ff_v){ // iterate for all the pairs of far field nodes
        BlockCluster_Y* pair_ = static_cast<BlockCluster_Y*>(pair);
        pair_->setKernel(kernel_);
        pair_->setMatrix(); // here because needed for each pair of nodes,
                            // cannot be moved to pre-processing
        pair_->setCVc();
    }
}

// post-processing: compute vector Vx*CVc for all far field xnodes and add it to vector f in the right place
template<>
void LowRankApp<BlockCluster_Y,Node_Y>::postProcess(std::vector<Node*> ff_v_x, Eigen::VectorXd& f)
{
    for(auto& xnode : ff_v_x){ // iterate for all the far field xnodes
        Eigen::VectorXd CVc = xnode->getCVc_Node();
        Eigen::MatrixXd  Vx = xnode->getV_Node();
        Eigen::VectorXd f_seg = Vx * CVc;
        for(int i=0; i<xnode->getPoints().size(); i++){
            f[xnode->getPoints()[i].getId()] += f_seg[i]; // add contribution of far field to ``f''
        }
    }
}

// compute far field contribution
template<>
void LowRankApp<BlockCluster_Y,Node_Y>::ff_contribution(std::vector<BlockCluster*> ff_v,
                                                        std::vector<Node*> ff_v_x, std::vector<Node*> ff_v_y,
                                                        const Eigen::VectorXd& c, Eigen::VectorXd& f)
{
    preProcess(ff_v_x, ff_v_y, c);
    blockProcess(ff_v);
    postProcess(ff_v_x, f);
}

// compute near-field contribution
template<>
void LowRankApp<BlockCluster_Y,Node_Y>::nf_contribution(std::vector<BlockNearF*> nf_v,
                                                        const Eigen::VectorXd& c, Eigen::VectorXd& f)
{
    for(auto& pair : nf_v){ // iterate for all the near field xnodes
        Node* xnode = pair->getXNode();
        Node* ynode = pair->getYNode();
        pair->setKernel(kernel_);
        pair->setMatrix();
        Eigen::MatrixXd C = pair->getMatrix();
        for(int i=0; i<xnode->getPoints().size(); i++){
            for(int j=0; j<ynode->getPoints().size(); j++){
                f(xnode->getPoints()[i].getId()) += C(i,j) * c(ynode->getPoints()[j].getId()); // add near field contribution to ``f''
            }
        }
    }
}

// approximate matrix-vector multiplication
template<>
Eigen::VectorXd LowRankApp<BlockCluster_Y,Node_Y>::mvProd(const Eigen::VectorXd& c)
{
    // compute Far and Near Field relationships between Nodes of the Cluster Tree
    HP_.setNearFar();
    // compute far field contribution
    Eigen::VectorXd f_approx = Eigen::VectorXd::Zero(c.size());
    ff_contribution(HP_.getFF(),
                    HP_.getFFxnds(), HP_.getFFynds(),
                    c, f_approx);
    // compute near-field contribution
    nf_contribution(HP_.getNF(),
                    c, f_approx);

    return f_approx;
}
