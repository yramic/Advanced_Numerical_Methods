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
LowRankApp<BlockCluster_Y,Node_Y>::LowRankApp(Kernel* kernel, const std::vector<Point>& GPoints, double eta, unsigned deg):
    kernel_(kernel), GPoints_(GPoints), HP_(GPoints,eta,deg), deg_(deg), nops_(0)
{ }

// pre-processing: initialize matrix V and vector Vc for all far field nodes
// do these steps once for each node, not every time the node appears in a pair
template<>
void LowRankApp<BlockCluster_Y,Node_Y>::preProcess(std::vector<Node*> ff_v_x, std::vector<Node*> ff_v_y, const Eigen::VectorXd& c)
{
    for(auto& xnode : ff_v_x){ // iterate for all the far field xnodes
        nops_ += xnode->setV();
    }
    for(auto& ynode : ff_v_y){ // iterate for all the far field ynodes
        Node_Y* ynode_ = static_cast<Node_Y*>(ynode);
        nops_ += ynode_->setV();
        nops_ += ynode_->setVc(c);
    }
}

// block-processing: compute vector CVc for all far field pairs and store it into xnode
// all vectors CVc of an xnode can already be summed together
template<>
void LowRankApp<BlockCluster_Y,Node_Y>::blockProcess(std::vector<BlockCluster*> ff_v)
{
    for(auto& pair : ff_v){ // iterate for all the pairs of far field nodes
        BlockCluster_Y* pair_ = static_cast<BlockCluster_Y*>(pair);
        nops_ += pair_->setMatrix(kernel_); // here because needed for each pair of nodes,
                                            // cannot be moved to pre-processing
        nops_ += pair_->setCVc();
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
        nops_ += Vx.rows()*Vx.cols();
        for(int i=0; i<xnode->getPoints().size(); i++){
            f[xnode->getPoints()[i].getId()] += f_seg[i]; // add contribution of far field to ``f''
        }
    }
}

// count far-field ynodes contributing to each row of the approximate low-rank matrix
template<>
void LowRankApp<BlockCluster_Y,Node_Y>::calc_numb_approx_per_row(std::vector<BlockCluster*> ff_v, Eigen::VectorXd& f_approx_ff_contr)
{
    for(auto& pair : ff_v){ // iterate for all the pairs of far field nodes
        Node* xnode = pair->getXNode();
        Node* ynode = pair->getYNode();
        for(int i=0; i<xnode->getPoints().size(); i++){
            f_approx_ff_contr(xnode->getPoints()[i].getId()) += ynode->getPoints().size();
        }
    }
}

// compute far field contribution
template<>
void LowRankApp<BlockCluster_Y,Node_Y>::ff_contribution(std::vector<BlockCluster*> ff_v,
                                                        std::vector<Node*> ff_v_x, std::vector<Node*> ff_v_y,
                                                        const Eigen::VectorXd& c, Eigen::VectorXd& f, Eigen::VectorXd& f_approx_ff_contr)
{
    preProcess(ff_v_x, ff_v_y, c);
    blockProcess(ff_v);
    postProcess(ff_v_x, f);
    calc_numb_approx_per_row(ff_v, f_approx_ff_contr);
}

// compute near-field contribution
template<>
void LowRankApp<BlockCluster_Y,Node_Y>::nf_contribution(std::vector<BlockNearF*> nf_v,
                                                        const Eigen::VectorXd& c, Eigen::VectorXd& f, Eigen::VectorXd& f_approx_nf_contr)
{
    for(auto& pair : nf_v){ // iterate for all the near field xnodes
        Node* xnode = pair->getXNode();
        Node* ynode = pair->getYNode();
        nops_ += pair->setMatrix(kernel_);
        Eigen::MatrixXd C = pair->getMatrix();
        for(int i=0; i<xnode->getPoints().size(); i++){
            for(int j=0; j<ynode->getPoints().size(); j++){
                f(xnode->getPoints()[i].getId()) += C(i,j) * c(ynode->getPoints()[j].getId()); // add near field contribution to ``f''
                ++f_approx_nf_contr(xnode->getPoints()[i].getId());
                nops_ += C.rows()*C.cols();
            }
        }
    }
}

// approximate matrix-vector multiplication
template<>
Eigen::VectorXd LowRankApp<BlockCluster_Y,Node_Y>::mvProd(const Eigen::VectorXd& c)
{
    nops_ = 0;

    // compute Far and Near Field relationships between Nodes of the Cluster Tree
    HP_.setNearFar();

    Eigen::VectorXd f_approx_ff_contr = Eigen::VectorXd::Zero(c.size());
    Eigen::VectorXd f_approx_nf_contr = Eigen::VectorXd::Zero(c.size());
    // compute far field contribution
    Eigen::VectorXd f_approx = Eigen::VectorXd::Zero(c.size());
    ff_contribution(HP_.getFF(),
                    HP_.getFFxnds(), HP_.getFFynds(),
                    c, f_approx, f_approx_ff_contr);
    std::cout << "Far Field Contribution for each row" << std::endl;
    std::cout << f_approx_ff_contr << std::endl;

    // compute near-field contribution
    nf_contribution(HP_.getNF(),
                    c, f_approx, f_approx_nf_contr);
    std::cout << "Near Field Contribution for each row" << std::endl;
    std::cout << f_approx_nf_contr << std::endl;

    // number of Near and Far Field Pairs
    int near = HP_.getNF().size(), far = HP_.getFF().size();
    std::cout << "Near Field Nodes: " << near << " Far Field Nodes: " << far << std::endl;
    std::cout << "Near Field Nodes: " << (double)near/(near+far)*100. << "% " << "Far Field Nodes: " << (double)far/(near+far)*100. << "%" << std::endl;

    std::cout << "Number of matrix operations performed for low-rank approximation: " << nops_ << std::endl;

    return f_approx;
}
