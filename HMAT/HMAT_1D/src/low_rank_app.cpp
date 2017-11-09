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
#include "../include/block_nearf.hpp"
#include "../include/ctree.hpp"
#include "../include/kernel.hpp"
#include "../include/node.hpp"
#include "../include/point.hpp"
#include <iostream>

// constructor
template<>
LowRankApp<BlockCluster,Node>::LowRankApp(Kernel* kernel, const std::vector<Point>& GPoints, double eta, unsigned deg):
    kernel_(kernel), GPoints_(GPoints), HP_(GPoints,eta,deg), deg_(deg)
{ }

// pre-processing: initialize matrix V and vector Vc for all far field nodes
// do these steps once for each node, not every time the node appears in a pair
template<>
void LowRankApp<BlockCluster,Node>::preProcess(std::vector<Node*> ff_v_x, std::vector<Node*> ff_v_y, const Eigen::VectorXd& c)
{
    for(auto& xnode : ff_v_x){ // iterate for all the far field xnodes
        xnode->setV();
    }
    for(auto& ynode : ff_v_y){ // iterate for all the far field ynodes
        ynode->setV();
        ynode->setVc(c);
    }
}

// block-processing: compute vector CVc for all far field pairs and store it into xnode
// all vectors CVc of an xnode can already be summed together
template<>
void LowRankApp<BlockCluster,Node>::blockProcess(std::vector<BlockCluster*> ff_v)
{
    for(auto& pair : ff_v){ // iterate for all the pairs of far field nodes
        pair->setMatrix(kernel_); // here because needed for each pair of nodes,
                                  // cannot be moved to pre-processing
        pair->setCVc();
    }
}

// post-processing: compute vector Vx*CVc for all far field xnodes and add it to vector f in the right place
template<>
void LowRankApp<BlockCluster,Node>::postProcess(std::vector<Node*> ff_v_x, Eigen::VectorXd& f)
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

// count far-field ynodes contributing to each row of the approximate low-rank matrix
void LowRankApp::calc_numb_approx_per_row(std::vector<BlockCluster*> ff_v, Eigen::VectorXd& f_approx_ff_contr)
{
    for(auto& pair : ff_v){ // iterate for all the pairs of far field nodes
        Node* xnode = pair->getXNode();
        Node* ynode = pair->getYNode();
        for(int i=0; i<xnode->getPPoints().size(); i++){
            f_approx_ff_contr(xnode->getPPoints()[i].getId()) += ynode->getPPoints().size();
        }
    }
}

// compute far field contribution
template<>
void LowRankApp<BlockCluster,Node>::ff_contribution(std::vector<BlockCluster*> ff_v,
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
void LowRankApp<BlockCluster,Node>::nf_contribution(std::vector<BlockNearF*> nf_v,
                                                    const Eigen::VectorXd& c, Eigen::VectorXd& f, Eigen::VectorXd& f_approx_nf_contr)
{
    for(auto& pair : nf_v){ // iterate for all the near field xnodes
        Node* xnode = pair->getXNode();
        Node* ynode = pair->getYNode();
        pair->setMatrix(kernel_);
        Eigen::MatrixXd C = pair->getMatrix();
        for(int i=0; i<xnode->getPoints().size(); i++){
            for(int j=0; j<ynode->getPoints().size(); j++){
                f(xnode->getPoints()[i].getId()) += C(i,j) * c(ynode->getPoints()[j].getId()); // add near field contribution to ``f''
                ++f_approx_nf_contr(xnode->getPPoints()[i].getId());
            }
        }
    }
}

// approximate matrix-vector multiplication
template<>
Eigen::VectorXd LowRankApp<BlockCluster,Node>::mvProd(const Eigen::VectorXd& c)
{
    // compute the Near and Far Field Pairs
    HP_.setNearFar();
    // number of Near and Far Field Pairs
    int near = HP_.getNF().size(), far = HP_.getFF().size();
    std::cout << "Near Field Nodes: " << near << " Far Field Nodes: " << far << std::endl;
    std::cout << "Near Field Nodes: " << (double)near/(near+far)*100. << "% " << "Far Field Nodes: " << (double)far/(near+far)*100. << "%" << std::endl;

    // compute far field contribution
    Eigen::VectorXd f_approx = Eigen::VectorXd::Zero(c.size());
    ff_contribution(HP_.getFF(),
                    HP_.getFFxnds(), HP_.getFFynds(),
                    c, f_approx);
    std::cout << "Far Field Contribution for each row" << std::endl;
    std::cout << f_approx_ff_contr << std::endl;

    // compute near-field contribution
    nf_contribution(HP_.getNF(),
                    c, f_approx);
    std::cout << "Near Field Contribution for each row" << std::endl;
    std::cout << f_approx_nf_contr << std::endl;

    return f_approx;
}
