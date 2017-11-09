#include "../include/low_rank_app.hpp"
#include "../include/block_cluster.hpp"
#include "../include/ctree.hpp"
#include "../include/kernel.hpp"
#include "../include/node.hpp"
#include <iostream>
#include <Eigen/Dense>

// constructor for solving the 2D problem
LowRankApp::LowRankApp(Kernel* kernel, const std::vector<Point> &pp, double eta, unsigned deg):
    kernel_(kernel), HP_(pp,eta,deg), deg_(deg)
{ }

// approximate matrix-vector multiplication
Eigen::VectorXd LowRankApp::mvProd(Eigen::VectorXd& c, double eta, unsigned deg)
{
    // Compute the Near and Far Field Pairs
    HP_.setNearFar();
    // Number of Near and Far Field Pairs
    int near = HP_.getNF().size(), far = HP_.getFF().size();
    std::cout << "Near Field Nodes: " << near << " Far Field Nodes: " << far << std::endl;
    std::cout << "Near Field Nodes: " << (double)near/(double)(near+far)*100 << "% " << "Far Field Nodes: " << (double)far/(double)(near+far)*100 << "%" << std::endl;
    //PPointsTree_.getRoot()->printree(0);          // printing the tree for testing
    Eigen::VectorXd f_approx_ff_contr = Eigen::VectorXd::Zero(c.size());
    Eigen::VectorXd f_approx_nf_contr = Eigen::VectorXd::Zero(c.size());
    // compute far field contribution
    Eigen::VectorXd f_approx = Eigen::VectorXd::Zero(c.size());
    ff_contribution(HP_.getFF(),
                    HP_.getFFxnds(), HP_.getFFynds(),
                    c, f_approx, f_approx_ff_contr);
    std::cout << "Far Field Contribution for each row" << std::endl << std::flush;
    std:: cout << f_approx_ff_contr << std::endl << std::flush;
    // compute near-field contribution
    nf_contribution(HP_.getNF(),
                    c, f_approx, f_approx_nf_contr);

    std::cout << "Near Field Contribution for each row" << std::endl << std::flush;
    std:: cout << f_approx_nf_contr << std::endl << std::flush;
    return f_approx;
}

// pre-processing: initialize matrix V and vector Vc for all far field nodes
// do these steps once for each node, not every time the node appears in a pair
void LowRankApp::preProcess(std::vector<Node*> ff_v_x, std::vector<Node*> ff_v_y, const Eigen::VectorXd& c)
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
void LowRankApp::blockProcess(std::vector<BlockCluster> ff_v)
{
    for(auto& pair : ff_v){ // iterate for all the pairs of far field nodes
        pair.setCVc(kernel_);
    }
}

// post-processing: compute vector Vx*CVc for all far field xnodes and add it to vector f in the right place
void LowRankApp::postProcess(std::vector<Node*> ff_v_x, Eigen::VectorXd& f)
{
    for(auto& xnode : ff_v_x){ // iterate for all the far field xnodes
        Eigen::VectorXd CVc = xnode->getCVc_Node();
        Eigen::MatrixXd  Vx = xnode->getV_node();
        Eigen::VectorXd f_seg = Vx * CVc;
        for(int i=0; i<xnode->getPPoints().size(); i++){
            f[xnode->getPPoints()[i].getId()] += f_seg[i]; // add contribution of far field to ``f''
        }
    }
}

// calculate the number of y points that are used for the far field contribution of each row of the product vector
void LowRankApp::calc_numb_approx_per_row(std::vector<BlockCluster> ff_v, Eigen::VectorXd& f_approx_ff_contr)
{
    for(auto& block : ff_v){ // iterate for all the far field xnodes
        Node* xnode = block.getXNode();
        Node* ynode = block.getYNode();
        for(int i=0; i<xnode->getPPoints().size(); i++){
            f_approx_ff_contr[xnode->getPPoints()[i].getId()] += ynode->getPPoints().size(); // add contribution of far field to ``f''
        }
    }
}

// compute far field contribution
void LowRankApp::ff_contribution(std::vector<BlockCluster> ff_v,
                                 std::vector<Node*> ff_v_x, std::vector<Node*> ff_v_y,
                                 const Eigen::VectorXd& c, Eigen::VectorXd& f, Eigen::VectorXd& f_approx_ff_contr)
{
    preProcess(ff_v_x, ff_v_y, c);
    blockProcess(ff_v);
    postProcess(ff_v_x, f);
    calc_numb_approx_per_row(ff_v, f_approx_ff_contr);
}

// compute near-field contribution
void LowRankApp::nf_contribution(std::vector<BlockNearF> nf_v,
                                 const Eigen::VectorXd& c, Eigen::VectorXd& f, Eigen::VectorXd& f_approx_nf_contr)
{
    for(auto& pair : nf_v){ // iterate for all the near field xnodes
        Node* xnode = pair.getXNode();
        Node* ynode = pair.getYNode();
        //pair.setKernel(kernel_);
        pair.setMatrix(kernel_);
        Eigen::MatrixXd C = pair.getMatrix();
        for(int i=0; i<xnode->getPPoints().size(); i++){
            for(int j=0; j<ynode->getPPoints().size(); j++){
                f(xnode->getPPoints()[i].getId()) += C(i,j) * c(ynode->getPPoints()[j].getId()); // add near field contribution to ``f''
                f_approx_nf_contr(xnode->getPPoints()[i].getId())++; // calculate the number of y points that are used for the near field contribution of each row of the product vector

            }
        }
    }
}
