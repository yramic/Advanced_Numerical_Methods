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
#include "../include/low_rank_app.hpp"
#include "../include/block_cluster.hpp"
#include "../include/ctree.hpp"
#include "../include/kernel.hpp"
#include "../include/node.hpp"
#include <Eigen/Dense>
#include <iostream>

// constructor for solving the 2D problem
LowRankApp::LowRankApp(Kernel* kernel, const std::vector<Segment> &pp, double eta, unsigned deg):
    kernel_(kernel), HP_(pp,eta,deg), deg_(deg), nops_(0)
{
    assert(deg < 32);
}

// approximate matrix-vector multiplication
Eigen::VectorXd LowRankApp::mvProd(const Eigen::VectorXd& c)
{
    nops_ = 0;

    // compute the Near and Far Field Pairs
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
    //PPointsTree_.getRoot()->printree(0); // printing the tree for testing

    std::cout << "Number of matrix operations performed for low-rank approximation: " << nops_ << std::endl;

    return f_approx;
}

// pre-processing: initialize matrix V and vector Vc for all far field nodes
// do these steps once for each node, not every time the node appears in a pair
void LowRankApp::preProcess(std::vector<Node*> ff_v_x, std::vector<Node*> ff_v_y, const Eigen::VectorXd& c)
{
    #pragma omp parallel
    {
        unsigned lsum = 0;
        #pragma omp for
        for(int i = 0; i < ff_v_x.size(); i++){
            Node* xnode = ff_v_x[i];
            lsum += xnode->setV();
        }
        #pragma omp atomic
        nops_ += lsum;

        lsum = 0;
        #pragma omp for
        for(int i = 0; i < ff_v_y.size(); i++){
            Node* ynode = ff_v_y[i];
            lsum += ynode->setV();
            lsum += ynode->setVc(c);
        }
        #pragma omp atomic
        nops_ += lsum;
    }
}

// block-processing: compute vector CVc for all far field pairs and store it into xnode
// all vectors CVc of an xnode can already be summed together
void LowRankApp::blockProcess(std::vector<BlockCluster*> ff_v)
{
    #pragma omp parallel
    {
        unsigned lsum = 0;
        #pragma omp for
        for(int i = 0; i < ff_v.size(); i++){
            BlockCluster* pair = ff_v[i];
            lsum += pair->setMatrix(kernel_); // here because needed for each pair of nodes,
                                              // cannot be moved to pre-processing
            lsum += pair->setCVc();
        }
        #pragma omp atomic
        nops_ += lsum;
    }
}

// post-processing: compute vector Vx*CVc for all far field xnodes and add it to vector f in the right place
void LowRankApp::postProcess(std::vector<Node*> ff_v_x, Eigen::VectorXd& f)
{
    #pragma omp parallel
    {
        unsigned lsum = 0;
        #pragma omp for
        for(int i = 0; i < ff_v_x.size(); i++){
            Node* xnode = ff_v_x[i];
            Eigen::VectorXd CVc = xnode->getCVc_Node();
            Eigen::MatrixXd  Vx = xnode->getV_node();
            Eigen::VectorXd f_seg = Vx * CVc;
            lsum += Vx.rows()*Vx.cols();
            for(int i=0; i<xnode->getSegments().size(); i++){
                f[xnode->getSegments()[i].getId()] += f_seg[i]; // add contribution of far field to ``f''
            }
        }
        #pragma omp atomic
        nops_ += lsum;
    }
}

// count far-field ynodes contributing to each row of the approximate low-rank matrix
void LowRankApp::calc_numb_approx_per_row(std::vector<BlockCluster*> ff_v, Eigen::VectorXd& f_approx_ff_contr)
{
    for(auto& pair : ff_v){ // iterate for all the pairs of far field nodes
        Node* xnode = pair->getXNode();
        Node* ynode = pair->getYNode();
        for(int i=0; i<xnode->getSegments().size(); i++){
            f_approx_ff_contr(xnode->getSegments()[i].getId()) += ynode->getSegments().size();
        }
    }
}

// compute far field contribution
void LowRankApp::ff_contribution(std::vector<BlockCluster*> ff_v,
                                 std::vector<Node*> ff_v_x, std::vector<Node*> ff_v_y,
                                 const Eigen::VectorXd& c, Eigen::VectorXd& f, Eigen::VectorXd& f_approx_ff_contr)
{
    preProcess(ff_v_x, ff_v_y, c);
    blockProcess(ff_v);
    postProcess(ff_v_x, f);
    calc_numb_approx_per_row(ff_v, f_approx_ff_contr);
}

// compute near-field contribution
void LowRankApp::nf_contribution(std::vector<BlockNearF*> nf_v,
                                 const Eigen::VectorXd& c, Eigen::VectorXd& f, Eigen::VectorXd& f_approx_nf_contr)
{
    for(auto& pair : nf_v){ // iterate for all the near field xnodes
        Node* xnode = pair->getXNode();
        Node* ynode = pair->getYNode();
        nops_ += pair->setMatrix(kernel_);
        Eigen::MatrixXd C = pair->getMatrix();
        for(int i=0; i<xnode->getSegments().size(); i++){
            for(int j=0; j<ynode->getSegments().size(); j++){
                f(xnode->getSegments()[i].getId()) += C(i,j) * c(ynode->getSegments()[j].getId()); // add near field contribution to ``f''
                ++f_approx_nf_contr(xnode->getSegments()[i].getId());
                nops_ += C.rows()*C.cols();
            }
        }
    }
}
