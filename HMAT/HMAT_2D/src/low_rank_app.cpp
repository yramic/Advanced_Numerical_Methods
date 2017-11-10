#include "../include/low_rank_app.hpp"
#include "../include/block_cluster.hpp"
#include "../include/ctree.hpp"
#include "../include/kernel.hpp"
#include "../include/node.hpp"
#include <iostream>
#include <Eigen/Dense>

// constructor for solving the 2D problem
LowRankApp::LowRankApp(Kernel* kernel, const std::vector<Point> &pp, double eta, unsigned deg):
    kernel_(kernel), HP_(pp,eta,deg), deg_(deg), nops_(0)
{ }

// approximate matrix-vector multiplication
Eigen::VectorXd LowRankApp::mvProd(Eigen::VectorXd& c)
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
    //omp_set_num_threads(8);
    //int num_threads = omp_get_num_threads();
    //std::cout << num_threads << std::endl;
#pragma omp parallel
{
    /*for(auto& xnode : ff_v_x){ // iterate for all the far field xnodes
        nops_ += xnode->setV();
    }*/
    unsigned lsum = 0;
    #pragma omp for
    for(int i = 0; i < ff_v_x.size(); i++){
        Node* xnode = ff_v_x[i];
//#pragma omp atomic
        lsum += xnode->setV();
    }
#pragma omp atomic
    nops_ += lsum;

    /*for(auto& ynode : ff_v_y){ // iterate for all the far field ynodes
        nops_ += ynode->setV();
        nops_ += ynode->setVc(c);
    }*/
    lsum = 0;
    #pragma omp for
    for(int i = 0; i < ff_v_y.size(); i++){
        Node* ynode = ff_v_y[i];
//#pragma omp atomic
        lsum += ynode->setV();
//#pragma omp atomic
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
    /*for(auto& pair : ff_v){ // iterate for all the pairs of far field nodes
        nops_ += pair->setMatrix(kernel_); // here because needed for each pair of nodes,
                                           // cannot be moved to pre-processing
        nops_ += pair->setCVc();
    }*/
#pragma omp parallel
    {
        unsigned lsum = 0;
#pragma omp for
        for(int i = 0; i < ff_v.size(); i++){
            BlockCluster* pair = ff_v[i];
//#pragma omp atomic
            lsum += pair->setMatrix(kernel_); // here because needed for each pair of nodes,
                                               // cannot be moved to pre-processing
//#pragma omp atomic
            lsum += pair->setCVc();
        }
#pragma omp atomic
        nops_ += lsum;
    }

}

// post-processing: compute vector Vx*CVc for all far field xnodes and add it to vector f in the right place
void LowRankApp::postProcess(std::vector<Node*> ff_v_x, Eigen::VectorXd& f)
{
    /*for(auto& xnode : ff_v_x){ // iterate for all the far field xnodes
        Eigen::VectorXd CVc = xnode->getCVc_Node();
        Eigen::MatrixXd  Vx = xnode->getV_node();
        Eigen::VectorXd f_seg = Vx * CVc;
        nops_ += Vx.rows()*Vx.cols();
        for(int i=0; i<xnode->getPPoints().size(); i++){
            f[xnode->getPPoints()[i].getId()] += f_seg[i]; // add contribution of far field to ``f''
        }
    }*/
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
            for(int i=0; i<xnode->getPPoints().size(); i++){
                f[xnode->getPPoints()[i].getId()] += f_seg[i]; // add contribution of far field to ``f''
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
        for(int i=0; i<xnode->getPPoints().size(); i++){
            f_approx_ff_contr(xnode->getPPoints()[i].getId()) += ynode->getPPoints().size();
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
        for(int i=0; i<xnode->getPPoints().size(); i++){
            for(int j=0; j<ynode->getPPoints().size(); j++){
                f(xnode->getPPoints()[i].getId()) += C(i,j) * c(ynode->getPPoints()[j].getId()); // add near field contribution to ``f''
                ++f_approx_nf_contr(xnode->getPPoints()[i].getId());
                nops_ += C.rows()*C.cols();
            }
        }
    }
}
