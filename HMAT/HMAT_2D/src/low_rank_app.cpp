#include "../include/low_rank_app.hpp"
#include "../include/block_cluster.hpp"
#include "../include/ctree.hpp"
#include "../include/kernel.hpp"
#include "../include/node.hpp"
#include <iostream>
#include <Eigen/Dense>

// constructor for solving the 2D problem
LowRankApp::LowRankApp(Kernel* kernel, const std::vector<Point>& pp, double eta, unsigned deg):
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
    //std:: cout << f_approx_ff_contr << std::endl << std::flush;
    // compute near-field contribution
    nf_contribution(HP_.getNF(),
                    c, f_approx, f_approx_nf_contr);

    std::cout << "Near Field Contribution for each row" << std::endl << std::flush;
    //std:: cout << f_approx_nf_contr << std::endl << std::flush;
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
void LowRankApp::blockProcess(std::vector<BlockCluster*> ff_v)
{
    for(auto& pair : ff_v){ // iterate for all the pairs of far field nodes
        pair->setMatrix(kernel_); // here because needed for each pair of nodes,
                                  // cannot be moved to pre-processing
        pair->setCVc();
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

/*// compute far field contribution
void LowRankApp::ff_contribution(Eigen::VectorXd& f, std::vector<std::pair<Node*,Node*>> ff_v, unsigned deg, Eigen::VectorXd& c, Eigen::VectorXd& f_approx_ff_contr)
{
    int n = ff_v.size();
    for(int i = 0; i<n; i++){   // iterate for all the pairs of far field nodes
        Node* xnode = ff_v[i].first;
        Node* ynode = ff_v[i].second;
        Eigen::VectorXd XVc = Eigen::VectorXd::Zero((deg +1)*(deg+1));                                          // auxiliary variable
        Eigen::VectorXd Vc((deg+1)*(deg+1));   // deg+1*deg+1                                                   // V*c restricted to the indices of **iter
        BlockCluster X_(xnode->getTkx(), xnode->getTky(), ynode->getTkx(), ynode->getTky(), deg, kernel_);    // calculation of matrix $X_{\sigma,\mu}$
        Eigen::MatrixXd X = X_.getMatrix();
        xnode->setV_node(deg);
        ynode->setV_node(deg);
        Eigen::MatrixXd Vy = ynode->getV_node();
        Eigen::MatrixXd Vm = Vy.transpose();
        int ny = ynode->getPPoints().size();
        Eigen::VectorXd c_seg(ny);
        for(int j=0; j<ny; j++){
            c_seg[j] = c[ynode->getPPoints()[j].getId()];
            for (int k = 0; k<xnode->getPPoints().size(); k++) {
                f_approx_ff_contr[xnode->getPPoints()[k].getId()]++;
            }
        }
        Vc = Vm * c_seg;
        XVc += X * Vc;
        Eigen::MatrixXd Vs = Eigen::MatrixXd::Zero(xnode->getPPoints().size(),(deg+1)*(deg+1));                    // $V_{\sigma}$
        Vs = xnode->getV_node();
        Eigen::VectorXd f_seg(xnode->getPPoints().size());                                                       // add contribution of far field to "f"
        f_seg = Vs * XVc;
        for (int j = 0; j<xnode->getPPoints().size(); j++) {
            f[xnode->getPPoints()[j].getId()] += f_seg[j];
        }
    }

}*/

// compute far field contribution
void LowRankApp::ff_contribution(std::vector<BlockCluster*> ff_v,
                                 std::vector<Node*> ff_v_x, std::vector<Node*> ff_v_y,
                                 const Eigen::VectorXd& c, Eigen::VectorXd& f, Eigen::VectorXd& f_approx_ff_contr)
{
    preProcess(ff_v_x, ff_v_y, c);
    blockProcess(ff_v);
    postProcess(ff_v_x, f);
}

// compute near-field contribution
void LowRankApp::nf_contribution(std::vector<BlockNearF*> nf_v,
                                 const Eigen::VectorXd& c, Eigen::VectorXd& f, Eigen::VectorXd& f_approx_nf_contr)
{
    for(auto& pair : nf_v){ // iterate for all the near field xnodes
        Node* xnode = pair->getXNode();
        Node* ynode = pair->getYNode();
        pair->setMatrix(kernel_);
        Eigen::MatrixXd C = pair->getMatrix();
        for(int i=0; i<xnode->getPPoints().size(); i++){
            for(int j=0; j<ynode->getPPoints().size(); j++){
                f(xnode->getPPoints()[i].getId()) += C(i,j) * c(ynode->getPPoints()[j].getId()); // add near field contribution to ``f''
            }
        }
    }
}
