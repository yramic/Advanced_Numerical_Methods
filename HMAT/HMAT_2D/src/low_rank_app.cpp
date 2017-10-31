#include "../include/low_rank_app.hpp"
#include "../include/block_cluster.hpp"
#include "../include/ctree.hpp"
#include "../include/kernel.hpp"
#include "../include/node.hpp"
#include <iostream>
#include <Eigen/Dense>

// constructor for solving the 2D problem
LowRankApp::LowRankApp(Kernel* kernel, const std::vector<Point> &pp, int n, double eta, unsigned deg):
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
    ff_contribution(f_approx, HP_.getFF(), deg, c, f_approx_ff_contr);
    // compute near-field contribution
    nf_contribution(f_approx, HP_.getNF(), c, f_approx_nf_contr);
    //std::cout << "Far Field Contribution for each row" << std::endl;
    //std:: cout << f_approx_ff_contr << std::endl;
    //std::cout << "Near Field Contribution for each row" << std::endl;
    //std:: cout << f_approx_nf_contr << std::endl;
    return f_approx;
}

// compute V-matrix of node
Eigen::MatrixXd LowRankApp::setV_node(Node* x, unsigned deg) //tt==PPointsTree??
{
    std::vector<Point> xpoints = x->getPPoints();
    int ppts = xpoints.size();
    Eigen::MatrixXd Vx = Eigen::MatrixXd::Constant(ppts, (deg+1)*(deg+1), 1);

    for(unsigned i=0; i<=ppts-1; ++i) { // calculation of Vx combined with Vy
        for(unsigned j1=0; j1<=deg; ++j1) {
            for(unsigned k1=0; k1<j1; ++k1) {
                for(unsigned j2=0; j2<=deg; ++j2) {
                    Vx(i,j1*(deg+1) + j2) *= (xpoints[i].getX() - x->getTkx()[k1]);
                }
            }
            // Skip "k1 == j1"
            for(unsigned k1=j1+1; k1<=deg; ++k1) {
                for(unsigned j2=0; j2<=deg; ++j2) {
                    Vx(i,j1*(deg+1) + j2) *= (xpoints[i].getX() - x->getTkx()[k1]);
                }
            }
            for(unsigned j2=0; j2<=deg; ++j2) {
                for(unsigned k2=0; k2<j2; ++k2) {
                    Vx(i,j1*(deg+1) + j2) *= (xpoints[i].getY() - x->getTky()[k2]);
                }
                // Skip "k2 == j2"
                for(unsigned k2=j2+1; k2<=deg; ++k2) {
                    Vx(i,j1*(deg+1) + j2) *= (xpoints[i].getY() - x->getTky()[k2]);
                }
                Vx(i,j1*(deg+1) + j2) *= x->getWkx()[j1] * x->getWky()[j2];
            }
        }
    }
    return Vx;

    // Alternate way of computing V matrix
    /*Eigen::MatrixXd VnodeX = Eigen::MatrixXd::Constant(ppts, (deg+1), 1);
    Eigen::MatrixXd VnodeY = Eigen::MatrixXd::Constant(ppts, (deg+1), 1);
    for(unsigned i=0; i<=ppts-1; ++i) {
        for(unsigned j=0; j<=deg; ++j) {
            for(unsigned k=0; k<j; ++k) {
                VnodeX(i,j) *= PPointsTree_[i].getX() - tkx[k];
            }
            // Skip "k == j"
            for(unsigned k=j+1; k<=deg; ++k) {
                VnodeX(i,j) *= PPointsTree_[i].getX() - tkx[k];
            }
            VnodeX(i,j) *= wkx(j);
        }
    }
    for(unsigned i=0; i<=ppts-1; ++i) {
        for(unsigned j=0; j<=deg; ++j) {
            for(unsigned k=0; k<j; ++k) {
                VnodeY(i,j) *= PPointsTree_[i].getY() - tky[k];
            }
            // Skip "k == j"
            for(unsigned k=j+1; k<=deg; ++k) {
                VnodeY(i,j) *= PPointsTree_[i].getY() - tky[k];
            }
            VnodeY(i,j) *= wky(j);
        }
    }*/

    /*Eigen::MatrixXd V_node_new(ppts, (deg+1)*(deg+1));
    for(unsigned i=0; i<=ppts-1; ++i) {
        for(unsigned j=0; j<=deg; ++j) {
            V_node_new.block(i, j*(deg+1), 1, deg+1) = VnodeX(i,j) * VnodeY.row(i);
        }
    }*/

    /*std::cout << "VnodeX" << std::endl;
    std::cout << VnodeX << std::endl;
    std::cout << "VnodeY" << std::endl;
    std::cout << VnodeY << std::endl;
    std::cout << "V_Node" << std::endl;
    std::cout << V_node_ << std::endl;
    std::cout << "V_Node_new" << std::endl;
    std::cout << V_node_new << std::endl;*/
}

// compute far field contribution
void LowRankApp::ff_contribution(Eigen::VectorXd& f, std::vector<std::pair<Node*,Node*>> ff_v, unsigned deg, Eigen::VectorXd& c, Eigen::VectorXd& f_approx_ff_contr)
{
    int n = ff_v.size();
    for(int i = 0; i<n; i++){   // iterate for all the pairs of far field nodes
        Node* xnode = ff_v[i].first;
        Node* ynode = ff_v[i].second;
        Eigen::VectorXd XVc = Eigen::VectorXd::Zero((deg +1)*(deg+1));                                          // auxiliary variable
        Eigen::VectorXd Vc((deg+1)*(deg+1));   // deg+1*deg+1                                                   // V*c restricted to the indices of **iter
        std::cout << xnode->getTkx() << std::endl;
        BlockCluster X_(xnode->getTkx(), xnode->getTky(), ynode->getTkx(), ynode->getTky(), deg, kernel_);    // calculation of matrix $X_{\sigma,\mu}$
        Eigen::MatrixXd X = X_.getMatrix();
        Eigen::MatrixXd Vy = setV_node(ynode, deg);
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
        Vs = setV_node(xnode, deg);
        Eigen::VectorXd f_seg(xnode->getPPoints().size());                                                       // add contribution of far field to "f"
        f_seg = Vs * XVc;
        for (int j = 0; j<xnode->getPPoints().size(); j++) {
            f[xnode->getPPoints()[j].getId()] += f_seg[j];
        }
    }
}


// compute near-field contribution
void LowRankApp::nf_contribution(Eigen::VectorXd& f, std::vector<std::pair<Node*,Node*>> nf_v, const Eigen::VectorXd& c, Eigen::VectorXd& f_approx_nf_contr)
{
    int n = nf_v.size();
    for(int i = 0; i<n; i++){
        Node* xnode = nf_v[i].first;
        Node* ynode = nf_v[i].second;
        for(int j=0; j<xnode->getPPoints().size(); j++){
            for(int k=0; k<ynode->getPPoints().size(); k++){
                f(xnode->getPPoints()[j].getId()) += (*kernel_)(xnode->getPPoints()[j].getX(), xnode->getPPoints()[j].getY(), ynode->getPPoints()[k].getX(), ynode->getPPoints()[k].getY()) * c(ynode->getPPoints()[k].getId()); // add near field contribution
            }
        }
    }
}
