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
    std::cout << "Far Field Contribution for each row" << std::endl << std::flush;
    //std:: cout << f_approx_ff_contr << std::endl << std::flush;
    // compute near-field contribution
    nf_contribution(f_approx, HP_.getNF(), c, f_approx_nf_contr);

    std::cout << "Near Field Contribution for each row" << std::endl << std::flush;
    //std:: cout << f_approx_nf_contr << std::endl << std::flush;
    return f_approx;
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
                f_approx_nf_contr(xnode->getPPoints()[j].getId())++;
            }
        }
    }
}
