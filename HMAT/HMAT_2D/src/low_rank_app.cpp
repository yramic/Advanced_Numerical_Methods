#include "../include/low_rank_app.hpp"
#include "../include/block_cluster.hpp"
#include "../include/ctree.hpp"
#include "../include/kernel.hpp"
#include "../include/node.hpp"
#include <iostream>
#include <Eigen/Dense>

// constructor for solving the 2D problem
LowRankApp::LowRankApp(Kernel* kernel, std::vector<Point> pp, int n):
    kernel_(kernel), PPointsTree_(pp)
{ }


// approximate matrix-vector multiplication
Eigen::VectorXd LowRankApp::mvProd(Eigen::VectorXd& c, double eta, unsigned deg)
{

    // compute V-matrices
    PPointsTree_.setV(deg);

    // add pointers to near and far field nodes of the tree
    int near = 0, far = 0;
    PPointsTree_.setNearFar(eta, PPointsTree_, near, far);     // for the 2D problem
    std::cout << "Near Field Nodes: " << near << " Far Field Nodes: " << far << std::endl;
    std::cout << "Near Field Nodes: " << (double)near/(double)(near+far)*100 << "% " << "Far Field Nodes: " << (double)far/(double)(near+far)*100 << "%" << std::endl;
    //PPointsTree_.getRoot()->printree(0);          // printing the tree for testing
    Eigen::VectorXd f_approx_ff_contr = Eigen::VectorXd::Zero(c.size());
    Eigen::VectorXd f_approx_nf_contr = Eigen::VectorXd::Zero(c.size());
    // compute far field contribution
    Eigen::VectorXd f_approx = Eigen::VectorXd::Zero(c.size());
    ff_contribution(f_approx, PPointsTree_.getRoot(), deg, c, f_approx_ff_contr);
    // compute near-field contribution
    nf_contribution(f_approx, PPointsTree_.getRoot(), c, f_approx_nf_contr);
    //std::cout << "Far Field Contribution for each row" << std::endl;
    //std:: cout << f_approx_ff_contr << std::endl;
    //std::cout << "Near Field Contribution for each row" << std::endl;
    //std:: cout << f_approx_nf_contr << std::endl;
    return f_approx;
}


// compute far field contribution
void LowRankApp::ff_contribution(Eigen::VectorXd& f, Node* tx, unsigned deg, Eigen::VectorXd& c, Eigen::VectorXd& f_approx_ff_contr)
{
    if(tx != NULL && (*tx).getPPoints().size()>1){
        if(!(*tx).getFarF().empty()){       // probably it isn´t needed
            Eigen::VectorXd XVc = Eigen::VectorXd::Zero((deg +1)*(deg+1));                                          // auxiliary variable
            std::vector<Node*> ffx = (*tx).getFarF();                                                               // far field of *tx_root
            for(std::vector<Node*>::iterator iter=ffx.begin(); iter!=ffx.end(); ++iter) {                           // iterating on the farfield nodes of the cluster tx
                Eigen::VectorXd Vc((deg+1)*(deg+1));   // deg+1*deg+1                                                      // V*c restricted to the indices of **iter
                BlockCluster X_(tx->getXl_b(), tx->getXr_b(), tx->getYl_b(), tx->getYr_b(), (*iter)->getXl_b(), (*iter)->getXr_b(), (*iter)->getYl_b(), (*iter)->getYr_b(), deg, kernel_);    // calculation of matrix $X_{\sigma,\mu}$
                Eigen::MatrixXd X = X_.getMatrix();                                                                 // matrix $X_{\sigma,\mu}$

                Eigen::MatrixXd Vm = (*iter)->getV_node().transpose();

                Eigen::VectorXd c_seg((*iter)->getPPoints().size());
                for (int i = 0; i<(*iter)->getPPoints().size(); i++) {
                    c_seg[i] = c[(*iter)->getPPoints()[i].getId()];
                    for (int i = 0; i<tx->getPPoints().size(); i++) {
                        f_approx_ff_contr[tx->getPPoints()[i].getId()]++;
                    }
                }
                Vc = Vm * c_seg;                                                                                    // computation of V*cm
                XVc += X * Vc;                                                                                      // add contribution of block **iter to "s"
            }
            Eigen::MatrixXd Vs = Eigen::MatrixXd::Zero(tx->getPPoints().size(),(deg+1)*(deg+1));                    // $V_{\sigma}$
            Vs = (*tx).getV_node();
            Eigen::VectorXd f_seg((*tx).getPPoints().size());                                                       // add contribution of far field to "f"
            f_seg = Vs * XVc;
            for (int i = 0; i<tx->getPPoints().size(); i++) {
                f[tx->getPPoints()[i].getId()] += f_seg[i];
            }
        }
        // add contribution of leaves of *tx
        ff_contribution(f, tx->getTl_Child(), deg, c, f_approx_ff_contr);
        ff_contribution(f, tx->getTr_Child(), deg, c, f_approx_ff_contr);
        ff_contribution(f, tx->getBl_Child(), deg, c, f_approx_ff_contr);
        ff_contribution(f, tx->getBr_Child(), deg, c, f_approx_ff_contr);
    }
}


// compute near-field contribution
void LowRankApp::nf_contribution(Eigen::VectorXd& f, Node* tx, const Eigen::VectorXd& c, Eigen::VectorXd& f_approx_nf_contr)
{
    if(tx != NULL) {
        // iterate over near field of *tx
            std::vector<Node*> nfx = (*tx).getNearF(); // near field of *tx
            for(std::vector<Node*>::iterator iter=nfx.begin(); iter!=nfx.end(); ++iter) {
                for(int j=0; j<(*tx).getPPoints().size(); j++){
                    for(int i=0; i<(*iter)->getPPoints().size(); i++){
                        std::vector<Point> t = (*tx).getPPoints();
                        f(t[j].getId()) += (*kernel_)(t[j].getX(),t[j].getY(),(*iter)->getPPoints()[i].getX(),(*iter)->getPPoints()[i].getY()) * c((*iter)->getPPoints()[i].getId());   // add near field contribution
                        f_approx_nf_contr[tx->getPPoints()[j].getId()]++;
                        //std::cout << t[j].getId() << " " << (*kernel_)(t[j].getX(),t[j].getY(),(*iter)->getPPoints()[(i)].getX(),(*iter)->getPPoints()[(i)].getY()) << " " << (*iter)->getPPoints()[(i)].getId() << std::endl;
                    }
                }
            }
        nf_contribution(f, tx->getTl_Child(), c, f_approx_nf_contr);
        nf_contribution(f, tx->getTr_Child(), c, f_approx_nf_contr);
        nf_contribution(f, tx->getBl_Child(), c, f_approx_nf_contr);
        nf_contribution(f, tx->getBr_Child(), c, f_approx_nf_contr);
    }
}
