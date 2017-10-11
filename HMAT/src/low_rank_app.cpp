#include "../include/low_rank_app.hpp"
#include "../include/block_cluster.hpp"
#include "../include/ctree.hpp"
#include "../include/kernel.hpp"
#include "../include/node.hpp"
#include <iostream>


// constructor
LowRankApp::LowRankApp(Kernel kernel, std::vector<Point> pp):
    kernel_(kernel), PPointsTree_(pp)
{ }

LowRankApp::LowRankApp(Kernel kernel, const Eigen::VectorXd& x, const Eigen::VectorXd& y):
    kernel_(kernel), Tx_(x), Ty_(y)
{ }

// approximate matrix-vector multiplication
Eigen::VectorXd LowRankApp::mvProd(const Eigen::VectorXd& c, double eta, unsigned deg)
{

    // compute V-matrices
    //Tx_.setV(deg);
    //Ty_.setV(deg);
    std::cout << "mvProd test1" << std::endl;
    PPointsTree_.setV(deg);
    std::cout << "mvProd test2" << std::endl;

    // compute V*c restricted to node indices of the tree
    //Ty_.setVc(c);
    PPointsTree_.setVc(c);
    std::cout << "mvProd test3" << std::endl;

    // add pointers to near and far field nodes of the tree
    //Tx_.setNearFar(eta, Ty_);
    PPointsTree_.setNearFar(eta, PPointsTree_);

    PPointsTree_.getRoot()->printree(0);
    // compute far field contribution
    Eigen::VectorXd f_approx = Eigen::VectorXd::Zero(c.size());
    //ff_contribution(f_approx, Tx_.getRoot(), deg);
    // compute near-field contribution
    //nf_contribution(f_approx, Tx_.getRoot(), c);

    return f_approx;
}


// compute far field contribution
void LowRankApp::ff_contribution(Eigen::VectorXd& f, Node* tx, unsigned deg)
{
    if((*tx).getLChild() != NULL) {

        Eigen::VectorXd XVc = Eigen::VectorXd::Zero(deg+1); // auxiliary variable
        int ixl = (*tx).getLInd(); // start index of cluster *tx
        int ixr = (*tx).getRInd(); // last  index of cluster *tx

        // iterate over far field of *tx_root
        std::vector<Node*> ffx = (*tx).getFarF(); // far field of *tx_root
        for(std::vector<Node*>::iterator iter=ffx.begin(); iter!=ffx.end(); ++iter) {

            int iyl = (**iter).getLInd(); // start index of current cluster in the far field
            int iyr = (**iter).getRInd(); // last  index of current cluster in the far field
            BlockCluster X_(Tx_.getVals()[ixl], Tx_.getVals()[ixr], Ty_.getVals()[iyl], Ty_.getVals()[iyr], deg, kernel_);
            Eigen::MatrixXd X = X_.getMatrix(); // matrix $X_{\sigma,\mu}$
            Eigen::VectorXd Vc = (**iter).getVc_node(); // V*c restricted to the indices of **iter
            XVc += X * Vc; // add contribution of block **iter to "s"
        };

        Eigen::MatrixXd Vx = (*tx).getV_node(); // $V_{\sigma}$
        f.segment(ixl, ixr-ixl+1) += Vx * XVc; // add contribution of far field to "f"

        Node* xl_c = (*tx).getLChild(); // pointer to left  child of *tx
        Node* xr_c = (*tx).getRChild(); // pointer to right child of *tx

        // add contribution of leaves of *tx
        ff_contribution(f, xl_c, deg);
        ff_contribution(f, xr_c, deg);
    }
}


// compute near-field contribution
void LowRankApp::nf_contribution(Eigen::VectorXd& f, Node* tx, const Eigen::VectorXd& c)
{
    if(tx != NULL) {

        unsigned ixl=(*tx).getLInd(); // start index of cluster *tx
        unsigned ixr=(*tx).getRInd(); // last  index of cluster *tx

        // iterate over near field of *tx
        std::vector<Node*> nfx = (*tx).getNearF(); // near field of *tx
        for(std::vector<Node*>::iterator iter=nfx.begin(); iter!=nfx.end(); ++iter) {
            unsigned iyl = (**iter).getLInd(); // start index of current cluster in the near field
            unsigned iyr = (**iter).getRInd(); // last  index of current cluster in the near field
            for(unsigned i=ixl; i<=ixr; ++i) {
                for(unsigned j=iyl; j<=iyr; ++j) {
                    f(i) += kernel_(Tx_.getVals()[i], Ty_.getVals()[j]) * c(j); // add near field contribution
                }
            }
        }

        Node* xl_c = (*tx).getLChild(); // pointer to left child of *tx
        if(xl_c != NULL) { // if *tx is not a leaf
            Node* xr_c = (*tx).getRChild(); // pointer to right child of *tx

            // add contribution of leaves of *tx
            nf_contribution(f, xl_c, c);
            nf_contribution(f, xr_c, c);
        }
    }
}
