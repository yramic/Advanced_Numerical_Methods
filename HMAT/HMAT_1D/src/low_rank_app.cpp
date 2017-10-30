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
#include "../include/ctree.hpp"
#include "../include/kernel.hpp"
#include "../include/node.hpp"
#include "../include/point.hpp"
#include <iostream>
//#define ver1
#define ver2
#ifdef ver2

// constructor
LowRankApp::LowRankApp(Kernel kernel, const std::vector<Point> &GPoints, double eta, unsigned deg):
    kernel_(kernel), GPoints_(GPoints), HP_(GPoints,eta,deg), deg_(deg)
{ }

// approximate matrix-vector multiplication
Eigen::VectorXd LowRankApp::mvProd(const Eigen::VectorXd& c)
{
    // compute Far and Near Field relationships between Nodes of the Cluster Tree
    HP_.setNearFar();
    // compute far field contribution
    Eigen::VectorXd f_approx = Eigen::VectorXd::Zero(c.size());
    ff_contribution(f_approx, HP_.getFF(), c);
    // compute near-field contribution
    nf_contribution(f_approx, HP_.getNF(), c);

    return f_approx;
}

// compute V-matrix of node
Eigen::MatrixXd LowRankApp::setV(Node* x)
{
    std::vector<Point> node_points = x->getPoints();
    int n = node_points.size();
    Eigen::VectorXd tk = x->getTK(); // Chebyshew interpolation nodes
    Eigen::VectorXd wk = x->getWK(); // weights of Lagrange polynomial

    Eigen::MatrixXd V_node = Eigen::MatrixXd::Constant(n, deg_+1, 1);
    // V-matrix computation
    for(int i = 0; i<n; i++){
        for(int j = 0; j<=deg_; j++){
            for(int k=0; k<j; k++){
                V_node(i,j) *= node_points[i].getX() - tk[k];
            }
            // Skip "k == j"
            for(int k=j+1; k<=deg_; ++k) {
                V_node(i,j) *= node_points[i].getX() - tk[k];
            }
            V_node(i,j) *= wk(j);
        }
    }
    return V_node;
}
// compute V*c restricted to node indices
Eigen::MatrixXd LowRankApp::setVc(Node* x, const Eigen::VectorXd& c)
{
    std::vector<Point> node_points = x->getPoints();
    int n = node_points.size();
    Eigen::VectorXd c_seg = Eigen::VectorXd::Zero(n);
    for(int i=0; i<n; i++){     // get only the part of vector c we need
        c_seg[i] = c(node_points[i].getId());
    }
    Eigen::MatrixXd Vx = setV(x);   // calculate the V matrix for the x node
    Eigen::MatrixXd Vc = Vx.transpose() * c_seg;    // Vc matrix calculation
    return Vc;
}
// compute far field contribution
void LowRankApp::ff_contribution(Eigen::VectorXd& f, std::vector<std::pair<Node*,Node*>> ff_v, const Eigen::VectorXd& c)
{
    int n = ff_v.size();
    for(int i = 0; i<n; i++){   // iterate for all the pairs of far field nodes
        Node* xnode = ff_v[i].first;
        Node* ynode = ff_v[i].second;
        Eigen::MatrixXd Vc = setVc(ynode, c);
        std::vector<Point> xp = xnode->getPoints();
        std::vector<Point> yp = ynode->getPoints();
        BlockCluster X_(xnode->getTK(), ynode->getTK(), deg_, kernel_);
        Eigen::MatrixXd X = X_.getMatrix(); // matrix $X_{\sigma,\mu}$
        Eigen::VectorXd XVc = Eigen::VectorXd::Zero(deg_+1);
        XVc += X * Vc;
        Eigen::MatrixXd Vx = setV(xnode);
        Eigen::VectorXd f_seg(xp.size());
        f_seg = Vx * XVc;
        for(int j=0; j<xnode->getPoints().size(); j++){
            f[xnode->getPoints()[j].getId()] += f_seg[j];   // add contribution of far field to "f"
        }
    }
}
// compute near-field contribution
void LowRankApp::nf_contribution(Eigen::VectorXd& f, std::vector<std::pair<Node*,Node*>> nf_v, const Eigen::VectorXd& c)
{
    int n = nf_v.size();
    for(int i = 0; i<n; i++){
        Node* xnode = nf_v[i].first;
        Node* ynode = nf_v[i].second;
        for(int j=0; j<xnode->getPoints().size(); j++){
            for(int k=0; k<ynode->getPoints().size(); k++){
                f(xnode->getPoints()[j].getId()) += kernel_(xnode->getPoints()[j].getX(), ynode->getPoints()[k].getX()) * c(ynode->getPoints()[k].getId()); // add near field contribution
            }
        }
    }
}

#endif

#ifdef ver1
// constructor
LowRankApp::LowRankApp(Kernel kernel, const Eigen::VectorXd& x, const Eigen::VectorXd& y):
    kernel_(kernel), Tx_(x), Ty_(y)
{ }


// approximate matrix-vector multiplication
Eigen::VectorXd LowRankApp::mvProd(const Eigen::VectorXd& c, double eta, unsigned deg)
{
    // compute V-matrices
    Tx_.setV(deg);
    Ty_.setV(deg);

    // compute V*c restricted to node indices of the tree
    Ty_.setVc(c);

    // add pointers to near and far field nodes of the tree
    Tx_.setNearFar(eta, Ty_);

    // compute far field contribution
    Eigen::VectorXd f_approx = Eigen::VectorXd::Zero(c.size());
    //ff_contribution(f_approx, Tx_.getRoot(), deg);
    // compute near-field contribution
    nf_contribution(f_approx, Tx_.getRoot(), c);

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
#endif
