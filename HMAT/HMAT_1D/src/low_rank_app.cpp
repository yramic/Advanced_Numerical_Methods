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

// compute far field contribution
void LowRankApp::ff_contribution(Eigen::VectorXd& f, std::vector<std::pair<Node*,Node*>> ff_v, const Eigen::VectorXd& c)
{
    int n = ff_v.size();
    for(int i = 0; i<n; i++){   // iterate for all the pairs of far field nodes
        Node* xnode = ff_v[i].first;
        Node* ynode = ff_v[i].second;
        xnode->setV();
        ynode->setV();
        ynode->setVc(c);
        Eigen::MatrixXd Vc = ynode->getVc_Node();
        std::vector<Point> xp = xnode->getPoints();
        std::vector<Point> yp = ynode->getPoints();
        BlockCluster X_(xnode->getTK(), ynode->getTK(), deg_, kernel_);
        Eigen::MatrixXd X = X_.getMatrix(); // matrix $X_{\sigma,\mu}$
        Eigen::VectorXd XVc = Eigen::VectorXd::Zero(deg_+1);
        XVc += X * Vc;
        Eigen::MatrixXd Vx = xnode->getV_Node();
        Eigen::VectorXd f_seg(xp.size());
        f_seg = Vx * XVc;
        for(int j=0; j<xnode->getPoints().size(); j++){
            f[xnode->getPoints()[j].getId()] += f_seg[j];   // add contribution of far field to ``f''
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
                f(xnode->getPoints()[j].getId()) += kernel_(xnode->getPoints()[j].getX(), ynode->getPoints()[k].getX()) * c(ynode->getPoints()[k].getId()); // add near field contribution to ``f''
            }
        }
    }
}
