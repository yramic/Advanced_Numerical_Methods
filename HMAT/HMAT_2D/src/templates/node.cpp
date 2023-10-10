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
#include "../../include/node.hpp"
#include "../../include/cheby.hpp"
#include "../../include/point.hpp"
#include <iostream>

#define equal_clusters

// actual  constructor: creates the root of the Cluster Tree and the recursivly creates the leaves
Node::Node(std::vector<Point> Points, unsigned deg):
    tl_child_(NULL), tr_child_(NULL), bl_child_(NULL), br_child_(NULL), deg_(deg), PPointsTree_(Points), CVc_node_(Eigen::VectorXd::Zero((deg+1)*(deg+1)))
{
    setSons();
}

// destructor
Node::~Node()
{
    if(tl_child_ != NULL) delete tl_child_;
    if(tr_child_ != NULL) delete tr_child_;
    if(bl_child_ != NULL) delete bl_child_;
    if(br_child_ != NULL) delete br_child_;
}

// calculate the rectangle defined by the points of the node
void Node::getRect()
{
    // TODO
}

// build tree recursively
void Node::setSons()
{
    if(!PPointsTree_.empty() && PPointsTree_.size()>1) { // if there are points in the PPointsTree vector of points then they are equaly divided into the node´s children
#ifdef equal_clusters
        // TODO
#endif
#ifdef inertia
        // TODO
#endif
        getRect(); // calculate the rectangle defined by the points of the node
    }
}

// compute V-matrix of node
unsigned Node::setV()
{
    unsigned ppts = PPointsTree_.size();

    // TODO

    return ppts * (deg_+2)*(deg_+1)/2; // return no. of 'operations' performed
}

// compute V*c restricted to node indices
unsigned Node::setVc(const Eigen::VectorXd& c)
{
    int n = PPointsTree_.size();
    Eigen::VectorXd c_seg = Eigen::VectorXd::Zero(n);
    for(int i=0; i<n; i++){ // get only the part of vector c needed
        c_seg[i] = c(PPointsTree_[i].getId());
    }
    Vc_node_ = V_node_.transpose() * c_seg; // Vc matrix calculation
    return V_node_.rows()*V_node_.cols(); // return no. of 'operations' performed
}

// update C*V*c restricted to node indices
unsigned Node::setCVc(const Eigen::VectorXd& CVc)
{
    CVc_node_ += CVc;
    return CVc_node_.size(); // return no. of 'operations' performed
}

void Node::printree(int n)
{
    std::cout << "Node " << n << std::endl;
    for (std::vector<Point>::iterator it=PPointsTree_.begin(); it!=PPointsTree_.end(); it++) {
        std::cout << it->getId() << " ";
    }
    std::cout << std::endl;
    if (tl_child_!=NULL) {
        std::cout << "tl " << n << " ";
        tl_child_->printree(n+1);
    }
    if (tr_child_!=NULL) {
        std::cout << "tr " << n << " ";
        tr_child_->printree(n+1);
    }
    if (bl_child_!=NULL) {
        std::cout << "bl " << n << " ";
        bl_child_->printree(n+1);
    }
    if (br_child_!=NULL) {
        std::cout << "br " << n << " ";
        br_child_->printree(n+1);
    }
}
