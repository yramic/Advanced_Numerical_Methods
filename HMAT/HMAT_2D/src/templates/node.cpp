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

#define inertia

// actual  constructor: creates the root of the Cluster Tree and the recursivly creates the leaves
Node::Node(std::vector<Point> Points, unsigned deg):
    tl_child_(NULL), tr_child_(NULL), bl_child_(NULL), br_child_(NULL), deg_(deg), PPointsTree_(Points), CVc_node_(Eigen::VectorXd::Zero((deg+1)*(deg+1)))
{ setSons(); }

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
    double maxX,minX,maxY,minY;
    maxX = PPointsTree_.begin()->getX();
    minX = PPointsTree_.begin()->getX();
    maxY = PPointsTree_.begin()->getY();
    minY = PPointsTree_.begin()->getY();
    for (std::vector<Point>::iterator it=PPointsTree_.begin()+1; it!=PPointsTree_.end(); it++){
        if (it->getX()>maxX) {
            maxX = it->getX();
        }
        else if (it->getX()<minX) {
            minX = it->getX();
        }
        if (it->getY()>maxY) {
            maxY = it->getY();
        }
        else if (it->getY()<minY) {
            minY = it->getY();
        }
        Cheby cbx(x1_b_, x2_b_, deg_);
        Cheby cby(y1_b_, y2_b_, deg_);
        tkx_ = cbx.getNodes(); // Chebyshew nodes for x axis
        wkx_ = cbx.getWghts(); // weights of Lagrange polynomial for x axis
        tky_ = cby.getNodes(); // Chebyshew nodes for y axis
        wky_ = cby.getWghts(); // weights of Lagrange polynomial for y axis
    }
    x1_ = minX;
    x2_ = maxX;
    y1_ = minY;
    y2_ = maxY;
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
        // fix for the rectangle if it is a segment
        if(std::abs(x1_-x2_)<10*std::numeric_limits<double>::epsilon()){
            x2_++;
        }
        if(std::abs(y1_-y2_)<10*std::numeric_limits<double>::epsilon()){
            y2_++;
        }
        Cheby cbx(x1_, x2_, deg_); // Chebvchev interpolation on the edges of the rectangle
        Cheby cby(y1_, y2_, deg_);
        tkx_ = cbx.getNodes(); // Chebyshew nodes for x axis
        wkx_ = cbx.getWghts(); // weights of Lagrange polynomial for x axis
        tky_ = cby.getNodes(); // Chebyshew nodes for y axis
        wky_ = cby.getWghts(); // weights of Lagrange polynomial for y axis
    }
}

// compute V-matrix of node
unsigned Node::setV()
{
    int ppts = PPointsTree_.size();

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
