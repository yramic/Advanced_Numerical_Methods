/***********************************************************************
 *                                                                     *
 * Code for Course "Advanced Numerical Methods for CSE"                *
 * (Prof. Dr. R. Hiptmair)                                             *
 * Author: Daniele Casati                                              *
 * Date: 11/2017                                                       *
 * (C) Seminar for Applied Mathematics, ETH Zurich                     *
 * This code can be freely used for non-commercial purposes as long    *
 * as this header is left intact.                                      *
 ***********************************************************************/
#include "../../include/node.hpp"
#include "../../include/cheby.hpp"
#include "../../include/segment.hpp"
extern "C" {
#include "../../../BEM/CppHilbert/Library/source/gaussQuadrature.h"
}
#include <iostream>

#define equal_clusters

// actual  constructor: creates the root of the Cluster Tree and the recursivly creates the leaves
Node::Node(std::vector<Segment> segments, unsigned deg):
    tl_child_(NULL), tr_child_(NULL), bl_child_(NULL), br_child_(NULL), deg_(deg), segments_(segments), CVc_node_(Eigen::VectorXd::Zero((deg+1)*(deg+1)))
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
    if(!segments_.empty() && segments_.size()>1) { // if there are points in the PPointsTree vector of points then they are equaly divided into the node´s children
#ifdef equal_clusters
        // TODO
#endif
#ifdef inertia
        // TODO
#endif
        getRect(); // calculate the rectangle defined by the points of the node
    }
}

// evaluate Lagrange polynomial
double Node::evalLagrange(unsigned j, double tk)
{
    double result_n = 1., result_d = 1.;
    for(unsigned k=0; k<j; ++k) {
        result_n *= tk - tkx_[k];
        result_d *= tkx_[j] - tkx_[k];
    }
    // Skip "k == j"
    for(unsigned k=j+1; k<=deg_; ++k) {
        result_n *= tk - tkx_[k];
        result_d *= tkx_[j] - tkx_[k];
    }
//  return result_n * wkx_(j);
    return result_n / result_d;
}

// compute V-matrix of node
unsigned Node::setV()
{
    /* SAM_LISTING_BEGIN_3 */
    unsigned order;
    if(deg_ < 2)
        order = 2;
    else if(deg_ < 4)
        order = 4;
    else if(deg_ < 8)
        order = 8;
    else if(deg_ < 16)
        order = 16;
    else
        order = 32;

    const double* gauss_point = getGaussPoints(order);
    const double* gauss_wht   = getGaussWeights(order);

    unsigned segs = segments_.size();

    // TODO

    return segs * order * (deg_+1)*(deg_+1); // return no. of 'operations' performed
}

// compute V*c restricted to node indices
unsigned Node::setVc(const Eigen::VectorXd& c)
{
    int n = segments_.size();
    Eigen::VectorXd c_seg = Eigen::VectorXd::Zero(n);
    for(int i=0; i<n; i++){ // get only the part of vector c needed
        c_seg[i] = c(segments_[i].getId());
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
    for (std::vector<Segment>::iterator it=segments_.begin(); it!=segments_.end(); it++) {
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
