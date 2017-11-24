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
#include <Eigen/SVD>
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
    /* SAM_LISTING_BEGIN_2 */
    double maxX,minX,maxY,minY;
    maxX = segments_.begin()->getA().x();
    minX = segments_.begin()->getA().x();
    maxY = segments_.begin()->getA().y();
    minY = segments_.begin()->getA().y();
    for (std::vector<Segment>::iterator it=segments_.begin()+1; it!=segments_.end(); it++){
        if (it->getA().x() > maxX) {
            maxX = it->getA().x();
        }
        else if (it->getA().x() < minX) {
            minX = it->getA().x();
        }
        if (it->getA().y() > maxY) {
            maxY = it->getA().y();
        }
        else if (it->getA().y() < minY) {
            minY = it->getA().y();
        }
        if (it->getB().x() > maxX) {
            maxX = it->getB().x();
        }
        else if (it->getB().x() < minX) {
            minX = it->getB().x();
        }
        if (it->getB().y() > maxY) {
            maxY = it->getB().y();
        }
        else if (it->getB().y() < minY) {
            minY = it->getB().y();
        }
    }
    x1_ = minX;
    x2_ = maxX;
    y1_ = minY;
    y2_ = maxY;
    Cheby cbx(x1_, x2_, deg_);
    Cheby cby(y1_, y2_, deg_);
    tkx_ = cbx.getNodes(); // Chebyshew nodes for x axis
    wkx_ = cbx.getWghts(); // weights of Lagrange polynomial for x axis
    tky_ = cby.getNodes(); // Chebyshew nodes for y axis
    wky_ = cby.getWghts(); // weights of Lagrange polynomial for y axis
    /* SAM_LISTING_END_2 */
}

// build tree recursively
void Node::setSons()
{
    if(!segments_.empty() && segments_.size()>1) { // if there are points in the PPointsTree vector of points then they are equaly divided into the node´s children
#ifdef equal_clusters
        /* SAM_LISTING_BEGIN_0 */
        auto checkX = [](Segment a, Segment b) -> bool { return a.getMean().x() < b.getMean().x(); };
        std::sort(segments_.begin(),segments_.end(),checkX);
        std::vector<Segment>::iterator it;
        it=segments_.begin()+(segments_.size()+1)/2; // set iterator in the middle of the vector of points of this node
        std::vector<Segment> l_segments,r_segments;  // now sort l\_points and r\_points based on their y-coordinates
        l_segments.assign(segments_.begin(), it);
        r_segments.assign(it, segments_.end());
        auto checkY = [](Segment a, Segment b) -> bool { return a.getMean().y() < b.getMean().y(); };
        std::sort(l_segments.begin(),l_segments.end(),checkY); // sort left and right vectors into top and bottom based on the y-coordinates
        std::sort(r_segments.begin(),r_segments.end(),checkY);
        std::vector<Segment> tl_segments, tr_segments, bl_segments, br_segments; // creation of vectors of points of child nodes
        it=l_segments.begin()+(l_segments.size()+1)/2;
        tl_segments.assign(it,l_segments.end());
        bl_segments.assign(l_segments.begin(),it);
        it=r_segments.begin()+(r_segments.size()+1)/2;
        tr_segments.assign(it,r_segments.end());
        br_segments.assign(r_segments.begin(),it);

        if (!tl_segments.empty()) tl_child_ = new Node(tl_segments, deg_); // recursive construction of the Cluster Tree levels below root
        if (!tr_segments.empty()) tr_child_ = new Node(tr_segments, deg_);
        if (!bl_segments.empty()) bl_child_ = new Node(bl_segments, deg_);
        if (!br_segments.empty()) br_child_ = new Node(br_segments, deg_);
        /* SAM_LISTING_END_0 */
#endif
#ifdef inertia
        /* SAM_LISTING_BEGIN_1 */
        double sumX, sumY;
        for(int i = 0; i < segments_.size(); i++){
            sumX += segments_[i].getMean().x();
            sumY += segments_[i].getMean().y();
        }
        double avgX = sumX/(double)segments_.size(), avgY = sumY/(double)segments_.size();
        Eigen::MatrixXd A(2, 2*segments_.size());
        for(unsigned i=0; i<segments_.size(); ++i){
            A(0,2*i) = segments_[i].getA().x() - avgX;
            A(1,2*i) = segments_[i].getA().y() - avgY;
            A(0,2*i+1) = segments_[i].getB().x() - avgX;
            A(1,2*i+1) = segments_[i].getB().y() - avgY;
        }
        Eigen::MatrixXd M = A * A.transpose();
        Eigen::JacobiSVD<Eigen::MatrixXd> svdOfM(M, Eigen::ComputeThinV); // 'M' is square, so ComputeFullV and ComputeThinV are the same
        Eigen::MatrixXd V = svdOfM.matrixV();
        auto y = [](double x, double x1, double y1, double avgX, double avgY) -> double { return avgY+(x-avgX)*y1/x1; };
        std::vector<Segment> top_segments, bottom_segments, left_segments, right_segments; // creation of vectors of points of child nodes
        for(unsigned i=0; i<segments_.size(); ++i){
            double y1 = y(segments_[i].getMean().x(),V(0,0),V(1,0),avgX,avgY); // y value of the line that is defined by the point {avgX,avgY} and the vector corresponding to the biggest eigen value
            double y2 = y(segments_[i].getMean().x(),V(0,1),V(1,1),avgX,avgY); // y value of the line that is defined by the point {avgX,avgY} and the vector corresponding to the second biggest eigen value
            if(y2<=segments_[i].getMean().y()){
                if(y1<=segments_[i].getMean().y()){
                    top_segments.push_back(segments_[i]);
                } else {
                    right_segments.push_back(segments_[i]);
                }
            } else {
                if(y1<=segments_[i].getMean().y()){
                    left_segments.push_back(segments_[i]);
                } else {
                    bottom_segments.push_back(segments_[i]);
                }
            }
        }

        if (!top_segments.empty()) tl_child_ = new Node(top_segments, deg_); // recursive construction of the Cluster Tree levels below root
        if (!right_segments.empty()) tr_child_ = new Node(right_segments, deg_);
        if (!bottom_segments.empty()) bl_child_ = new Node(bottom_segments, deg_);
        if (!left_segments.empty()) br_child_ = new Node(left_segments, deg_);
        /* SAM_LISTING_END_1 */
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
    V_node_ = Eigen::MatrixXd::Zero(segs, (deg_+1)*(deg_+1));
    for(unsigned i=0; i<segs; ++i) {

        for(unsigned j1=0; j1<=deg_; ++j1) {
            for(unsigned j2=0; j2<=deg_; ++j2) {
                double sum = 0.;
                for(unsigned k=0; k<order; ++k) {

                    /* transformation of quadrature nodes from [-1,1] to [a,b] */
                    Eigen::Vector2d tk = 0.5 * (segments_[i].getB() - segments_[i].getA()) * gauss_point[k] + 0.5 * (segments_[i].getB() + segments_[i].getA());
//                           double wk = 0.5 * (segments_[i].getB() - segments_[i].getA()).norm() * gauss_wht[k];

                    sum += gauss_wht[k] * evalLagrange(j1, tk.x()) * evalLagrange(j2, tk.y());
                }

                sum *= 0.5 * (segments_[i].getB() - segments_[i].getA()).norm();
                V_node_(i, j1*(deg_+1) + j2) = sum;
            }
        }
    }
    /* SAM_LISTING_END_3 */
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
