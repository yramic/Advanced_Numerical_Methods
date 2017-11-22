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
#include "../include/node.hpp"
#include "../include/cheby.hpp"
#include "../include/segment.hpp"
#include <iostream>

#define inertia

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
    if(!segments_.empty() && segments_.size()>1) { // if there are points in the PPointsTree vector of points then they are equaly divided into the node´s children
#ifdef equal_clusters
        /* SAM_LISTING_BEGIN_0 */
        auto checkX = [](Point a, Point b) -> bool { return a.getX()<b.getX(); };
        std::sort(PPointsTree_.begin(),PPointsTree_.end(),checkX);
        std::vector<Point>::iterator it;
        it=PPointsTree_.begin()+(PPointsTree_.size()+1)/2; // set iterator in the middle of the vector of points of this node
                                                           // alternative: sort points by x-coordinates and split them in half
        std::vector<Point> l_points,r_points;              // now sort l\_points and r\_points based on their y-coordinates
        l_points.assign(PPointsTree_.begin(), it);
        r_points.assign(it, PPointsTree_.end());
        auto checkY = [](Point a, Point b) -> bool { return a.getY()<b.getY(); };
        std::sort(l_points.begin(),l_points.end(),checkY); // sort left and right vectors into top and bottom based on the y-coordinates
        std::sort(r_points.begin(),r_points.end(),checkY);
        std::vector<Point> tl_PPoints, tr_PPoints, bl_PPoints, br_PPoints; // creation of vectors of points of child nodes
        it=l_points.begin()+(l_points.size()+1)/2;
        tl_PPoints.assign(it,l_points.end());
        bl_PPoints.assign(l_points.begin(),it);
        it=r_points.begin()+(r_points.size()+1)/2;
        tr_PPoints.assign(it,r_points.end());
        br_PPoints.assign(r_points.begin(),it);

        if (!tl_PPoints.empty()) tl_child_ = new Node(tl_PPoints, deg_); // recursive construction of the Cluster Tree levels below root
        if (!tr_PPoints.empty()) tr_child_ = new Node(tr_PPoints, deg_);
        if (!bl_PPoints.empty()) bl_child_ = new Node(bl_PPoints, deg_);
        if (!br_PPoints.empty()) br_child_ = new Node(br_PPoints, deg_);
        /* SAM_LISTING_END_0 */
#endif
#ifdef inertia
        /* SAM_LISTING_BEGIN_1 */
        Eigen::MatrixXd A(segments_.size(),2);
        for(int i = 0; i < segments_.size(); i++){
            A(i,0) = segments_[i].getX();
            A(i,1) = segments_[i].getY();
        }
        Eigen::BDCSVD<Eigen::MatrixXd> svdOfA(A,Eigen::DecompositionOptions::ComputeEigenvectors | Eigen::DecompositionOptions::ComputeFullV);
        Eigen::MatrixXd V = svdOfA.matrixV();
        double sumX, sumY;
        for(int i = 0; i < segments_.size(); i++){
            sumX += segments_[i].getX();
            sumY += segments_[i].getY();
        }
        double avgX = sumX/(double)segments_.size(), avgY = sumY/(double)segments_.size();
        auto y = [](double x, double x1, double y1, double avgX, double avgY) -> double {return avgY+(x-avgX)*y1/x1; };
        std::vector<Segment> top_PPoints, bottom_PPoints, left_PPoints, right_PPoints; // creation of vectors of points of child nodes
        for(int i = 0; i < segments_.size(); i++){
            double y1 = y(segments_[i].getX(),V(0,0),V(1,0),avgX,avgY); // y value of the line that is defined by the point {avgX,avgY} and the vector corresponding to the biggest eigen value
            double y2 = y(segments_[i].getX(),V(0,1),V(1,1),avgX,avgY); // y value of the line that is defined by the point {avgX,avgY} and the vector corresponding to the second biggest eigen value
            if(y2<=segments_[i].getY()){
                if(y1<=segments_[i].getY()){
                    top_PPoints.push_back(segments_[i]);
                } else {
                    right_PPoints.push_back(segments_[i]);
                }
            } else {
                if(y1<=segments_[i].getY()){
                    left_PPoints.push_back(segments_[i]);
                } else {
                    bottom_PPoints.push_back(segments_[i]);
                }
            }
        }

        if (!top_PPoints.empty()) tl_child_ = new Node(top_PPoints, deg_); // recursive construction of the Cluster Tree levels below root
        if (!right_PPoints.empty()) tr_child_ = new Node(right_PPoints, deg_);
        if (!bottom_PPoints.empty()) bl_child_ = new Node(bottom_PPoints, deg_);
        if (!left_PPoints.empty()) br_child_ = new Node(left_PPoints, deg_);
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

// compute V-matrix of node
unsigned Node::setV()
{
    Eigen::MatrixXd VnodeX = Eigen::MatrixXd::Ones(segments_.size(), (deg+1));
    Eigen::MatrixXd VnodeY = Eigen::MatrixXd::Ones(segments_.size(), (deg+1));
    for(unsigned i=0; i<=ppts-1; ++i) {
        for(unsigned j=0; j<=deg; ++j) {
            for(unsigned k=0; k<j; ++k) {
                VnodeX(i,j) *= PPointsTree_[i].getX() - tkx_[k];
            }
            // Skip "k == j"
            for(unsigned k=j+1; k<=deg; ++k) {
                VnodeX(i,j) *= PPointsTree_[i].getX() - tkx_[k];
            }
            VnodeX(i,j) *= wkx_(j);
        }
    }
    for(unsigned i=0; i<=ppts-1; ++i) {
        for(unsigned j=0; j<=deg; ++j) {
            for(unsigned k=0; k<j; ++k) {
                VnodeY(i,j) *= PPointsTree_[i].getY() - tky_[k];
            }
            // Skip "k == j"
            for(unsigned k=j+1; k<=deg; ++k) {
                VnodeY(i,j) *= PPointsTree_[i].getY() - tky_[k];
            }
            VnodeY(i,j) *= wky_(j);
        }
    }

    Eigen::MatrixXd V_node_new(ppts, (deg+1)*(deg+1));
    for(unsigned i=0; i<=ppts-1; ++i) {
        for(unsigned j=0; j<=deg; ++j) {
            V_node_new.block(i, j*(deg+1), 1, deg+1) = VnodeX(i,j) * VnodeY.row(i);
        }
    }
    V_node_ = V_node_new;





    return ppts * (deg_+2)*(deg_+1)/2; // return no. of 'operations' performed
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
