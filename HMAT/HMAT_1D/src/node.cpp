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
#include "../include/node.hpp"
#include "../include/cheby.hpp"
#include <Eigen/Dense>
#include <iostream>

// actual constructor: adds a tree below the node if left_index != right_index
Node::Node(std::vector<Point> points, int& id, unsigned deg):
    l_child_(NULL), r_child_(NULL), node_points_(points), nodeId_(id), deg_(deg), CVc_node_(Eigen::VectorXd::Zero(deg+1))
{
    setLeaves(id);
}

// destructor
Node::~Node()
{
    if((l_child_ == NULL) && (r_child_ == NULL))
        std::cout << "leaves destroyed" << std::endl;
    if(l_child_ != NULL) delete l_child_;
    if(r_child_ != NULL) delete r_child_;
}

// build tree recursively
void Node::setLeaves(int& id)
{
    if(node_points_.size() > 1) {
        std::vector<Point>::iterator it;                    // we assume that the points are in ascending order
        it=node_points_.begin()+(node_points_.size()+1)/2;  // set iterator in the middle of the vector of points in this node
        std::vector<Point> lc,rc;                           // left and right child´s vectors
        lc.assign(node_points_.begin(), it);                // division of node´s points into it´s childs
        rc.assign(it, node_points_.end());
        id++;                                               // increase the identification number of the nodes of the tree
        l_child_ = new Node(lc,id,deg_);                    // recursivly do the same for the left child
        id++;                                               // increase the identification number of the nodes of the tree
        r_child_ = new Node(rc,id,deg_);                    // recursivly do the same for the right child
        double xmin = node_points_.front().getX();          // find min coordinate of the axis
        double xmax = node_points_.back().getX();           // find max coordinate of the axis
        Cheby cb(xmin, xmax, deg_);                         // Checyshev interpolation
        tk_ = cb.getNodes();                                // Chebyshev interpolation nodes
        wk_ = cb.getWghts();                                // Weights of Lagrange polynomial
    }
}

// compute V-matrix of node
unsigned Node::setV()
{
    int n = node_points_.size();
    V_node_ = Eigen::MatrixXd::Constant(n, deg_+1, 1);
    // V-matrix computation
    for(int i = 0; i<n; i++){
        for(int j = 0; j<=deg_; j++){
            for(int k=0; k<j; k++){
                V_node_(i,j) *= node_points_[i].getX() - tk_[k];
            }
            // Skip "k == j"
            for(int k=j+1; k<=deg_; ++k) {
                V_node_(i,j) *= node_points_[i].getX() - tk_[k];
            }
            V_node_(i,j) *= wk_[j];
        }
    }
    return n * (deg_+2)*(deg_+1)/2; // return no. of 'operations' performed
}

// compute V*c restricted to node indices
unsigned Node::setVc(const Eigen::VectorXd& c)
{
    int n = node_points_.size();
    Eigen::VectorXd c_seg = Eigen::VectorXd::Zero(n);
    for(int i=0; i<n; i++){ // get only the part of vector c needed
        c_seg[i] = c(node_points_[i].getId());
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
