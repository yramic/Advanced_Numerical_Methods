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
//#define ver1
#define ver2
#ifdef ver2
// actual constructor: adds a tree below the node if left_index != right_index
Node::Node(std::vector<Point> points, int& id, unsigned deg):
  l_child_(NULL), r_child_(NULL), node_points_(points), nodeId_(id), deg_(deg)
{ setLeaves(id); }


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
    it=node_points_.begin()+(node_points_.size()+1)/2;  // set  iterator in the middle of the vector of points in this node
    std::vector<Point> lc,rc;                           // left and right child´s vectors
    lc.assign(node_points_.begin(), it);                // division of node´s points into it´s childs
    rc.assign(it,node_points_.end());
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

#endif
#ifdef ver1
// actual constructor: adds a tree below the node if left_index != right_index
Node::Node(unsigned l_ind, unsigned r_ind):
  l_child_(NULL), r_child_(NULL), l_ind_(l_ind), r_ind_(r_ind), near_f_(0), far_f_(0)
{ setLeaves(); }


// destructor
Node::~Node()
{
  if((l_child_ == NULL) && (r_child_ == NULL))
    std::cout << "leaves destroyed" << std::endl;
  if(l_child_ != NULL) delete l_child_;
  if(r_child_ != NULL) delete r_child_;
}

// build tree recursively
void Node::setLeaves()
{
  if(r_ind_ - l_ind_ > 0) {
    l_child_ = new Node(l_ind_, (l_ind_+r_ind_)/2);
    r_child_ = new Node((l_ind_+r_ind_)/2 + 1, r_ind_);
  }
}

// compute V-matrix of node
void Node::setV_node(const Eigen::VectorXd &x, unsigned deg)
{
  if(r_ind_ - l_ind_ > 0) {
    double xmin = x[l_ind_]; // left bound of node interval
    double xmax = x[r_ind_]; // right bound of node interval
    Cheby cb(xmin, xmax, deg);
    Eigen::VectorXd tk = cb.getNodes(); // Chebyshew interpolation nodes
    Eigen::VectorXd wk = cb.getWghts(); // weights of Lagrange polynomial
    V_node_ = Eigen::MatrixXd::Constant(r_ind_-l_ind_+1, deg+1, 1);

    for(unsigned i=0; i<=r_ind_-l_ind_; ++i) { // loop: collocation  points of node
      for(unsigned j=0; j<=deg; ++j) { // loop: Lagrange polynomials
    for(unsigned k=0; k<j; ++k) { // evaluation loop
      V_node_(i,j) *= x[i+l_ind_] - tk[k];
    }
    // Skip "k == j"
    for(unsigned k=j+1; k<=deg; ++k) {
      V_node_(i,j) *= x[i+l_ind_] - tk[k];
    }
    V_node_(i,j) *= wk(j);
      } // end for j=
    } // end for i =
  }
}


// compute V*c restricted to node indices
void Node::setVc_node(const Eigen::VectorXd& c)
{
  if(r_ind_ - l_ind_ > 0) {
    Eigen::VectorXd c_seg = c.segment(l_ind_, r_ind_-l_ind_+1);
    Vc_node_ = V_node_.transpose() * c_seg;
  }
}
#endif
