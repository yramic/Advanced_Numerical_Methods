#include "../include/node.hpp"
#include "../include/cheby.hpp"
#include <Eigen/Dense>
#include <iostream>


// actual constructor: adds a tree below the node if left_index != right_index
Node::Node(unsigned left_index, unsigned right_index):
    l_child_(NULL), r_child_(NULL), l_ind_ (left_index), r_ind_(right_index), near_f_(0), far_f_(0)
{
    setLeaves();
}


// destructor
Node::~Node()
{
    if((l_child_ == NULL) && (r_child_ == NULL))
        std::cout << "leaves destroyed" << std::endl;
    if(l_child_ != NULL)
        delete l_child_;
    if(r_child_ != NULL)
        delete r_child_;
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

        double xmin = x[l_ind_];
        double xmax = x[r_ind_];
        Cheby cb(xmin, xmax, deg);
        Eigen::VectorXd tk = cb.getNodes(); // Chebyshew nodes
        Eigen::VectorXd wk = cb.getWghts(); // weights of Lagrange polynomial
        V_node_ = Eigen::MatrixXd::Constant(r_ind_-l_ind_+1, deg+1, 1);

        for(unsigned i=0; i<=r_ind_-l_ind_; ++i) {
            for(unsigned j=0; j<=deg; ++j) {
                for(unsigned k=0; k<j; ++k) {
                    V_node_(i,j) = V_node_(i,j) * (x[i+l_ind_]-tk[k]);
                }
                for(unsigned k=j+1; k<=deg; ++k) {
                    V_node_(i,j) = V_node_(i,j) * (x[i+l_ind_]-tk[k]);
                }
                V_node_(i,j) *= wk(j);
            }
        }
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
