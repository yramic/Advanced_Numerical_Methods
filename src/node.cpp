#include "../include/node.hpp"
#include "../include/cheby.hpp"
#include <iostream>


// destructor of the class
node::~node()
{
  if ((l_child_ == NULL) && (r_child_ == NULL))
     { std::cout << "leaf destroyed" << std::endl; }
   if (l_child_ != NULL) delete l_child_;
   if (r_child_ != NULL) delete r_child_;
}


// deconstructor functions
void node::destroy_tree()
{
  destroy_tree(l_child_);
  destroy_tree(r_child_);
}


void node::destroy_tree(node *leaf)
{
  if (leaf!=NULL)	{
    destroy_tree(leaf->l_child_);
    destroy_tree(leaf->r_child_);
    delete leaf;	
  };	
}


// build tree recursively	
void node::add_leaf()
{
  if (r_ind_-l_ind_>0){
    l_child_=new node(l_ind_,(l_ind_+r_ind_)/2);
    r_child_=new node((l_ind_+r_ind_)/2+1,r_ind_);
    //(*lchild).add_leaf();
    //(*rchild).add_leaf();
  }	
}


// construct V-matrix
void node::fill_V(const Eigen::VectorXd &x, int deg)
{
  if(r_ind_-l_ind_>0){
    double xmin=x[l_ind_];
    double xmax=x[r_ind_];
    cheby cp(xmin, xmax,deg);
    const std::vector<double>& t_k=cp.cheb_pts(); // Chebyshew points
    const Eigen::VectorXd& omega=cp.get_omega(); // weights for Lagrange polynomials
    V_=Eigen::MatrixXd::Constant(r_ind_-l_ind_+1,deg+1,1);
		
    for(unsigned int i=0;i<=r_ind_-l_ind_;++i){
      for(unsigned int j=0;j<=deg;++j){
	for(unsigned int k=0;k<j;++k){
      V_(i,j)=V_(i,j)*(x[i+l_ind_]-t_k[k]);};
    for(unsigned int k=j+1;k<=deg;++k){
      V_(i,j)=V_(i,j)*(x[i+l_ind_]-t_k[k]);};
    V_(i,j)=V_(i,j)*omega(j);
      };
    };
  };
};


// multiplication of c with V restricted to indices of node
void node::Vc(const Eigen::VectorXd& c)
{
  if(r_ind_-l_ind_>0){
    Eigen::VectorXd c_seg=c.segment(l_ind_,r_ind_-l_ind_+1);
    c_node_=V_.transpose()*c_seg;
  };
}
