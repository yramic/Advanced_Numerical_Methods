// node class

#include "../include/node.hpp"
#include "../include/cheby.hpp"
#include <iostream>

// deconstructor of the class
node::~node() {
  if ((lchild == NULL) && (rchild == NULL)) 
     { std::cout << "leaf destroyed" << std::endl; }
   if (lchild != NULL) delete lchild; 
   if (rchild != NULL) delete rchild;  
};
// deconstructor functions
void node::destroy_tree(){
  destroy_tree(lchild);
  destroy_tree(rchild);	
};
void node::destroy_tree(node *leaf){
  if (leaf!=NULL)	{
    destroy_tree(leaf->lchild);
    destroy_tree(leaf->rchild);
    delete leaf;	
  };	
};
// build tree recursively	
void node::add_leaf(){
  if (r_ind-l_ind>0){
    lchild=new node(l_ind,(l_ind+r_ind)/2);
    rchild=new node((l_ind+r_ind)/2+1,r_ind);
    //(*lchild).add_leaf();
    //(*rchild).add_leaf();
  }	
}; 

// construct V-matrix
void node::fill_V(const std::vector<double>& Xvals, int degree_p){
  if(r_ind-l_ind>0){
    double xmin=Xvals[l_ind]; 
    double xmax=Xvals[r_ind];
    cheby cp(xmin, xmax,degree_p);
    const std::vector<double>& t_k=cp.cheb_pts(); // Chebyshew points
    const Eigen::VectorXd& omega=cp.get_omega(); // weights for Lagrange polynomials
    V=Eigen::MatrixXd::Constant(r_ind-l_ind+1,degree_p+1,1);
		
    for(unsigned int i=0;i<=r_ind-l_ind;++i){
      for(unsigned int j=0;j<=degree_p;++j){
	for(unsigned int k=0;k<j;++k){
	  V(i,j)=V(i,j)*(Xvals[i+l_ind]-t_k[k]);};
	for(unsigned int k=j+1;k<=degree_p;++k){
	  V(i,j)=V(i,j)*(Xvals[i+l_ind]-t_k[k]);};
	V(i,j)=V(i,j)*omega(j);
      };
    };
  };
};
// multiplication of c with V restricted to indices of node
void node::V_c(const Eigen::VectorXd& cvector){
  if(r_ind-l_ind>0){	
    Eigen::VectorXd c_seg=cvector.segment(l_ind,r_ind-l_ind+1);
    c_node=V.transpose()*c_seg;
  };
};


