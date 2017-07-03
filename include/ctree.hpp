// header of the binary tree class
#ifndef CTREE_H
#define CTREE_H

#include <vector>
#include "Eigen/Dense"
#include "node.hpp"

class ctree{
 public:
  // constructors
  // default constructor
 ctree(): root(NULL), X(0){ }
  // constructor building a tree
  ctree(std::vector<double> Xvec);

  //deconstructor	
  virtual ~ctree(){//if (root != NULL) delete root;
  };

  // member access functions
  const std::vector<double>& get_X() const{return X;}	
  // returns the vector X
  node* get_root(){return root;}
  // returns a pointer to the root of the ctree	

  // public computation member functions
  void add_V(int polynom_degree){V_recursion(root,polynom_degree);}		
  // build the V-matrices for the nodes of the tree (contain evaluations of Chebyshew polynomials at the corresponding points)
  void vc_mult(const Eigen::VectorXd& c_vec){c_recursion(root,c_vec);}	
  // multiplication V*c restricted to the indices of each node of the tree, PRE: V already constructed
  void near_far(double eta, ctree Ty){divide_tree(root,Ty.root,eta, Ty);}
  // add lists of pointers to near and far field nodes to each node of the tree	
  void make_fflist(std::vector<double>& xlist, std::vector<double>& ylist){rec_fflist(xlist,ylist, root);} 
  // just for testing, makes lists with the boundaries of the bounding boxes

 private:
  // private member functions
  void V_recursion(node* cluster, int polynom_degree);
  // needed for member fcn "add_V(...)"
  void c_recursion(node* cluster,const Eigen::VectorXd& c_vec);
  // needed for member fcn "vc_mult(...)"
  void divide_tree(node* xnode, node* ynode, double eta, ctree Ty);
  // needed for member fcn "near_far(...)"
  void rec_fflist(std::vector<double>& xlist, std::vector<double>& ylist, node* xnode);
  // needed for member fcn "make_fflist(...)"
  // members:
  node* root; // pointer to the root of the tree
  const std::vector<double> X; // vector associated to ctree
  friend class node;
};
#endif 

