// header of the node class
#ifndef NODE_H
#define NODE_H

#include <vector>
#include "Eigen/Dense"
#include <math.h>

class node{
  typedef std::vector<node*> vector_type;
 public:
  
  // constructor(s)
  //default constructor
 node() : lchild(NULL), rchild(NULL), l_ind (0), r_ind (0), V(), c_node(), near_f(0), far_f(0) {} 
  // single node constructor, adds a tree below the node if left_index!=right_index
  node (int left_index, int right_index) : lchild(NULL), rchild(NULL),l_ind (left_index), r_ind (right_index), V(), c_node(), near_f(0), far_f(0) {add_leaf();}
  
  // tree node constructor
 node(node* left_child, node* right_child, int left_index, int right_index) : lchild(left_child), rchild(right_child), l_ind(left_index), r_ind(right_index), V(),c_node(), near_f(0), far_f(0) {}

  // tree node constructor with V matrix
  node(node* left_child, node* right_child, int left_index, int right_index, Eigen::MatrixXd Vmatrix) : lchild(left_child), rchild(right_child), l_ind(left_index), r_ind(right_index), V(Vmatrix),c_node(), near_f(0), far_f(0){}
  
  // deconstructor
  virtual ~node();
  
  // member access functions:
  int left_ind() const {return l_ind;}			
  // returns index of leftmost point in cluster of node
  int right_ind() const {return r_ind;}			
  // returns index of rightmost point in cluster of node
  node* leftchild() const {return lchild;}		
  // returns a pointer to the left child of the node	
  node* rightchild() const {return rchild;}		     
  // returns a pointer to the right child of the node
  const Eigen::MatrixXd& Vv() const {return V;}	
  // returns the matrix $V_{\sigma}$, where $\sigma$ denotes the cluster
  const Eigen::VectorXd& get_Vc() const {return c_node;}
  // returns V*c restricted to the indices of the cluster
  vector_type get_nearf() const {return near_f;}	
  // returns a list of pointers to the nodes belonging to near field of the node
  vector_type get_farf() const {return far_f;}		
  // returns a list of pointers to the nodes belonging to far field of the node
  
  // set indices, children
  void set_leftindex(int left){l_ind=left;}		
  // set lind to "left"
  void set_rightindex(int right){r_ind=right;}	
  // set rind to "right"
  void set_leftchild(node* leftkid){lchild=leftkid;}	
  // set lchild to "leftkid"
  void set_rightchild(node* rightkid){rchild=rightkid;}	
  // set rchild to "rightkid"
  
  // member functions
  void add_leaf();			
  // build tree (recursively), used in constructors of class ctree
  void fill_V(const std::vector<double>& Xvals, int degree_p); // to build V-matrices
  void V_c(const Eigen::VectorXd& cvector); 
  // V*c restricted to node indices
  void push_nearf(node* near_node){near_f.push_back(near_node);} 	
  // add node pointer to near field list
  void push_farf(node* far_node){far_f.push_back(far_node);}
  // add pointer to node to far field list
  
 private:
  void destroy_tree(node *leaf);
  void destroy_tree(); 
  // members
  node* lchild; // left child of node
  node* rchild;	// right child of node
  int l_ind;    // smallest index in cluster of node
  int r_ind;	// largest index in cluster of node  
  Eigen::MatrixXd V; // $V_{\sigma}$ where $\sigma$ is the cluster of node
  Eigen::VectorXd c_node; // $V_{\sigma}* c_{|\sigma}$ where $\sigma$ is the cluster of node
  vector_type near_f; // list with pointers to nodes of the near field of the node
  vector_type far_f; // list with pointers to nodes of the far field of the node
  friend class ctree; // our friend and partner:-)

 };

#endif 
