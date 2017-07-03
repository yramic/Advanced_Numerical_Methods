#include "../include/ctree.hpp"
#include "../include/comparison.hpp"


cTree::cTree():
    root(NULL), x_(0)
{ }

// constructor with X only
cTree::cTree(std::vector<double> Xvec) : root(NULL), x_(Xvec){
  int N=Xvec.size();
  if(N>1){ // building the tree
    root=new node(0,N-1); 
    // root is a node, leafs are added 
  }
  else{
    root=new node(NULL,NULL,0,0);
    // if the vector X has less than two elements the tree consists of a single node
  };
};

// construct V- matrices (contain evaluations of Chebyshew polynomials at the corresponding points)
void cTree::V_recursion(node* cluster, int polynom_degree){
  if((*cluster).lchild!=NULL){
    (*cluster).fill_V(x_,polynom_degree);
    // build V-matrix for *cluster
    V_recursion((*cluster).lchild,polynom_degree);
    // recursively call the function for the left child of *cluster
    V_recursion((*cluster).rchild,polynom_degree);
    // same for the right child
  };
};

// for the clusters in the tree, do multiplication V*c restricted to cluster, PRE: V already constructed
void cTree::c_recursion(node* cluster, const Eigen::VectorXd& c_vec){
  if((*cluster).lchild!=NULL){
    (*cluster).V_c(c_vec);	
    // do the multiplication V*c restricted to the indices beloning to the cluster of *cluster
    c_recursion((*cluster).lchild,c_vec);
    // call the function for the left child of *cluster	
    c_recursion((*cluster).rchild,c_vec);
    // same for the right child of *cluster
  };
};

// build near and far field, add lists of far field and near field nodes of Ty to the nodes of ctree
void cTree::divide_tree(node* xnode, node* ynode, double eta, cTree Ty){
  // if *xnode or *ynode is a leaf, we add it to the near field
  if((*xnode).lchild==NULL || (*ynode).lchild==NULL){
    (*xnode).push_nearf(ynode);
  }
  else{	// if the cluster corresponding to *xnode and *ynode is admissible, we add ynode to the far field list of *xnode.	
    if(is_admissible(x_[(*xnode).l_ind],x_[(*xnode).r_ind],(Ty.getVals())[(*ynode).l_ind],(Ty.getVals())[(*ynode).r_ind],eta)){
      (*xnode).push_farf(ynode);
    }
    else {// we consider the children of *xnode and *ynode and check whether their clusters are admissible
      divide_tree((*xnode).lchild,(*ynode).lchild, eta, Ty);
      divide_tree((*xnode).rchild,(*ynode).lchild, eta, Ty);
      divide_tree((*xnode).lchild,(*ynode).rchild, eta, Ty);
      divide_tree((*xnode).rchild,(*ynode).rchild, eta, Ty);
    };
  };
};

// makes two lists with the x- and y-coordinates of boundaries of the bounding boxes of the clusters, odd entries of the lists are coordinates of the left boundaries and even entries are coordinates of the right boundaries of the bounding boxes
void cTree::rec_fflist(std::vector<double>& xlist, std::vector<double>& ylist, node* xnode){
  std::vector<node*> ffx=(*xnode).get_farf(); 
  // ffx is a list of pointers to the nodes of the far field	
  int ixl=(*xnode).left_ind();
  // ixl is the left index of *xnode		
  int ixr=(*xnode).right_ind();
  // ixr is the right index of *xnode
  for(std::vector<node*>::iterator iter=ffx.begin(); iter!=ffx.end();++iter){	
    int iyl=(**iter).left_ind(); 
    // start index of current cluster in the far field
    int iyr=(**iter).right_ind(); 
    // last index of current cluster in the far field
    xlist.push_back(x_[ixl]);  // add the coordinate to the list
    xlist.push_back(x_[ixr]); // gaensefuessli
    ylist.push_back(x_[iyl]); // gaensefuessli
    ylist.push_back(x_[iyr]); // nomal gaensefuessli:-)
  };
  if ((*xnode).lchild!=0){ // *xnode is not a leaf, so we do the same for its children
    rec_fflist(xlist,ylist, (*xnode).lchild);
    rec_fflist(xlist,ylist, (*xnode).rchild);
  };
};
