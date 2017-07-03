#include "../include/ctree.hpp"
#include "../include/is_admissible.hpp"


// constructor with X only
cTree::cTree(const Eigen::VectorXd& x):
    root_(NULL), x_(x)
{
    unsigned n = x.size();
    if(n > 1) { // build the tree
        root_ = new node(0, n-1); // root is a node, leafs are added
    } else {
        root_ = new node(NULL,NULL,0,0);
        // if vector x has less than 2 elements, the tree consists of a single node
    }
}


// build V-matrices for nodes of the tree (contain evaluations of Chebyshew polynomials at corresponding points)
void cTree::V_recursion(node* cluster, unsigned deg)
{
    if((*cluster).lchild != NULL) {
        // build V-matrix for *cluster
        (*cluster).fill_V(x_, deg);
        // recursively call the function for the left  child of *cluster
        V_recursion((*cluster).lchild, deg);
        // recursively call the function for the right child of *cluster
        V_recursion((*cluster).rchild, deg);
    }
}


// product V*c restricted to the indices of each node of the tree. IMPO.: V should already be constructed!
void cTree::c_recursion(node* cluster, const Eigen::VectorXd& c)
{
    if((*cluster).lchild != NULL) {
        // do the multiplication V*c restricted to the indices beloning to the cluster of *cluster
        (*cluster).V_c(c);
        // call the function for the left child of *cluster
        c_recursion((*cluster).lchild, c);
        // same for the right child of *cluster
        c_recursion((*cluster).rchild, c);
    }
}


// add lists of pointers to near and far field nodes to each node of the tree
void cTree::divide_tree(node* xnode, node* ynode, double eta, cTree Ty)
{
    // if *xnode or *ynode is a leaf, we add it to the near field
    if((*xnode).lchild == NULL || (*ynode).lchild == NULL) {

        (*xnode).push_nearf(ynode);

    } else { // if the cluster corresponding to *xnode and *ynode is admissible, we add ynode to the far field list of *xnode

        if(is_admissible(x_[(*xnode).l_ind], x_[(*xnode).r_ind], (Ty.getVals())[(*ynode).l_ind], (Ty.getVals())[(*ynode).r_ind], eta)) {

            (*xnode).push_farf(ynode);

        } else { // we consider the children of *xnode and *ynode and check whether their clusters are admissible
            divide_tree((*xnode).lchild, (*ynode).lchild, eta, Ty);
            divide_tree((*xnode).rchild, (*ynode).lchild, eta, Ty);
            divide_tree((*xnode).lchild, (*ynode).rchild, eta, Ty);
            divide_tree((*xnode).rchild, (*ynode).rchild, eta, Ty);
        }
    }
}


// makes two lists with the x- and y-coordinates of boundaries of the bounding boxes of the clusters:
// odd entries of the lists are coordinates of the left boundaries and even entries are coordinates of the right boundaries of the bounding boxes
void cTree::rec_fflist(node* cluster, std::vector<double>& xlist, std::vector<double>& ylist)
{
    // list of pointers to the nodes of the far field
    std::vector<node*> ffx=(*cluster).get_farf();
    // left  index of *xnode
    unsigned ixl = (*cluster).left_ind();
    // right index of *xnode
    unsigned ixr = (*cluster).right_ind();

    for(std::vector<node*>::iterator iter=ffx.begin(); iter!=ffx.end(); ++iter) {

        // start index of current cluster in the far field
        unsigned iyl = (**iter).left_ind();
        // last  index of current cluster in the far field
        unsigned iyr = (**iter).right_ind();

        xlist.push_back(x_[ixl]); // add the coordinate to the list
        xlist.push_back(x_[ixr]); // gaensefuessli
        ylist.push_back(x_[iyl]); // gaensefuessli
        ylist.push_back(x_[iyr]); // normal gaensefuessli :-)
    }
    if((*cluster).lchild!=0) { // *xnode is not a leaf, so we do the same for its children
        rec_fflist((*cluster).lchild, xlist, ylist);
        rec_fflist((*cluster).rchild, xlist, ylist);
    }
}
