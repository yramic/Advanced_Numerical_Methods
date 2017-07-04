#include "../include/ctree.hpp"
#include "../include/is_admissible.hpp"
#include <Eigen/Dense>


// constructor with X only
cTree::cTree(const Eigen::VectorXd& x):
    root_(NULL), x_(x)
{
    unsigned n = x.size();
    if(n > 1) { // build the tree
        root_ = new Node(0, n-1); // root is a node, leafs are added
    } else {
        root_ = new Node(0, 0);
        // if vector x has less than 2 elements, the tree consists of a single node
    }
}


// build V-matrices for nodes of the tree (contain evaluations of Chebyshew polynomials at corresponding points)
void cTree::V_recursion(Node* cluster, unsigned deg)
{
    if((*cluster).l_child_ != NULL) {
        // build V-matrix for *cluster
        (*cluster).setV(x_, deg);
        // recursively call the function for the left  child of *cluster
        V_recursion((*cluster).l_child_, deg);
        // recursively call the function for the right child of *cluster
        V_recursion((*cluster).r_child_, deg);
    }
}


// product V*c restricted to the indices of each node of the tree. IMPO.: V should already be constructed!
void cTree::c_recursion(Node* cluster, const Eigen::VectorXd& c)
{
    if((*cluster).l_child_ != NULL) {
        // do the multiplication V*c restricted to the indices beloning to the cluster of *cluster
        (*cluster).setVc_node(c);
        // call the function for the left child of *cluster
        c_recursion((*cluster).l_child_, c);
        // same for the right child of *cluster
        c_recursion((*cluster).r_child_, c);
    }
}


// add lists of pointers to near and far field nodes to each node of the tree
void cTree::divide_tree(Node* xnode, Node* ynode, double eta, cTree Ty)
{
    // if *xnode or *ynode is a leaf, we add it to the near field
    if((*xnode).l_child_ == NULL || (*ynode).l_child_ == NULL) {

        (*xnode).push2NearF(ynode);

    } else { // if the cluster corresponding to *xnode and *ynode is admissible, we add ynode to the far field list of *xnode

        if(is_admissible(x_[(*xnode).l_ind_], x_[(*xnode).r_ind_], (Ty.getVals())[(*ynode).l_ind_], (Ty.getVals())[(*ynode).r_ind_], eta)) {

            (*xnode).push2FarF(ynode);

        } else { // we consider the children of *xnode and *ynode and check whether their clusters are admissible
            divide_tree((*xnode).l_child_, (*ynode).l_child_, eta, Ty);
            divide_tree((*xnode).r_child_, (*ynode).l_child_, eta, Ty);
            divide_tree((*xnode).l_child_, (*ynode).r_child_, eta, Ty);
            divide_tree((*xnode).r_child_, (*ynode).r_child_, eta, Ty);
        }
    }
}


// makes two lists with the x- and y-coordinates of boundaries of the bounding boxes of the clusters:
// odd entries of the lists are coordinates of the left boundaries and even entries are coordinates of the right boundaries of the bounding boxes
void cTree::rec_fflist(Node* cluster, std::vector<double>& xlist, std::vector<double>& ylist)
{
    // list of pointers to the nodes of the far field
    std::vector<Node*> ffx=(*cluster).getFarF();
    // left  index of *xnode
    unsigned ixl = (*cluster).getLInd();
    // right index of *xnode
    unsigned ixr = (*cluster).getRInd();

    for(std::vector<Node*>::iterator iter=ffx.begin(); iter!=ffx.end(); ++iter) {

        // start index of current cluster in the far field
        unsigned iyl = (**iter).left_ind();
        // last  index of current cluster in the far field
        unsigned iyr = (**iter).right_ind();

        xlist.push_back(x_[ixl]); // add the coordinate to the list
        xlist.push_back(x_[ixr]); // gaensefuessli
        ylist.push_back(x_[iyl]); // gaensefuessli
        ylist.push_back(x_[iyr]); // normal gaensefuessli :-)
    }
    if((*cluster).l_child_!=0) { // *xnode is not a leaf, so we do the same for its children
        rec_fflist((*cluster).l_child_, xlist, ylist);
        rec_fflist((*cluster).r_child_, xlist, ylist);
    }
}
