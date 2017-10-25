/*!
  * \author Ioannis Magkanaris, Daniele Casati
  * \date 11/2017
  * \mainpage Low Rank Approximation for BEM
  */
#include "../include/ctree.hpp"
#include "../include/point.hpp"
#include "../include/is_admissible.hpp"
#include "../include/node.hpp"
#include <Eigen/Dense>
#include <vector>
#include <iostream>


// actual constructor
cTree::cTree(const std::vector<Point> PPointsTree):
    root_(NULL), PPointsTree_(PPointsTree)
{
    unsigned n = PPointsTree.size();
    if(n > 0) { // build the tree
        root_ = new Node(PPointsTree); // root is a node, leaves are added
    }
}

cTree::cTree(const Eigen::VectorXd& x):
    root_(NULL), grid_(x)
{
    unsigned n = x.size();
    if(n > 1) { // build the tree
        root_ = new Node(0, n-1); // root is a node, leaves are added
    } else {
        root_ = new Node(0, 0);
        // if "x" has less than 2 elements, the tree consists of a single node
    }
}


// compute V-matrices for nodes of the tree (contains evaluations of Chebyshew polynomials at corresponding points)
void cTree::setV_recursion(Node* cluster, unsigned deg)
{
    //if((*cluster).l_child_ != NULL) {
    if(!(*cluster).PPointsTree_.empty()) {
        // compute V-matrix for *cluster
        (*cluster).setV_node(PPointsTree_, deg);    // passing PPointsTree_ here may is not needed

        // recursively call the function for the childs of *cluster
        if((*cluster).tl_child_ != NULL) setV_recursion((*cluster).tl_child_, deg);
        if((*cluster).tr_child_ != NULL) setV_recursion((*cluster).tr_child_, deg);
        if((*cluster).bl_child_ != NULL) setV_recursion((*cluster).bl_child_, deg);
        if((*cluster).br_child_ != NULL) setV_recursion((*cluster).br_child_, deg);
    }
}


// compute V*c restricted to node indices of the tree
void cTree::setVc_recursion(Node* cluster, const Eigen::VectorXd& c)
{
    if(!(*cluster).PPointsTree_.empty()) {
        // compute V*c restricted to node indices beloning to *cluster
        (*cluster).setVc_node(c);

        // recursively call the function for the childs of *cluster
        if((*cluster).tl_child_ != NULL) setVc_recursion((*cluster).tl_child_, c);
        if((*cluster).tr_child_ != NULL) setVc_recursion((*cluster).tr_child_, c);
        if((*cluster).bl_child_ != NULL) setVc_recursion((*cluster).bl_child_, c);
        if((*cluster).br_child_ != NULL) setVc_recursion((*cluster).br_child_, c);
    }
}


// add pointers to near and far field nodes of the tree
void cTree::setNearFar_recursion(Node* xnode, Node* ynode, double eta, cTree &Ty,int& near, int& far)
{
    // if *xnode or *ynode don`t exist, we have got nothing to compare
    if(xnode == NULL || ynode == NULL) {
        return;
    }
    // admissibility for 4D Matrices
    AdmissibilityH adm;
    if((*ynode).PPointsTree_.size()<=1 || (*xnode).PPointsTree_.size()<=1){
            (*ynode).near_f_.push_back(xnode);
            near++;
    }
    else {
        // if the cluster corresponding to *xnode and *ynode is admissible, we add them to the far field list of each one
        if(adm.is_admissible(xnode, ynode, eta)) {
            (*xnode).far_f_.push_back(ynode);
            far++;
        } else {    // else we consider all the different combinations of the children of *xnode and *ynode and check whether their clusters are admissible
            setNearFar_recursion((*xnode).tl_child_, (*ynode).tl_child_, eta, Ty, near, far);
            setNearFar_recursion((*xnode).tl_child_, (*ynode).tr_child_, eta, Ty, near, far);
            setNearFar_recursion((*xnode).tl_child_, (*ynode).bl_child_, eta, Ty, near, far);
            setNearFar_recursion((*xnode).tl_child_, (*ynode).br_child_, eta, Ty, near, far);
            setNearFar_recursion((*xnode).tr_child_, (*ynode).tl_child_, eta, Ty, near, far);
            setNearFar_recursion((*xnode).tr_child_, (*ynode).tr_child_, eta, Ty, near, far);
            setNearFar_recursion((*xnode).tr_child_, (*ynode).bl_child_, eta, Ty, near, far);
            setNearFar_recursion((*xnode).tr_child_, (*ynode).br_child_, eta, Ty, near, far);
            setNearFar_recursion((*xnode).bl_child_, (*ynode).tl_child_, eta, Ty, near, far);
            setNearFar_recursion((*xnode).bl_child_, (*ynode).tr_child_, eta, Ty, near, far);
            setNearFar_recursion((*xnode).bl_child_, (*ynode).bl_child_, eta, Ty, near, far);
            setNearFar_recursion((*xnode).bl_child_, (*ynode).br_child_, eta, Ty, near, far);
            setNearFar_recursion((*xnode).br_child_, (*ynode).tl_child_, eta, Ty, near, far);
            setNearFar_recursion((*xnode).br_child_, (*ynode).tr_child_, eta, Ty, near, far);
            setNearFar_recursion((*xnode).br_child_, (*ynode).bl_child_, eta, Ty, near, far);
            setNearFar_recursion((*xnode).br_child_, (*ynode).br_child_, eta, Ty, near, far);
        }
    }

}

// add pointers to near and far field nodes of the tree for debugging
void cTree::setNearFar_recursion(Node* xnode, Node* ynode, double eta, cTree &Ty, Eigen::MatrixXd& cmatrix)
{
    // if *xnode or *ynode don`t exist, we have got nothing to compare
    if(xnode == NULL || ynode == NULL) {
        return;
    }
    // admissibility for 4D Matrices
    AdmissibilityH adm;
    if((*ynode).PPointsTree_.size()<=1 || (*xnode).PPointsTree_.size()<=1){
        (*ynode).near_f_.push_back(xnode);
        cmatrix((*ynode).getNodeID(),(*xnode).getNodeID())++;
    }
    else {
        // if the cluster corresponding to *xnode and *ynode is admissible, we add them to the far field list of each one
        if(adm.is_admissible(xnode, ynode, eta)) {
            (*xnode).far_f_.push_back(ynode);
            cmatrix((*xnode).getNodeID(),(*ynode).getNodeID())++;
        } else {    // else we consider all the different combinations of the children of *xnode and *ynode and check whether their clusters are admissible
            setNearFar_recursion((*xnode).tl_child_, (*ynode).tl_child_, eta, Ty, cmatrix);
            setNearFar_recursion((*xnode).tl_child_, (*ynode).tr_child_, eta, Ty, cmatrix);
            setNearFar_recursion((*xnode).tl_child_, (*ynode).bl_child_, eta, Ty, cmatrix);
            setNearFar_recursion((*xnode).tl_child_, (*ynode).br_child_, eta, Ty, cmatrix);
            setNearFar_recursion((*xnode).tr_child_, (*ynode).tl_child_, eta, Ty, cmatrix);
            setNearFar_recursion((*xnode).tr_child_, (*ynode).tr_child_, eta, Ty, cmatrix);
            setNearFar_recursion((*xnode).tr_child_, (*ynode).bl_child_, eta, Ty, cmatrix);
            setNearFar_recursion((*xnode).tr_child_, (*ynode).br_child_, eta, Ty, cmatrix);
            setNearFar_recursion((*xnode).bl_child_, (*ynode).tl_child_, eta, Ty, cmatrix);
            setNearFar_recursion((*xnode).bl_child_, (*ynode).tr_child_, eta, Ty, cmatrix);
            setNearFar_recursion((*xnode).bl_child_, (*ynode).bl_child_, eta, Ty, cmatrix);
            setNearFar_recursion((*xnode).bl_child_, (*ynode).br_child_, eta, Ty, cmatrix);
            setNearFar_recursion((*xnode).br_child_, (*ynode).tl_child_, eta, Ty, cmatrix);
            setNearFar_recursion((*xnode).br_child_, (*ynode).tr_child_, eta, Ty, cmatrix);
            setNearFar_recursion((*xnode).br_child_, (*ynode).bl_child_, eta, Ty, cmatrix);
            setNearFar_recursion((*xnode).br_child_, (*ynode).br_child_, eta, Ty, cmatrix);
        }
    }

}
// make two lists with the x- and y-coordinates of boundaries of the bounding boxes of the clusters:
// odd entries of the lists are coordinates of the left boundaries // even entries are coordinates of the right boundaries
void cTree::setLists_recursion(Node* cluster, std::vector<double>& xlist, std::vector<double>& ylist)
{
    // list of pointers to nodes of the far field
    std::vector<Node*> ffx = (*cluster).getFarF();
    // left  index of *xnode
    unsigned ixl = (*cluster).getLInd();
    // right index of *xnode
    unsigned ixr = (*cluster).getRInd();

    for(std::vector<Node*>::iterator iter=ffx.begin(); iter!=ffx.end(); ++iter) {

        // start index of current cluster in the far field
        unsigned iyl = (**iter).l_ind_;
        // last  index of current cluster in the far field
        unsigned iyr = (**iter).r_ind_;

        // add the coordinates to the lists
        xlist.push_back(grid_[ixl]);
        xlist.push_back(grid_[ixr]);
        ylist.push_back(grid_[iyl]);
        ylist.push_back(grid_[iyr]);
    }
    if((*cluster).l_child_ != 0) { // *cluster is not a leaf, so we do the same for its children
        setLists_recursion((*cluster).l_child_, xlist, ylist);
        setLists_recursion((*cluster).r_child_, xlist, ylist);
    }
}

