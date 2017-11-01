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
#include "../include/ctree.hpp"
#include "../include/is_admissible.hpp"
#include "../include/node.hpp"
#include <Eigen/Dense>
#include <vector>
//#define ver1
#define ver2

#ifdef ver2
cTree::cTree(const std::vector<Point>& GPoints, unsigned deg):
    root_(NULL)
{
    unsigned n = GPoints.size();
    if(n > 1) { // build the tree
        int node_id = 0;
        root_ = new Node(GPoints, node_id, deg); // root is a node, leaves are added
    }
}
/*
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
*/
#endif
#ifdef ver1
// actual constructor
cTree::cTree(const Eigen::VectorXd& x):  root_(NULL), grid_(x)
{
    unsigned n = x.size();
    if(n > 1) { // build the tree
        root_ = new Node(0, n-1); // root is a node, leaves are added
    } else {
        root_ = new Node(0, 0);
        // if "x" has less than 2 elements, the tree consists of a single node
    }
}


// compute V-matrices for nodes of the tree
// (involves evaluations of Chebyshew polynomials at corresponding points)
void cTree::setV_recursion(Node* cluster, unsigned deg)
{
    if((*cluster).l_child_ != NULL) {
        // compute V-matrix for *cluster
        (*cluster).setV_node(grid_, deg);
        // recursively call the function for the left  child of *cluster
        setV_recursion((*cluster).l_child_, deg);
        // recursively call the function for the right child of *cluster
        setV_recursion((*cluster).r_child_, deg);
    }
}


// compute V*c restricted to node indices of the tree
void cTree::setVc_recursion(Node* cluster, const Eigen::VectorXd& c)
{
    if((*cluster).l_child_ != NULL) {
        // compute V*c restricted to node indices beloning to *cluster
        (*cluster).setVc_node(c);
        // recursively call the function for the left  child of *cluster
        setVc_recursion((*cluster).l_child_, c);
        // recursively call the function for the right child of *cluster
        setVc_recursion((*cluster).r_child_, c);
    }
}


// add pointers to near and far field nodes of the tree
void cTree::setNearFar_recursion(Node* xnode, Node* ynode, double eta, cTree Ty)
{
    // if *xnode or *ynode is a leaf, we add it to the near field
    if((*xnode).l_child_ == NULL || (*ynode).l_child_ == NULL) {

        (*xnode).near_f_.push_back(ynode);

    } else {

        // if the cluster corresponding to *xnode and *ynode is admissible, we add *ynode to the far field list of *xnode
        if(is_admissible(grid_[(*xnode).l_ind_], grid_[(*xnode).r_ind_], (Ty.getVals())[(*ynode).l_ind_], (Ty.getVals())[(*ynode).r_ind_], eta)) {
            // the line above checks the admissibility condition (eta)

            (*xnode).far_f_.push_back(ynode);

        } else { // else we consider the children of *xnode and *ynode and check whether their clusters are admissible
            setNearFar_recursion((*xnode).l_child_, (*ynode).l_child_, eta, Ty);
            setNearFar_recursion((*xnode).r_child_, (*ynode).l_child_, eta, Ty);
            setNearFar_recursion((*xnode).l_child_, (*ynode).r_child_, eta, Ty);
            setNearFar_recursion((*xnode).r_child_, (*ynode).r_child_, eta, Ty);
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
#endif
