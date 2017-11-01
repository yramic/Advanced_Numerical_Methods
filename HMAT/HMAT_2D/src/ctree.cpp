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
cTree::cTree(const std::vector<Point> &PPoints, unsigned deg):
    root_(NULL), PPointsTree_(PPoints)
{
    unsigned n = PPoints.size();
    if(n > 0) { // build the tree
        root_ = new Node(PPoints,deg); // root is a node, leaves are added
    }
}

// make two lists with the x- and y-coordinates of boundaries of the bounding boxes of the clusters:
// odd entries of the lists are coordinates of the left boundaries // even entries are coordinates of the right boundaries
/*void cTree::setLists_recursion(Node* cluster, std::vector<double>& xlist, std::vector<double>& ylist)
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
}*/

