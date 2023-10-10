/***********************************************************************
 *                                                                     *
 * Code for Course "Advanced Numerical Methods for CSE"                *
 * (Prof. Dr. R. Hiptmair)                                             *
 * Author: Daniele Casati                                              *
 * Date: 11/2017                                                       *
 * (C) Seminar for Applied Mathematics, ETH Zurich                     *
 * This code can be freely used for non-commercial purposes as long    *
 * as this header is left intact.                                      *
 ***********************************************************************/
#include "../../include/hierarchical_partition.hpp"
#include "../../include/block_nearf.hpp"
#include "../../include/cheby.hpp"
#include "../../include/ctree.hpp"
#include "../../include/is_admissible.hpp"
#include "../../include/node.hpp"
#include <Eigen/Dense>
#include <vector>

// Constructor: creates the cluster tree and
HierarchicalPartitioning::HierarchicalPartitioning(const std::vector<Segment>& segments, double eta, unsigned deg):
    Tx_(segments,deg), Ty_(Tx_), eta_(eta)
{}

// compute the Far and Near Field pairs
void HierarchicalPartitioning::setNearFar()
{
    setNearFar_recursion(Tx_.getRoot(), Ty_.getRoot(), eta_);
    auto checkpointers = [](Node* x, Node* y) -> bool { return x<y; };
    std::sort(FarField_xnds_.begin(), FarField_xnds_.end(), checkpointers);
    FarField_xnds_.erase(std::unique(FarField_xnds_.begin(), FarField_xnds_.end()), FarField_xnds_.end());
    std::sort(FarField_ynds_.begin(), FarField_ynds_.end(), checkpointers);
    FarField_ynds_.erase(std::unique(FarField_ynds_.begin(), FarField_ynds_.end()), FarField_ynds_.end());
}

// add pointers to near and far field nodes of the tree
void HierarchicalPartitioning::setNearFar_recursion(Node* xnode, Node* ynode, double eta)
{
    // if *xnode or *ynode don`t exist, we have got nothing to compare
    if(xnode == NULL || ynode == NULL) {
        return;
    }
    // admissibility for 4D Matrices
    AdmissibilityH adm;

    // TODO
}

// return the bounding box corresponding to index i of the far-field vector
std::pair<std::pair<std::pair<double,double>,std::pair<double,double>>,std::pair<std::pair<double,double>,std::pair<double,double>> > HierarchicalPartitioning::getBB(int i)
{
    Node* xnode = FarField_[i]->getXNode();
    Node* ynode = FarField_[i]->getYNode();
    return {{{xnode->getX1(),xnode->getY1()},{xnode->getX2(),xnode->getY2()}},{{ynode->getX1(),ynode->getY1()},{ynode->getX2(),ynode->getY2()}}};
}
