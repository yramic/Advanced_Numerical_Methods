/***********************************************************************
 *                                                                     *
 * Code for Course "Advanced Numerical Methods for CSE"                *
 * (Prof. Dr. R. Hiptmair)                                             *
 * Author: Ioannis Magkanaris                                          *
 * Date: 11/2017                                                       *
 * (C) Seminar for Applied Mathematics, ETH Zurich                     *
 * This code can be freely used for non-commercial purposes as long    *
 * as this header is left intact.                                      *
 ***********************************************************************/
#include "../include/hierarchical_partition.hpp"
#include "../include/block_nearf.hpp"
#include "../include/cheby.hpp"
#include "../include/ctree.hpp"
#include "../include/is_admissible.hpp"
#include "../include/node.hpp"
#include <Eigen/Dense>
#include <vector>
#include <iostream>

// Constructor: creates the cluster tree and
HierarchicalPartitioning::HierarchicalPartitioning(const std::vector<Point>& GPoints, double eta, unsigned deg):
    Tx_(GPoints,deg), Ty_(Tx_), eta_(eta)
{}
// compute the Far and Near Field pairs
void HierarchicalPartitioning::setNearFar()
{
    setNearFar_recursion(Tx_.getRoot(), Ty_.getRoot(), eta_, Ty_);
    auto checkpointers = [](Node* x, Node* y) -> bool { return x<y; };
    std::sort(FarField_xnds_.begin(), FarField_xnds_.end(), checkpointers);
    FarField_xnds_.erase(std::unique(FarField_xnds_.begin(), FarField_xnds_.end()), FarField_xnds_.end());
    std::sort(FarField_ynds_.begin(), FarField_ynds_.end(), checkpointers);
    FarField_ynds_.erase(std::unique(FarField_ynds_.begin(), FarField_ynds_.end()), FarField_ynds_.end());
}
// compute the Far and Near Field pairs
void HierarchicalPartitioning::setNearFar(Eigen::MatrixXd& cmatrix)
{
    setNearFar_recursion(Tx_.getRoot(), Ty_.getRoot(), eta_, Ty_, cmatrix);
    auto checkpointers = [](Node* x, Node* y) -> bool { return x<y; };
    std::sort(FarField_xnds_.begin(), FarField_xnds_.end(), checkpointers);
    FarField_xnds_.erase(std::unique(FarField_xnds_.begin(), FarField_xnds_.end()), FarField_xnds_.end());
    std::sort(FarField_ynds_.begin(), FarField_ynds_.end(), checkpointers);
    FarField_ynds_.erase(std::unique(FarField_ynds_.begin(), FarField_ynds_.end()), FarField_ynds_.end());
}
// add pointers to near and far field nodes of the tree
void HierarchicalPartitioning::setNearFar_recursion(Node* xnode, Node* ynode, double eta, cTree &Ty)
{
    // if *xnode or *ynode don`t exist, we have got nothing to compare
    if(xnode == NULL || ynode == NULL) {
        return;
    }
    // admissibility for 4D Matrices
    AdmissibilityH adm;
    if((*ynode).getPPoints().size()<=1 || (*xnode).getPPoints().size()<=1){
            //(*ynode).near_f_.push_back(xnode);
            NearField_.push_back(new BlockNearF(xnode,ynode));
    }
    else {
        // if the cluster corresponding to *xnode and *ynode is admissible, we add them to the far field list of each one
        if(adm.is_admissible(xnode, ynode, eta)) {
            //(*xnode).far_f_.push_back(ynode);
            FarField_.push_back(new BlockCluster(xnode,ynode));
            FarField_xnds_.push_back(xnode); FarField_ynds_.push_back(ynode);
        } else {    // else we consider all the different combinations of the children of *xnode and *ynode and check whether their clusters are admissible
            setNearFar_recursion((*xnode).getTl_Child(), (*ynode).getTl_Child(), eta, Ty);
            setNearFar_recursion((*xnode).getTl_Child(), (*ynode).getTr_Child(), eta, Ty);
            setNearFar_recursion((*xnode).getTl_Child(), (*ynode).getBl_Child(), eta, Ty);
            setNearFar_recursion((*xnode).getTl_Child(), (*ynode).getBr_Child(), eta, Ty);
            setNearFar_recursion((*xnode).getTr_Child(), (*ynode).getTl_Child(), eta, Ty);
            setNearFar_recursion((*xnode).getTr_Child(), (*ynode).getTr_Child(), eta, Ty);
            setNearFar_recursion((*xnode).getTr_Child(), (*ynode).getBl_Child(), eta, Ty);
            setNearFar_recursion((*xnode).getTr_Child(), (*ynode).getBr_Child(), eta, Ty);
            setNearFar_recursion((*xnode).getBl_Child(), (*ynode).getTl_Child(), eta, Ty);
            setNearFar_recursion((*xnode).getBl_Child(), (*ynode).getTr_Child(), eta, Ty);
            setNearFar_recursion((*xnode).getBl_Child(), (*ynode).getBl_Child(), eta, Ty);
            setNearFar_recursion((*xnode).getBl_Child(), (*ynode).getBr_Child(), eta, Ty);
            setNearFar_recursion((*xnode).getBr_Child(), (*ynode).getTl_Child(), eta, Ty);
            setNearFar_recursion((*xnode).getBr_Child(), (*ynode).getTr_Child(), eta, Ty);
            setNearFar_recursion((*xnode).getBr_Child(), (*ynode).getBl_Child(), eta, Ty);
            setNearFar_recursion((*xnode).getBr_Child(), (*ynode).getBr_Child(), eta, Ty);
        }
    }

}

// add pointers to near and far field nodes of the tree for debugging
void HierarchicalPartitioning::setNearFar_recursion(Node* xnode, Node* ynode, double eta, cTree &Ty, Eigen::MatrixXd& cmatrix)
{
    // if *xnode or *ynode don`t exist, we have got nothing to compare
    if(xnode == NULL || ynode == NULL) {
        return;
    }
    // admissibility for 4D Matrices
    AdmissibilityH adm;
    if((*ynode).getPPoints().size()<=1 || (*xnode).getPPoints().size()<=1){
        //(*ynode).near_f_.push_back(xnode);
        NearField_.push_back(new BlockNearF(xnode,ynode));
        cmatrix((*ynode).getNodeID(),(*xnode).getNodeID())++;
    }
    else {
        // if the cluster corresponding to *xnode and *ynode is admissible, we add them to the far field list of each one
        if(adm.is_admissible(xnode, ynode, eta)) {
            //(*xnode).far_f_.push_back(ynode);
            FarField_.push_back(new BlockCluster(xnode,ynode));
            cmatrix((*xnode).getNodeID(),(*ynode).getNodeID())++;
        } else {    // else we consider all the different combinations of the children of *xnode and *ynode and check whether their clusters are admissible
            setNearFar_recursion((*xnode).getTl_Child(), (*ynode).getTl_Child(), eta, Ty, cmatrix);
            setNearFar_recursion((*xnode).getTl_Child(), (*ynode).getTr_Child(), eta, Ty, cmatrix);
            setNearFar_recursion((*xnode).getTl_Child(), (*ynode).getBl_Child(), eta, Ty, cmatrix);
            setNearFar_recursion((*xnode).getTl_Child(), (*ynode).getBr_Child(), eta, Ty, cmatrix);
            setNearFar_recursion((*xnode).getTr_Child(), (*ynode).getTl_Child(), eta, Ty, cmatrix);
            setNearFar_recursion((*xnode).getTr_Child(), (*ynode).getTr_Child(), eta, Ty, cmatrix);
            setNearFar_recursion((*xnode).getTr_Child(), (*ynode).getBl_Child(), eta, Ty, cmatrix);
            setNearFar_recursion((*xnode).getTr_Child(), (*ynode).getBr_Child(), eta, Ty, cmatrix);
            setNearFar_recursion((*xnode).getBl_Child(), (*ynode).getTl_Child(), eta, Ty, cmatrix);
            setNearFar_recursion((*xnode).getBl_Child(), (*ynode).getTr_Child(), eta, Ty, cmatrix);
            setNearFar_recursion((*xnode).getBl_Child(), (*ynode).getBl_Child(), eta, Ty, cmatrix);
            setNearFar_recursion((*xnode).getBl_Child(), (*ynode).getBr_Child(), eta, Ty, cmatrix);
            setNearFar_recursion((*xnode).getBr_Child(), (*ynode).getTl_Child(), eta, Ty, cmatrix);
            setNearFar_recursion((*xnode).getBr_Child(), (*ynode).getTr_Child(), eta, Ty, cmatrix);
            setNearFar_recursion((*xnode).getBr_Child(), (*ynode).getBl_Child(), eta, Ty, cmatrix);
            setNearFar_recursion((*xnode).getBr_Child(), (*ynode).getBr_Child(), eta, Ty, cmatrix);
        }
    }
}

// return the bounding box corresponding to index i of the far-field vector
std::pair<std::pair<std::pair<double,double>,std::pair<double,double>>,std::pair<std::pair<double,double>,std::pair<double,double>> > HierarchicalPartitioning::getBB(int i)
{
    Node* xnode = FarField_[i]->getXNode();
    Node* ynode = FarField_[i]->getYNode();
    return {{{xnode->getX1(),xnode->getY1()},{xnode->getX2(),xnode->getY2()}},{{ynode->getX1(),ynode->getY1()},{ynode->getX2(),ynode->getY2()}}};
}
