#include "../include/ctree.hpp"
#include "../include/is_admissible.hpp"
#include "../include/node.hpp"
#include "../include/cheby.hpp"
#include "../include/hierarchical_partition.hpp"
#include <Eigen/Dense>
#include <vector>
#include <iostream>

// Constructor: creates the cluster tree and
HierarchicalPartitioning::HierarchicalPartitioning(const std::vector<Point> &GPoints, double eta, unsigned deg):
    Tx_(GPoints,deg), Ty_(Tx_), eta_(eta), deg_(deg)
{}
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
            NearField_.push_back(std::make_pair(xnode,ynode));
    }
    else {
        // if the cluster corresponding to *xnode and *ynode is admissible, we add them to the far field list of each one
        if(adm.is_admissible(xnode, ynode, eta)) {
            //(*xnode).far_f_.push_back(ynode);
            FarField_.push_back(std::make_pair(xnode,ynode));
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
        NearField_.push_back(std::make_pair(xnode,ynode));
        cmatrix((*ynode).getNodeID(),(*xnode).getNodeID())++;
    }
    else {
        // if the cluster corresponding to *xnode and *ynode is admissible, we add them to the far field list of each one
        if(adm.is_admissible(xnode, ynode, eta)) {
            //(*xnode).far_f_.push_back(ynode);
            FarField_.push_back(std::make_pair(xnode,ynode));
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
