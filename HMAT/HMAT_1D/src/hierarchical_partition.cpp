#include "../include/ctree.hpp"
#include "../include/is_admissible.hpp"
#include "../include/node.hpp"
#include "../include/cheby.hpp"
#include "../include/hierarchical_partition.hpp"
#include <Eigen/Dense>
#include <vector>
#include <iostream>

// add pairs of pointers to Node of the Cluster Tree in the Near and Far Field Vectors of pairs
void HierarchicalPartitioning::setNearFar_recursion(Node* xnode, Node* ynode, double eta, cTree Ty)
{
    // if *xnode or *ynode is a leaf, we add the pair (*xnode,*ynode) to the near field vector
    if((*xnode).getLChild() == NULL || (*ynode).getRChild() == NULL) {
        NearField_.push_back(std::make_pair(xnode,ynode));
    } else {
        double xl = (*xnode).getPoints().front().getX();
        double xr = (*xnode).getPoints().back().getX();
        double yl = (*ynode).getPoints().front().getX();
        double yr = (*ynode).getPoints().back().getX();
        // if the cluster corresponding to *xnode and *ynode is admissible, we add the pair (*xnode,*ynode) to the far field vector
        if(is_admissible(xl,xr,yl,yr, eta)) {
            // the line above checks the admissibility condition (eta)
            FarField_.push_back(std::make_pair(xnode,ynode));
        } else { // else we consider the children of *xnode and *ynode and check whether their clusters are admissible
            setNearFar_recursion((*xnode).getLChild(), (*ynode).getLChild(), eta, Ty);
            setNearFar_recursion((*xnode).getRChild(), (*ynode).getLChild(), eta, Ty);
            setNearFar_recursion((*xnode).getLChild(), (*ynode).getRChild(), eta, Ty);
            setNearFar_recursion((*xnode).getRChild(), (*ynode).getRChild(), eta, Ty);
        }
    }
}
