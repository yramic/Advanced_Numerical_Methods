#include "../include/hierarchical_partition.hpp"
#include "../include/block_cluster.hpp"
#include "../include/cheby.hpp"
#include "../include/ctree.hpp"
#include "../include/is_admissible.hpp"
#include "../include/node.hpp"
#include "../include/point.hpp"
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
            FarField_.push_back(BlockCluster(xnode,ynode));
        } else { // else we consider the children of *xnode and *ynode and check whether their clusters are admissible
            setNearFar_recursion((*xnode).getLChild(), (*ynode).getLChild(), eta, Ty);
            setNearFar_recursion((*xnode).getRChild(), (*ynode).getLChild(), eta, Ty);
            setNearFar_recursion((*xnode).getLChild(), (*ynode).getRChild(), eta, Ty);
            setNearFar_recursion((*xnode).getRChild(), (*ynode).getRChild(), eta, Ty);
        }
    }
}

// return the bounding box corresponding to index i of the far-field vector
std::pair<std::pair<double,double>,std::pair<double,double> > HierarchicalPartitioning::getBB(int i)
{
    Node* xnode = FarField_[i].getXNode();
    std::vector<Point> xpts = xnode->getPoints();
    double xmin = xpts[0].getX();
    double xmax = xpts[0].getX();
    for(int j=1; j<xpts.size(); j++){
        if(xpts[j].getX() < xmin){
            xmin = xpts[j].getX();
        }
        if(xpts[j].getX() > xmax){
            xmax = xpts[j].getX();
        }
    }
    Node* ynode = FarField_[i].getYNode();
    std::vector<Point> ypts = ynode->getPoints();
    double ymin = ypts[0].getX();
    double ymax = ypts[0].getX();
    for(int j=1; j<ypts.size(); j++){
        if(ypts[j].getX() < ymin){
            ymin = ypts[j].getX();
        }
        if(ypts[j].getX() > ymax){
            ymax = ypts[j].getX();
        }
    }
    return {{xmin,ymin},{xmax,ymax}};
}
