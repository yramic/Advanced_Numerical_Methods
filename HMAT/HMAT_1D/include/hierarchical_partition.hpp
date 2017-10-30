#ifndef HIERARCHICAL_PARTITION_HPP
#define HIERARCHICAL_PARTITION_HPP
//#define ver1
#define ver2
#include "ctree.hpp"
#include "kernel.hpp"
#include "point.hpp"
#ifdef ver2
class HierarchicalPartitioning
{
public:
    // Constructor: creates the cluster tree and
    HierarchicalPartitioning(const std::vector<Point> &GPoints, double eta, unsigned deg):
        Tx_(GPoints,deg), Ty_(Tx_), eta_(eta), deg_(deg)
    {}
    // return Far Field vector
    std::vector<std::pair<Node*,Node*>> getFF(){
        return FarField_;
    }
    // return Near Field vector
    std::vector<std::pair<Node*,Node*>> getNF(){
        return NearField_;
    }
    // compute the Near and Far Field vectors
    void setNearFar() { setNearFar_recursion(Tx_.getRoot(), Ty_.getRoot(), eta_, Ty_); }
private:
    // needed for "setNearFar(...)"
    void setNearFar_recursion(Node* xnode, Node* ynode, double eta, cTree Ty);
    cTree Tx_, Ty_; // Cluster trees for comparison
    std::vector<std::pair<Node*,Node*>> FarField_, NearField_;  // Vectors for Near and Far Field vectors
    unsigned deg_;  // degree of interpolation
    double eta_;    // eta-admissibility constant
};
#endif
#endif // HIERARCHICAL_PARTITION_HPP
