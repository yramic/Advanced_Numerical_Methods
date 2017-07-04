#ifndef CTREE_HPP
#define CTREE_HPP

#include "node.hpp"
#include <Eigen/Dense>


class cTree
{
public:

    // default constructor
    cTree():
        root_(NULL), x_(0)
    { }
    // constructor building a tree
    cTree(const Eigen::VectorXd& x);

    // destructor
    virtual ~cTree() { }

    // returns the whole vector x_
    Eigen::VectorXd getVals() const {
        return x_;
    }
    // returns a pointer to root of cTree
    Node* get_root() {
        return root_;
    }

    // build V-matrices for nodes of the tree (contain evaluations of Chebyshew polynomials at corresponding points)
    void add_V(unsigned deg) {
        V_recursion(root_, deg);
    }
    // product V*c restricted to the indices of each node of the tree. IMPO.: V should already be constructed!
    void vc_mult(const Eigen::VectorXd& c) {
        c_recursion(root_, c);
    }
    // add lists of pointers to near and far field nodes to each node of the tree
    void near_far(double eta, cTree Ty) {
        divide_tree(root_, Ty.root_, eta, Ty);
    }
    // just for testing, makes lists with the boundaries of the bounding boxes
    void make_fflist(std::vector<double>& xlist, std::vector<double>& ylist) {
        rec_fflist(root_, xlist, ylist);
    }

private:

    // needed for "add_V(...)"
    void V_recursion(Node* cluster, unsigned deg);

    // needed for "vc_mult(...)"
    void c_recursion(Node* cluster, const Eigen::VectorXd& c);

    // needed for "near_far(...)"
    void divide_tree(Node* xnode, Node* ynode, double eta, cTree Ty);

    // needed for "make_fflist(...)"
    void rec_fflist(Node* cluster, std::vector<double>& xlist, std::vector<double>& ylist);

    Node* root_; // pointer to root of cTree
    const Eigen::VectorXd x_; // vector associated to cTree
    friend class Node;
};

#endif // CTREE_HPP
