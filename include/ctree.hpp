#ifndef CTREE_HPP
#define CTREE_HPP

#include <Eigen/Dense>
#include <vector>
#include "node.hpp"


class cTree
{
public:

    // default constructor
    cTree();
    // constructor building a tree
    cTree(std::vector<double> Xvec);

    // destructor
    virtual ~cTree() { }

    // returns the whole vector x_
    Eigen::VectorXd getVals() const {
        return x_;
    }
    // returns a pointer to root of cTree
    node* get_root() {
        return root;
    }

    // build V-matrices for nodes of the tree (contain evaluations of Chebyshew polynomials at corresponding points)
    void add_V(unsigned deg) {
        V_recursion(root, deg);
    }
    // product V*c restricted to the indices of each node of the tree, PRE: V already constructed
    void vc_mult(const Eigen::VectorXd& c) {
        c_recursion(root, c);
    }
    // add lists of pointers to near and far field nodes to each node of the tree
    void near_far(double eta, cTree Ty) {
        divide_tree(root, Ty.root, eta, Ty);
    }
    // just for testing, makes lists with the boundaries of the bounding boxes
    void make_fflist(std::vector<double>& xlist, std::vector<double>& ylist){
        rec_fflist(xlist, ylist, root);
    }

private:

     // needed for "add_V(...)"
    void V_recursion(node* cluster, int polynom_degree);

    // needed for "vc_mult(...)"
    void c_recursion(node* cluster, const Eigen::VectorXd& c_vec);

    // needed for "near_far(...)"
    void divide_tree(node* xnode, node* ynode, double eta, cTree Ty);

    // needed for "make_fflist(...)"
    void rec_fflist(std::vector<double>& xlist, std::vector<double>& ylist, node* xnode);

    friend class node;
    node* root; // pointer to root of the tree
    const Eigen::VectorXd x_; // vector associated to cTree
};

#endif // CTREE_HPP
