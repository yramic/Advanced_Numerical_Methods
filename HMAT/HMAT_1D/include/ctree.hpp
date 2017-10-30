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
#ifndef CTREE_HPP
#define CTREE_HPP
//#define ver1
#define ver2
#include <Eigen/Dense>
#include <vector>
#include "point.hpp"

// forward declaration to avoid cross-referencing
class Node;

#ifdef ver2
/**
* \brief Cluster tree class
*/
class cTree
{
public:

    /**
    * \brief Constructors
    */
    // default constructor
    cTree(): root_(NULL) { }
    // actual  constructor
    cTree(const std::vector<Point>& x, unsigned deg);
    // copy the root of a tree
    cTree(const cTree& T): root_(T.getRoot()) {}
    /**
    * \brief Getters
    */
    // return a pointer to node-root of "cTree"
    Node* getRoot() const { return root_; }
    /**
    * \brief Setters
    */
    // make lists with boundaries of bounding boxes (just for testing)
    void setLists(std::vector<double>& xlist, std::vector<double>& ylist) {
        setLists_recursion(root_, xlist, ylist);
    }

private:

    /**
    * \brief Recursions
    */
    // needed for "setLists(...)"
    void setLists_recursion(Node* cluster, std::vector<double>& xlist, std::vector<double>& ylist);

    Node* root_; // pointer to node-root of "cTree"
    friend class Node;
};
#endif
#ifdef ver1
    /**
    * \brief Cluster tree class
    */
    class cTree
    {
    public:

        /**
        * \brief Constructors
        */
        // default constructor
        cTree(): root_(NULL) { }
        // actual  constructor
        cTree(const Eigen::VectorXd& x);

        /**
        * \brief Getters
        */
        // return a pointer to node-root of "cTree"
        Node* getRoot() const { return root_; }
        // return the whole vector "x_"
        Eigen::VectorXd getVals() const { return grid_; }

        /**
        * \brief Setters
        */
        // compute V-matrices for nodes of the tree (contains evaluations of Chebyshew polynomials at corresponding points)
        void setV(unsigned deg) { setV_recursion(root_, deg); }
        // compute V*c restricted to node indices of the tree
        void setVc(const Eigen::VectorXd& c) { setVc_recursion(root_, c); }
        // add pointers to near and far field nodes of the tree
        void setNearFar(double eta, cTree Ty) {
            setNearFar_recursion(root_, Ty.root_, eta, Ty);
        }
        // make lists with boundaries of bounding boxes (just for testing)
        void setLists(std::vector<double>& xlist, std::vector<double>& ylist) {
            setLists_recursion(root_, xlist, ylist);
        }

    private:

        /**
        * \brief Recursions
        */
        // needed for "setV(...)"
        void setV_recursion(Node* cluster, unsigned deg);
        // needed for "setVc(...)"
        void setVc_recursion(Node* cluster, const Eigen::VectorXd& c);
        // needed for "setNearFar(...)"
        void setNearFar_recursion(Node* xnode, Node* ynode, double eta, cTree Ty);
        // needed for "setLists(...)"
        void setLists_recursion(Node* cluster, std::vector<double>& xlist, std::vector<double>& ylist);

        Node* root_; // pointer to node-root of "cTree"
        const Eigen::VectorXd grid_; // vector associated to "cTree"
        friend class Node;
    };
#endif
#endif // CTREE_HPP
