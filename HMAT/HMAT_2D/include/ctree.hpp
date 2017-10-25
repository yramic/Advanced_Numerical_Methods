#ifndef CTREE_HPP
#define CTREE_HPP

#include <Eigen/Dense>
#include <vector>
#include "point.hpp"


// forward declaration to avoid cross-referencing
class Node;


/*!
* \brief Cluster tree class
*/
class cTree
{
public:

    /*!
    * \brief Default Constructor
    */
    // default constructor
    cTree():
        root_(NULL)
    { }
    // actual  constructor
    /*!
    * \brief Actual Constructor for 2D
    */
    cTree(const Eigen::VectorXd& x);
    /*!
    * \brief Actual Constructor for 4D
    */
    cTree(const std::vector<Point> PPointsTree);

    /*!
    * \brief return a pointer to node-root of "cTree"
    */
    // return a pointer to node-root of "cTree"
    Node* getRoot() const {
        return root_;
    }
    /*!
    * \brief return the whole vector "x_"
    */
    // return the whole vector "x_"
    Eigen::VectorXd getVals() const {
        return grid_;
    }

    /*!
    * \brief compute V-matrices for nodes of the tree (contains evaluations of Chebyshew polynomials at corresponding points)
    */
    // compute V-matrices for nodes of the tree (contains evaluations of Chebyshew polynomials at corresponding points)
    void setV(unsigned deg) {
        setV_recursion(root_, deg);
    }
    /*!
    * \brief compute V*c restricted to node indices of the tree
    */
    // compute V*c restricted to node indices of the tree
    void setVc(const Eigen::VectorXd& c) {
        setVc_recursion(root_, c);
    }
    /*!
    * \brief add pointers to near and far field nodes of the tree
    */
    // add pointers to near and far field nodes of the tree
    void setNearFar(double eta, cTree &Ty,int& near, int& far) {
        setNearFar_recursion(root_, Ty.root_, eta, Ty, near, far);
    }
    /*!
    * \brief add pointers to near and far field nodes of the tree for debugging
    */
    // add pointers to near and far field nodes of the tree debugging
    void setNearFar(double eta, cTree &Ty, Eigen::MatrixXd& cmatrix) {
        setNearFar_recursion(root_, Ty.root_, eta, Ty, cmatrix);
    }
    /*!
    * \brief make lists with boundaries of bounding boxes (just for testing)
    */
    // make lists with boundaries of bounding boxes (just for testing)
    void setLists(std::vector<double>& xlist, std::vector<double>& ylist) {
        setLists_recursion(root_, xlist, ylist);
    }

private:

    /*!
    * \brief needed for "setV(...)"
    */
    // needed for "setV(...)"
    void setV_recursion(Node* cluster, unsigned deg);
    /*!
    * \brief needed for "setVc(...)"
    */
    // needed for "setVc(...)"
    void setVc_recursion(Node* cluster, const Eigen::VectorXd& c);
    /*!
    * \brief needed for "setNearFar(...)"
    */
    // needed for "setNearFar(...)"
    void setNearFar_recursion(Node* xnode, Node* ynode, double eta, cTree &Ty, int& near, int& far);
    /*!
    * \brief needed for "setNearFar(...)" for debugging
    */
    // needed for "setNearFar(...)" for debugging
    void setNearFar_recursion(Node* xnode, Node* ynode, double eta, cTree &Ty, Eigen::MatrixXd& cmatrix);
    /*!
    * \brief needed for "setLists(...)"
    */
    // needed for "setLists(...)"
    void setLists_recursion(Node* cluster, std::vector<double>& xlist, std::vector<double>& ylist);

    Node* root_; //!< pointer to node-root of "cTree"
    const Eigen::VectorXd grid_; //!< vector associated to "cTree"
    const std::vector<Point> PPointsTree_;  //!< vector that has all the Polygon Points
    friend class Node;
};

#endif // CTREE_HPP
