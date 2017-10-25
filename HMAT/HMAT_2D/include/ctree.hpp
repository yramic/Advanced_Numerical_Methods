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
    cTree():
        root_(NULL)
    { }

    /*!
    * \brief Actual Constructor for 2D
    * \param PPointsTree Vector of Polygon Points
    */
    cTree(const std::vector<Point> PPointsTree);

    /*!
    * \brief Return a pointer to node-root of "cTree"
    */
    Node* getRoot() const {
        return root_;
    }

    /*!
    * \brief Compute V-matrices for nodes of the tree (contains evaluations of Chebyshew polynomials at corresponding points)
    * \param deg Degree of interpolation
    */
    void setV(unsigned deg) {
        setV_recursion(root_, deg);
    }

    /*!
    * \brief Add pointers to near and far field nodes of the tree
    * \param eta eta variable of admissibility
    * \param Ty Reference tree for recursion
    * \param near Number of Near Field relationships between Nodes
    * \param far Number of Far Field relationships between Nodes
    */
    void setNearFar(double eta, cTree &Ty, int& near, int& far) {
        setNearFar_recursion(root_, Ty.root_, eta, Ty, near, far);
    }

    /*!
    * \brief Add pointers to near and far field nodes of the tree for debugging
    * \param eta eta variable of admissibility
    * \param Ty Reference tree for recursion
    * \param cmatrix Matrix for saving calculations used for checking
    */
    void setNearFar(double eta, cTree &Ty, Eigen::MatrixXd& cmatrix) {
        setNearFar_recursion(root_, Ty.root_, eta, Ty, cmatrix);
    }
    // I should delete this function maybe also
    /*!
    * \brief Make lists with boundaries of bounding boxes (just for testing)
    */
    void setLists(std::vector<double>& xlist, std::vector<double>& ylist) {
        setLists_recursion(root_, xlist, ylist);
    }

private:

    /*!
    * \brief Needed for "setV(...)"
    * \param cluster Node that contains the cluster
    * \param deg Degree of interpolation
    */
    void setV_recursion(Node* cluster, unsigned deg);

    /*!
    * \brief Needed for "setNearFar(...)"
    * \param xnode First node for checking
    * \param ynode Second node for checking
    * \param eta eta admissibility variable
    * \param Ty Tree for reference in the recursion
    * \param near Number of Near Field relationships between Nodes
    * \param far Number of Far Field relationships between Nodes
    */
    void setNearFar_recursion(Node* xnode, Node* ynode, double eta, cTree &Ty, int& near, int& far);

    /*!
    * \brief Needed for "setNearFar(...)" for debugging
    * \param xnode First node for checking
    * \param ynode Second node for checking
    * \param eta eta admissibility variable
    * \param Ty Tree for reference in the recursion
    * \param cmatrix Matrix for saving calculations used for checking
    */
    void setNearFar_recursion(Node* xnode, Node* ynode, double eta, cTree &Ty, Eigen::MatrixXd& cmatrix);

    /*!
    * \brief needed for "setLists(...)"
    */
    void setLists_recursion(Node* cluster, std::vector<double>& xlist, std::vector<double>& ylist);

    Node* root_; //!< pointer to node-root of "cTree"
    const std::vector<Point> PPointsTree_;  //!< vector that has all the Polygon Points
    friend class Node;
};

#endif // CTREE_HPP
