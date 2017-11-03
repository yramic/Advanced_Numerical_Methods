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
    cTree(): root_(NULL) {}
    /*!
    * \brief Actual Constructor for 2D
    * \param PPointsTree Vector of Polygon Points
    */
    cTree(const std::vector<Point> &PPointsTree, unsigned deg);
    /*!
    * \brief Copy Constructor for 2D
    * \param T The cTree to copy
    */
    cTree(const cTree& T): root_(T.getRoot()) {}
    /*!
    * \brief Return a pointer to node-root of "cTree"
    */
    Node* getRoot() const {
        return root_;
    }

private:
    Node* root_; //!< pointer to node-root of "cTree"
    const std::vector<Point> PPointsTree_;  //!< vector that has all the Polygon Points
    friend class Node;
};

#endif // CTREE_HPP
