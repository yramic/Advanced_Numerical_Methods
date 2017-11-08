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
#ifndef NODE_Y_HPP
#define NODE_Y_HPP

#include "../node.hpp"
#include <Eigen/Dense>

// forward declaration to avoid cross-referencing
template<typename NODE>
class cTree;

/**
* \brief Node_Y of a cluster tree ("cTree" class) for uni-directional interpolation
*/
class Node_Y : public Node
{
    using Node::Node; // C++11 inheritance of constructors

public:
    /*!
     * \brief Actual constructor
     * \param points Vector of points that is contained in the node
     * \param id ID of the node in the cTree
     * \param deg Degree of interpolation
     */
    Node_Y(std::vector<Point> points, int& id, unsigned deg);
    /*!
     * \brief Destructor
     */
    virtual ~Node_Y();
    /*!
     * \brief build tree recursively
     * \param id ID of node
     */
    void setLeaves(int& id);
    /*!
     * \brief compute fake V-matrix of ynode with uni-directional interpolation
     */
    void setV();
};

#endif // NODE_Y_HPP
