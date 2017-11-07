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
#ifndef BLOCK_NEARF_HPP
#define BLOCK_NEARF_HPP

#include <Eigen/Dense>
#include "kernel.hpp"
#include "node.hpp"

/*!
 * \brief Helper class to compute near field block matrices
 */
class BlockNearF
{
public:

    /*!
     * \brief Constructor
     * \param ndx XNode
     * \param ndy YNode
     */
    BlockNearF(Node* xnode, Node* ynode);

    /*!
     * \brief Constructor
     * \param ndx XNode
     * \param ndy YNode
     * \param G Kernel Function
     */
    BlockNearF(Node* xnode, Node* ynode, Kernel* G);

    /*!
     * \brief return near field block matrix
     */
    Eigen::MatrixXd getMatrix() const {
        return C_;
    }

    /*!
     * \brief compute near field block matrix
     */
    void setMatrix(Kernel* G);

    /*!
     * \brief return pointer to xnode
     */
    Node* getXNode() {
        return pair_.first;
    }

    /*!
     * \brief return pointer to ynode
     */
    Node* getYNode() {
        return pair_.second;
    }

private:
    //Kernel G_; //!< kernel
    Eigen::MatrixXd C_; //!< near field block matrix
    std::pair<Node*,Node*> pair_;
};
#endif // BLOCK_NEARF_HPP
