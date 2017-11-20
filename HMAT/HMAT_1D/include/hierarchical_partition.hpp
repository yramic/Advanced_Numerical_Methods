/***********************************************************************
 *                                                                     *
 * Code for Course "Advanced Numerical Methods for CSE"                *
 * (Prof. Dr. R. Hiptmair)                                             *
 * Author: Daniele Casati                                              *
 * Date: 11/2017                                                       *
 * (C) Seminar for Applied Mathematics, ETH Zurich                     *
 * This code can be freely used for non-commercial purposes as long    *
 * as this header is left intact.                                      *
 ***********************************************************************/
#ifndef HIERARCHICAL_PARTITION_HPP
#define HIERARCHICAL_PARTITION_HPP

#include "block_cluster.hpp"
#include "block_nearf.hpp"
#include "ctree.hpp"
#include "kernel.hpp"
#include "point.hpp"

/*!
 * \brief Master class for doing the Hierarchical Partitioning of the Cluster Tree
 */
template<typename BLOCK_CLUSTER = BlockCluster,
         typename NODE_Y = Node>
class HierarchicalPartitioning
{
public:
    /*!
     * \brief Constructor for the Hierarchical Partitioning Class
     * \param GPoints Grid Points
     * \param eta eta-admissibility constant
     * \param deg Degree of itnerpolation
     */
    HierarchicalPartitioning(const std::vector<Point>& GPoints, double eta, unsigned deg):
        Tx_(GPoints,deg), Ty_(GPoints,deg), eta_(eta)
    {}
    /*!
     * \brief Return the Far Field pairs vector
     */
    std::vector<BlockCluster*> getFF() {
        return FarField_;
    }
    /*!
     * \brief Return the Far Field vector of unique pointers to xnodes
     */
    std::vector<Node*> getFFxnds() {
        return FarField_xnds_;
    }
    /*!
     * \brief Return the Far Field vector of unique pointers to ynodes
     */
    std::vector<Node*> getFFynds() {
        return FarField_ynds_;
    }
    /*!
     * \brief Return the Near Field pairs vector
     */
    std::vector<BlockNearF*> getNF() {
        return NearField_;
    }
    /*!
     * \brief Compute the Far and Near Field pairs
     */
    void setNearFar();
    /*!
     * \brief Return the bounding box corresponding to index i of the far-field vector
     */
    std::pair<std::pair<double,double>,
              std::pair<double,double> > getBB(int i);
private:
    /*!
     * \brief Needed for "setNearFar(...)"
     * \param xnode First node for checking
     * \param ynode Second node for checking
     * \param eta eta admissibility variable
     * \param Ty Tree for reference in the recursion
     */
    void setNearFar_recursion(Node* xnode, Node* ynode, double eta, cTree<NODE_Y> Ty);
    cTree<Node>   Tx_; //!< Cluster tree of xnodes
    cTree<NODE_Y> Ty_; //!< Cluster tree of ynodes
    std::vector<BlockCluster*> FarField_;      //!< Vector for Far Field
    std::vector<Node*>         FarField_xnds_; //!< Vector for Far Field XNodes
    std::vector<Node*>         FarField_ynds_; //!< Vector for Far Field YNodes
    std::vector<BlockNearF*>  NearField_;      //!< Vector for Near Field
    double eta_;  //!< eta-admissibility constant
};
#endif // HIERARCHICAL_PARTITION_HPP
