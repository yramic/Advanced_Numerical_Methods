#ifndef HIERARCHICAL_PARTITION_HPP
#define HIERARCHICAL_PARTITION_HPP

#include "ctree.hpp"
#include "kernel.hpp"
#include "point.hpp"

/*!
 * \brief Master class for doing the Hierarchical Partitioning of the Cluster Tree
 */
class HierarchicalPartitioning
{
public:
    /*!
     * \brief Constructor for the Hierarchical Partitioning Class
     * \param GPoints Grid Points
     * \param eta eta-admissibility constant
     * \param deg Degree of itnerpolation
     */
    HierarchicalPartitioning(const std::vector<Point> &GPoints, double eta, unsigned deg):
        Tx_(GPoints,deg), Ty_(Tx_), eta_(eta), deg_(deg)
    {}
    /*!
     * \brief Return the Far Field pairs vector
     */
    std::vector<std::pair<Node*,Node*>> getFF(){
        return FarField_;
    }
    /*!
     * \brief Return the Near Field pairs vector
     */
    std::vector<std::pair<Node*,Node*>> getNF(){
        return NearField_;
    }
    /*!
     * \brief Compute the Far and Near Field pairs
     */
    void setNearFar() { setNearFar_recursion(Tx_.getRoot(), Ty_.getRoot(), eta_, Ty_); }
private:
    /*!
     * \brief Needed for "setNearFar(...)"
     * \param xnode First node for checking
     * \param ynode Second node for checking
     * \param eta eta admissibility variable
     * \param Ty Tree for reference in the recursion
     */
    void setNearFar_recursion(Node* xnode, Node* ynode, double eta, cTree Ty);
    cTree Tx_, Ty_; //!< Cluster trees for comparison
    std::vector<std::pair<Node*,Node*>> FarField_, NearField_;  //!< Vectors for Near and Far Field vectors
    unsigned deg_;  //!< degree of interpolation
    double eta_;    //!< eta-admissibility constant
};
#endif // HIERARCHICAL_PARTITION_HPP
