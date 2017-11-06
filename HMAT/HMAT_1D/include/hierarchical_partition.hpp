#ifndef HIERARCHICAL_PARTITION_HPP
#define HIERARCHICAL_PARTITION_HPP

#include "block_cluster.hpp"
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
    HierarchicalPartitioning(const std::vector<Point>& GPoints, double eta, unsigned deg):
        Tx_(GPoints,deg), Ty_(Tx_), eta_(eta)
    {}
    /*!
     * \brief Return the Far Field pairs vector
     */
    std::vector<BlockCluster> getFF() {
        return FarField_;
    }
    /*!
     * \brief Return the Near Field pairs vector
     */
    std::vector<std::pair<Node*,Node*> > getNF() {
        return NearField_;
    }
    /*!
     * \brief Compute the Far and Near Field pairs
     */
    void setNearFar() { setNearFar_recursion(Tx_.getRoot(), Ty_.getRoot(), eta_, Ty_); }
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
    void setNearFar_recursion(Node* xnode, Node* ynode, double eta, cTree Ty);
    cTree Tx_, Ty_; //!< Cluster trees for comparison
    std::vector<BlockCluster> FarField_; //!< Vector for Far Field
    std::vector<std::pair<Node*,Node*> > NearField_; //!< Vector for Near Field
    double eta_;  //!< eta-admissibility constant
};
#endif // HIERARCHICAL_PARTITION_HPP
