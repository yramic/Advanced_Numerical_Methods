#ifndef HIERARCHICAL_PARTITION_HPP
#define HIERARCHICAL_PARTITION_HPP

#include "ctree.hpp"
#include "kernel.hpp"
#include "point.hpp"
#include "block_cluster.hpp"
#include "block_nearf.hpp"
#include <iostream>

/*!
 * \brief Master class for doing the Hierarchical Partitioning of the Cluster Tree
 */
class HierarchicalPartitioning
{
public:
    /*!
    * \brief Constructor for the Hierarchical Partitioning Class
    * \param PPoints Polygon Points
    * \param eta eta variable of admissibility
    * \param deg Degree of itnerpolation
    */
    HierarchicalPartitioning(const std::vector<Point> &PPoints, double eta, unsigned deg);
    /*!
     * \brief Return the Far Field pairs vector
     */
    std::vector<BlockCluster*> getFF() {
        return FarField_;
    }
    /*!
     * \brief Return the Far Field vector of unique pointer to xnodes
     */
    std::vector<Node*> getFFxnds() {
        return FarField_xnds_;
    }
    /*!
     * \brief Return the Far Field vector of unique pointer to ynodes
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
    * \brief Return the Cluster Tree
    */
    cTree getTx(){
        return Tx_;
    }

    /*!
    * \brief Compute the Far and Near Field pairs
    */
    void setNearFar();
    /*!
    * \brief Compute the Far and Near Field pairs and creates the cmatrix(matrix that contains which nodes of the cTree were checked) for debugging
    */
    void setNearFar(Eigen::MatrixXd& cmatrix);
    /*!
     * \brief Return the bounding box corresponding to index i of the far-field vector
     */
    std::pair<std::pair<std::pair<double,double>,std::pair<double,double>>,
              std::pair<std::pair<double,double>,std::pair<double,double>> > getBB(int i);
private:
    /*!
    * \brief Needed for "setNearFar(...)"
    * \param xnode First node for checking
    * \param ynode Second node for checking
    * \param eta eta admissibility variable
    * \param Ty Tree for reference in the recursion
    */
    void setNearFar_recursion(Node* xnode, Node* ynode, double eta, cTree &Ty);
    /*!
    * \brief Needed for "setNearFar(...)" for debugging
    * \param xnode First node for checking
    * \param ynode Second node for checking
    * \param eta eta admissibility variable
    * \param Ty Tree for reference in the recursion
    * \param cmatrix Matrix for saving calculations used for checking
    */
    void setNearFar_recursion(Node* xnode, Node* ynode, double eta, cTree &Ty, Eigen::MatrixXd& cmatrix);

    cTree Tx_, Ty_; //!< Cluster trees for comparison
    //std::vector<std::pair<Node*,Node*>> FarField_, NearField_;  //!< Vectors for Near and Far Field Pairs
    std::vector<BlockCluster*> FarField_;      //!< Vector for Far Field
    std::vector<Node*>         FarField_xnds_; //!< Vector for Far Field XNodes
    std::vector<Node*>         FarField_ynds_; //!< Vector for Far Field YNodes
    //std::vector<std::pair<Node*,Node*> > NearField_; //!< Vector for Near Field
    std::vector<BlockNearF*>  NearField_;      //!< Vector for Near Field
    double eta_;    //!< eta-admissibility constant
};

#endif // HIERARCHICAL_PARTITION_HPP
