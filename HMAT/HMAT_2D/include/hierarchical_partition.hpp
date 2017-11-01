#ifndef HIERARCHICAL_PARTITION_HPP
#define HIERARCHICAL_PARTITION_HPP

#include "ctree.hpp"
#include "kernel.hpp"
#include "point.hpp"
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
    * \brief Return the Cluster Tree
    */
    cTree getTx(){
        return Tx_;
    }

    /*!
    * \brief Compute the Far and Near Field pairs
    */
    void setNearFar() {
        setNearFar_recursion(Tx_.getRoot(), Ty_.getRoot(), eta_, Ty_);
    }
    /*!
    * \brief Compute the Far and Near Field pairs and creates the cmatrix(matrix that contains which nodes of the cTree were checked) for debugging
    */
    void setNearFar(Eigen::MatrixXd& cmatrix) {
        setNearFar_recursion(Tx_.getRoot(), Ty_.getRoot(), eta_, Ty_, cmatrix);
    }

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
    std::vector<std::pair<Node*,Node*>> FarField_, NearField_;  //!< Vectors for Near and Far Field Pairs
    unsigned deg_;  //!< degree of interpolation
    double eta_;    //!< eta-admissibility constant
};
#endif // HIERARCHICAL_PARTITION_HPP
