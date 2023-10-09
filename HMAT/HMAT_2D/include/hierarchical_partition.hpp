/***********************************************************************
 *                                                                     *
 * Code for Course "Advanced Numerical Methods for CSE"                *
 * (Prof. Dr. R. Hiptmair)                                             *
 * Author: Ioannis Magkanaris                                          *
 * Date: 11/2017                                                       *
 * (C) Seminar for Applied Mathematics, ETH Zurich                     *
 * This code can be freely used for non-commercial purposes as long    *
 * as this header is left intact.                                      *
 ***********************************************************************/
#ifndef HIERARCHICAL_PARTITION_HPP
#define HIERARCHICAL_PARTITION_HPP

#include <iostream>

#include "block_cluster.hpp"
#include "block_nearf.hpp"
#include "ctree.hpp"
#include "kernel.hpp"
#include "point.hpp"

/*!
 * \brief Master class for doing the Hierarchical Partitioning of the Cluster
 * Tree
 */
class HierarchicalPartitioning {
 public:
  /*!
   * \brief Constructor for the Hierarchical Partitioning Class
   * \param PPoints Polygon Points
   * \param eta eta variable of admissibility
   * \param deg Degree of itnerpolation
   */
  HierarchicalPartitioning(const std::vector<Point>& PPoints, double eta,
                           unsigned deg);
  /*!
   * \brief Return the Far Field pairs vector
   */
  std::vector<BlockCluster*> getFF() { return FarField_; }
  /*!
   * \brief Return the Far Field vector of unique pointer to xnodes
   */
  std::vector<Node*> getFFxnds() { return FarField_xnds_; }
  /*!
   * \brief Return the Far Field vector of unique pointer to ynodes
   */
  std::vector<Node*> getFFynds() { return FarField_ynds_; }
  /*!
   * \brief Return the Near Field pairs vector
   */
  std::vector<BlockNearF*> getNF() { return NearField_; }
  /*!
   * \brief Return the Cluster Tree
   */
  cTree getTx() { return Tx_; }
  /*!
   * \brief Compute the Far and Near Field pairs
   */
  void setNearFar();
  /*!
   * \brief Compute the Far and Near Field pairs and creates the cmatrix(matrix
   * that contains which nodes of the cTree were checked) for debugging
   */
  void setNearFar(Eigen::MatrixXd& cmatrix);
  /*!
   * \brief Return the bounding box corresponding to index i of the far-field
   * vector
   */
  std::pair<std::pair<std::pair<double, double>, std::pair<double, double>>,
            std::pair<std::pair<double, double>, std::pair<double, double>>>
  getBB(int i);

 private:
  /*!
   * \brief Needed for "setNearFar(...)"
   * \param xnode First node for checking
   * \param ynode Second node for checking
   * \param eta eta admissibility variable
   */
  void setNearFar_recursion(Node* xnode, Node* ynode, double eta);
  /*!
   * \brief Needed for "setNearFar(...)" for debugging
   * \param xnode First node for checking
   * \param ynode Second node for checking
   * \param eta eta admissibility variable
   * \param cmatrix Matrix for saving calculations used for checking
   */
  void setNearFar_recursion(Node* xnode, Node* ynode, double eta,
                            Eigen::MatrixXd& cmatrix);

  cTree Tx_, Ty_;                        //!< Cluster trees for comparison
  std::vector<BlockCluster*> FarField_;  //!< Vector for Far Field
  std::vector<Node*> FarField_xnds_;     //!< Vector for Far Field XNodes
  std::vector<Node*> FarField_ynds_;     //!< Vector for Far Field YNodes
  std::vector<BlockNearF*> NearField_;   //!< Vector for Near Field
  double eta_;                           //!< eta-admissibility constant
};

#endif  // HIERARCHICAL_PARTITION_HPP
