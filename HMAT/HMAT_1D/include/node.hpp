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
#ifndef NODE_HPP
#define NODE_HPP

#include <Eigen/Dense>
#include <iostream>
#include <vector>

#include "point.hpp"

// forward declaration to avoid cross-referencing
template <typename NODE>
class cTree;

/**
* \brief Node of a cluster tree ("cTree" class)
*/
class Node {
  typedef std::vector<Node*> vector_t;

 public:
  /*!
     * \brief Default constructor
     */
  Node() : l_child_(NULL), r_child_(NULL), node_points_(std::vector<Point>()) {}
  /*!
     * \brief Actual constructor
     * \param points Vector of points that is contained in the node
     * \param id ID of the node in the cTree
     * \param deg Degree of interpolation
     */
  Node(std::vector<Point> points, int& id, unsigned deg);
  /*!
     * \brief Destructor
     */
  virtual ~Node() {
    if ((l_child_ == NULL) && (r_child_ == NULL))
      std::cout << "leaves destroyed" << std::endl;
    if (l_child_ != NULL) delete l_child_;
    if (r_child_ != NULL) delete r_child_;
  }
  /*!
     * \brief return vector of node´s points
     */
  std::vector<Point> getPoints() const { return node_points_; }
  /*!
     * \brief return a pointer to the left  child of the node
     */
  Node* getLChild() const { return l_child_; }
  /*!
     * \brief return a pointer to the right child of the node
     */
  Node* getRChild() const { return r_child_; }
  /*!
     * \brief return the Chebyshev nodes of interpolation for the segment of points of this node
     */
  Eigen::VectorXd getTK() const { return tk_; }
  /*!
     * \brief return this node´s Weights of Lagrange polynomial
     */
  Eigen::VectorXd getWK() const { return wk_; }
  /*!
     * \brief return the V matrix
     */
  Eigen::MatrixXd getV_Node() { return V_node_; }
  /*!
     * \brief return the Vc vector
     */
  Eigen::VectorXd getVc_Node() { return Vc_node_; }
  /*!
     * \brief return the CVc vector
     */
  Eigen::VectorXd getCVc_Node() { return CVc_node_; }
  /*!
     * \brief return node´s ID
     */
  int getId() const { return nodeId_; }
  /*!
     * \brief build tree recursively
     * \param id ID of node
     */
  void setSons(int& id);
  /*!
     * \brief compute V matrix
     * \return no. of 'operations' performed
     */
  unsigned setV();
  /*!
     * \brief compute Vc vector
     * \param c vector for multiplication with the matrix
     * \return no. of 'operations' performed
     */
  unsigned setVc(const Eigen::VectorXd& c);
  /*!
     * \brief update CVc vector
     * \param c vector for multiplication with the matrix
     * \return no. of 'operations' performed
     */
  unsigned setCVc(const Eigen::VectorXd& CVc);
  /*!
     * \brief reset CVc vector to zero
     */
  void resetCVc();

 protected:
  Node* l_child_;                   //!< left  child of node
  Node* r_child_;                   //!< right child of node
  std::vector<Point> node_points_;  //!< points of the node
  int nodeId_;                      //!< ID of the node in the Cluster Tree
  unsigned deg_;                    //!< degree of interpolation
  Eigen::VectorXd tk_;              //!< Chebyshev interpolation nodes
  Eigen::VectorXd wk_;              //!< Weights of Lagrange polynomial
  Eigen::MatrixXd V_node_;          //!< V  matrix
  Eigen::VectorXd Vc_node_;         //!< Vc vector
  Eigen::VectorXd CVc_node_;        //!< CVc vector
  template <typename NODE>
  friend class cTree;
};

#endif  // NODE_HPP
