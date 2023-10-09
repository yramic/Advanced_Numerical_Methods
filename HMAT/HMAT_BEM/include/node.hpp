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
#include <vector>

#include "segment.hpp"

// forward declaration to avoid cross-referencing
class cTree;

/**
 * \brief Node of a cluster tree ("cTree" class)
 */
class Node {
  typedef std::vector<Node*> vector_t;

 public:
  /*!
   * \brief Default Constructor
   */
  Node()
      : tl_child_(NULL),
        tr_child_(NULL),
        bl_child_(NULL),
        br_child_(NULL),
        segments_(std::vector<Segment>()) {}
  /*!
   * \brief Constructor for the 2D problem
   * \details Actual Constructor: creates the root of the Cluster Tree and then
   * recursivly creates the leaves \param segments Vector of segments \param deg
   * Degree of interpolation
   */
  Node(const std::vector<Segment> segments, unsigned deg);
  /*!
   * \brief Constructor for the 2D problem
   * \details Actual Constructor: creates the leaves of the Cluster Tree
   * \param segments Vector of segments
   * \param x1 x coordinate of left edge of cluster
   * \param x2 x coordinate of right edge of cluster
   * \param y1 y coordinate of bottom edge of cluster
   * \param y2 y coordinate of top edge of cluster
   * \param deg Degree of interpolation
   */
  Node(const std::vector<Segment> segments, double x1, double x2, double y1,
       double y2, unsigned deg);
  /*!
   * \brief Default Destructor
   */
  virtual ~Node();
  /*!
   * \brief return a pointer to the top left child of the node
   */
  Node* getTl_Child() const { return tl_child_; }
  /*!
   * \brief return a pointer to the top right child of the node
   */
  Node* getTr_Child() const { return tr_child_; }
  /*!
   * \brief return a pointer to the bottom left  child of the node
   */
  Node* getBl_Child() const { return bl_child_; }
  /*!
   * \brief return a pointer to the bottom right child of the node
   */
  Node* getBr_Child() const { return br_child_; }
  /*!
   * \brief return Bounding Box Xl coordinate
   */
  double getX1() const { return x1_; }
  /*!
   * \brief return Bounding Box Xr coordinate
   */
  double getX2() const { return x2_; }
  /*!
   * \brief return Bounding Box Yl coordinate
   */
  double getY1() const { return y1_; }
  /*!
   * \brief return Bounding Box Yr coordinate
   */
  double getY2() const { return y2_; }
  /*!
   * \brief return the vector of segments of this node
   */
  std::vector<Segment> getSegments() const { return segments_; }
  /*!
   * \brief evaluate Lagrange polynomial
   * \param j Index of polynomial
   * \param tk x or y coordinate for eval
   */
  double evalLagrange(unsigned j, double tk);
  /*!
   * \brief return the V matrix of this node
   */
  Eigen::MatrixXd getV_Node() const { return V_node_; }
  /*!
   * \brief return the Vc vector
   */
  Eigen::VectorXd getVc_Node() { return Vc_node_; }
  /*!
   * \brief return the CVc vector
   */
  Eigen::VectorXd getCVc_Node() { return CVc_node_; }
  /*!
   * \brief return the id of this node
   */
  int getNodeID() { return nodeId_; }
  /*!
   * \brief return Chebyshev Nodes for the x axis of the Bounding Box
   */
  Eigen::VectorXd getTkx() const { return tkx_; }
  /*!
   * \brief return Chebyshev Nodes for the y axis of the Bounding Box
   */
  Eigen::VectorXd getTky() const { return tky_; }
  /*!
   * \brief return weights of Lagrange polynomial for the x axis of the bounding
   * box
   */
  Eigen::VectorXd getWkx() const { return wkx_; }
  /*!
   * \brief return weights of Lagrange polynomial for the y axis of the bounding
   * box
   */
  Eigen::VectorXd getWky() const { return wky_; }
  /*!
   * \brief return rectangle of the cluster
   */
  void getRect();
  /*!
   * \brief Build tree recursively
   */
  void setSons();
  /*!
   * \brief Compute V matrix
   * \param deg Degree of interpolation
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
   * \brief Function to print the cluster tree for debugging
   * \param n Level of the node in the cluster tree(0 = root)
   */
  void printree(int n);
  /*!
   * \brief Setter for X1 coordinate of the bounding box
   */
  void setX1_b(double x1) { x1_ = x1; }
  /*!
   * \brief Setter for X1 coordinate of the bounding box
   */
  void setX2_b(double x2) { x2_ = x2; }
  /*!
   * \brief Setter for X2 coordinate of the bounding box
   */
  void setY1_b(double y1) { y1_ = y1; }
  /*!
   * \brief Setter for Y2 coordinate of the bounding box
   */
  void setY2_b(double y2) { y2_ = y2; }
  /*!
   * \brief Setter for id of the node
   */
  void setNodeId(double n) { nodeId_ = n; }

 private:
  Node* tl_child_;                 //!< top left child of node
  Node* tr_child_;                 //!< top right child of node
  Node* bl_child_;                 //!< bottom left child of node
  Node* br_child_;                 //!< bottom right child of node
  std::vector<Segment> segments_;  //!< vector of segments
  double x1_, x2_, y1_, y2_;       //!< bounding box coordinates
  unsigned deg_;                   //!< degree of interpolation
  Eigen::VectorXd tkx_, tky_;      //!< Chebyshev Nodes
  Eigen::VectorXd wkx_, wky_;      //!< Lagrange polynomial weights
  Eigen::MatrixXd V_node_;         //!< V matrix
  Eigen::VectorXd Vc_node_;        //!< Vc vector
  Eigen::VectorXd CVc_node_;       //!< CVc vector
  int nodeId_;                     //!< id of the node in the cluster tree
  friend class cTree;
};

#endif  // NODE_HPP
