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
#ifndef NODE_HPP
#define NODE_HPP

#include <Eigen/Dense>
#include <vector>


// forward declaration to avoid cross-referencing
class cTree;


/**
* \brief Node of a cluster tree ("cTree" class)
*/
class Node
{
    typedef std::vector<Node*> vector_t;

public:
    /**
    * \brief Constructors
    */
    // default constructor
    Node():
        l_child_(NULL), r_child_(NULL), l_ind_(0), r_ind_(0), near_f_(0), far_f_(0)
    { }
    // actual  constructor: adds a tree below the node if left_index != right_index
    Node(unsigned l_ind, unsigned r_ind);

    // destructor
    virtual ~Node();

    /**
    * \brief Read access methods
    */
    // return index of leftmost  point in cluster of node
    unsigned getLInd() const { return l_ind_;  }
    // return index of rightmost point in cluster of node
    unsigned getRInd() const {  return r_ind_; }
    // return a pointer to the left  child of the node
    Node* getLChild() const { return l_child_; }
    // return a pointer to the right child of the node
    Node* getRChild() const {  return r_child_; }
    // return matrix $V_{\sigma}$, where $\sigma$ denotes the cluster
    Eigen::MatrixXd getV_node() const { return V_node_; }
    // return V*c restricted to node indices of the cluster
    Eigen::VectorXd getVc_node() const { return Vc_node_; }
    // return a list of pointers to the nodes belonging to the near field of the node
    vector_t getNearF() const { return near_f_; }
    // return a list of pointers to the nodes belonging to the  far field of the node
    vector_t getFarF() const { return far_f_;  }

    /**
    * \brief Setters
    */
    // build tree recursively
    void setLeaves();
    // compute V-matrix of cluster
    void setV_node( const Eigen::VectorXd& x, unsigned deg);
    // compute V*c restricted to node indices of the cluster
    void setVc_node(const Eigen::VectorXd& c);
  
private:
    Node* l_child_;  // left  child of node
    Node* r_child_;  // right child of node
    unsigned l_ind_; // smallest index in cluster of node
    unsigned r_ind_; // largest  index in cluster of node
    Eigen::MatrixXd V_node_;  // $V_{\sigma}$, where $\sigma$ is the cluster of node
    Eigen::VectorXd Vc_node_; // $V_{\sigma}* c_{|\sigma}$, where $\sigma$ is the cluster of node
    vector_t near_f_; // list with pointers to nodes of the near field of the node
    vector_t  far_f_; // list with pointers to nodes of the  far field of the node
    friend class cTree;
 };

#endif // NODE_HPP
