#ifndef NODE_HPP
#define NODE_HPP

#include "ctree.hpp"
#include <Eigen/Dense>
#include <vector>


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
    // actual constructor: adds a tree below the node if left_index != right_index
    Node(unsigned left_index, unsigned right_index);

    // destructor
    virtual ~Node();

    /**
    * \brief Getters
    */
    // return index of leftmost  point in cluster of node
    unsigned getLInd() const {
        return l_ind_;
    }
    // return index of rightmost point in cluster of node
    unsigned getRInd() const {
        return r_ind_;
    }
    // return a pointer to the left  child of the node
    Node* getLChild() const {
        return l_child_;
    }
    // return a pointer to the right child of the node
    Node* getRChild() const {
        return r_child_;
    }
    // return matrix $V_{\sigma}$, where $\sigma$ denotes the cluster
    Eigen::MatrixXd getV() const {
        return V_;
    }
    // return a list of pointers to the nodes belonging to the near field of the node
    vector_t getNearF() const {
        return near_f_;
    }
    // return a list of pointers to the nodes belonging to the  far field of the node
    vector_t getFarF() const {
        return far_f_;
    }

    /**
    * \brief Setters
    */
    // build tree recursively
    void setLeaves();
    // compute V-matrix of node
    void setV(const Eigen::VectorXd& x, unsigned deg);
    // compute V*c restricted to node indices
    void setVc_node(const Eigen::VectorXd& c);
    // add node pointer to near field list
    void push2NearF(Node* near_node){
        near_f_.push_back(near_node);
    }
    // add node pointer to  far field list
    void push2FarF(Node* far_node){
        far_f_.push_back(far_node);
    }
  
private:

    Node* l_child_; // left  child of node
    Node* r_child_;	// right child of node
    unsigned l_ind_; // smallest index in cluster of node
    unsigned r_ind_; // largest  index in cluster of node
    Eigen::MatrixXd V_; // $V_{\sigma}$, where $\sigma$ is the cluster of node
    Eigen::VectorXd Vc_node_; // $V_{\sigma}* c_{|\sigma}$, where $\sigma$ is the cluster of node
    vector_t near_f_; // list with pointers to nodes of the near field of the node
    vector_t  far_f_; // list with pointers to nodes of the  far field of the node
    friend class ctree;
 };

#endif // NODE_HPP
