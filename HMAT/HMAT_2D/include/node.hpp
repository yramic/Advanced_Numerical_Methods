#ifndef NODE_HPP
#define NODE_HPP

#include <Eigen/Dense>
#include <vector>
#include "point.hpp"

// forward declaration to avoid cross-referencing
class cTree;


/**
* \brief Node of a cluster tree ("cTree" class)
*/
class Node
{
    typedef std::vector<Node*> vector_t;

public:

    /*!
    * \brief Default Constructor
    */
    // default constructor
    Node():
        tl_child_(NULL), tr_child_(NULL), bl_child_(NULL), br_child_(NULL), PPointsTree_(std::vector<Point>()), near_f_(vector_t()), far_f_(vector_t())
    { }
    // actual  constructor: adds a tree below the node if left_index != right_index for the 2d problem(grid)
    /*!
    * \brief Constructor for the 2D problem
    * \details Actual  Constructor: adds a tree below the node if left_index != right_index for the 2d problem(grid)
    */
    Node(unsigned l_ind, unsigned r_ind);
    // actual  constructor: creates the root of the Cluster Tree and the recursivly creates the leaves
    /*!
    * \brief Constructor for the 4D problem
    * \details Actual  Constructor: creates the root of the Cluster Tree and then recursivly creates the leaves
    */
    Node(const std::vector<Point> PPointsTree);
    // recursive constructor for the leaves of the Cluster Tree
    /*!
    * \brief Constructor for the 4D problem
    * \details Actual  Constructor: creates the leaves of the Cluster Tree
    */
    Node(const std::vector<Point> PPointsTree, double x1, double x2, double y1, double y2);
    // destructor
    /*!
    * \brief Default Destructor
    */
    virtual ~Node();

    /*!
    * \brief return index of leftmost  point in cluster of node
    */
    // return index of leftmost  point in cluster of node
    unsigned getLInd() const {
        return l_ind_;
    }
    /*!
    * \brief return index of rightmost point in cluster of node
    */
    // return index of rightmost point in cluster of node
    unsigned getRInd() const {
        return r_ind_;
    }
    /*!
    * \brief return a pointer to the left  child of the node
    */
    // return a pointer to the left  child of the node
    Node* getLChild() const {
        return l_child_;
    }
    /*!
    * \brief return a pointer to the right child of the node
    */
    // return a pointer to the right child of the node
    Node* getRChild() const {
        return r_child_;
    }
    /*!
    * \brief return a pointer to the top left child of the node
    */
    // return a pointer to the top left child of the node
    Node* getTl_Child() const {
        return tl_child_;
    }
    /*!
    * \brief return a pointer to the top right child of the node
    */
    // return a pointer to the top right child of the node
    Node* getTr_Child() const {
        return tr_child_;
    }
    /*!
    * \brief return a pointer to the bottom left  child of the node
    */
    // return a pointer to the bottom left  child of the node
    Node* getBl_Child() const {
        return bl_child_;
    }
    /*!
    * \brief return a pointer to the bottom right child of the node
    */
    // return a pointer to the bottom right child of the node
    Node* getBr_Child() const {
        return br_child_;
    }
    /*!
    * \brief return matrix $V_{\sigma}$, where $\sigma$ denotes the cluster
    */
    // return matrix $V_{\sigma}$, where $\sigma$ denotes the cluster
    Eigen::MatrixXd getV_node() const {
        return V_node_;
    }
    /*!
    * \brief return V*c restricted to node indices of the cluster
    */
    // return V*c restricted to node indices of the cluster
    Eigen::VectorXd getVc_node() const {
        return Vc_node_;
    }
    /*!
    * \brief return a list of pointers to the nodes belonging to the near field of the node
    */
    // return a list of pointers to the nodes belonging to the near field of the node
    vector_t getNearF() const {
        return near_f_;
    }
    /*!
    * \brief return a list of pointers to the nodes belonging to the  far field of the node
    */
    // return a list of pointers to the nodes belonging to the  far field of the node
    vector_t getFarF() const {
        return far_f_;
    }
    /*!
    * \brief return Bounding Box Xl coordinate
    */
    // return Bounding Box x1
    double getXl_b() const {
        return x1_b_;
    }
    /*!
    * \brief return Bounding Box Xr coordinate
    */
    // return Bounding Box x2
    double getXr_b() const {
        return x2_b_;
    }
    /*!
    * \brief return Bounding Box Yl coordinate
    */
    // return Bounding Box y1
    double getYl_b() const {
        return y1_b_;
    }
    /*!
    * \brief return Bounding Box Yr coordinate
    */
    // return Bounding Box y2
    double getYr_b() const {
        return y2_b_;
    }
    /*!
    * \brief return the vector of points of this node
    */
    std::vector<Point> getPPoints() const {
        return PPointsTree_;
    }
    /*!
    * \brief return the id of this node
    */
    int getNodeID(){
        return nodeId_;
    }

    /*!
    * \brief Build tree recursively
    */
    // build tree recursively
    void setLeaves();
    /*!
    * \brief Function for splitting the square space in 4 smaller square spaces
    */
    void setLeaves(double x1, double x2, double y1, double y2);
    // compute V-matrix of cluster
    //void setV_node( const Eigen::VectorXd& x, unsigned deg);
    /*!
    * \brief Compute V-matrix of cluster
    */
    void setV_node(const std::vector<Point>& t, unsigned deg);
    /*!
    * \brief Compute V*c restricted to node indices of the cluster
    */
    // compute V*c restricted to node indices of the cluster
    void setVc_node(const Eigen::VectorXd& c);
    /*!
    * \brief Function to print the cluster tree for debugging
    */
    void printree(int n);
    /*!
    * \brief Setter for X1 coordinate of the bounding box
    */
    void setX1_b(double x1) {
        x1_b_=x1;
    }
    /*!
    * \brief Setter for X1 coordinate of the bounding box
    */
    void setX2_b(double x2) {
        x2_b_=x2;
    }
    /*!
    * \brief Setter for X2 coordinate of the bounding box
    */
    void setY1_b(double y1) {
        y1_b_=y1;
    }
    /*!
    * \brief Setter for Y2 coordinate of the bounding box
    */
    void setY2_b(double y2) {
        y2_b_=y2;
    }
    /*!
    * \brief Setter for id of the node
    */
    void setNodeId(double n){
        nodeId_ = n;
    }

private:

    Node* l_child_; //!< left  child of node
    Node* r_child_;	//!< right child of node
    Node* tl_child_;  //!< top left child of node
    Node* tr_child_;  //!< top right child of node
    Node* bl_child_;  //!< bottom left child of node
    Node* br_child_;  //!< bottom right child of node
    /* maybe it´s not needed in node*/std::vector<Point> PPointsTree_; //!< vector of median points of the polygone's edges
    double x1_,x2_,y1_,y2_; //!< cluster coordinates
    double x1_b_,x2_b_,y1_b_,y2_b_; //!< bounding box coordinates
    unsigned l_ind_; //!< smallest index in cluster of node
    unsigned r_ind_; //!< largest  index in cluster of node
    Eigen::MatrixXd V_node_; //!< $V_{\sigma}$, where $\sigma$ is the cluster of node
    Eigen::VectorXd Vc_node_; //!< $V_{\sigma}* c_{|\sigma}$, where $\sigma$ is the cluster of node
    vector_t near_f_; //!< list with pointers to nodes of the near field of the node
    vector_t  far_f_; //!< list with pointers to nodes of the  far field of the node
    int nodeId_;    //!< id of the node in the cluster tree used for debugging
    friend class cTree;
 };

#endif // NODE_HPP
