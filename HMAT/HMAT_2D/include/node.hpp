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
    Node():
        tl_child_(NULL), tr_child_(NULL), bl_child_(NULL), br_child_(NULL), PPointsTree_(std::vector<Point>()), near_f_(vector_t()), far_f_(vector_t())
    { }

    /*!
     * \brief Constructor for the 2D problem
     * \details Actual  Constructor: creates the root of the Cluster Tree and then recursivly creates the leaves
     * \param PPointsTree Vector of Polygon points
     */
    Node(const std::vector<Point> PPointsTree);

    /*!
     * \brief Constructor for the 2D problem
     * \details Actual  Constructor: creates the leaves of the Cluster Tree
     * \param PPointsTree Vector of Polygon points
     * \param x1 x coordinate of left edge of cluster
     * \param x2 x coordinate of right edge of cluster
     * \param y1 y coordinate of bottom edge of cluster
     * \param y2 y coordinate of top edge of cluster
     */
    Node(const std::vector<Point> PPointsTree, double x1, double x2, double y1, double y2);

    /*!
     * \brief Default Destructor
     */
    virtual ~Node();

    /*!
     * \brief return a pointer to the top left child of the node
     */
    Node* getTl_Child() const {
        return tl_child_;
    }

    /*!
     * \brief return a pointer to the top right child of the node
     */
    Node* getTr_Child() const {
        return tr_child_;
    }

    /*!
     * \brief return a pointer to the bottom left  child of the node
     */
    Node* getBl_Child() const {
        return bl_child_;
    }

    /*!
    * \brief return a pointer to the bottom right child of the node
    */
    Node* getBr_Child() const {
        return br_child_;
    }

    /*!
     * \brief return matrix \f$V_{\sigma}\f$, where \f$\sigma\f$ denotes the cluster
     */
    Eigen::MatrixXd getV_node() const {
        return V_node_;
    }

    /*!
     * \brief return a list of pointers to the nodes belonging to the near field of the node
     */
    vector_t getNearF() const {
        return near_f_;
    }

    /*!
     * \brief return a list of pointers to the nodes belonging to the  far field of the node
     */
    vector_t getFarF() const {
        return far_f_;
    }

    /*!
     * \brief return Bounding Box Xl coordinate
     */
    double getXl_b() const {
        return x1_b_;
    }

    /*!
     * \brief return Bounding Box Xr coordinate
     */
    double getXr_b() const {
        return x2_b_;
    }

    /*!
     * \brief return Bounding Box Yl coordinate
     */
    double getYl_b() const {
        return y1_b_;
    }

    /*!
     * \brief return Bounding Box Yr coordinate
     */
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
    void setLeaves();

    /*!
     * \brief Function for splitting the square space in 4 smaller square spaces
     * \param x1 x coordinate of left edge of cluster
     * \param x2 x coordinate of right edge of cluster
     * \param y1 y coordinate of bottom edge of cluster
     * \param y2 y coordinate of top edge of cluster
     */
    void setLeaves(double x1, double x2, double y1, double y2);

    /*!
     * \brief Compute V-matrix of cluster
     * \param t Vector of points that this node´s cluster contains
     */
    void setV_node(const std::vector<Point>& t, unsigned deg);

    /*!
     * \brief Function to print the cluster tree for debugging
     * \param n Level of the node in the cluster tree(0 = root)
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
    Node* tl_child_;  //!< top left child of node
    Node* tr_child_;  //!< top right child of node
    Node* bl_child_;  //!< bottom left child of node
    Node* br_child_;  //!< bottom right child of node
    std::vector<Point> PPointsTree_; //!< vector of median points of the polygone's edges
    double x1_,x2_,y1_,y2_; //!< cluster coordinates
    double x1_b_,x2_b_,y1_b_,y2_b_; //!< bounding box coordinates
    Eigen::MatrixXd V_node_; //!< \f$V_{\sigma}\f$, where \f$\sigma\f$ is the cluster of node
    vector_t near_f_; //!< list with pointers to nodes of the near field of the node
    vector_t  far_f_; //!< list with pointers to nodes of the  far field of the node
    int nodeId_;    //!< id of the node in the cluster tree used for debugging
    friend class cTree;
 };

#endif // NODE_HPP
