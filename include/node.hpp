#ifndef NODE_HPP
#define NODE_HPP


#include <Eigen/Dense>
#include <cmath>
#include <vector>


class node
{
    typedef std::vector<node*> vector_type;

public:
  
    // default constructor
    node():
        l_child_(NULL), r_child_(NULL), l_ind_ (0), r_ind_ (0), V_(), c_node_(), near_f(0), far_f(0)
    { }
    // single node constructor: adds a tree below the node if left_index != right_index
    node(unsigned left_index, unsigned right_index):
        l_child_(NULL), r_child_(NULL),l_ind_ (left_index), r_ind_ (right_index), V_(), c_node_(), near_f(0), far_f(0)
    {
        add_leaf();
    }

    // tree node constructor
    node(node* left_child, node* right_child, unsigned left_index, unsigned right_index):
        l_child_(left_child), r_child_(right_child), l_ind_(left_index), r_ind_(right_index), V_(),c_node_(), near_f(0), far_f(0)
    { }

    // tree node constructor with V matrix
    node(node* left_child, node* right_child, unsigned left_index, unsigned right_index, Eigen::MatrixXd Vmatrix):
        l_child_(left_child), r_child_(right_child), l_ind_(left_index), r_ind_(right_index), V_(Vmatrix),c_node_(), near_f(0), far_f(0)
    { }

    // destructor
    virtual ~node();

    // returns index of leftmost  point in cluster of node
    unsigned left_ind() const {
        return l_ind_;
    }
    // returns index of rightmost point in cluster of node
    unsigned right_ind() const {
        return r_ind_;
    }
    // returns a pointer to the left  child of the node
    node* leftchild() const {
        return l_child_;
    }
    // returns a pointer to the right child of the node
    node* rightchild() const {
        return r_child_;
    }
    // returns the matrix $V_{\sigma}$, where $\sigma$ denotes the cluster
    Eigen::MatrixXd getVv() const {
        return V_;
    }
//    // returns V*c restricted to the indices of the cluster
//    Eigen::VectorXd getVc() const {
//        return c_node_;
//    }
    // returns a list of pointers to the nodes belonging to near field of the node
    vector_type get_nearf() const {
        return near_f;
    }
    // returns a list of pointers to the nodes belonging to far field of the node
    vector_type get_farf() const {
        return far_f;
    }

//    // set l_ind_ to "left"
//    void set_leftindex(unsigned left){
//        l_ind_ = left;
//    }
//    // set r_ind_ to "right"
//    void set_rightindex(unsigned right) {
//        r_ind_ = right;
//    }
//    // set l_child_ to "leftkid"
//    void set_leftchild(node* leftkid) {
//        l_child_ = leftkid;
//    }
//    // set r_child_ to "rightkid"
//    void set_rightchild(node* rightkid) {
//        r_child_ = rightkid;
//    }

    // build tree (recursively), used in constructors of class ctree
    void add_leaf();
    // build V-matrices
    void fill_V(const Eigen::VectorXd& x, int deg);
    // V*c restricted to node indices
    void Vc(const Eigen::VectorXd& c);
    // add node pointer to near field list
    void push_nearf(node* near_node){near_f.push_back(near_node);}
    // add pointer to node to far field list
    void push_farf(node* far_node){far_f.push_back(far_node);}
  
private:

    void destroy_tree(node *leaf);
    void destroy_tree();

    node* l_child_; // left child of node
    node* r_child_;	// right child of node
    unsigned l_ind_; // smallest index in cluster of node
    unsigned r_ind_; // largest index in cluster of node
    Eigen::MatrixXd V_; // $V_{\sigma}$ where $\sigma$ is the cluster of node
    Eigen::VectorXd c_node_; // $V_{\sigma}* c_{|\sigma}$ where $\sigma$ is the cluster of node
    vector_type near_f; // list with pointers to nodes of the near field of the node
    vector_type far_f; // list with pointers to nodes of the far field of the node
    friend class ctree; // our friend and partner :-)
 };

#endif // NODE_HPP
