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
#ifndef LOW_RANK_APP_HPP
#define LOW_RANK_APP_HPP
//#define ver1
#define ver2
#include "ctree.hpp"
#include "kernel.hpp"
#include "point.hpp"
#include "hierarchical_partition.hpp"
#include <Eigen/Dense>

#ifdef ver2
/**
* \brief Master class for low-rank approximation
*/
class LowRankApp
{
public:

    /**
    * \brief Constructor
    */
    LowRankApp(Kernel kernel, const std::vector<Point> &GPoints, double eta, unsigned deg);

    // approximate matrix-vector multiplication
    Eigen::VectorXd mvProd(const Eigen::VectorXd& c);

private:
    // compute V-matrix of cluster
    Eigen::MatrixXd setV( Node* x);
    // compute V*c restricted to node indices
    Eigen::MatrixXd setVc(Node* x, const Eigen::VectorXd& c);
    // compute far  field contribution
    void ff_contribution(Eigen::VectorXd& f, std::vector<std::pair<Node*,Node*>> ff_v, const Eigen::VectorXd &c);
    // compute near field contribution
    void nf_contribution(Eigen::VectorXd& f, std::vector<std::pair<Node*,Node*>> nf_v, const Eigen::VectorXd& c);

    unsigned deg_;  // degree of interpolation
    Kernel kernel_; // kernel
    HierarchicalPartitioning HP_;   // Hierarchical Partiotion class for constructing the tree and calculate near and far field nodes
    std::vector<Point> GPoints_;    // Vector of points of the axis
};
#endif
#ifdef ver1
    /**
    * \brief Master class for low-rank approximation
    */
    class LowRankApp
    {
    public:

        /**
        * \brief Constructor
        */
        LowRankApp(Kernel kernel, const Eigen::VectorXd& x, const Eigen::VectorXd& y);

        // approximate matrix-vector multiplication
        Eigen::VectorXd mvProd(const Eigen::VectorXd& c, double eta, unsigned deg);
    private:
        // compute far  field contribution
        void ff_contribution(Eigen::VectorXd& f, Node* tx, unsigned deg);
        // compute near field contribution
        void nf_contribution(Eigen::VectorXd& f, Node* tx, const Eigen::VectorXd& c);

        Kernel kernel_; // kernel
        cTree Tx_; // cluster tree of x-values
        cTree Ty_; // cluster tree of y-values
    };
#endif
#endif // LOW_RANK_APP_HPP
