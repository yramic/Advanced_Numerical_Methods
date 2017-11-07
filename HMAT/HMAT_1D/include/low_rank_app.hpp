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

#include <Eigen/Dense>
#include "block_cluster.hpp"
#include "block_nearf.hpp"
#include "hierarchical_partition.hpp"
#include "kernel.hpp"
#include "point.hpp"

/**
* \brief Master class for low-rank approximation (Far and Near Field distribution computation)
*/
class LowRankApp
{
public:

    /*!
     * \brief Constructor for 1D Low Rank Approximation
     * \param kernel Kernel used for the matrix multiplication
     * \param Gpoints Vector of points in space
     * \param eta eta-admissibility constant
     * \param deg Degree of interpolation
     */
    LowRankApp(Kernel kernel, const std::vector<Point> &GPoints, double eta, unsigned deg);

    /*!
     * \brief Approximate matrix-vector multiplication
     * \param c Vector c
     */
    Eigen::VectorXd mvProd(const Eigen::VectorXd& c);

private:
    /*!
     * \brief Pre-processing: initialize matrix V and vector Vc for all far field nodes
     * \param ff_v_x Vector of Far Field XNodes
     * \param ff_v_y Vector of Far Field YNodes
     * \param c Vector c to multiply
     */
    void preProcess(std::vector<Node*> ff_v_x, std::vector<Node*> ff_v_y, const Eigen::VectorXd& c);

    /*!
     * \brief Block-processing: compute vector CVc for all far field pairs and store it into xnode
     * \param ff_v Vector of Far Field Pairs
     */
    void blockProcess(std::vector<BlockCluster> ff_v);

    /*!
     * \brief Post-processing: compute vector Vx*CVc for all far field xnodes and add it to vector f in the right place
     * \param ff_v_x Vector of Far Field XNodes
     * \param f Output product vector
     */
    void postProcess(std::vector<Node*> ff_v_x, Eigen::VectorXd& f);

    /*!
     * \brief Compute far field contribution
     * \param ff_v Vector of Far Field Pairs
     * \param ff_v Vector of Far Field unique xnodes
     * \param ff_v Vector of Far Field unique ynodes
     * \param c Vector c to multiply
     * \param f Output product vector
     */
    void ff_contribution(std::vector<BlockCluster> ff_v,
                         std::vector<Node*> ff_v_x, std::vector<Node*> ff_v_y,
                         const Eigen::VectorXd& c, Eigen::VectorXd& f);
    /*!
     * \brief Compute near field contribution
     * \param nf_v Vector of Near Field pairs
     * \param c Vector c to multiply
     * \param f Output product vector
     */
    void nf_contribution(std::vector<BlockNearF> nf_v,
                         const Eigen::VectorXd& c, Eigen::VectorXd& f);

    unsigned  deg_; //!< degree of interpolation
    Kernel kernel_; //!< kernel
    HierarchicalPartitioning HP_; //!< Hierarchical Partiotion class for constructing the tree and calculate near and far field nodes
    std::vector<Point> GPoints_;  //!< Vector of points of the axis
};
#endif // LOW_RANK_APP_HPP
