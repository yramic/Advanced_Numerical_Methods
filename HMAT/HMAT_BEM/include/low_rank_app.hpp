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
#ifndef LOW_RANK_APP_HPP
#define LOW_RANK_APP_HPP

#include "ctree.hpp"
#include "hierarchical_partition.hpp"
#include "kernel.hpp"
#include <Eigen/Dense>

/*!
* \brief Master class for low-rank approximation (Far and Near Field distribution computation)
*/
class LowRankApp
{
public:
    /*!
     * \brief Constructor for 2D Low Rank Approximation
     * \param kernel Kernel used for the matrix multiplication
     * \param pp Vector of points in space
     * \param eta eta variable of admissibility
     * \param deg Degree of interpolation
     */
    LowRankApp(Kernel* kernel,const std::vector<Segment> &pp, double eta, unsigned deg);
    /*!
     * \brief Approximate matrix-vector multiplication
     * \param c Vector c
     */
    Eigen::VectorXd mvProd(Eigen::VectorXd& c);
    /*!
     * \brief Count number of corresponding far field points for each row of the product vector
     * \param ff_v Vector of BlockClusters
     * \param f_approx_ff_contr Vector for saving the number of corresponding far field points for each row of the product vector
     */
    void calc_numb_approx_per_row(std::vector<BlockCluster*> ff_v, Eigen::VectorXd& f_approx_ff_contr);

private:
    /*!
     * \brief Compute far field contribution
     * \param ff_v Vector of Far Field Pairs
     * \param ff_v Vector of Far Field unique xnodes
     * \param ff_v Vector of Far Field unique ynodes
     * \param c Vector c to multiply
     * \param f Output product vector
     * \param f_aprox_ff_contr Number of far field contributions
     */
    void ff_contribution(std::vector<BlockCluster*> ff_v,
                         std::vector<Node*> ff_v_x, std::vector<Node*> ff_v_y,
                         const Eigen::VectorXd& c, Eigen::VectorXd& f, Eigen::VectorXd& f_aprox_ff_contr);
    /*!
     * \brief Compute near field contribution
     * \param nf_v Vector of Near Field pairs
     * \param c Vector c to multiply
     * \param f Output product vector
     * \param f_aprox_nf_contr Number of near field contributions
     */
    void nf_contribution(std::vector<BlockNearF*> nf_v,
                         const Eigen::VectorXd& c, Eigen::VectorXd& f, Eigen::VectorXd& f_aprox_nf_contr);
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
    void blockProcess(std::vector<BlockCluster*> ff_v);
    /*!
     * \brief Post-processing: compute vector Vx*CVc for all far field xnodes and add it to vector f in the right place
     * \param ff_v_x Vector of Far Field XNodes
     * \param f Output product vector
     */
    void postProcess(std::vector<Node*> ff_v_x, Eigen::VectorXd& f);

    Kernel* kernel_; //!< pointer for kernel
    HierarchicalPartitioning HP_; //!< Hierarchical Partiotion class for constructing the tree and calculate near and far field nodes
    unsigned   deg_; //!< degree of interpolation
    unsigned  nops_; //!< number of 'operations' performed
};

#endif // LOW_RANK_APP_HPP
