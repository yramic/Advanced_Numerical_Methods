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

#include "ctree.hpp"
#include "kernel.hpp"
#include "point.hpp"
#include "hierarchical_partition.hpp"
#include <Eigen/Dense>

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
     * \brief Compute far field contribution
     * \param f Output product vector
     * \param ff_v Vector of Far Field Pairs
     * \param c Vector c to multiply
     */
    void ff_contribution(Eigen::VectorXd& f, std::vector<std::pair<Node*,Node*>> ff_v, const Eigen::VectorXd &c);
    /*!
     * \brief Compute near field contribution
     * \param f Output product vector
     * \param nf_v Vector of Near Field pairs
     * \param c Vector c to multiply
     */
    void nf_contribution(Eigen::VectorXd& f, std::vector<std::pair<Node*,Node*>> nf_v, const Eigen::VectorXd& c);

    unsigned  deg_; //!< degree of interpolation
    Kernel kernel_; //!< kernel
    HierarchicalPartitioning HP_; //!< Hierarchical Partiotion class for constructing the tree and calculate near and far field nodes
    std::vector<Point> GPoints_;  //!< Vector of points of the axis
};
#endif // LOW_RANK_APP_HPP
