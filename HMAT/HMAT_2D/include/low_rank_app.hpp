#ifndef LOW_RANK_APP_HPP
#define LOW_RANK_APP_HPP

#include "ctree.hpp"
#include "kernel.hpp"
#include "hierarchical_partition.hpp"
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
     * \param n Number of points
     * \param eta eta variable of admissibility
     * \param deg Degree of interpolation
     */
    LowRankApp(Kernel* kernel,const std::vector<Point> &pp, int n, double eta, unsigned deg);
    /*!
     * \brief Approximate matrix-vector multiplication
     * \param c Vector c
     * \param eta eta variable of admissibility
     * \param deg Degree of itnerpolation
     */
    Eigen::VectorXd mvProd(Eigen::VectorXd &c, double eta, unsigned deg);

private:
    /*!
     * \brief Compute far field contribution
     * \param f Output product vector
     * \param ff_v Vector of Far Field Pairs
     * \param deg Degree of interpolation
     * \param c Vector c to multiply
     * \param f_aprox_ff_contr Number of far field contributions
     */
    void ff_contribution(Eigen::VectorXd& f, std::vector<std::pair<Node*,Node*>> ff_v, unsigned deg, Eigen::VectorXd& c, Eigen::VectorXd& f_aprox_ff_contr);
    /*!
     * \brief Compute near field contribution
     * \param f Output product vector
     * \param nf_v Vector of Near Field pairs
     * \param c Vector c to multiply
     * \param f_aprox_nf_contr Number of near field contributions
     */
    void nf_contribution(Eigen::VectorXd& f, std::vector<std::pair<Node*,Node*>> nf_v, const Eigen::VectorXd& c, Eigen::VectorXd& f_aprox_nf_contr);

    Kernel* kernel_;    //!< pointer for kernel
    HierarchicalPartitioning HP_;   //!< Hierarchical Partiotion class for constructing the tree and calculate near and far field nodes
    unsigned deg_;  //!< degree of interpolation
};

#endif // LOW_RANK_APP_HPP
