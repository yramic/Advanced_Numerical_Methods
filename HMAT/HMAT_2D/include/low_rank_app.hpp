#ifndef LOW_RANK_APP_HPP
#define LOW_RANK_APP_HPP

#include "ctree.hpp"
#include "kernel.hpp"
#include <Eigen/Dense>


/*!
* \brief Master class for low-rank approximation
*/
class LowRankApp
{
public:

    /*!
     * \brief Constructor for 2D Low Rank Approximation
     * \param kernel Kernel used for the matrix multiplication
     * \param pp Vector of points in space
     * \param n Number of points
     */
    LowRankApp(Kernel* kernel, std::vector<Point> pp, int n);

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
     * \param tx Node to process
     * \param deg Degree of interpolation
     * \param c Vector c to multiply
     * \param f_aprox_ff_contr Number of far field contributions
     */
    void ff_contribution(Eigen::VectorXd& f, Node* tx, unsigned deg, Eigen::VectorXd& c, Eigen::VectorXd& f_aprox_ff_contr);

    /*!
     * \brief Compute near field contribution
     * \param f Output product vector
     * \param tx Node to process
     * \param c Vector c to multiply
     * \param f_aprox_nf_contr Number of near field contributions
     */
    void nf_contribution(Eigen::VectorXd& f, Node* tx, const Eigen::VectorXd& c, Eigen::VectorXd& f_aprox_nf_contr);

    Kernel* kernel_;    //!< pointer for kernel
    cTree PPointsTree_; //!< cluster tree of ppoints
};

#endif // LOW_RANK_APP_HPP
