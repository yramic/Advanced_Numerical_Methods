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
    //void ff_contribution(Eigen::VectorXd& f, std::vector<std::pair<Node*,Node*>> ff_v, unsigned deg, Eigen::VectorXd& c, Eigen::VectorXd& f_aprox_ff_contr);
    /*!
     * \brief Compute near field contribution
     * \param f Output product vector
     * \param nf_v Vector of Near Field pairs
     * \param c Vector c to multiply
     * \param f_aprox_nf_contr Number of near field contributions
     */
    //void nf_contribution(Eigen::VectorXd& f, std::vector<std::pair<Node*,Node*>> nf_v, const Eigen::VectorXd& c, Eigen::VectorXd& f_aprox_nf_contr);

    /*!
     * \brief Compute far field contribution
     * \param ff_v Vector of Far Field Pairs
     * \param ff_v Vector of Far Field unique xnodes
     * \param ff_v Vector of Far Field unique ynodes
     * \param c Vector c to multiply
     * \param f Output product vector
     * \param f_aprox_ff_contr Number of far field contributions
     */
    void ff_contribution(std::vector<BlockCluster> ff_v,
                         std::vector<Node *> ff_v_x, std::vector<Node *> ff_v_y,
                         const Eigen::VectorXd& c, Eigen::VectorXd& f, Eigen::VectorXd& f_aprox_ff_contr);
    /*!
     * \brief Compute near field contribution
     * \param nf_v Vector of Near Field pairs
     * \param c Vector c to multiply
     * \param f Output product vector
     * \param f_aprox_nf_contr Number of near field contributions
     */
    void nf_contribution(std::vector<std::pair<Node*,Node*> > nf_v,
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
    void blockProcess(std::vector<BlockCluster> ff_v);

    /*!
     * \brief Post-processing: compute vector Vx*CVc for all far field xnodes and add it to vector f in the right place
     * \param ff_v_x Vector of Far Field XNodes
     * \param f Output product vector
     */
    void postProcess(std::vector<Node*> ff_v_x, Eigen::VectorXd& f);

    Kernel* kernel_;    //!< pointer for kernel
    HierarchicalPartitioning HP_;   //!< Hierarchical Partiotion class for constructing the tree and calculate near and far field nodes
    unsigned deg_;  //!< degree of interpolation
};

#endif // LOW_RANK_APP_HPP
