#ifndef LOW_RANK_APP_HPP
#define LOW_RANK_APP_HPP

#include <Eigen/Dense>
#include <vector>
#include "BC.hpp"
#include "ctree.hpp"
#include "node.hpp"


class LowRankApp
{
public:

    LowRankApp( Kernel kernel, const Eigen::VectorXd& x, const Eigen::VectorXd& y );

    // Approximate matrix-vector multiplication
    Eigen::VectorXd mvProd( const Eigen::VectorXd& c, double eta, unsigned deg );

private:

    // Compute far-field contribution
    void ff_contribution( Eigen::VectorXd& f, node* tx, int deg );

    // Compute near-field contribution
    void nf_contribution( Eigen::VectorXd& f, node* tx, const Eigen::VectorXd& c_ );

    Kernel kernel_;
    cTree Tx_;
    cTree Ty_;
};

#endif // LOW_RANK_APP_HPP
