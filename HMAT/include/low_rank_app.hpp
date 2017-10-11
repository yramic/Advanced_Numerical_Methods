#ifndef LOW_RANK_APP_HPP
#define LOW_RANK_APP_HPP

#include "ctree.hpp"
#include "kernel.hpp"
#include <Eigen/Dense>


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
    LowRankApp(Kernel kernel, std::vector<Point> pp);

    // approximate matrix-vector multiplication
    Eigen::VectorXd mvProd(const Eigen::VectorXd& c, double eta, unsigned deg);

private:

    // compute far  field contribution
    void ff_contribution(Eigen::VectorXd& f, Node* tx, unsigned deg);
    // compute near field contribution
    void nf_contribution(Eigen::VectorXd& f, Node* tx, const Eigen::VectorXd& c);

    Kernel kernel_; // kernel
    cTree PPointsTree_; // cluster tree of ppoints
    cTree Tx_; // cluster tree of x-values
    cTree Ty_; // cluster tree of y-values

};

#endif // LOW_RANK_APP_HPP
