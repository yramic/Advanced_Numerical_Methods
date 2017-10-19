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
    * \brief Constructors
    */
    LowRankApp(Kernel2D kernel, const Eigen::VectorXd& x, const Eigen::VectorXd& y, int n);
    //LowRankApp(Kernel4D* kernel, std::vector<Point> pp, int n);
    //LowRankApp(PolynomialKernel* kernel, std::vector<Point> pp, int n);
    //LowRankApp(ConstantKernel* kernel, std::vector<Point> pp, int n);
    LowRankApp(Kernel* kernel, std::vector<Point> pp, int n);

    // approximate matrix-vector multiplication
    Eigen::VectorXd mvProd(Eigen::VectorXd &c, double eta, unsigned deg);

private:

    // compute far  field contribution
    void ff_contribution(Eigen::VectorXd& f, Node* tx, unsigned deg, Eigen::VectorXd& c);
    // compute near field contribution
    void nf_contribution(Eigen::VectorXd& f, Node* tx, const Eigen::VectorXd& c);

    Kernel2D kernel2d_; // kernel for 2D
    //Kernel4D kernel4d_; // kernel for 4D
    //PolynomialKernel PolynomialKernel_; // Polynomial Kernel
    //ConstantKernel CKernel_;    // Constant Kernel
    Kernel* kernel_;
    cTree PPointsTree_; // cluster tree of ppoints
    cTree Tx_; // cluster tree of x-values
    cTree Ty_; // cluster tree of y-values
    Eigen::MatrixXd HM_;
};

#endif // LOW_RANK_APP_HPP
