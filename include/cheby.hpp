#ifndef CHEBY_HPP
#define CHEBY_HPP

#include <Eigen/Dense>
#include <cmath>
#include <vector>


class Cheby
{
public:

    // default constructor
    Cheby():
        xl_(0), xr_(0), deg_(0), tk_(0), omega_(Eigen::VectorXd::Zero(0))
    { }
    // meaningful constructor
    Cheby(double x_l, double x_r, unsigned degree_pol);

    // returns the left  boundary xl of the domain
    double left_bd() const {
        return xl_;
    }
    // returns the right boundary xr of the domain
    double right_bd() const {
        return xr_;
    }
    // returns the degree deg_pol of the interpolating polynomials
    unsigned get_degree() const {
        return deg_;
    }
    // returns the Chebyshew nodes on domain [xl,xr]
    Eigen::VectorXd cheb_pts() const {
        return tk_;
    }
    // returns the weights of the Lagrange polynomials
    Eigen::VectorXd get_omega() const {
        return omega_;
    }

    // set xl_ to "lbd"
    void set_leftbd(double lbd) {
        xl_ = lbd;
    }
    // set xr_ to "rbd"
    void set_rightbd(double rbd) {
        xr_ = rbd;
    }
    // set deg_ to "deg"
    void set_deg(unsigned deg) {
        deg_ = deg;
    }
    // compute Chebyshew nodes on domain [xl,xr]
    void set_chebpt();
    // PRE: Chebyshew points already computed, compute omega, the weights of the Lagrange polynomials
    void set_omega();

private:

    double xl_; // left boundary of the domain, on which we construct the Chebyshew points
    double xr_; // right boundary of the domain, on which we construct the Chebyshew points
    unsigned deg_; // degree of polynomials, $\in \mathbb{N}_{>0}$ please
    Eigen::VectorXd tk_; // Chebyshew points
    Eigen::VectorXd omega_;  // weights of Lagrange polynomials
};

#endif // CHEBY_HPP
