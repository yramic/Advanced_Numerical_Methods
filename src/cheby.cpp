#include "../include/cheby.hpp"
#include <Eigen/Dense>


//  meaningful constructor
Cheby::Cheby(double x_l, double x_r, unsigned degree_pol):
    xl_(x_l), xr_(x_r), deg_(degree_pol), tk_(0), omega_(Eigen::VectorXd::Zero(degree_pol+1))
{
    double val = 0;
    for(unsigned j=0; j<=deg_; ++j) { // calculate Chebyshew nodes
        val = (xl_+xr_)/2 + (xr_-xl_)/2 * cos((2.*j+1)/(2*(deg_+1.))*M_PI);
        tk_.push_back(val);
    }
    set_omega();
}


// define Chebyshev nodes
void Cheby::set_chebpt()
{
    double val = 0;
    for(unsigned j=0; j<=deg_; ++j) {
        val = (xl_+xr_)/2 + (xr_-xl_)/2 * cos((2.*j+1)/(2.*(deg_+1.))*M_PI);
        tk_.push_back(val);
    }
}


// calculate omega (factors for Langrange polynomials), PRE: Chebyshew points already computed
void Cheby::set_omega()
{
    omega_ = Eigen::VectorXd::Zero(deg_+1);
    for(unsigned j=0; j<=deg_; ++j) {
        double hc = 1;
        for(unsigned k=0; k<j; ++k)
            hc *= tk_[j]-tk_[k];
        for(unsigned k=j+1; k<=deg_; ++k)
            hc *= tk_[j]-tk_[k];
        omega_(j) = 1./hc;
    }
}
