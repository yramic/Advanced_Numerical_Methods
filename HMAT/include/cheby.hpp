#ifndef CHEBY_HPP
#define CHEBY_HPP

#include <Eigen/Dense>


/**
* \brief Compute Chebyshew nodes "tk_" and weights "wk_"
* on domain [xl,xr] for a certain polynomial degree
*/
class Cheby
{
public:

    /**
    * \brief Constructor
    */
    Cheby(double xl, double xr, unsigned deg);

    /**
    * \brief Getters
    */
    // return Chebyshew nodes on domain [xl,xr]
    Eigen::VectorXd getNodes() const {
        return tk_;
    }
    // return weights of Lagrange polynomial
    Eigen::VectorXd getWghts() const {
        return wk_;
    }

    /**
    * \brief Setters
    */
    // compute Chebyshew nodes on domain [xl,xr]
    void setNodes();
    // compute weights of Lagrange polynomial
    void setWghts();

private:

    double xl_; // left  boundary of domain on which we compute the Chebyshew nodes
    double xr_; // right boundary of domain on which we compute the Chebyshew nodes
    unsigned deg_; // degree of Lagrange polynomial
    Eigen::VectorXd tk_; // Chebyshew nodes
    Eigen::VectorXd wk_; // weights of Lagrange polynomial
};

#endif // CHEBY_HPP