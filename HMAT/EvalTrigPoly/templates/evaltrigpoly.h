/**
 * @file evaltrigpoly.h
 * @brief NPDE homework EvalTrigPoly code
 * @author Ralf Hiptmair
 * @date August 2023
 * @copyright Developed at SAM, ETH Zurich
 */

#include <kernmatllrapprox.h>

#include <Eigen/Dense>
#include <iomanip>
#include <iostream>
#include <vector>

namespace EvalTrigPoly {

// Evaluation of trigonometric polynomial in equidistant points
Eigen::VectorXcd evalTrigPolyEquid(const Eigen::VectorXcd &coeffs);

// Direct evaluation of a trigonometric polynomial
Eigen::VectorXcd evalTrigPoly(const Eigen::VectorXcd &gamma,
                              const Eigen::VectorXd &x);

// Detecting problematic evaluation points
std::vector<unsigned int> checkX(const Eigen::VectorXcd &gamma,
                                 const Eigen::VectorXd &x, double tau);

// Approximate fast evaluation of a trigonometric polynomial
Eigen::VectorXcd evalTrigPolyApprox(const Eigen::VectorXcd &gamma,
                                    const Eigen::VectorXd &x, unsigned int q);

}  // namespace EvalTrigPoly
