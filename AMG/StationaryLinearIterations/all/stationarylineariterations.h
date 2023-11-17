/**
 * @file stationarylineariterations.h
 * @brief ADVNCSE homework StationaryLinearIterations code
 * @author D. Casati, Bob Schreiner
 * @date October 2018
 * @copyright Developed at SAM, ETH Zurich
 */

#ifndef STATIONARYLINEARITERATIONS_H
#define STATIONARYLINEARITERATIONS_H

#include <Eigen/Dense>
#include <Eigen/Sparse>
#include <iomanip>
#include <iostream>

namespace StationaryLinearIterations {
/** @brief Initialization of Poisson matrix as sparse matrix */
Eigen::SparseMatrix<double> poissonMatrix(unsigned int n);

/** @brief Gauss-Seidel iterative solver for a sparse linear system of equations
 */
void gaussSeidel(const Eigen::SparseMatrix<double> &A,
                 const Eigen::VectorXd &phi, Eigen::VectorXd &mu,
                 int maxItr = 1000, double TOL = 1.0E-06);

/* @brief Determines the asymptotic convergence rate of the Gauss-Seidel
 * iterative method using power iteration,
 * for the matrix \$\mathbf{X} = \mathbf{A} + c \mathbf{I}_{N}\$.
 * Here \$\mathbf{A}\$ is the Poisson matrix on a unit square domain.
 * @param n matrix size
 * @param c linear combination coefficient
 * @return asymptotic convergence rate
 */
double gaussSeidelRate(unsigned int n, double c, double TOL = 1.0E-03);

}  // namespace StationaryLinearIterations

#endif  //STATIONARYLINEARITERATIONS_H