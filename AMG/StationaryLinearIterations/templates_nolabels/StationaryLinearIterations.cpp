//// 
//// Copyright (C) 2017 SAM (D-MATH) @ ETH Zurich
//// Author(s): curzuato < > 
//// Contributors:  dcasati 
//// This file is part of the AdvNumCSE repository.
////
#include <Eigen/Dense>
#include <Eigen/Sparse>
#include <unsupported/Eigen/FFT>
#include <cmath>
#include <iomanip>
#include <iostream>
#include <ctime>

using namespace Eigen;
using namespace std;


/* @brief Poisson matrix for a uniform triangulation on a unit square domain
 * \param N matrix size
 * \\return sparse Poisson matrix
 */
SparseMatrix<double> poissonMatrix(const int n)
{
    // TODO: Set-up the Poisson matrix
}


/* @brief The Gauss-Seidel iterative method
 * to solve sparse linear system of equations
 * \param A Sparse matrix, $\IR^{N \times N}$
 * \param phi right-hand-side vector, $\IR^{N}$
 * \param mu initial guess and output solution, $\IR^{N \times N}$
 * \param TOL iteration error tolerance for termination criteria
 */
void gaussSeidel(const SparseMatrix<double> &A, const VectorXd &phi, VectorXd &mu, double TOL = 1.0E-06)
{
    // TODO: Implement the Gauss-Seidel iterative method to solve a sparse linear system of equations
}


double comp_lmax_gaussSeidel(const SparseMatrix<double> &X, double TOL=1.0E-03)
{
    int N = X.rows();
    
    VectorXd v = MatrixXd::Random(N,1);
    double lambda_new = 0;
    double lambda_old = 1;
    
    int itr=0;
    do
    {
        lambda_old = lambda_new;
        v /= v.norm();
        v -= X.triangularView<Lower>().solve(X*v); // E = I - M*A; M = inv(tril(A))
        lambda_new = v.norm();
    }
    while (fabs(lambda_new-lambda_old) > TOL*lambda_new);
    
    return lambda_new;
}


/* @brief Determines the asymptotic convergence rate of the Gauss-Seidel
 * iterative method using power iteration, 
 * for the matrix $\mathbf{X} = \mathbf{A} + c \mathbf{I}_{N}$.
 * Here $\mathbf{A}$ is the Poisson matrix on a unit square domain.
 * \param N matrix size
 * \param c linear combination coefficient
 * \\return asymptotic convergence rate
 */
double gaussSeidelRate(const int n, double c, double TOL=1.0E-03)
{
    // TODO: Implement power iteration to compute the asymptotic convergence rate
    // of the Gauss-Seidel iterative method for the matrix $\mathbf{X}$.
}




int main() {

    // TODO: Run tests

}
