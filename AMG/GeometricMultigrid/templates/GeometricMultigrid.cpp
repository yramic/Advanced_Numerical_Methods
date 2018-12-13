#include <Eigen/Dense>
#include <Eigen/Sparse>
#include <Eigen/SparseLU>
#include <cmath>
#include <ctime>
#include <iomanip>
#include <iostream>
#include <unsupported/Eigen/FFT>
#include <cassert>


using namespace Eigen;
using namespace std;


/* SAM_LISTING_END_0 */


/* @brief Galerkin matrix for FE discretisation of $-\nabla$
 * for a uniform triangulation on an equilateral triangle
 * \param l refinement level
 * \\return sparse Galerkin matrix
 */
/* SAM_LISTING_BEGIN_1 */
SparseMatrix<double> genGalerkinMat(const int l) {
    // TODO: Implement the Galerkin matrix generator
}
/* SAM_LISTING_END_1 */


/* @brief Prolongation matrix $\VP_{\ell-1, \ell}$
 * \param l refinement level
 * \\return sparse Prolongation matrix
 */
/* SAM_LISTING_BEGIN_2 */
SparseMatrix<double> genProlongationMat(const int l) 
{
    // TODO: Implement the prolongation matrix generator
}
/* SAM_LISTING_END_2 */




/* @brief Multi-grid method
 * to solve sparse linear system of equations
 * \param phi right-hand-side vector, $\IR^{N}$
 * \param mu initial guess and output solution, $\IR^{N \times N}$
 * \param max_n_steps maximum number of steps
 * \param TOL error tolerance for termination criteria
 */
/* SAM_LISTING_BEGIN_4 */
void multiGridIter(const VectorXd &phi, VectorXd &mu, int l, int max_n_steps, double TOL=1E-04) {
    // TODO: Implement multi-grid iteration
}
/* SAM_LISTING_END_4 */


int main() {
    
    int l = 2;
    double TOL = 1E-04;
    
    SparseMatrix<double> A = genGalerkinMat(l);
    cout << "\n\nGalerkin matrix for l="<< l << ":\n\n" << A << endl;
    
    SparseMatrix<double> P = genProlongationMat(l);
    cout << "\n\nProlongation matrix for l="<< l << ":\n\n" << P << endl;
    
    unsigned n = std::pow(2, l);
    unsigned N = 0.5 * (n + 2) * (n + 1); // Number of nodes
    VectorXd phi = VectorXd::Zero(N);
    VectorXd mu = VectorXd::Random(N);
    
    int max_nsteps = 1000;
    multiGridIter(phi, mu, l, max_nsteps, TOL);
    //cout << "\n" << mu.transpose() << endl;
    
}


// End of file
