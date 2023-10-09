#include <Eigen/Dense>
#include <Eigen/Sparse>
#include <Eigen/SparseLU>
#include <cassert>
#include <cmath>
#include <ctime>
#include <iomanip>
#include <iostream>
#include <unsupported/Eigen/FFT>

using namespace Eigen;
using namespace std;

#if SOLUTION
/* @brief Generate the nodes
 * for a uniform triangulation on an equilateral triangle
 * \param l refinement level
 * \\return nodes (co-ordindates)
 */
MatrixXd genNodes(const int l) {
  int n = pow(2, l);
  int N = (n + 2) * (n + 1) / 2;
  double h = 1.0 / n;

  MatrixXd nodes(2, N);

  int count = 0;
  for (int k = 0; k <= n; k++) {
    for (int j = 0; j <= (n - k); j++) {
      nodes(0, count) = -0.5 + 0.5 * k * h + j * h;
      nodes(1, count) = 0.5 * k * h * sqrt(3);
      count++;
    }
  }

  return nodes;
}

/* @brief Generate the cell connectivity (anti-clockwise) matrix
 * for a uniform triangulation on an equilateral triangle
 * \param l refinement level
 * \\return connectivity matrix
 */
/* SAM_LISTING_BEGIN_0 */
MatrixXi genConnectivity(const int l) {
  int n = pow(2, l);
  int Ne = n * n;  // Number of cells/elements
  MatrixXi cells(3, Ne);

  int shift = 0;
  int count = 0;
  for (int k = 0; k <= n; k++) {
    int nn_local = (n + 1) - k;
    int j_start = shift + 1;  // numbering starts from 1
    int j_end = shift + n - k;

    for (int j = j_start; j <= j_end; j++) {
      {  // cell type I
        int v1 = j;
        int v2 = j + 1;
        int v3 = j + nn_local;

        cells(0, count) = v1;
        cells(1, count) = v2;
        cells(2, count) = v3;

        count++;
      }

      if (j < j_end) {  // cell type II
        int v1 = j + 1;
        int v2 = j + nn_local + 1;
        int v3 = j + nn_local;

        cells(0, count) = v1;
        cells(1, count) = v2;
        cells(2, count) = v3;

        count++;
      }
    }
    shift += nn_local;
  }

  return cells;
}
#endif
/* SAM_LISTING_END_0 */

/* @brief Galerkin matrix for FE discretisation of $-\nabla$
 * for a uniform triangulation on an equilateral triangle
 * \param l refinement level
 * \\return sparse Galerkin matrix
 */
/* SAM_LISTING_BEGIN_1 */
SparseMatrix<double> genGalerkinMat(const int l) {
#if SOLUTION
  // matrix for an arbitrary triangle $K$ of $\mathcal{M}_{\ell}$
  MatrixXd A_el(3, 3);
  A_el << 2, -1, -1, -1, 2, -1, -1, -1, 2;
  A_el *= 0.5 / std::sqrt(3);

  // system dimensions
  unsigned n = std::pow(2, l);
  unsigned N = 0.5 * (n + 2) * (n + 1);  // Number of nodes
  unsigned Ne = n * n;                   // Number of triangular elements

  // define vector of triplets
  typedef Triplet<double> triplet;
  std::vector<triplet> coo_data;

  // get the connectivity/mapping matrix
  MatrixXi map = genConnectivity(l);

  // assemble local element matrix entries to the global matrix
  for (unsigned elId = 0; elId < Ne; ++elId)
    for (unsigned i = 0; i < 3; ++i)
      for (unsigned j = 0; j < 3; ++j) {
        int rowId = map(i, elId) - 1;
        int colId = map(j, elId) - 1;
        coo_data.push_back(triplet(rowId, colId, A_el(i, j)));
      }

  // initializing the sparse Galerkin Matrix
  SparseMatrix<double> A(N, N);
  A.setFromTriplets(coo_data.begin(), coo_data.end());

  return A;
#else  // TEMPLATE
  // TODO: Implement the Galerkin matrix generator
#endif
}
/* SAM_LISTING_END_1 */

/* @brief Prolongation matrix $\VP_{\ell-1, \ell}$
 * \param l refinement level
 * \\return sparse Prolongation matrix
 */
/* SAM_LISTING_BEGIN_2 */
SparseMatrix<double> genProlongationMat(const int l) {
#if SOLUTION
  unsigned nh = std::pow(2, l);      // Cells on edge of triangle for fine mesh
  unsigned nH = std::pow(2, l - 1);  // Cells on edge for coarse mesh
  unsigned Nh = (nh + 2) * (nh + 1) / 2;  // Nodes in fine mesh
  unsigned NH = (nH + 2) * (nH + 1) / 2;  // Nodes in coarse mesh
  unsigned nodeid = 0;

  // define vector of triplets
  typedef Triplet<double> triplet;
  std::vector<triplet> coo_data;

  // Loop over the nodes of fine mesh
  for (unsigned k = 0; k < nh + 1;
       ++k)  // k represents the cell row from the base of the triangle
  {
    for (unsigned i = 0; i < nh + 1 - k;
         ++i)  // i represents the local cell index in the cell row number k
    {
      unsigned idh =
          (nh + 1) * k - (k - 1) * k / 2 + i;  // node id for fine mesh

      // coinciding points on even cell rows
      if (k % 2 == 0 && i % 2 == 0) {
        unsigned kc = k / 2;  // corresponding cell row value for coarse mesh
        unsigned ic = i / 2;  // corresponding cell index value for coarse mesh
        unsigned idH =
            (nH + 1) * kc - (kc - 1) * kc / 2 + ic;  // node id for coarse mesh
        coo_data.push_back(triplet(idh, idH, 1));
      }

      // for the points which are at the center of the edges
      // parallel to x axis, of elements in coarse mesh
      if (k % 2 == 0 && i % 2 == 1) {
        unsigned kc = k / 2;  // cell row value for coarse mesh
        unsigned icl =
            (i - 1) / 2;  // cell index value of left point for coarse mesh
        unsigned icr =
            (i + 1) / 2;  // cell index value for right point coarse mesh

        unsigned idH = (nH + 1) * kc - (kc - 1) * kc / 2 + icl;
        coo_data.push_back(triplet(idh, idH, 0.5));

        idH = (nH + 1) * kc - (kc - 1) * kc / 2 + icr;
        coo_data.push_back(triplet(idh, idH, 0.5));
      }

      // for the points which are at the center of the left
      // edges of elements in coarse mesh
      if (k % 2 == 1 && i % 2 == 0) {
        unsigned kct =
            (k + 1) / 2;  // cell row value for the upper row for coarse mesh
        unsigned kcb =
            (k - 1) / 2;  // cell row value for the lower row for coarse mesh
        unsigned ic = i / 2;  // cell index value for coarse mesh

        unsigned idH = (nH + 1) * kct - (kct - 1) * kct / 2 + ic;
        coo_data.push_back(triplet(idh, idH, 0.5));

        idH = (nH + 1) * kcb - (kcb - 1) * kcb / 2 + ic;
        coo_data.push_back(triplet(idh, idH, 0.5));
      }

      // for the points which are at the center of the right
      // edges of elements in coarse mesh
      if (k % 2 == 1 && i % 2 == 1) {
        unsigned kct =
            (k + 1) / 2;  // cell row value for the upper row for coarse mesh
        unsigned kcb =
            (k - 1) / 2;  // cell row value for the lower row for coarse mesh
        unsigned ic =
            (i - 1) / 2;  // cell index value for upper row of coarse mesh

        unsigned idH = (nH + 1) * kct - (kct - 1) * kct / 2 + ic;
        coo_data.push_back(triplet(idh, idH, 0.5));

        ic = (i + 1) / 2;  // cell index value for lower row of coarse mesh
        idH = (nH + 1) * kcb - (kcb - 1) * kcb / 2 + ic;
        coo_data.push_back(triplet(idh, idH, 0.5));
      }
    }
  }

  // initializing the prolongation matrix with appropriate size
  SparseMatrix<double> P(Nh, NH);
  P.setFromTriplets(coo_data.begin(), coo_data.end());

  return P;
#else  // TEMPLATE
  // TODO: Implement the prolongation matrix generator
#endif
}
/* SAM_LISTING_END_2 */

#if SOLUTION
/* @brief The Gauss-Seidel iterative method
 * to solve sparse linear system of equations
 * \param A Sparse matrix, $\IR^{N \times N}$
 * \param phi right-hand-side vector, $\IR^{N}$
 * \param mu initial guess and output solution, $\IR^{N \times N}$
 * \param maxItr maximum number of iterations for termination criteria
 * \param TOL iteration error tolerance for termination criteria
 */
void gaussSeidel(const SparseMatrix<double> &A, const VectorXd &phi,
                 VectorXd &mu, int maxItr, double TOL = 1.0E-06) {
  VectorXd delta(A.rows());
  for (int i = 0; i < maxItr; i++) {
    delta = A.triangularView<Lower>().solve(phi - A * mu);
    mu += delta;
    // cout << i << "\t" << delta.norm() << endl;
    if (delta.norm() <= TOL * mu.norm()) break;
  }
}

/* @brief Augment the singular system matrix
 * \param A Sparse matrix, $\IR^{N \times N}$
 * \\return Dense augmented matrix
 */
/* SAM_LISTING_BEGIN_3 */
MatrixXd augmentA(const SparseMatrix<double> &A) {
  unsigned N = A.rows();
  MatrixXd fullA(A);
  MatrixXd augA = MatrixXd::Zero(N + 1, N + 1);

  augA.block(0, 0, N, N) = fullA;
  augA.block(N, 0, 1, N) = VectorXd::Constant(N, 1).transpose();
  augA.block(0, N, N, 1) = VectorXd::Constant(N, 1);

  return augA;
}
/* SAM_LISTING_END_3 */
#endif

/* @brief Multi-grid method
 * to solve sparse linear system of equations
 * \param phi right-hand-side vector, $\IR^{N}$
 * \param mu initial guess and output solution, $\IR^{N \times N}$
 * \param max_n_steps maximum number of steps
 * \param TOL error tolerance for termination criteria
 */
/* SAM_LISTING_BEGIN_4 */
void multiGridIter(const VectorXd &phi, VectorXd &mu, int l, int max_n_steps,
                   double TOL = 1E-04) {
#if SOLUTION
  // get the Galerkin matrix
  SparseMatrix<double> A = genGalerkinMat(l);

  // direct solve using augmentation
  if (l == 0) {
    MatrixXd Aaug = augmentA(A);  // augmented A matrix
    VectorXd rhsaug(4);
    rhsaug << phi, 0;  // augmented rhs vector

    // solve the augmented system
    VectorXd muaug = Aaug.lu().solve(rhsaug);
    assert((Aaug * muaug - rhsaug).norm() < TOL);
    mu = muaug.segment(0, 3);
  } else {
    // get the prolongation matrix P
    SparseMatrix<double> P = genProlongationMat(l);
    for (unsigned nsteps = 0; nsteps < max_n_steps; ++nsteps) {
      VectorXd mu_old = mu;
      gaussSeidel(A, phi, mu, 1, TOL);  // One gauss seidel step

      VectorXd rhoh = phi - A * mu;
      VectorXd rhoH = P.transpose() * rhoh;
      VectorXd vH = VectorXd::Zero(rhoH.rows());
      multiGridIter(rhoH, vH, l - 1, TOL, max_n_steps);
      mu += P * vH;

      if ((mu - mu_old).norm() < TOL * mu.norm()) break;
    }
  }
#else  // TEMPLATE
  // TODO: Implement multi-grid iteration
#endif
}
/* SAM_LISTING_END_4 */

/* @brief Computes the residual for the Multi-grid method with max_n_steps=1
 * \param L mesh level
 * \param max_itr maximum number of multi-grid iterations
 * \param TOL error tolerance for termination criteria
 */
/* SAM_LISTING_BEGIN_5 */
VectorXd test_multigrid_residual(const int L, int max_itr,
                                 double TOL = 1.0E-08) {
  VectorXd rho_norm(max_itr);

  unsigned n = std::pow(2, L);
  unsigned N = 0.5 * (n + 2) * (n + 1);  // Number of nodes

  SparseMatrix<double> A = genGalerkinMat(L);
  VectorXd phi = VectorXd::Zero(N);
  VectorXd mu = VectorXd::Random(N);

  unsigned int count = 0;
  for (int i = 0; i < max_itr; i++) {
    multiGridIter(phi, mu, L, 1);
    rho_norm(i) = (phi - A * mu).norm();

    count++;
    if (rho_norm(i) <= TOL) break;
  }

  return rho_norm.head(count);
}
/* SAM_LISTING_END_5 */

/* @brief Computes the asymptotic convergence rates of the Multi-grid method
 * by measuring the quotient of the residuals from successive iterations.
 */
/* SAM_LISTING_BEGIN_6 */
void test_multigrid_convergence() {
  int max_itr = 100;
  cout << "L\trate" << endl;
  for (int L = 2; L <= 8; L++) {
    VectorXd rho_norm = test_multigrid_residual(L, max_itr);
    VectorXd rates(rho_norm.size() - 1);
    for (int i = 0; i < rates.size(); i++)
      rates(i) = rho_norm(i + 1) / rho_norm(i);
    cout << L << "\t" << rates.tail(1) << endl;
  }
}
/* SAM_LISTING_END_6 */

int main() {
  int l = 2;
  double TOL = 1E-04;

  SparseMatrix<double> A = genGalerkinMat(l);
  cout << "\n\nGalerkin matrix for l=" << l << ":\n\n" << A << endl;

  SparseMatrix<double> P = genProlongationMat(l);
  cout << "\n\nProlongation matrix for l=" << l << ":\n\n" << P << endl;

  unsigned n = std::pow(2, l);
  unsigned N = 0.5 * (n + 2) * (n + 1);  // Number of nodes
  VectorXd phi = VectorXd::Zero(N);
  VectorXd mu = VectorXd::Random(N);

  int max_nsteps = 1000;
  multiGridIter(phi, mu, l, max_nsteps, TOL);
  // cout << "\n" << mu.transpose() << endl;

#if SOLUTION
  cout << "Multi-grid convergence test:" << endl;
  test_multigrid_convergence();
#endif
}

// End of file
