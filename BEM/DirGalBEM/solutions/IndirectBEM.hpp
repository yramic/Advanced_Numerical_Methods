#ifndef INDIRECT_BEM_HPP
#define INDIRECT_BEM_HPP

#include <Eigen/Dense>
#include <cmath>
#include <unsupported/Eigen/IterativeSolvers>
// CppHilbert includes
#include "source/BoundaryMesh.hpp"
#include "source/buildM.hpp"
#include "source/evaluateK.hpp"
#include "source/evaluateV.hpp"

namespace IndirectFirstKind {

/*
 * @brief Build and solve indirect first kind BIE arising from Dirichlet
 *        Laplace problem (Interior BVP).
 * \param[in] mesh
 * \param[in] g Dirichlet data
 * \returns coefficient vector of \f$ \mathcal{S}^{-1}_0(\mathcal{G}\f$
 *          corresponding to the jump of Neumann trace of the BEM solution.
 */
/* SAM_LISTING_BEGIN_0 */
template <typename FUNC>
Eigen::VectorXd solveDirichlet(const BoundaryMesh& mesh, const FUNC& g) {
  // 1. Assemble bilinear form of \cop{$\biov$}, see \lectref{par:iddir}
  Eigen::MatrixXd V;
  computeV(V, mesh, 1e-05);
  // 2. Assemble right hand side
  // \cop{$(\Fg,\psi)\mapsto\int\limits_{\Gamma}\Fg(\Bx)\,\psi(\Bx)\,\mathrm{d}S(\Bx)$},
  // see \lectref{eq:iddirVv} Compute the mass matrix, Galerkin matrix for
  // \cop{$\Ltwo[\Gamma]$}-inner product (using different test and trial space
  // for discretization).
  Eigen::SparseMatrix<double> M01(mesh.numElements(), mesh.numElements());
  computeM01(M01, mesh);
  Eigen::MatrixXd M = Eigen::MatrixXd(M01);
  // - Compute coefficient vector for g (in $\mathcal{S}^{0}_1(\mathcal{G})$
  //   (we do this by interpolation).
  Eigen::VectorXd G(mesh.numVertices());
  for (int i = 0; i < mesh.numVertices(); i++) {
    G(i) = g(mesh.getVertex(i));
  }
  // - Construct RHS
  Eigen::VectorXd RHS = (M * G).eval();

  // 3. Solve system
  Eigen::VectorXd sol = V.lu().solve(RHS);

  return sol;
}
/* SAM_LISTING_END_0 */

/*
 * @brief Reconstructs and evaluates solution from density psi using Single
 *        Layer Potential
 * \param[in] X Evaluation point
 * \param[in] psi Coefficient vector corresponding to density
 * \param[in] mesh
 */
/* SAM_LISTING_BEGIN_1 */
double reconstructSolution(const Eigen::Vector2d& X, const Eigen::VectorXd& phi,
                           const BoundaryMesh& mesh) {
  Eigen::VectorXd SLphi_x(1);
  evaluateV(SLphi_x, mesh, phi, X.transpose(), 1e-05);

  return SLphi_x(0);
}
/* SAM_LISTING_END_1 */

}  // namespace IndirectFirstKind

namespace IndirectSecondKind {

/*
 * @brief Build and solve indirect second kind BIE arising from Dirichlet
 *        Laplace problem (Interior BVP).
 * \param[in] mesh
 * \param[in] g Dirichlet data
 * \returns coefficient vector of \f$ \mathcal{S}^{-1}_0(\mathcal{G}\f$
 *          corresponding to jump of the Neumann trace of the BEM solution.
 */
/* SAM_LISTING_BEGIN_2 */
template <typename FUNC>
Eigen::VectorXd solveDirichlet(const BoundaryMesh& mesh, const FUNC& g) {
  // 1. Assemble bilinear form (see page 31 on tablet's notes)
  // - Compute K
  Eigen::MatrixXd K;
  computeK00(K, mesh, 1e-05);
  // - Compute Mass matrix for
  // $\mathcal{S}^{-1}_0(\mathcal{G})/\mathcal{S}^{-1}_0(\mathcal{G})$
  Eigen::SparseMatrix<double> M00aux(mesh.numElements(), mesh.numElements());
  computeM00(M00aux, mesh);
  Eigen::MatrixXd M = Eigen::MatrixXd(M00aux);
  Eigen::MatrixXd LHS = (-0.5 * M + K).eval();

  // 2. Assemble right hand side using <g, psi>
  // - Compute coefficient vector for g (in $\mathcal{S}^{0}_1(\mathcal{G}$)
  Eigen::VectorXd G(mesh.numVertices());
  for (int i = 0; i < mesh.numVertices(); i++) {
    G(i) = g(mesh.getVertex(i));
  }
  // - Compute Mass matrix for
  // $\mathcal{S}^{-1}_0(\mathcal{G})/\mathcal{S}^{0}_1(\mathcal{G})$
  Eigen::SparseMatrix<double> M01(mesh.numElements(), mesh.numVertices());
  computeM01(M01, mesh);
  // - Construct RHS
  Eigen::VectorXd RHS = Eigen::MatrixXd(M01) * G;

  // 3. Solve system
  Eigen::VectorXd sol = LHS.lu().solve(RHS);

  return sol;
}
/* SAM_LISTING_END_2 */

/*
 * @brief Reconstructs and evaluates solution from density v using Double
 *        Layer Potential
 * \param[in] X Evaluation point
 * \param[in] v Coefficient vector corresponding to density
 * \param[in] mesh
 */
/* SAM_LISTING_BEGIN_3 */
double reconstructSolution(const Eigen::Vector2d& X, const Eigen::VectorXd& f,
                           const BoundaryMesh& mesh) {
  Eigen::VectorXd DLf_x(1);
  evaluateK(DLf_x, mesh, f, X.transpose(), 1e-05);

  return DLf_x(0);
}
/* SAM_LISTING_END_3 */

}  // namespace IndirectSecondKind

#endif
