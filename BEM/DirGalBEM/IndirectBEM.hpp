#ifndef INDIRECT_BEM_HPP
#define INDIRECT_BEM_HPP

#include <cmath>
#include <Eigen/Dense>
#include <unsupported/Eigen/IterativeSolvers>
// CppHilbert includes
#include "../CppHilbert/Library/source/BoundaryMesh.hpp"
#include "../CppHilbert/Library/source/evaluateV.hpp"
#include "../CppHilbert/Library/source/evaluateK.hpp"
#include "../CppHilbert/Library/source/buildM.hpp"

namespace IndirectFirstKind{

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
  Eigen::VectorXd solveDirichlet(const BoundaryMesh& mesh, const FUNC& g){
    #if SOLUTION
    // 1. Assemble bilinear form of V (see par 1.3.132)
    Eigen::MatrixXd V;
    computeV(V, mesh, 1e-05);
    // 2. Assemble right hand side using <g, psi> (see 1.3.134)
    // - Compute Mass Matrix
    Eigen::SparseMatrix<double> M01(mesh.numElements(), mesh.numVertices());
    computeM01(M01, mesh);
    Eigen::MatrixXd M; M = Eigen::MatrixXd(M01);
    // - Compute coefficient vector for g (in \f$\mathcal{S}^{0}_1(\mathcal{G}\f$)
    //   (we do this by interpolation).
    Eigen::VectorXd G(mesh.numVertices());
    for(int i=0; i<mesh.numVertices(); i++){
      G(i) = g(mesh.getVertex(i));
    }
    // - Construct RHS
    Eigen::VectorXd RHS = (M*G).eval();
    
    // 3. Solve system
    Eigen::VectorXd sol = V.lu().solve(RHS);

    #else // TEMPLATE
    // TODO: ASSEMBLE AND SOLVE INDIRECT FIRST-KIND BIE
    Eigen::VectorXd sol(mesh.numElements());
    #endif // TEMPLATE

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
			     const BoundaryMesh& mesh){
    Eigen::VectorXd SLphi_x(1);
    #if SOLUTION
    evaluateV(SLphi_x, mesh, phi, X.transpose(), 1e-05);
    
    #else // TEMPLATE
    // TODO: USE SINGLE LAYER POTENTIAL TO EVALUATE U(X)
    #endif // TEMPLATE
    return SLphi_x(0);
  }
  /* SAM_LISTING_END_1 */

}  // end namespace Indirect1stKind



namespace IndirectSecondKind{

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
  Eigen::VectorXd solveDirichlet(const BoundaryMesh& mesh, const FUNC& g){
    #if SOLUTION
    // 1. Assemble bilinear form (see page 31 on tablet's notes)
    // - Compute K
    Eigen::MatrixXd K;
    computeK(K, mesh, 1e-05);
    // - Compute Mass matrix for p.w.c/p.w.l
    Eigen::SparseMatrix<double> M01aux(mesh.numElements(), mesh.numVertices());
    computeM01(M01aux, mesh);
    Eigen::MatrixXd M; M = Eigen::MatrixXd(M01aux);
    Eigen::MatrixXd LHS = (-0.5*M + K).eval();
    
    // 2. Assemble right hand side using <g, psi> 
    // - Compute coefficient vector for g (in \f$\mathcal{S}^{0}_1(\mathcal{G}\f$)
    //   (we do this by interpolation).
    Eigen::VectorXd G(mesh.numVertices());
    for(int i=0; i<mesh.numVertices(); i++){
      G(i) = g(mesh.getVertex(i));
    }
    // - Construct RHS
    Eigen::VectorXd RHS = M*G;

    // 3. Solve system    
    Eigen::VectorXd sol = LHS.lu().solve(RHS);

    #else // TEMPLATE
    // TODO: ASSEMBLE AND SOLVE INDIRECT SECOND-KIND BIE
    Eigen::VectorXd sol(mesh.numVertices());
    #endif // TEMPLATE

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
			     const BoundaryMesh& mesh){
    Eigen::VectorXd DLf_x(1);
    #if SOLUTION
    evaluateK(DLf_x, mesh, f, X.transpose(), 1e-05);
    
    #else // TEMPLATE
    // TODO: USE DOUBLE LAYER POTENTIAL TO EVALUATE U(X)
    #endif // TEMPLATE
    return DLf_x(0);
  }
  /* SAM_LISTING_END_3 */

} // end namespace Indirect2ndkind


#endif
