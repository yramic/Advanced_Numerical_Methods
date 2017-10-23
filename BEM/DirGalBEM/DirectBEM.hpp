#ifndef DIRECT_BEM_HPP
#define DIRECT_BEM_HPP

#include <cmath>
#include <iostream>
#include <iomanip>
#include <Eigen/Dense>
#include <unsupported/Eigen/IterativeSolvers>
// CppHilbert includes
#include "../CppHilbert/Library/source/BoundaryMesh.hpp"
#include "../CppHilbert/Library/source/buildV.hpp"
#include "../CppHilbert/Library/source/buildW.hpp"
#include "../CppHilbert/Library/source/buildK.hpp"
#include "../CppHilbert/Library/source/buildM.hpp"
#include "../CppHilbert/Library/source/evaluateW.hpp"


namespace DirectFirstKind{

  /* 
   * @brief Build and solve direct first kind BIE arising from Dirichlet Laplace
   *        problem (Interior BVP).
   * \param[in] mesh 
   * \param[in] g Dirichlet data. Should take a 2d-vector and return a double.
   * \returns coefficient vector of \f$\mathcal{S}^{-1}_0(\mathcal{G}\f$ 
   *          corresponding to BEM solution (TNu).
   */
  /* SAM_LISTING_BEGIN_0 */
  template <typename FUNC>
  Eigen::VectorXd solveDirichlet(const BoundaryMesh& mesh, const FUNC& g){
    #if SOLUTION
    // 1. Assemble bilinear form of V as in (1.3.107)
    Eigen::MatrixXd V; computeV(V, mesh, 1e-05);
    // 2. Assemble right hand side using <(1/2Id + K)g, psi> as in (1.3.107)
    // - Compute K
    Eigen::MatrixXd K; computeK00(K, mesh, 1e-05);
    // - Compute Mass Matrix
    Eigen::SparseMatrix<double> M0(mesh.numElements(), mesh.numElements());
    computeM00(M0, mesh);
    Eigen::MatrixXd M = Eigen::MatrixXd(M0);
    // - Compute coefficient vector for g (in $\mathcal{S}^{0}_1(\mathcal{G}$)
    //   (we do this by interpolation).
    Eigen::VectorXd G(mesh.numElements());
    for(int i=0; i<mesh.numElements(); i++){
      Eigen::Vector2d a,b;
      std::tie(a,b) = mesh.getElementVertices(i);
      G(i) = g(0.5*(a+b).eval());
    }
    // - Put all pieces together and construct RHS
    Eigen::VectorXd RHS = ((0.5*M+K)*G).eval();
    
    // 3. Solve system
    Eigen::VectorXd sol = V.lu().solve(RHS);

    #else // TEMPLATE
    // TODO: ASSEMBLE AND SOLVE DIRECT FIRST-KIND BIE
    Eigen::VectorXd sol(mesh.numElements());
    #endif // TEMPLATE

    return sol;
  }
  /* SAM_LISTING_END_0 */

} // end namespace Direct1stKind



namespace DirectSecondKind{

   /* 
   * @brief Build and solve direct second kind BIE arising from Dirichlet Laplace 
   *        problem (Interior BVP). 
   * \param[in] mesh 
   * \param[in] g Dirichlet data. Should take a 2d-vector and return a double.
   * \returns coefficient vector of \f$ \mathcal{S}^{-1}_0(\mathcal{G}\f$ 
   *          corresponding to BEM solution (TNu)
   */
  /* SAM_LISTING_BEGIN_1 */
  template <typename FUNC>
  Eigen::VectorXd solveDirichlet(const BoundaryMesh& mesh, const FUNC& g){
    #if SOLUTION
    // 1. Assemble bilinear form as in (1.3.122)
    // - Compute K
    Eigen::MatrixXd K; computeK00(K, mesh, 1e-05);
    // - Compute Mass matrix for p.w.c/p.w.l
    Eigen::SparseMatrix<double> M00aux(mesh.numElements(), mesh.numElements());
    computeM00(M00aux, mesh);
    Eigen::MatrixXd M0 = Eigen::MatrixXd(M00aux);
    Eigen::MatrixXd LHS = (0.5*M0 - K).transpose();
    
    // 2. Assemble right hand side using bilinear form of W as in (1.3.122)
    Eigen::MatrixXd W; computeW(W, mesh, 1e-05);
    // - Compute coefficient vector for g (in $\mathcal{S}^{0}_1(\mathcal{G}$)
    //   (we do this by interpolation).
    Eigen::VectorXd G(mesh.numVertices());
    for(int i=0; i<mesh.numVertices(); i++){
      G(i) = g(mesh.getVertex(i));
    }
    // - Put all pieces together and construct RHS
    Eigen::VectorXd RHS = W*G; // funnyly this works although it's wrong
    // this one doesn't
    /*
    Eigen::VectorXd G(mesh.numElements());
    Eigen::MatrixXd X(mesh.numElements(),2), n(mesh.numElements(),2);
    for(int i=0; i<mesh.numElements(); i++){
      Eigen::Vector2d a,b;
      std::tie(a,b) = mesh.getElementVertices(i);
      n.row(i) = unitNormal(a,b);
      X.row(i) = (a+b)/2.;      
      G(i) = g(mesh.getVertex(i));
    }
    Eigen::VectorXd WG_x;
    evaluateW(WG_x, mesh, G, X, n, 1e-05);
    Eigen::VectorXd RHS = M0*WG_x;
    */
    
    // 3. Solve system
    Eigen::VectorXd sol = LHS.lu().solve(RHS);
    
    #else // TEMPLATE
    // TODO: ASSEMBLE AND SOLVE DIRECT SECOND-KIND BIE
    Eigen::VectorXd sol(mesh.numElements());
    #endif // TEMPLATE

    return sol;
  }
  /* SAM_LISTING_END_1 */

} // end namespace Direct2ndKind

#endif
