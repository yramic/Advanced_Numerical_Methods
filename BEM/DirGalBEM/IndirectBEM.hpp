#ifndef INDIRECT_BEM_HPP
#define INDIRECT_BEM_HPP

#include <cmath>
#include <Eigen/Dense>
#include <unsupported/Eigen/IterativeSolvers>

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
  template <typename FUNC>
  Eigen::VectorXd solveDirichlet(const BoundaryMesh& mesh, const FUNC& g){
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

    return sol;
  }


  /* 
   * @brief Reconstructs and evaluates solution from density psi using Single 
   *        Layer Potential
   * \param[in] X Evaluation point
   * \param[in] psi Coefficient vector corresponding to density
   * \param[in] mesh 
   */
  double reconstructSolution(const Eigen::Vector2d& X, const Eigen::VectorXd& psi,
			     const BoundaryMesh& mesh){
    Eigen::VectorXd SLpsi_x;
    evaluateV(SLpsi_x, mesh, psi, X.transpose(), 1e-05);
    return SLpsi_x(0);
  }

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
  template <typename FUNC>
  Eigen::VectorXd solveDirichlet(const BoundaryMesh& mesh, const FUNC& g){
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

    return sol;
  }

  
  /* 
   * @brief Reconstructs and evaluates solution from density v using Double 
   *        Layer Potential
   * \param[in] X Evaluation point
   * \param[in] v Coefficient vector corresponding to density
   * \param[in] mesh 
   */
  double reconstructSolution(const Eigen::Vector2d& X, const Eigen::VectorXd& v,
			     const BoundaryMesh& mesh){
    Eigen::VectorXd DLv_x;
    evaluateK(DLv_x, mesh, v, X.transpose(), 1e-05);
    return DLv_x(0);
  }

} // end namespace Indirect2ndkind


#endif
