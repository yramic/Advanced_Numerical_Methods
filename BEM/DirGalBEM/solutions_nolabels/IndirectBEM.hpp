//// 
//// Copyright (C) 2017 SAM (D-MATH) @ ETH Zurich
//// Author(s): curzuato < > 
//// Contributors:  dcasati 
//// This file is part of the AdvNumCSE repository.
////
#ifndef INDIRECT_BEM_HPP
#define INDIRECT_BEM_HPP

#include <cmath>
#include <Eigen/Dense>
#include <unsupported/Eigen/IterativeSolvers>
// CppHilbert includes
#include "source/BoundaryMesh.hpp"
#include "source/evaluateV.hpp"
#include "source/evaluateK.hpp"
#include "source/buildM.hpp"

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
    Eigen::MatrixXd V; computeV(V, mesh, 1e-05);
    // 2. Assemble right hand side using <g, psi> (see 1.3.134)
    // - Compute Mass Matrix
    Eigen::SparseMatrix<double> M00(mesh.numElements(), mesh.numElements());
    computeM00(M00, mesh);
    Eigen::MatrixXd M = Eigen::MatrixXd(M00);
    // - Compute coefficient vector for g (in $\mathcal{S}^{-1}_0(\mathcal{G}$)
    Eigen::VectorXd G(mesh.numElements());
    for(int i=0; i<mesh.numElements(); i++){
      Eigen::Vector2d a,b; std::tie(a,b) = mesh.getElementVertices(i);
      G(i) = g(0.5*(a+b).eval());
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
  double reconstructSolution(const Eigen::Vector2d& X, const Eigen::VectorXd& phi,
			     const BoundaryMesh& mesh){
    Eigen::VectorXd SLphi_x(1);
    evaluateV(SLphi_x, mesh, phi, X.transpose(), 1e-05);
    
    return SLphi_x(0);
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
    Eigen::MatrixXd K; computeK00(K, mesh, 1e-05);
    // - Compute Mass matrix for p.w.c/p.w.c
    Eigen::SparseMatrix<double> M00aux(mesh.numElements(), mesh.numElements());
    computeM00(M00aux, mesh);
    Eigen::MatrixXd M = Eigen::MatrixXd(M00aux);
    Eigen::MatrixXd LHS = (-0.5*M + K).eval();
    
    // 2. Assemble right hand side using <g, psi> 
    // - Compute coefficient vector for g (in $\mathcal{S}^{-1}_0(\mathcal{G}$)
    Eigen::VectorXd G(mesh.numElements());
    for(int i=0; i<mesh.numElements(); i++){
      Eigen::Vector2d a,b;
      std::tie(a,b) = mesh.getElementVertices(i);
      G(i) = g(0.5*(a+b).eval());
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
  double reconstructSolution(const Eigen::Vector2d& X, const Eigen::VectorXd& f,
			     const BoundaryMesh& mesh){
    Eigen::VectorXd DLf_x(1);
    evaluateK(DLf_x, mesh, f, X.transpose(), 1e-05);
    
    return DLf_x(0);
  }

} // end namespace Indirect2ndkind


#endif
