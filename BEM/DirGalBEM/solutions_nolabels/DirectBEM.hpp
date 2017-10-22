//// 
//// Copyright (C) 2017 SAM (D-MATH) @ ETH Zurich
//// Author(s): curzuato < > 
//// Contributors:  dcasati 
//// This file is part of the AdvNumCSE repository.
////
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


namespace DirectFirstKind{

  /* 
   * @brief Build and solve direct first kind BIE arising from Dirichlet Laplace
   *        problem (Interior BVP).
   * \param[in] mesh 
   * \param[in] g Dirichlet data. Should take a 2d-vector and return a double.
   * \returns coefficient vector of \f$\mathcal{S}^{-1}_0(\mathcal{G}\f$ 
   *          corresponding to BEM solution (TNu).
   */
  template <typename FUNC>
  Eigen::VectorXd solveDirichlet(const BoundaryMesh& mesh, const FUNC& g){
    // 1. Assemble bilinear form of V as in (1.3.107)
    Eigen::MatrixXd V;
    computeV(V, mesh, 1e-05);
    // 2. Assemble right hand side using <(1/2Id + K)g, psi> as in (1.3.107)
    // - Compute K
    Eigen::MatrixXd K;
    computeK(K, mesh, 1e-05);
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
    // - Put all pieces together and construct RHS
    Eigen::VectorXd RHS = ((0.5*M+K)*G).eval();
    
    // 3. Solve system
    Eigen::VectorXd sol = V.lu().solve(RHS);


    return sol;
  }

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
  template <typename FUNC>
  Eigen::VectorXd solveDirichlet(const BoundaryMesh& mesh, const FUNC& g){
    // 1. Assemble bilinear form as in (1.3.122)
    // - Compute K
    Eigen::MatrixXd K;
    computeK(K, mesh, 1e-05);
    // - Compute Mass matrix for p.w.c/p.w.l
    Eigen::SparseMatrix<double> M01aux(mesh.numElements(), mesh.numVertices());
    computeM01(M01aux, mesh);
    Eigen::MatrixXd M01(mesh.numElements(), mesh.numVertices());
    M01 = Eigen::MatrixXd(M01aux);
    Eigen::MatrixXd LHS = (0.5*M01 - K).transpose();
    
    // 2. Assemble right hand side using bilinear form of W as in (1.3.122)
    Eigen::MatrixXd W;
    computeW(W, mesh, 1e-05);
    // - Compute coefficient vector for g (in \f$\mathcal{S}^{0}_1(\mathcal{G}\f$)
    //   (we do this by interpolation).
    Eigen::VectorXd G(mesh.numVertices());
    for(int i=0; i<mesh.numVertices(); i++){
      G(i) = g(mesh.getVertex(i));
    }
    // - Put all pieces together and construct RHS
    Eigen::VectorXd RHS = W*G;

    // 3. Solve system
    Eigen::VectorXd sol = LHS.lu().solve(RHS);
    

    return sol;
  }

} // end namespace Direct2ndKind

#endif
