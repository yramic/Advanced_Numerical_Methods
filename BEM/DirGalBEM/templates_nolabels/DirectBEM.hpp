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
    // TODO: ASSEMBLE AND SOLVE DIRECT FIRST-KIND BIE
    Eigen::VectorXd sol(mesh.numElements());

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
    // TODO: ASSEMBLE AND SOLVE DIRECT SECOND-KIND BIE
    Eigen::VectorXd sol(mesh.numElements());

    return sol;
  }

} // end namespace Direct2ndKind

#endif