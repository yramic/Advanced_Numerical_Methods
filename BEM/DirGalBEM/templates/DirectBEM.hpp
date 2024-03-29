#ifndef DIRECT_BEM_HPP
#define DIRECT_BEM_HPP

#include <cmath>
#include <iostream>
#include <iomanip>
#include <Eigen/Dense>
#include <unsupported/Eigen/IterativeSolvers>
// CppHilbert includes
#include "source/BoundaryMesh.hpp"
#include "source/buildV.hpp"
#include "source/buildW.hpp"
#include "source/buildK.hpp"
#include "source/buildM.hpp"


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
    // TODO: ASSEMBLE AND SOLVE DIRECT FIRST-KIND BIE
    Eigen::VectorXd sol(mesh.numElements());

    return sol;
  }
  /* SAM_LISTING_END_0 */

} // end namespace Direct1stKind



namespace DirectSecondKind{

   /* 
   * @brief Build and solve direct second kind BIE arising from Dirichlet Laplace 
   *        problem (Interior BVP) with different test and trial space (non-stable)
   * \param[in] mesh 
   * \param[in] g Dirichlet data. Should take a 2d-vector and return a double.
   * \returns coefficient vector of \f$ \mathcal{S}^{-1}_0(\mathcal{G}\f$ 
   *          corresponding to BEM solution (TNu)
   */
  /* SAM_LISTING_BEGIN_1 */
  template <typename FUNC>
  Eigen::VectorXd solveDirichlet(const BoundaryMesh& mesh, const FUNC& g){
    // TODO: ASSEMBLE AND SOLVE DIRECT SECOND-KIND BIE
    Eigen::VectorXd sol(mesh.numElements());

    return sol;
  }
  /* SAM_LISTING_END_1 */

} // end namespace Direct2ndKind

#endif
