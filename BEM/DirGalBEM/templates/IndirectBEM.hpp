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
    // TODO: ASSEMBLE AND SOLVE INDIRECT FIRST-KIND BIE
    Eigen::VectorXd sol(mesh.numElements());

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
    // TODO: USE SINGLE LAYER POTENTIAL TO EVALUATE U(X)
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
    // TODO: ASSEMBLE AND SOLVE INDIRECT SECOND-KIND BIE
    Eigen::VectorXd sol(mesh.numVertices());

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
    // TODO: USE DOUBLE LAYER POTENTIAL TO EVALUATE U(X)
    return DLf_x(0);
  }
  /* SAM_LISTING_END_3 */

} // end namespace Indirect2ndkind


#endif
