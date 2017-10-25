#ifndef MESH_GEN_HPP
#define MESH_GEN_HPP

#include <cmath>
#include <Eigen/Dense>
// CppHilbert includes
#include "source/BoundaryMesh.hpp"

/*
 * @brief Set up an uniform BoundaryMesh for a Square [0, 0.5]^2
 * \param[in] N Number of panels on BoundaryMesh
 */
/* SAM_LISTING_BEGIN_0 */
BoundaryMesh createMiniSquareMesh(const int& N){
  Eigen::MatrixXd coordinates(N,2);
  Eigen::MatrixXi elements(N,2);
  // TODO: CREATE MESH FOR [0, 0.5]^2
  
  BoundaryMesh squareMesh(coordinates, elements);

  return squareMesh;
}
/* SAM_LISTING_END_0 */


/*
 * @brief Set up an uniform BoundaryMesh using polygonal approximation for a 
 *        given parametrized curve. The parametrization should be given on the 
 *        interval I=[-1,1].
 * \tparam[in] gamma Function describing parametrized curve
 * \param[in] N Number of panels on BoundaryMesh
 */
/* SAM_LISTING_BEGIN_1 */
template<typename PARAM>
BoundaryMesh createMeshwithGamma(const PARAM& gamma, const int& N){
  Eigen::MatrixXd coordinates(N,2);
  Eigen::MatrixXi elements(N,2);
  // TODO: CREATE MESH FOR GAMMA([-1,1])
  
  BoundaryMesh gammaMesh(coordinates, elements);

  return gammaMesh;
}
/* SAM_LISTING_END_1 */

#endif
