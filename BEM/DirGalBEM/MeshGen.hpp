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
  #if SOLUTION
  // N should be divisible by 4
  assert(N%4 ==0);
  // N should be larger or equal than 4
  assert(N>=4);

  // Create variable indicating number of points on each side of square
  int NperSide = N/4;
  
  // 1. Create boundary vertices for unit square [0,0.5]x[0,0.5].
  // - Begin by creating auxiliary vector discretizing interval [0,0.5]  
  Eigen::VectorXd auxLinSp; auxLinSp.setLinSpaced(NperSide+1, 0, 0.5);
  // - Create auxiliary vectors with ones and zeros
  Eigen::VectorXd auxZero(NperSide), auxOne(NperSide);
  auxZero.setZero(); auxOne.setOnes();

  // - Create Matrix containing the coordinates of the mesh vertices
  Eigen::MatrixXd coordinates(N,2);
  // Use auxiliary vectors to fill coordinate matrix. Remember we need to
  // generate the points counter-clockwise.
  coordinates << auxLinSp.segment(0,NperSide),             auxZero, // bottom line
              auxOne*0.5,             auxLinSp.segment(0,NperSide), // right line
              (auxLinSp.reverse()).segment(0,NperSide), 0.5*auxOne, // upper line
              auxZero,    (auxLinSp.reverse()).segment(0,NperSide); // left line

  // 2. Create matrix specifying the indices of the vertices of each element
  Eigen::MatrixXi elements(N,2);
  elements <<  Eigen::VectorXi::LinSpaced(N, 0, N-1),
               Eigen::VectorXi::LinSpaced(N, 1, N);
  elements(N-1,1)=0; // correct last entry

  // 3. Construct BoundaryMesh from the computed coordinates and elements.
  BoundaryMesh squareMesh(coordinates, elements);

  #else // TEMPLATE
  Eigen::MatrixXd coordinates(N,2);
  Eigen::MatrixXi elements(N,2);
  // TODO: CREATE MESH FOR [0, 0.5]^2
  
  BoundaryMesh squareMesh(coordinates, elements);
  #endif // TEMPLATE

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
  #if SOLUTION
  // 1. Create boundary vertices for gamma([-1,1])
  // - Create auxiliary vector for interval I
  Eigen::VectorXd Xi; Xi.setLinSpaced(N+1, -1, 1);
  // - Create Matrix containing the coordinates of the mesh vertices
  Eigen::MatrixXd coordinates(N,2);
  // - Compute the coordinates and fill the matrix
  for(int i=0; i<N; i++){
    coordinates.row(i) = gamma(Xi(i)).transpose();
  }

  // 2. Create matrix specifying the indices of the vertices of each element
  Eigen::MatrixXi elements(N,2);
  elements <<  Eigen::VectorXi::LinSpaced(N, 0, N-1),
               Eigen::VectorXi::LinSpaced(N, 1, N);
  elements(N-1,1)=0; // correct last entry

  // 3. Construct BoundaryMesh from the computed coordinates and elements.
  BoundaryMesh gammaMesh(coordinates, elements);

  #else // TEMPLATE
  Eigen::MatrixXd coordinates(N,2);
  Eigen::MatrixXi elements(N,2);
  // TODO: CREATE MESH FOR GAMMA([-1,1])
  
  BoundaryMesh gammaMesh(coordinates, elements);
  #endif // TEMPLATE

  return gammaMesh;
}
/* SAM_LISTING_END_1 */

#endif
