#include <Eigen/Sparse>
#include <cmath>
#include <stdexcept>
#include <iostream>
#include <iomanip>

#include "../Library/meshReader.hpp"
#include "../Library/source/geometry.hpp"
#include "../Library/source/buildW.hpp"
#include "../Library/source/buildV.hpp"
#include "../Library/source/buildK.hpp"
#include "../Library/source/buildM.hpp"

//! Sparse Matrix type. Makes using this type easier.
typedef Eigen::SparseMatrix<double> SparseMatrix;

//! Vector type
typedef Eigen::VectorXd Vector;





int main(int, char**) {

  Eigen::Vector2d a,b,c;
  a << 1.0, 4.5;
  b << 2., 2;
  c << 5.1, 3.;
  std::cout << distancePointToSegment(a,b,c) << std::endl;

  Eigen::MatrixXi els;
  Eigen::MatrixXd coords;
  readMesh("Lshape", els, coords);

  
  Eigen::MatrixXd W;
  computeW(W, coords, els, 0.001);

  const std::string fnameW =  "W_Lshape_Hilbert.dat";
  std::ofstream outW( fnameW.c_str() );
  outW << std::setprecision(18) << W; 
  outW.close( );

  Eigen::MatrixXd V;
  computeV(V, coords, els, 0.001);

  const std::string fnameV =  "V_Lshape_Hilbert.dat";
  std::ofstream outV( fnameV.c_str() );
  outV << std::setprecision(18) << V; 
  outV.close( );

  Eigen::MatrixXd K;
  computeK(K, coords, els, 0.001);

  const std::string fnameK =  "K_Lshape_Hilbert.dat";
  std::ofstream outK( fnameK.c_str() );
  outK << std::setprecision(18) << K; 
  outK.close( );

  Eigen::SparseMatrix<double> M(coords.rows(), coords.rows());
  computeM11(M, coords, els);
  std::cout<< M << std::endl;
  
  std::cout << " Done " << std::endl;

  // loaded them in matlab and compared with Hilbert original implementation.
  // the norm of the difference is zero :)
  
  
  
}




