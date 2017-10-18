#include <iostream>
#include <fstream>
#include <istream> 
#include <cmath>
#include <Eigen/Dense>

#include "MeshGen.hpp"
#include "DirectBEM.hpp"
#include "../CppHilbert/Library/source/buildM.hpp"

double phi(const Eigen::Vector2d & X, const Eigen::Vector2d & a,
	   const Eigen::Vector2d & b){
  Eigen::Vector2d n;
  n << b[1]-a[1], a[0]-b[0];
  n /= (b-a).norm();

  Eigen::Vector2d grad;
  double t = std::atan2(X(1),X(0));
  grad<< 2/3.*std::pow(X.norm(), -1/3.)*(cos(t)*cos(2*t/3.)+sin(t)*sin(2*t/3.)),
         2/3.*std::pow(X.norm(), -1/3.)*(sin(t)*cos(2*t/3.)-cos(t)*sin(2*t/3.));

  return grad.dot(n);
}

Eigen::VectorXd evalPhi(const BoundaryMesh& mesh){
  Eigen::MatrixXi elems = mesh.getMeshElements();
  Eigen::VectorXd phival(mesh.numElements());
  for(int k=0; k<mesh.numElements(); k++){
    int aidx = mesh.getElementVertex(k,0);
    int bidx = mesh.getElementVertex(k,1);
    const Eigen::Vector2d& a = mesh.getVertex(aidx);
    const Eigen::Vector2d& b = mesh.getVertex(bidx);
    phival(k) = phi((a+b).eval()/2.,a,b);
  }
  return phival;
}

int main() {

  auto squareMesh = createUnitSquareMesh(8);

  std::function<Eigen::Vector2d(const double&)> gamma = [](const double& t){
    Eigen::Vector2d res;
    res << cos(M_PI*t) + 0.65*cos(2*M_PI*t), 1.5*sin(M_PI*t);
    return res;
  };

  auto gammaMesh = createMeshwithGamma(gamma, 20);
  gammaMesh.writeMeshToFile("gamma");

  std::cout << "done playing with meshes" << std::endl;
  // TEST DIRECT APPROACH FOR LSHAPE DOMAIN
  BoundaryMesh LshapeMesh("miniLshape");

  std::function<double(const Eigen::Vector2d&)> g = [](const Eigen::Vector2d& X){
    double t = std::atan2(X(1),X(0));
    return std::pow(X.norm(),2./3)*cos(2*t/3.);
  };

  std::cout << "solving 1st kind" << std::endl;
  Eigen::VectorXd sol1 = DirectFirstKind::solveDirichlet(LshapeMesh, g);
     std::cout << "solving 2nd kind" << std::endl;
  Eigen::VectorXd sol2 = DirectSecondKind::solveDirichlet(LshapeMesh, g);
  std::cout << "sol1 : " << sol1.transpose() << std::endl  << std::endl;
  std::cout << "sol2 : " << sol2.transpose() << std::endl;

  Eigen::VectorXd solex = evalPhi(LshapeMesh);
  std::cout << "sol ex : " << solex.transpose() << std::endl;
    
  return 0;

}
