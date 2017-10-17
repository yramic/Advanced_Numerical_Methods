#include <iostream>
#include <fstream>
#include <istream> 
#include <cmath>
#include <Eigen/Dense>

#include "MeshGen.hpp"


int main() {

  createUnitSquareMesh(4);

  createUnitSquareMesh(8);

  std::function<Eigen::Vector2d(const double&)> gamma = [](const double& t){
    Eigen::Vector2d res;
    res << cos(M_PI*t) + 0.65*cos(2*M_PI*t), 1.5*sin(M_PI*t);
    return res;
  };

  auto gammaMesh = createUnitSquareMesh(gamma, 20);
  gammaMesh.writeMeshToFile("gamma");

    
  return 0;

}
