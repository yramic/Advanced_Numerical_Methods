#include <iostream>
#include <fstream>
#include <istream>
#include <string> 
#include <cmath>
#include <Eigen/Dense>
// CppHilbert includes
#include "../CppHilbert/Library/source/buildM.hpp"
#include "../CppHilbert/Library/source/geometry.hpp"
// Own includes
#include "MeshGen.hpp"
#include "DirectBEM.hpp"
#include "IndirectBEM.hpp"


//------------------------------------------------------------------------------
// DATA FOR SQUARE AND "KITE-SHAPED DOMAIN"
/*
 * @brief Dirichlet data for the Laplace Dirichlet problem over these two domains
 */
double g(const Eigen::Vector2d& X){
  return sin(X(0)-X(1))*sinh(X(0)+X(1));
}


/*
 * @brief Evaluate on the point X in the element [a,b] the Neumann trace of the 
 *        exact solution for the Laplace Dirichlet problem over these two domains 
 *        (for the dirichlet data given above).
 */
/* SAM_LISTING_BEGIN_0 */
double TNu(const Eigen::Vector2d & X, const Eigen::Vector2d & a,
	   const Eigen::Vector2d & b){
  Eigen::Vector2d n = unitNormal(a,b);
  Eigen::Vector2d grad;
  
  grad<< cos(X(0)-X(1))*sinh(X(0)+X(1)) + sin(X(0)-X(1))*cosh(X(0)+X(1)),
        -cos(X(0)-X(1))*sinh(X(0)+X(1)) + sin(X(0)-X(1))*cosh(X(0)+X(1));
  

  return grad.dot(n);
}
/* SAM_LISTING_END_0 */


//------------------------------------------------------------------------------
/*
 * @brief Compute coefficient vector of \f$\mathcal{S}^{-1}_0(\mathcal{G}\f$ 
 *        corresponding to the Neumann trace of the exact solution for the Laplace 
 *        Dirichlet problem over the corresponding domain (given through the mesh)
 *
 * \tparam TNFUNC Function taking a point X and the end-points of the segment 
 *                [a,b] and returning the evaluation of the Neumann trace of the 
 *                exact solution on the domain of interest.
 * \param[in] mesh BoundaryMesh object corresponding to the boundary of the domain 
 *                 of interest.
 */
/* SAM_LISTING_BEGIN_1 */
template <typename TNFUNC>
Eigen::VectorXd ComputeTNu(const TNFUNC& tnu, const BoundaryMesh& mesh){
  Eigen::MatrixXi elems = mesh.getMeshElements();
  Eigen::VectorXd tnuval(mesh.numElements());
  for(int k=0; k<mesh.numElements(); k++){
    int aidx = mesh.getElementVertex(k,0);
    int bidx = mesh.getElementVertex(k,1);
    const Eigen::Vector2d& a = mesh.getVertex(aidx);
    const Eigen::Vector2d& b = mesh.getVertex(bidx);
    tnuval(k) = tnu((a+b).eval()/2.,a,b);
  }
  return tnuval;
}
/* SAM_LISTING_END_1 */


//------------------------------------------------------------------------------
int main() {

  auto squareMesh = createMiniSquareMesh(16);

  std::function<Eigen::Vector2d(const double&)> gamma = [](const double& t){
    Eigen::Vector2d res;
    res << 0.25*cos(M_PI*t) + 0.1625*cos(2*M_PI*t), 0.375*sin(M_PI*t);
    return res;
  };

  auto gammaMesh = createMeshwithGamma(gamma, 8);
  gammaMesh.writeMeshToFile("gamma");

  std::cout << "Done playing with meshes" << std::endl;

  //--------------------------------------------------------
  // RESULTS FOR KITE DOMAIN
  //--------------------------------------------------------
  std::cout << "RESULTS FOR KITE DOMAIN" << std::endl;
  /* SAM_LISTING_BEGIN_2 */
  int Nl = 7; //Number of levels
  Eigen::VectorXd errorD1(Nl), errorD2(Nl), errorD1L2(Nl), errorD2L2(Nl);
  Eigen::VectorXd errorI1(Nl), errorI2(Nl);
  Eigen::VectorXi Nall(7); Nall<< 50,100,200,400,800,1600,3200;
  Eigen::Vector2d X({0.,0.3});
  for(int k=0; k<Nl; k++){
    int N = Nall(k);
    std::cout << "Using N = " << N << " elements" << std::endl;
    auto mesh = createMeshwithGamma(gamma, N);

    Eigen::MatrixXd V; computeV(V, mesh, 1e-05);
    Eigen::SparseMatrix<double> M00(mesh.numElements(), mesh.numElements());
    computeM00(M00, mesh);
    Eigen::VectorXd solex = ComputeTNu(TNu, mesh);
    
    std::cout << "Solving 1st kind Direct. ";
    Eigen::VectorXd sol1 = DirectFirstKind::solveDirichlet(mesh, g);
    errorD1(k) = sqrt((sol1-solex).transpose()*V*(sol1-solex));
    errorD1L2(k) = sqrt((sol1-solex).transpose()*M00*(sol1-solex)); 
    std::cout << "Obtained error " << errorD1(k) << std::endl;
    
    std::cout << "Solving 2nd kind Direct. " ;
    Eigen::VectorXd sol2 = DirectSecondKind::solveDirichlet(mesh, g);
    errorD2(k) = sqrt((sol2-solex).transpose()*V*(sol2-solex));
    errorD2L2(k) = sqrt((sol2-solex).transpose()*M00*(sol2-solex));
    std::cout << "Obtained error " << errorD2(k) << std::endl;    

    std::cout << "Solving 1st kind Indirect : ";
    Eigen::VectorXd sol1i = IndirectFirstKind::solveDirichlet(mesh, g);
    double solEval1i = IndirectFirstKind::reconstructSolution(X, sol1i, mesh);
    errorI1(k) = fabs(solEval1i - g(X) );
    std::cout << "Obtained error " << errorI1(k) << std::endl;
    
    std::cout << "Solving 2nd kind Indirect : " ;
    Eigen::VectorXd sol2i = IndirectSecondKind::solveDirichlet(mesh, g);
    double solEval2i = IndirectSecondKind::reconstructSolution(X, sol2i, mesh);
    errorI2(k) = fabs(solEval2i - g(X) );
    std::cout << "Obtained error " << errorI2(k) << std::endl;
    
  }
  /* SAM_LISTING_END_2 */

  
  // OUTPUT ERRORS
  {
  std::ofstream out_error("DBEM1stK_errors.txt");
  out_error << std::setprecision(18) << errorD1; 
  out_error.close( );
  }
  {
  std::ofstream out_error("DBEM2ndK_errors.txt");
  out_error << std::setprecision(18) << errorD2; 
  out_error.close( );
  }
    {
  std::ofstream out_error("DBEM1stK_L2errors.txt");
  out_error << std::setprecision(18) << errorD1L2; 
  out_error.close( );
  }
  {
  std::ofstream out_error("DBEM2ndK_L2errors.txt");
  out_error << std::setprecision(18) << errorD2L2; 
  out_error.close( );
  }
  {
  std::ofstream out_error("IBEM1stK_errors.txt");
  out_error << std::setprecision(18) << errorI1; 
  out_error.close( );
  }
  {
  std::ofstream out_error("IBEM2ndK_errors.txt");
  out_error << std::setprecision(18) << errorI2; 
  out_error.close( );
  }
  std::ofstream out_N("BEM_N.txt");
  out_N << Nall.segment(0,Nl); 
  out_N.close( );
  
    
  return 0;

}
