//header of the class used to find the Chebyshew nodes and the weights of the Lagrange polynomials
#ifndef CHEBY_H
#define CHEBY_H

#include <vector>
#include "Eigen/Dense"
#include <math.h>
#include "cheby.hpp"


class cheby{
 public:
  // constructors
  // default constructor (very silly)
  cheby() : xl(0), xr(0), deg_pol(0), t_k(0), omega(Eigen::VectorXd::Zero(0)){}
  // another constructor
  cheby(double x_l, double x_r, int degree_pol, std::vector<double> chebpt): xl(x_l), xr(x_r), deg_pol(degree_pol), t_k(chebpt), omega(Eigen::VectorXd::Zero(degree_pol+1)) { set_omega();}
  // meaningful constructor
  cheby(double x_l, double x_r, int degree_pol);
  
  // deconstructor
  virtual ~cheby(){};
	
  // member access fcns
  double left_bd() const{ return xl;}		     
  // returns the left boundary xl of the domain
  double right_bd() const{ return xr;}
  // returns the right boundary xr of the domain
  int get_degree() const{return deg_pol;}
  // returns the degree deg_pol of the interpolating polynomials
  const std::vector<double>& cheb_pts() const{ return t_k;}
  // returns the Chebyshew nodes on domain [xl,xr]
  const Eigen::VectorXd& get_omega() const{ return omega;}	
  // returns the weights of the Lagrange polynomials

  // set members
  void set_leftbd(double lbd) {xl=lbd;}
  // set xl to "lbd"
  void set_rightbd(double rbd){xr=rbd;}
  // set xr to "rbd"			
  void set_deg(int degp){deg_pol=degp;}
  // set deg_pol to "degp"
  void set_chebpt();
  // compute Chebyshew nodes on domain [xl,xr]
  void set_omega(); 
  // PRE: Chebyshew points already computed, compute omega, the weights of the Lagrange polynomials

 private:
  double xl; // left boundary of the domain, on which we construct the Chebyshew points
  double xr; // right boundary of the domain, on which we construct the Chebyshew points
  int deg_pol; // degree of polynomials, $\in \mathbb{N}_{>0}$ please
  std::vector<double> t_k; // Chebyshew points
  Eigen::VectorXd omega;  // weights of Lagrange polynomials
};
#endif 
