// header for block cluster class, actually only needed to build the matrices
// $X_{\sigma,\mu}$
#ifndef BC_H
#define BC_H

#include <vector>
#include "Eigen/Dense"
#include "node.hpp"
#include "cheby.hpp"
#include "kernel.hpp"

class BC {
 public:
  // constructors
  // default constructors, creates empty instance with no information
 BC():XX(Eigen::MatrixXd::Zero(0,0)), x_l(0), x_r(0), y_l(0), y_r(0), deg(0), G_const(1) { }
  // constructor taking box arguments and the degree of Chebychev interpolant
  BC(double xl, double xr, double yl, double yr, int degree,double K=1):
  XX(Eigen::MatrixXd::Zero(degree+1,degree+1)), x_l(xl), x_r(xr), y_l(yl), y_r(yr), deg(degree), G_const(K)
    { fill_X(); }

  virtual ~BC() {};   //destructor
  
  // access member function
  // returns the matrix $X_{\sigma,\mu}$, where $\sigma$ and $\mu$ denote the clusters  
  const Eigen::MatrixXd& get_matrix() const { return XX; }
  
  // computation member functions
  void fill_X();	// construct matrix $X_{\sigma,\mu}$
 private:
  Eigen::MatrixXd XX; 	// matrix $X_{\sigma,\mu}$, where $\sigma$ and $\mu$ 
                        // denote the corresponding clusters
  double x_l;		// left boundary of bounding box of *xcluster		
  double x_r;		// right boundary of bounding box of *xcluster		
  double y_l;		// left boundary of bounding box of *ycluster		
  double y_r;		// right boundary of bounding box of *ycluster		
  int deg;		// degree of interpolation polynomials
  double G_const;	// constant for the kernel function, see file "kernel.hpp"
};
#endif 
