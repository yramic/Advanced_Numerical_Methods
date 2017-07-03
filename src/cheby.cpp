// class to find the Chebyshew nodes and the weights of the Lagrange polynomials

#include "../include/cheby.hpp"


//  meaningful constructor
cheby::cheby(double x_l, double x_r, int degree_pol)
	: xl(x_l), xr(x_r), deg_pol(degree_pol), t_k(0), omega(Eigen::VectorXd::Zero(degree_pol+1)){	
  double val=0;
  for(int j=0; j<=deg_pol;++j){ // calculate Chebyshew nodes
    val=(xl+xr)/2+(xr-xl)/2*cos((2.*j+1)/(2*(deg_pol+1.))*M_PI);
    t_k.push_back(val);	
  };
  set_omega();
};

// define Chebyshev nodes
void cheby::set_chebpt(){
  double val=0;
  for(int j=0; j<=deg_pol;++j){
    val=(xl+xr)/2+(xr-xl)/2*cos((2.*j+1)/(2.*(deg_pol+1.))*M_PI);
    t_k.push_back(val);	
  };
};

// calculate omega (factors for Langrange polynomials), PRE: Chebyshew points already computed
void cheby::set_omega(){
  omega=Eigen::VectorXd::Zero(deg_pol+1);
  for(int j=0;j<=deg_pol;++j){
    double hc=1;
    for(int k=0;k<j;++k)
      hc=hc*(t_k[j]-t_k[k]);
    for(int k=j+1;k<=deg_pol;++k)
      hc=hc*(t_k[j]-t_k[k]);
    omega(j)=1./hc;
  };
};
			




	
