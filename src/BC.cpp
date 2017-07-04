//block cluster class, actually only needed to build the matrices $X_{\sigma,\mu}$
#include "../include/BC.hpp"

// construct XX matrix
void BC::fill_X() {
  Cheby C_ptx(x_l, x_r, deg);		
  Cheby C_pty(y_l, y_r, deg);
  // Chebyshew nodes in interval [x_l,x_r]
  const std::vector<double>& cheb_ptx = C_ptx.cheb_pts();
  // Chebyshew nodes in interval [y_l,y_r]
  const std::vector<double>& cheb_pty = C_pty.cheb_pts();
  // kernel function with constant factor "G_const", 
  // see file "kernel.hpp" for definition
  Kernel func(G_const); 
  for(int i=0; i<=deg;++i)
    for(int j=0; j<=deg;++j)
      XX(i,j)=func(cheb_ptx[i],cheb_pty[j]);
};
