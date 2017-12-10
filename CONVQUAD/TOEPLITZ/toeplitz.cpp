/***********************************************************************
 *                                                                     *
 * Code for Course "Advanced Numerical Methods for CSE"                *
 * (Prof. Dr. R. Hiptmair)                                             * 
 * Author: R.H.                                                        *
 * Date: Nov 18, 2017                                                  *
 * (C) Seminar for Applied Mathematics, ETH Zurich                     *
 * This code can be freely used for non-commercial purposes as long    *
 * as this header is left intact.                                      *
 ***********************************************************************/
#include <iostream>
#include <Eigen/Dense>
#include <unsupported/Eigen/FFT>

using namespace std;
using namespace Eigen;

/** @brief DFT based periodic convolution */
/* SAM_LISTING_BEGIN_0 */
VectorXcd pconvfft(const VectorXcd& u, const VectorXcd& x) {
  FFT<double> fft; // Object implementing FFT algorithm
  VectorXcd tmp = ( fft.fwd(u) ).cwiseProduct( fft.fwd(x) );
  return fft.inv(tmp);
}
/* SAM_LISTING_END_0 */

/** @brief Construction of Toeplitz matrix 
    Filling of a dense matrix with Teoplitz structure from generating sequence */
/* SAM_LISTING_BEGIN_1 */
template <typename VECTOR>
MatrixXcd toeplitz(int m,int n,const VECTOR &u) {
  const size_t l_seq=u.size();
  if (l_seq != (m+n-1)) throw(runtime_error("sequence length mismatch"));
  MatrixXcd T(m,n);
  for (int l=0;l<l_seq;++l) {
    for(int i=max(m-1-l,0),j=max(0,l-m+1);(i<m) && (j<n);++i,++j) T(i,j) = u[l];
  }
  return(T);
}
/* SAM_LISTING_END_1 */

int main(int, char**) {
  VectorXd u(9); u << 1,2,3,4,5,6,7,8,9;
  cout << toeplitz(4,6,u) << endl;
  exit(0);
}
  
