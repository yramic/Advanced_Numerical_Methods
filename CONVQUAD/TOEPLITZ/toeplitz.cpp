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

/** @brief DFT based periodic convolution */
/* SAM_LISTING_BEGIN_0 */
VectorXcd pconvfft(const VectorXcd& u, const VectorXcd& x) {
  FFT<double> fft; // Object implementing FFT algorithm
  VectorXcd tmp = ( fft.fwd(u) ).cwiseProduct( fft.fwd(x) );
  return fft.inv(tmp);
}
/* SAM_LISTING_END_0 */

/** @brief Simple inefficient implementation of periodic convolution */
/* SAM_LISTING_BEGIN_2 */
VectorXcd pconv(const VectorXcd& u, const VectorXcd& x) {
  using idx_t = VectorXcd::Index; // may be unsigned !
  const idx_t n = x.size(); 
  VectorXcd z = VectorXcd::Zero(n);
  // Need signed indices when differences are formed
  for (long k = 0; k < n; ++k) {
    for (long j = 0; j < n; ++j) {
      long ind = (k - j < 0 ? n + k - j : k - j);
      z(k) += u(ind)*x(j);
    }}
  return z;
}
/* SAM_LISTING_END_2 */

/** @brief implementation of discrete convolution by reduction to periodic convolution.
    Note that this function returns a vector twice as long as the input vectors. */
/* SAM_LISTING_BEGIN_3 */
VectorXcd myconv(const VectorXcd& h, const VectorXcd& x) {
  const long n = h.size();
  // Zero padding
  VectorXcd hp(2*n - 1), xp(2*n - 1);
  hp << h, VectorXcd::Zero(n - 1);
  xp << x, VectorXcd::Zero(n - 1);
  // Periodic discrete convolution of length \Blue{$2n-1$}. 
  return pconv(hp, xp);
}
/* SAM_LISTING_END_3 */

/** @brief Straightforward implementation of discrete convolution.
    Note that this function returns a vector twice as long as the input vectors. */
/* SAM_LISTING_BEGIN_4 */
VectorXcd seqconv(const VectorXcd& h, const VectorXcd& x) {
  int n = h.size();
  if (n != x.size()) throw(runtime_error("sequence length mismatch"));
  VectorXcd z = VectorXcd::Zero(2*n-1);
  for (int j=0; j < 2*n-1; ++j) {
    for(int l=max(j-n+1,0);l<min(j-1,n);++l) z[j] += h[j-l]*x[l];
  }
  return z;
}
/* SAM_LISTING_END_4 */

int main(int, char**) {
  // First test: output Toeplitz matrix
  VectorXd u(9); u << 1,2,3,4,5,6,7,8,9;
  cout << toeplitz(4,6,u) << endl;
  // Second test: discrete convolution
  const VectorXcd x = VectorXcd::Random(9);
  const VectorXcd h = VectorXcd::Random(9);
  cout << "Std convolution h*x = " << seqconv(h,x) << endl;
  cout << "Fast convolution h*x = " << myconv(h,x) << endl;
  exit(0);
}
  
