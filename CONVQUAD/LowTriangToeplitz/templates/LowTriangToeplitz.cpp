#include <Eigen/Dense>
#include <unsupported/Eigen/FFT>
#include <cmath>
#include <iostream>

using namespace Eigen;
using namespace std;


VectorXcd pconvfft(const VectorXcd& u, const VectorXcd& x)
{
    FFT<double> fft;
    VectorXcd tmp = ( fft.fwd(u) ).cwiseProduct( fft.fwd(x) );
    return fft.inv(tmp);
}


VectorXcd myconv(const VectorXcd& h, const VectorXcd& x) {
  const long n = h.size();
  // Zero padding, cf. \eqref{eq:zeropad}
  VectorXcd hp(2*n - 1), xp(2*n - 1);
  hp << h, VectorXcd::Zero(n - 1);
  xp << x, VectorXcd::Zero(n - 1);
  // Periodic discrete convolution of length \Blue{$2n-1$}, \cref{cpp:pconffft}
  return pconvfft(hp, xp);
}


/* @brief Multiply two lower triangular Toeplitz matrices
 * \param f Vector of entries of first lower triangular Toeplitz matrix
 * \param g Vector of entries of second lower triangular Toeplitz matrix
 * \\return Vector of entries of output lower triangular Toeplitz matrix
 */
/* SAM_LISTING_BEGIN_0 */
VectorXcd ltpmult(const VectorXcd& f, const VectorXcd& g)
{
    assert(f.size() == g.size() &&
           "f and g vectors must have the same length!");

    // TODO: Multiply two lower triangular Toeplitz matrices
}
/* SAM_LISTING_END_0 */


VectorXcd toepmult(const VectorXcd& c, const VectorXcd& r,
                          const VectorXcd& x)
{
    assert(c.size() == x.size() &&
           r.size() == x.size() &&
           "c, r, x have different lengths!");
    size_t n = c.size();

    VectorXcd cr_tmp = c;
    cr_tmp.conservativeResize(2*n); cr_tmp.tail(n) = VectorXcd::Zero(n);
    cr_tmp.tail(n-1) = r.tail(n-1).reverse();

    VectorXcd  x_tmp = x;
    x_tmp.conservativeResize(2*n);   x_tmp.tail(n) = VectorXcd::Zero(n);

    VectorXcd y = pconvfft(cr_tmp, x_tmp);
    y.conservativeResize(n);

    return y;
}


/* @brief Solve a linear problem involving a lower triangular Toeplitz matrix
 * \param f Vector of entries of lower triangular Toeplitz matrix
 * \param y Right-hand side of linear problem
 * \\return Solution of linear problem
 */
/* SAM_LISTING_BEGIN_1 */
VectorXcd ltp_solve(const VectorXcd& f, const VectorXcd& y)
{
    assert(f.size() == y.size() &&
           "f and y vectors must have the same length!");
    assert(abs(f(0)) > 1e-10 &&
           "Lower triangular Toeplitz matrix must be invertible!");
    assert(log2(f.size()) == floor(log2(f.size())) &&
           "Size of f must be a power of 2!");

    // TODO: Solve a linear problem involving a lower triangular Toeplitz matrix
}
/* SAM_LISTING_END_1 */


MatrixXcd toeplitz(const VectorXcd& c, const VectorXcd& r)
{
    if(c(0) != r(0)) {
        cerr << "First entries of c and r are different!" <<
        endl << "We assign the first entry of c to the diagonal" << endl;
    }

    // Initialization
    size_t m = c.size();
    size_t n = r.size();
    MatrixXcd T(m, n);

    for(int i=0; i<n; ++i) {
        T.col(i).tail(m-i) = c.head(m-i);
    }
    for(int i=0; i<m; ++i) {
        T.row(i).tail(n-i-1) = r.segment(1,n-i-1);
    } // Do not reassign the diagonal!

    return T;
}


int main() {

    // Initialization
    size_t n = 4;
    VectorXcd c1(n), c2(n), r1(n), r2(n), y(n);
    c1 << 1, 2, 3, 4;
    r1 << 1, 0, 0, 0;
    c2 << 5, 6, 7, 8;
    r2 << 5, 0, 0, 0;
    y  << 9,10,11,12;
    MatrixXcd T1 = toeplitz(c1, r1);
    MatrixXcd T2 = toeplitz(c2, r2);

    cout << "Check that ltpmult is correct" << endl;
    VectorXcd c1c2 = ltpmult(c1, c2);
    MatrixXcd T1T2 = T1*T2;
    cout << "Error = " << (c1c2 - T1T2.col(0)).norm() << endl;

    cout << "Check that ltp_solve is correct" << endl;
    VectorXcd u_rec = ltp_solve(c1, y);
    VectorXcd u_sol = T1.triangularView<Lower>().solve(y);
    cout << "Error = " << (u_rec - u_sol).norm() << endl;
}
