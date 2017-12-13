#include <Eigen/Dense>
#include <unsupported/Eigen/FFT>
#include <cmath>
#include <iostream>


Eigen::VectorXcd pconvfft(const Eigen::VectorXcd& u, const Eigen::VectorXcd& x)
{
    Eigen::FFT<double> fft;
    Eigen::VectorXcd tmp = ( fft.fwd(u) ).cwiseProduct( fft.fwd(x) );
    return fft.inv(tmp);
}


Eigen::VectorXcd myconv(const Eigen::VectorXcd& h, const Eigen::VectorXcd& x) {
  const long n = h.size();
  // Zero padding, cf. \eqref{eq:zeropad}
  Eigen::VectorXcd hp(2*n - 1), xp(2*n - 1);
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
Eigen::VectorXcd ltpmult(const Eigen::VectorXcd& f, const Eigen::VectorXcd& g)
{
    assert(f.size() == g.size() &&
           "f and g vectors must have the same length!");

    size_t n = f.size();
    return myconv(f, g).head(n);
}
/* SAM_LISTING_END_0 */


Eigen::VectorXcd toepmult(const Eigen::VectorXcd& c, const Eigen::VectorXcd& r,
                          const Eigen::VectorXcd& x)
{
    assert(c.size() == x.size() &&
           r.size() == x.size() &&
           "c, r, x have different lengths!");
    size_t n = c.size();

    Eigen::VectorXcd cr_tmp = c;
    cr_tmp.conservativeResize(2*n); cr_tmp.tail(n) = Eigen::VectorXcd::Zero(n);
    cr_tmp.tail(n-1).real() = r.tail(n-1).reverse();

    Eigen::VectorXcd  x_tmp = x;
    x_tmp.conservativeResize(2*n);   x_tmp.tail(n) = Eigen::VectorXcd::Zero(n);

    Eigen::VectorXcd y = pconvfft(cr_tmp, x_tmp);
    y.conservativeResize(n);

    return y;
}


/* @brief Solve a linear problem involving a lower triangular Toeplitz matrix
 * \param f Vector of entries of lower triangular Toeplitz matrix
 * \param y Right-hand side of linear problem
 * \\return Solution of linear problem
 */
/* SAM_LISTING_BEGIN_1 */
Eigen::VectorXcd ltp_solve(const Eigen::VectorXcd& f, const Eigen::VectorXcd& y)
{
    assert(f.size() == y.size() &&
           "f and y vectors must have the same length!");
    assert(f(0) != 0 &&
           "Lower triangular Toeplitz matrix must be invertible!");
    assert(std::log2(f.size) = std::floor(std::log2(f.size)) &&
           "Size of f must be a power of 2!");

    size_t n = f.size();
    if(n == 1) {
        return y.cwiseQuotient(f);
    }

    Eigen::VectorXcd u_head = ltp_solve(f.head(n/2), y.head(n/2));
    Eigen::VectorXcd t = y.tail(n/2) - toepmult(f.tail(n/2), f.segment(1,n/2).reverse(), u_head);
    Eigen::VectorXcd u_tail = ltp_solve(f.tail(n/2), y.tail(n/2));

    Eigen::VectorXcd u(n); u << u_head, u_tail;
    return u;
}
/* SAM_LISTING_END_1 */


Eigen::MatrixXcd toeplitz(const Eigen::VectorXcd& c, const Eigen::VectorXcd& r)
{
    if(c(0) != r(0)) {
        std::cerr << "First entries of c and r are different!" <<
        std::endl << "We assign the first entry of c to the diagonal" << std::endl;
    }

    // Initialization
    size_t m = c.size();
    size_t n = r.size();
    Eigen::MatrixXcd T(m, n);

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
    size_t n = 3;
    Eigen::VectorXcd c1(n), c2(n), r1(n), r2(n), y(n);
    c1 << 1, 2, 3;
    r1 << 1, 0, 0;
    c2 << 4, 5, 6;
    r2 << 4, 0, 0;
    y  << 7, 8, 9;
    Eigen::MatrixXcd T1 = toeplitz(c1, r1);
    Eigen::MatrixXcd T2 = toeplitz(c2, r2);

    std::cout << "Check that ltpmult is correct" << std::endl;
    Eigen::VectorXcd c1c2 = ltpmult(c1, c2);
    Eigen::MatrixXcd T1T2 = T1*T2;
    std::cout << "Error = " << (c1c2 - T1T2.col(0)).norm() << std::endl;

    std::cout << "Check that ltp_solve is correct" << std::endl;
    Eigen::VectorXcd u_rec = ltp_solve(c1, y);
    Eigen::VectorXcd u_sol = T1.triangularView<Eigen::Lower>().solve(y);
    std::cout << "Error = " << (u_rec - u_sol).norm() << std::endl;
}
