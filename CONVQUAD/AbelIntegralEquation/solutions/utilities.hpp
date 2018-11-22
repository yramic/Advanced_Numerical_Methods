#include <iostream>
#include <cmath>
#include <Eigen/Dense>
#include <unsupported/Eigen/FFT>

using namespace Eigen;
using namespace std;


/* @brief Create a real lower triangular Toeplitz matrix
 * \param c real Vector of entries of first column of the Toeplitz matrix
 * \param r real Vector of entries of first row of the Toeplitz matrix
 * \param x real Vector
 * \\return Lower triangular Toeplitz matrix
 */
MatrixXd toeplitz_triangular(const VectorXd& c)
{
    size_t n = c.size();
    MatrixXd T = MatrixXd::Zero(n, n);
    for(int i=0; i<n; ++i) {
        T.col(i).tail(n-i) = c.head(n-i);
    }
    return T;
}


/* @brief Generate a Toeplitz matrix
 * \param c real Vector of entries of first column of the Toeplitz matrix
 * \param r real Vector of entries of first row of the Toeplitz matrix
 * \\return Toeplitz matrix
 */
MatrixXd toeplitz(const VectorXd& c, const VectorXd& r)
{
    if(c(0) != r(0)) {
        cerr << "First entries of c and r are different!" <<
        endl << "We assign the first entry of c to the diagonal" << endl;
    }

    // Initialization
    size_t m = c.size();
    size_t n = r.size();
    MatrixXd T(m, n);
    
    for(int i=0; i<n; ++i) {
        T.col(i).tail(m-i) = c.head(m-i);
    }
    for(int i=0; i<m; ++i) {
        T.row(i).tail(n-i-1) = r.segment(1,n-i-1);
    } // Do not reassign the diagonal!

    return T;
}


/* @brief Multiply a Circulant matrix with a vector, using FFT
 * \param u complex Generating vector for the Circulant matrix
 * \param x complex Vector
 * \\return Circulant(u)*x
 */
VectorXcd pconvfft(const VectorXcd& u, const VectorXcd& x)
{
    FFT<double> fft;
    VectorXcd tmp = ( fft.fwd(u) ).cwiseProduct( fft.fwd(x) );
    return fft.inv(tmp);
}


/* @brief Compute discrete convolution using FFT
 * \param h complex Vector
 * \param x complex Vector
 * \\return Convolution(h, x)
 */
VectorXcd myconv(const VectorXcd& h, const VectorXcd& x) {
  const long n = h.size();
  // Zero padding, cf. \eqref{eq:zeropad}
  VectorXcd hp(2*n - 1), xp(2*n - 1);
  hp << h, VectorXcd::Zero(n - 1);
  xp << x, VectorXcd::Zero(n - 1);
  // Periodic discrete convolution of length \Blue{$2n-1$}, \cref{cpp:pconffft}
  return pconvfft(hp, xp);
}


/* @brief Multiply a Toeplitz matrix with a vector, uses pconvfft
 * \param c real Vector of entries of first column of the Toeplitz matrix
 * \param r real Vector of entries of first row of the Toeplitz matrix
 * \param x real Vector
 * \\return toeplitz(c,r)*x
 */
VectorXd toepMatVecMult(const VectorXd& c, const VectorXd& r, const VectorXd& x)
{
    assert(c.size() == x.size() &&
           r.size() == x.size() &&
           "c, r, x have different lengths!");
    
    size_t n = c.size();
    VectorXd cr_tmp(2*n), x_tmp(2*n);
    
    cr_tmp.head(n) = c;
    cr_tmp.tail(n) = VectorXd::Zero(n);
    cr_tmp.tail(n-1) = r.tail(n-1).reverse();
    
    x_tmp.head(n) = x; 
    x_tmp.tail(n) = VectorXd::Zero(n);
    
    VectorXcd y = pconvfft(cr_tmp, x_tmp);
    
    return y.head(n).real();
}


/* @brief Multiply two lower triangular Toeplitz matrices
 * \param f real Vector of entries of first lower triangular Toeplitz matrix
 * \param g real Vector of entries of second lower triangular Toeplitz matrix
 * \\return real Vector of entries of output lower triangular Toeplitz matrix
 */
VectorXd ltpMult(const VectorXd& f, const VectorXd& g)
{
    assert(f.size() == g.size() &&
           "f and g vectors must have the same length!");

    size_t n = f.size();
    return toepMatVecMult(f, VectorXd::Zero(n), g);
}


/* @brief Solve a linear problem involving a lower triangular Toeplitz matrix
 * \param f real Vector of entries of lower triangular Toeplitz matrix
 * \param y real Right-hand side of linear problem
 * \\return Solution of linear problem
 */
VectorXd ltpSolve(const VectorXd& f, const VectorXd& y)
{
    assert(f.size() == y.size() &&
           "f and y vectors must have the same length!");
    assert(abs(f(0)) > 1e-10 &&
           "Lower triangular Toeplitz matrix must be invertible!");
    assert(log2(f.size()) == floor(log2(f.size())) &&
           "Size of f must be a power of 2!");

    size_t n = f.size();
    if(n == 1) {
        return y.cwiseQuotient(f);
    }
    
    VectorXd u_head = ltpSolve(f.head(n/2), y.head(n/2));
    VectorXd t = y.tail(n/2) - toepMatVecMult(f.tail(n/2), f.segment(1,n/2).reverse(), u_head);
    VectorXd u_tail = ltpSolve(f.head(n/2), t);
    
    VectorXd u(n); u << u_head, u_tail;
    return u;
}



// Convolution approximation routines, use Convolution Quadrature

// Implicit-Euler
template<typename FUNC>
VectorXd cq_ieul_abel_test(const FUNC& u, size_t N)
{
    VectorXd w(N+1); w(0) = 1.;
    for(int l=1; l<N+1; ++l) {
        w(l) = w(l-1) * (l - 0.5) / l; // denominator is factorial
    }
    w *= sqrt(M_PI/N);
    
    // Solve the convolution quadrature:
    
    VectorXd grid = VectorXd::LinSpaced(N+1,0.,1.);
    VectorXd u_N(N+1);
    for(int i=0; i<N+1; ++i) {
        u_N(i) = u(grid(i));
    }
    
    VectorXd y = ltpMult(w, u_N);
    
    return y;
}

// BDF-2
template<typename FUNC>
VectorXd cq_bdf2_abel_test(const FUNC& u, size_t N)
{
    VectorXd w1(N+1); w1(0) = 1.;
    for(int l=1; l<N+1; ++l) {
        w1(l) = w1(l-1) * (l - 0.5) / l; // denominator is factorial
    }
    
    VectorXd w2 = w1;
    for(int l=1; l<N+1; ++l) {
        w2(l) /= pow(3,l);
    }
    
    VectorXd w = myconv(w1, w2).head(N+1).real();
    w *= sqrt(M_PI/N) * sqrt(2./3.);
    
    // Solve the convolution quadrature:
    
    VectorXd grid = VectorXd::LinSpaced(N+1,0.,1.);
    VectorXd u_N(N+1);
    for(int i=0; i<N+1; ++i) {
        u_N(i) = u(grid(i));
    }
    
    VectorXd y = ltpMult(w, u_N);
    
    return y;
}



// Convergence tests for Convolution approximation routines

// Implicit-Euler
void test_approx_eq_ieul_abel()
{
    auto u = [](double t) { return 2./M_PI*sqrt(t); };
    auto y = [](double t) { return t; };
        
    cout << "\n\nConvolution Quadrature, Implicit-Euler, Function approximation test\n" << endl;
    for(int N=16; N<=4906; N*=2) {
        
        VectorXd grid = VectorXd::LinSpaced(N+1,0.,1.);
        VectorXd y_ex(N+1);
        for(int i=0; i<N+1; ++i) {
            y_ex(i) = y(grid(i));
        }

        VectorXd y_app = cq_ieul_abel_test(u, N);
        VectorXd diff  = y_ex - y_app;
        double err_max = diff.cwiseAbs().maxCoeff();
        cout <<   "N = " << N << setw(15)
             << "Max = "
             << scientific << setprecision(3)
             << err_max << endl;
    }
}

// BDF-2
void test_approx_eq_bdf2_abel()
{
    auto u = [](double t) { return 2./M_PI*sqrt(t); };
    auto y = [](double t) { return t; };
        
    cout << "\n\nConvolution Quadrature, BDF-2, Function approximation test\n" << endl;
    for(int N=16; N<=4906; N*=2) {
        
        VectorXd grid = VectorXd::LinSpaced(N+1,0.,1.);
        VectorXd y_ex(N+1);
        for(int i=0; i<N+1; ++i) {
            y_ex(i) = y(grid(i));
        }

        VectorXd y_app = cq_bdf2_abel_test(u, N);
        VectorXd diff  = y_ex - y_app;
        double err_max = diff.cwiseAbs().maxCoeff();
        cout <<   "N = " << N << setw(15)
             << "Max = "
             << scientific << setprecision(3)
             << err_max << endl;
    }
}


// End of file
