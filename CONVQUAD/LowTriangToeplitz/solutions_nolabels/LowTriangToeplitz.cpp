//// 
//// Copyright (C) 2017 SAM (D-MATH) @ ETH Zurich
//// Author(s): curzuato < > 
//// Contributors:  dcasati 
//// This file is part of the AdvNumCSE repository.
////
#include <Eigen/Dense>
#include <unsupported/Eigen/FFT>
#include <cmath>
#include <iostream>
#include <ctime>

using namespace Eigen;
using namespace std;


/* @brief Generate a Toeplitz matrix
 * \param c Vector of entries of first column of the Toeplitz matrix
 * \param r Vector of entries of first row of the Toeplitz matrix
 * \\return Toeplitz matrix
 */
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


/* @brief Multiply a Circulant matrix with a vector, using FFT
 * \param u Generating vector for the Circulant matrix
 * \param x Vector
 * \\return Circulant(u)*x
 */
VectorXcd pconvfft(const VectorXcd& u, const VectorXcd& x)
{
    FFT<double> fft;
    VectorXcd tmp = ( fft.fwd(u) ).cwiseProduct( fft.fwd(x) );
    return fft.inv(tmp);
}


/* @brief Multiply a Toeplitz matrix with a vector, uses pconvfft
 * \param c Vector of entries of first column of the Toeplitz matrix
 * \param r Vector of entries of first row of the Toeplitz matrix
 * \param x Vector
 * \\return toeplitz(c,r)*x
 */
VectorXcd toepMatVecMult(const VectorXcd& c, const VectorXcd& r, const VectorXcd& x)
{
    assert(c.size() == x.size() &&
           r.size() == x.size() &&
           "c, r, x have different lengths!");
    
    size_t n = c.size();
    VectorXcd cr_tmp(2*n), x_tmp(2*n);
    
    cr_tmp.head(n) = c;
    cr_tmp.tail(n) = VectorXcd::Zero(n);
    cr_tmp.tail(n-1) = r.tail(n-1).reverse();
    
    x_tmp.head(n) = x; 
    x_tmp.tail(n) = VectorXcd::Zero(n);
    
    VectorXcd y = pconvfft(cr_tmp, x_tmp);
    
    return y.head(n);
}


/* @brief Multiply two lower triangular Toeplitz matrices
 * \param f Vector of entries of first lower triangular Toeplitz matrix
 * \param g Vector of entries of second lower triangular Toeplitz matrix
 * \\return Vector of entries of output lower triangular Toeplitz matrix
 */
VectorXcd ltpMult(const VectorXcd& f, const VectorXcd& g)
{
    assert(f.size() == g.size() &&
           "f and g vectors must have the same length!");

    size_t n = f.size();
    return toepMatVecMult(f, VectorXcd::Zero(n), g);
}


/* @brief Solve a linear problem involving a lower triangular Toeplitz matrix
 * \param f Vector of entries of lower triangular Toeplitz matrix
 * \param y Right-hand side of linear problem
 * \\return Solution of linear problem
 */
VectorXcd ltpSolve(const VectorXcd& f, const VectorXcd& y)
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
    
    VectorXcd u_head = ltpSolve(f.head(n/2), y.head(n/2));
    VectorXcd t = y.tail(n/2) - toepMatVecMult(f.tail(n/2), f.segment(1,n/2).reverse(), u_head);
    VectorXcd u_tail = ltpSolve(f.head(n/2), t);
    
    VectorXcd u(n); u << u_head, u_tail;
    return u;
}


// check accuracy ltpMult
void test_accuracy_ltpMult() {

    size_t n = 4;
    VectorXcd c1(n), c2(n), r1(n), r2(n), y(n);
    c1 << 1, 2, 3, 4;
    r1 << 1, 0, 0, 0;
    c2 << 5, 6, 7, 8;
    r2 << 5, 0, 0, 0;
    MatrixXcd T1 = toeplitz(c1, r1);
    MatrixXcd T2 = toeplitz(c2, r2);
    
    cout << "\nCheck that ltpMult is correct" << endl;
    VectorXcd c1c2 = ltpMult(c1, c2);
    MatrixXcd T1T2 = T1*T2;
    cout << "Error = " << (c1c2 - T1T2.col(0)).norm() << endl;
}


// measure time ltpMult
void time_measure_ltpMult() {
    
    size_t nl = 12;
    size_t n_start = 4;
    size_t n_end = n_start*pow(2,nl-1);
    int num_repititions = 6;
    
    VectorXd error(nl), et_slow(nl), et_fast(nl);
    clock_t start_time, end_time;
    double et_sum;
    
    cout << "\nMatrix size, start: " << n_start << endl;
    cout << "Matrix size, end: " << n_end << endl;
    cout << "Number of matrices: " << nl << "\n" << endl;
    
    size_t n = n_start;
    for (int l=0; l<nl; l++) {
        
        VectorXcd c(n), r(n), v(n);
        c = VectorXcd::Random(n);
        v = VectorXcd::Constant(n, 1.0);
        r.setZero();
        r(0) = c(0);
        
        MatrixXcd T = toeplitz(c, r);
        
        et_sum = 0;
        VectorXcd T_mult_v;
        for (int k=0; k<num_repititions; k++) {
            start_time = clock();
            T_mult_v = T*v;
            end_time = clock();
            if (k>0)
                et_sum+= double(end_time-start_time)/CLOCKS_PER_SEC;
        }
        et_slow(l) = et_sum/(num_repititions-1);
        
        et_sum = 0;
        VectorXcd c_conv_v;
        for (int k=0; k<num_repititions; k++) {
            start_time = clock();
            c_conv_v = ltpMult(c, v);
            end_time = clock();
            if (k>0)
                et_sum+= double(end_time-start_time)/CLOCKS_PER_SEC;
        }
        et_fast(l) = et_sum/(num_repititions-1);
        
        error(l) = (c_conv_v - T_mult_v.col(0)).norm();
        cout << l << "\t" << n << "\t" << error(l) << "\t" << et_slow(l) << "\t" << et_fast(l) << endl;
        
        n*=2;
    }

}


// check accuracy ltpSolve
void test_accuracy_ltpSolve() {

    size_t n = 4;
    VectorXcd c(n), r(n), y(n);
    c << 1, 2, 3, 4;
    r.setZero(); r(0) = c(0);
    y  << 9,10,11,12;
    MatrixXcd T = toeplitz(c, r);
    
    cout << "\nCheck that ltpSolve is correct" << endl;
    VectorXcd u_rec = ltpSolve(c, y);
    VectorXcd u_sol = T.triangularView<Lower>().solve(y);
    cout << "Error = " << (u_rec - u_sol).norm() << endl;
}


// measure time ltpSolve
void time_measure_ltpSolve() {
    
    size_t nl = 13;
    size_t n_start = 4;
    size_t n_end = n_start*pow(2,nl-1);
    int num_repititions = 6;
    
    VectorXd error(nl), et_slow(nl), et_fast(nl);
    clock_t start_time, end_time;
    double et_sum;
    
    cout << "\nMatrix size, start: " << n_start << endl;
    cout << "Matrix size, end: " << n_end << endl;
    cout << "Number of matrices: " << nl << "\n" << endl;
    
    size_t n = n_start;
    for (int l=0; l<nl; l++) {
        
        VectorXcd c(n), r(n), v(n);
        for (int i=0; i<n; i++) {
            c(i) = i+1;
        }
        v = VectorXcd::Constant(n, 1.0);
        r.setZero();
        r(0) = c(0);
        
        MatrixXcd T = toeplitz(c, r);
        VectorXcd T_mult_v = ltpMult(c, v);
        
        et_sum = 0;
        VectorXcd u_sol;
        for (int k=0; k<num_repititions; k++) {
            start_time = clock();
            u_sol = T.triangularView<Lower>().solve(T_mult_v);
            end_time = clock();
            if (k>0)
                et_sum+= double(end_time-start_time)/CLOCKS_PER_SEC;
        }
        et_slow(l) = et_sum/(num_repititions-1);
        
        et_sum = 0;
        VectorXcd u_rec;
        for (int k=0; k<num_repititions; k++) {
            start_time = clock();
            u_rec = ltpSolve(c, T_mult_v);
            end_time = clock();
            if (k>0)
                et_sum+= double(end_time-start_time)/CLOCKS_PER_SEC;
        }
        et_fast(l) = et_sum/(num_repititions-1);
        
        error(l) = (u_sol - u_rec).norm()/(T_mult_v).norm();
        cout << l << "\t" << n << "\t" << error(l) << "\t" << et_slow(l) << "\t" << et_fast(l) << endl;
        
        n*=2;
    }

}


int main() {
    
    test_accuracy_ltpMult();
    time_measure_ltpMult();
    
    test_accuracy_ltpSolve();
    time_measure_ltpSolve();

}

// End of file
