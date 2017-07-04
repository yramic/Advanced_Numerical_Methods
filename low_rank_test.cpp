#include "include/kernel.hpp"
#include "include/low_rank_app.hpp"

#include <Eigen/Dense>
#include <cmath>
#include <ctime>
#include <iostream>

int main() {

    // Input

    std::cout << "Enter gridsize:" << std::endl;
    unsigned n; std::cin >> n;
    Eigen::VectorXd x = Eigen::VectorXd::LinSpaced(n, 0., (n-1)/n);
    Eigen::VectorXd c = Eigen::VectorXd::Random(n);

    std::cout << "Enter admissibility constant:" << std::endl;
    double eta; std::cin >> eta;

    std::cout << "Enter degree of interpolating polynomials:" << std::endl;
    double d; std::cin >> d;

    Kernel G(1.);

    // Compute exact matrix-vector product

    time_t start1; time(&start1);

    Eigen::MatrixXd M(n,n);
    for(int i=0; i<n; ++i)
        for(int j=0; j<n; ++j)
            M(i,j) = G(x[i], x[j]);
    Eigen::VectorXd f_exact = M * c;

    time_t end1; time(&end1);
    double time_diff1 = std::difftime(end1, start1);

    // Compute approximated matrix-vector product, given admissibility constant 'eta'

    time_t start2; time(&start2);

    LowRankApp lra(G, x, x);
    Eigen::VectorXd f_approx = lra.mvProd(c, eta, d);

//    Eigen::VectorXd c = Eigen::VectorXd::Zero(n);
//    LowRankApp lra(G, x, x);
//    Eigen::MatrixXd M_approx = Eigen::MatrixXd::Zero(n,n);
//    Eigen::VectorXd err_inf(d+1), err_2(d+1);

//    for(unsigned deg=0; deg<=d; ++deg) {
//        for(unsigned i=0; i<nx; ++i) {

//            x(i) = 1.;
//            M_approx.col(i) = lra.mvProd(c, eta, d);
//            x(i) = 0.;
//        }

//        err_inf(deg) = (M - M_approx).lpNorm<Eigen::Infinity>();
//        err_2(deg)   = (M - M_approx).lpNorm<2>();
//    }

    time_t end2; time(&end2);
    double time_diff2 = std::difftime(end2, start2);

    // Compute approximation error

    Eigen::VectorXd diff = f_exact - f_approx;

    std::cout << "Approximation error (l-inf norm): " << diff.lpNorm<Eigen::Infinity>() << std::endl
              << "Approximation error (l-2 norm): "   << diff.lpNorm<2>() << std::endl
              << "Time needed for exact multiplication: " << time_diff1 << std::endl
              << "Time needed for approximate multiplication: " << time_diff2 << std::endl;
}
