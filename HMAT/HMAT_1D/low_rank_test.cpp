#include "include/kernel.hpp"
#include "include/low_rank_app.hpp"
#include "include/point.hpp"

#include <Eigen/Dense>
#include <chrono>
#include <cmath>
#include <iostream>


int main() {

    // Input

    std::cout << "Enter gridsize:" << std::endl;
    unsigned n; std::cin >> n;
//    unsigned n = 1000;

    // grid points initialization
    Eigen::VectorXd grid = Eigen::VectorXd::LinSpaced(n, 0., (n-1.)/n);
    Eigen::VectorXd    c = Eigen::VectorXd::Random(n);

//    std::cout << "Enter admissibility constant:" << std::endl;
//    double eta; std::cin >> eta;
    double eta = 0.5;

//    std::cout << "Enter degree of interpolating polynomials:" << std::endl;
//    unsigned d; std::cin >> d;
    unsigned d = 3;

    KernelLog G(1.); // Kernel initilization

    std::vector<Point> GPoints; // initalizing Grid Points properties
    GPoints.reserve(n);
    int k = 0;
    for(int i=0; i<n; ++i){
        Point p;
        p.setId(k);
        p.setX(grid[i]);
        k++;
        GPoints.push_back(p);
    }

    // Compute exact matrix-vector product

    auto start1 = std::chrono::system_clock::now();

    Eigen::MatrixXd M(n,n);
    for(int i=0; i<n; ++i)
        for(int j=0; j<n; ++j)
            M(i,j) = G(GPoints[i].getX(), GPoints[j].getX());
    Eigen::VectorXd f_exact = M * c;

    auto end1 = std::chrono::system_clock::now();
    std::chrono::duration<double> time_diff1 = end1 - start1;

    // Compute approximated matrix-vector product, given admissibility constant 'eta'

    auto start2 = std::chrono::system_clock::now();

    LowRankApp<> HMat(&G, GPoints, eta, d);
    Eigen::VectorXd f_approx = HMat.mvProd(c);

    auto end2 = std::chrono::system_clock::now();
    std::chrono::duration<double> time_diff2 = end2 - start2;

    // Compute approximation error

    Eigen::VectorXd diff = f_exact - f_approx;

    std::cout << "Approximation error (l-inf norm): " << diff.lpNorm<Eigen::Infinity>() << std::endl
              << "Approximation error (l-2 norm): "   << diff.lpNorm<2>() << std::endl
              << "Relative Approximation error (l-2 norm): "   << diff.lpNorm<2>()/f_exact.lpNorm<2>() << std::endl
              << "Time needed for exact multiplication: "       << time_diff1.count() << " s" << std::endl
              << "Time needed for approximate multiplication: " << time_diff2.count() << " s" << std::endl;
}
