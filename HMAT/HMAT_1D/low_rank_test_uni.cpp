#include "../HMAT_1D/include/kernel.hpp"
#include "../HMAT_1D/include/low_rank_app.hpp"
#include "../HMAT_1D/include/point.hpp"
#include "../HMAT_1D/include/uni-direct/block_cluster_Y.hpp"
#include "../HMAT_1D/include/uni-direct/node_Y.hpp"

#include <Eigen/Dense>
#include <chrono>
#include <cmath>
#include <fstream>
#include <iomanip>
#include <iostream>


int main() {

    // Input

    std::cout << "Enter gridsize:" << std::endl;
    unsigned n; std::cin >> n;
//    unsigned n = 1000;

//    std::ofstream myfile;
//    myfile.open("test_hmat_1d_time_uni.txt");
//    for(unsigned n : {10, 50, 100, 500, 1000, 5000, 10000, 20000}) {
//    std::ofstream myfile;
//    myfile.open("test_hmat_1d_error_uni.txt");
//    for(unsigned n : {10, 20, 30, 40, 50, 60, 70, 80, 90, 100}) {

    // grid points initialization
    Eigen::VectorXd grid = Eigen::VectorXd::LinSpaced(n, 0., 1.);
    Eigen::VectorXd    c = Eigen::VectorXd::Random(n);

    std::cout << "Enter admissibility constant:" << std::endl;
    double eta; std::cin >> eta;
//    double eta = 0.5;

    std::cout << "Enter degree of interpolating polynomials:" << std::endl;
    unsigned d; std::cin >> d;
//    unsigned d = 2;

    KernelLog G(100.); // Kernel initialization

    std::vector<Point> GPoints; // initializing Grid Points properties
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

    LowRankApp<BlockCluster_Y, Node_Y> HMat(&G, GPoints, eta, d);
    Eigen::VectorXd f_approx = HMat.mvProd(c);

    auto end2 = std::chrono::system_clock::now();
    std::chrono::duration<double> time_diff2 = end2 - start2;

    std::cout << "Number of matrix operations performed for exact matrix: " << n*n << std::endl;

    // Compute approximation error

    Eigen::VectorXd diff = f_exact - f_approx;

    std::cout << "Approximation error (l-inf norm): " << diff.lpNorm<Eigen::Infinity>() << std::endl
              << "Approximation error (l-2 norm): "   << diff.lpNorm<2>() << std::endl
              << "Relative Approximation error (l-2 norm): "    << diff.lpNorm<2>()/f_exact.lpNorm<2>() << std::endl
              << "Time needed for exact multiplication: "       << time_diff1.count() << " s" << std::endl
              << "Time needed for approximate multiplication: " << time_diff2.count() << " s" << std::endl;

//    myfile << "time, " << n << ", " << std::setprecision(10) << time_diff1.count() - time_diff2.count() << std::endl;

//    Eigen::MatrixXd Mtilde(n,n);
//    for(int i=0; i<n; ++i) {
//        Mtilde.col(i) = HMat.mvProd(Eigen::VectorXd::Unit(n,i));
//    }
//    Eigen::MatrixXd diff_M = M - Mtilde;

//    myfile << "error_Frobenius, " << n << ", " << std::setprecision(10) << diff_M.norm()/n << std::endl;
//    myfile << "error_max, "       << n << ", " << std::setprecision(10) << diff_M.cwiseAbs().maxCoeff() << std::endl;
//    }
}
