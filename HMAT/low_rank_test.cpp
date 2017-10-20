#include "include/kernel.hpp"
#include "include/low_rank_app.hpp"
#include "include/global_interpolation_app.hpp"
#include "include/point.hpp"

#include <Eigen/Dense>
#include <cmath>
#include <ctime>
#include <iostream>
#include <chrono>

int main() {

    // Input

    //std::cout << "Enter number of edges:" << std::endl;
    //unsigned n; std::cin >> n;

    unsigned n=10;
    //Eigen::VectorXd grid = Eigen::VectorXd::LinSpaced(n, 0., (n-1.)/n);
    std::vector<Point> PPoints; // initalizing Polygon Points properties
    PPoints.reserve(n);

    // x,y,v vextors used for testing(n=16)
    //std::vector<int> x = {14,77,11,38,81,54,7,48,17,30}, y = {62,45,46,87,37,20,72,21,30,90}, v = {59,7,23,71,34,50,48,54,30,60};
    //std::vector<int> x = {24,22,73,63,14,17,39,99,83,41,40,4,83,30,65,23}, y = {24,55,87,91,30,1,9,28,67,85,10,84,57,72,86,56}, v = {66,63,27,46,83,58,46,8,95,57,2,79,34,21,64,95};
    //std::vector<int> x = {18,89,39,39,26,53,4,44}, y = {63,49,66,27,46,4,54,59}, v = {1,2,3,4,5,6,7,8}; // vector with problem with setnearfar
    //std::vector<int> x = {80,27,17,80,23,86,99,89,70,43}, y = {65,74,12,88,44,56,77,31,41,88}, v = {1,2,3,4,5,6,7,8,9,10};
    std::srand(std::time(0));   // initializing points properties randomly
    int tx,ty;
    int f=1;
    for (int i=0; i<n; i++){
        Point p;
        p.setId(i);
        p.setV(std::rand()%100);    // values 0-100
        tx = std::rand()%100;
        ty = std::rand()%100;
        while (f) {             // checking if a point exists 2 times
            f=0;
            for (std::vector<Point>::iterator it=PPoints.begin(); it!=PPoints.end(); it++) {
                if (it->getX() == tx && it->getY() == ty) {
                    f = 1;
                }
            }
            tx = std::rand()%100;
            ty = std::rand()%100;
        }
        p.setX(tx);
        p.setY(ty);
        //p.setX(x[i]);             // for testing
        //p.setY(y[i]);             // for testing
        //p.setV(v[i]);             // for testing
        PPoints.push_back(p);
    }
    Eigen::VectorXd    c = Eigen::VectorXd::Random(n);  // random initialization of vector c

    // vector c initialization for testing
    //Eigen::VectorXd c(n);
    //c << 0.168182, -0.835255, 0.661581, -0.827464, 0.0385175, 0.559586, 0.745203, -0.624835, 0.273339, -0.327815, -0.528608, -0.320443, -0.731551, -0.600541, 0.474063, 0.98031;

    //std::cout << "c" << std::endl;
    //std::cout << c << std::endl;

    // printing the Polygon Points
    //for (std::vector<Point>::iterator it=PPoints.begin(); it!=PPoints.end(); it++)
    //    std::cout << ' ' << it->getX() << ' ' << it->getY() << ' ' << it->getId() << ' ' << it->getV() << std::endl << std::flush;

    //std::cout << "Enter admissibility constant:" << std::endl;
    //double eta; std::cin >> eta;

    // initializing eta for testing
    double eta = 2;

    //std::cout << "Enter degree of interpolating polynomials:" << std::endl;
    //unsigned d; std::cin >> d;

    // initializing degree of interpolation for testing
    //unsigned d=1;
    std::vector<unsigned> t = {1,2,3,5,7,10,20,40,65,70,80,90,100};
    for(auto d : t){
    //Kernel G(1.);             // default kernel for 2d problem
    Kernel4D G;                 // initialization of kernel for 4d problem -1/(2*pi)*log||x-y||
    PolynomialKernel P;         // Polynomial Kernel initialization
    ConstantKernel CK(1);       // Constant Kernel initialization
    GlobalSmoothKernel gskernel;// Global Smooth Periodic Kernel initialization
    GaussKernel gkernel;        // Gauss Kernel initialization

    // Compute exact matrix-vector product

    auto start1 = std::chrono::system_clock::now();

    Eigen::MatrixXd M(n,n);
    for(int i=0; i<n; ++i)
        for(int j=0; j<n; ++j)
            M(i,j) = gskernel(PPoints[i].getX(), PPoints[i].getY(), PPoints[j].getX(), PPoints[j].getY());
    Eigen::VectorXd f_exact = M * c;


    auto end1 = std::chrono::system_clock::now();
    std::chrono::duration<double> time_diff1 = end1 - start1;

    // Compute approximated matrix-vector product, given admissibility constant 'eta'

    auto start2 = std::chrono::system_clock::now();

    //LowRankApp lra(G, grid, grid);    // initialization of low rank approximation for 2D
    LowRankApp lra(&gskernel, PPoints, n);         // initialization of low rank approximation for BEM approx for matrix multiplication

    Eigen::VectorXd f_approx = lra.mvProd(c, eta, d);   // calculation of the low rank approximation

    auto end2 = std::chrono::system_clock::now();
    std::chrono::duration<double> time_diff2 = end2 - start2;

    // debugging code
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


    std::cout << std::endl << "Interpolation with polynomial kernel" << std::endl;
    std::cout << "Number of points: " << n << std::endl;
    std::cout << "Degree: " << d << std::endl;
    // Compute approximation error
    std::cout << std::endl << "Local Interpolation" << std::endl;

    Eigen::VectorXd diff = f_exact - f_approx;

    // printing for testing
    std::cout << "f_exact   f_approx    diff" << std::endl;
    for(int i=0; i<n; i++){
        std::cout << f_exact(i) << "    " << f_approx(i) << "   " << diff(i) << std::endl;
    }
    std::cout << "Approximation error (l-inf norm): " << diff.lpNorm<Eigen::Infinity>() << std::endl
              << "Approximation error (l-2 norm): "   << diff.lpNorm<2>() << std::endl
              << "Time needed for exact multiplication: "       << time_diff1.count() << " s" << std::endl
              << "Time needed for approximate multiplication: " << time_diff2.count() << " s" << std::endl;

    // Global Interpolation
    std::cout << std::endl << "Global Interpolation" << std::endl;
    std::cout << "Global Smooth Kernel (cos)" << std::endl;

    // Compute approximated matrix-vector product, given admissibility constant 'eta'
    auto start3 = std::chrono::system_clock::now();

    GlobalInterpolationApp gip_gskernel(&gskernel, PPoints, n);
    Eigen::VectorXd f_g_approx = gip_gskernel.mvProd(c,d);

    auto end3 = std::chrono::system_clock::now();
    std::chrono::duration<double> time_diff3 = end3 - start3;

    // Compute exact matrix-vector product
    Eigen::MatrixXd MG(n,n);
    auto start4 = std::chrono::system_clock::now();

    for(int i=0; i<n; ++i)
        for(int j=0; j<n; ++j)
            MG(i,j) = gskernel(PPoints[i].getX(), PPoints[i].getY(), PPoints[j].getX(), PPoints[j].getY());
    Eigen::VectorXd f_g_exact = MG * c;
    auto end4 = std::chrono::system_clock::now();
    std::chrono::duration<double> time_diff4 = end4 - start4;

    // Compute approximation error
    Eigen::VectorXd diff_g = f_g_exact - f_g_approx;

    // printing for testing
    std::cout << "f_g_exact   f_g_approx    diff_g" << std::endl;
    for(int i=0; i<n; i++){
        std::cout << f_g_exact(i) << "    " << f_g_approx(i) << "   " << diff_g(i) << std::endl;
    }
    std::cout << "Approximation error (l-inf norm): " << diff_g.lpNorm<Eigen::Infinity>() << std::endl
              << "Approximation error (l-2 norm): "   << diff_g.lpNorm<2>() << std::endl
              << "Time needed for exact multiplication: "       << time_diff4.count() << " s" << std::endl
              << "Time needed for approximate multiplication: " << time_diff3.count() << " s" << std::endl;


    std::cout << "Gauss" << std::endl;

    // Compute approximated matrix-vector product, given admissibility constant 'eta'

    GlobalInterpolationApp gip_gkernel(&gkernel, PPoints, n);
    Eigen::VectorXd f_gg_approx = gip_gkernel.mvProd(c,d);

    // Compute exact matrix-vector product

    Eigen::MatrixXd MGG(n,n);
    for(int i=0; i<n; ++i)
        for(int j=0; j<n; ++j)
            MGG(i,j) = gkernel(PPoints[i].getX(), PPoints[i].getY(), PPoints[j].getX(), PPoints[j].getY());
    Eigen::VectorXd f_gg_exact = MGG * c;

    // Compute approximation error
    Eigen::VectorXd diff_gg = f_gg_exact - f_gg_approx;

    // printing for testing
    std::cout << "f_gg_exact   f_gg_approx    diff_gg" << std::endl;
    for(int i=0; i<n; i++){
        std::cout << f_gg_exact(i) << "    " << f_gg_approx(i) << "   " << diff_gg(i) << std::endl;
    }
    std::cout << "Approximation error (l-inf norm): " << diff_gg.lpNorm<Eigen::Infinity>() << std::endl
              << "Approximation error (l-2 norm): "   << diff_gg.lpNorm<2>() << std::endl;

    }
}
