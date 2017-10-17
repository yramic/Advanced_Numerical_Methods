#include "include/kernel.hpp"
#include "include/low_rank_app.hpp"
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

    unsigned n=50;  // initalizing n for testing

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
        tx = std::rand()%1000;
        ty = std::rand()%1000;
        while (f) {             // checking if a point exists 2 times
            f=0;
            for (std::vector<Point>::iterator it=PPoints.begin(); it!=PPoints.end(); it++) {
                if (it->getX() == tx && it->getY() == ty) {
                    f = 1;
                }
            }
            tx = std::rand()%1000;
            ty = std::rand()%1000;
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

    std::cout << "c" << std::endl;
    std::cout << c << std::endl;

    // printing the Polygon Points
    //for (std::vector<Point>::iterator it=PPoints.begin(); it!=PPoints.end(); it++)
    //    std::cout << ' ' << it->getX() << ' ' << it->getY() << ' ' << it->getId() << ' ' << it->getV() << std::endl << std::flush;
    for (std::vector<Point>::iterator it=PPoints.begin(); it!=PPoints.end(); it++)
        std::cout << "(" << it->getX() << "," << it->getY() << ") ";
    //std::cout << "Enter admissibility constant:" << std::endl;
    //double eta; std::cin >> eta;

    // initializing eta for testing
    double eta = 2;

    //std::cout << "Enter degree of interpolating polynomials:" << std::endl;
    //unsigned d; std::cin >> d;

    // initializing degree of interpolation for testing
    unsigned d=3;
    std::cout << "degree: " << d << std::endl;
    //Kernel G(1.); // default kernel for 2d problem
    Kernel4D G;     // kernel for 4d problem -1/(2*pi)*log||x-y||
    PolynomialKernel P;
    // Compute exact matrix-vector product

    auto start1 = std::chrono::system_clock::now();

    Eigen::MatrixXd M(n,n);
    for(int i=0; i<n; ++i)
        for(int j=0; j<n; ++j)
            //M(i,j) = G(grid[i], grid[j]);
            //if(i!=j)
                M(i,j) = P(PPoints[i].getX(), PPoints[i].getY(), PPoints[j].getX(), PPoints[j].getY());    // the diagonal of the matrix will be zero because the 4D kernel between the same point is inf
    Eigen::VectorXd f_exact = M * c;


    auto end1 = std::chrono::system_clock::now();
    std::chrono::duration<double> time_diff1 = end1 - start1;

    // printing for debbuging
    std::cout << "M  " << std::endl;
    std::cout << M << std::endl;
    //std::cout << "f_exact " << std::endl;
    //std::cout << f_exact << std::endl;

    // Compute approximated matrix-vector product, given admissibility constant 'eta'

    auto start2 = std::chrono::system_clock::now();

    //LowRankApp lra(G, grid, grid);    // initialization of low rank approximation for 2D
    LowRankApp lra(P, PPoints, n);         // initialization of low rank approximation for BEM approx for matrix multiplication

    Eigen::VectorXd f_approx = lra.mvProd(c, eta, d);   // calculation of the low rank approximation

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

    auto end2 = std::chrono::system_clock::now();
    std::chrono::duration<double> time_diff2 = end2 - start2;

    // printing for testing
    //std::cout << "f_approx" << std::endl;
    //std::cout << f_approx << std::endl;

    // Compute approximation error

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
    std::cout << "enddeg" << std::endl;

    std::cout << "degree+1: " << d+2 << std::endl;
    LowRankApp lra2(P, PPoints, n);         // initialization of low rank approximation for BEM approx for matrix multiplication

    Eigen::VectorXd f_approx2 = lra2.mvProd(c, eta, d+2);   // calculation of the low rank approximation
     Eigen::VectorXd diff2 = f_exact - f_approx2;
    std::cout << "f_exact   f_approx    diff" << std::endl;
    for(int i=0; i<n; i++){
        std::cout << f_exact(i) << "    " << f_approx2(i) << "   " << diff2(i) << std::endl;
    }
    std::cout << "Approximation error (l-inf norm): " << diff2.lpNorm<Eigen::Infinity>() << std::endl
              << "Approximation error (l-2 norm): "   << diff2.lpNorm<2>() << std::endl
              << "Time needed for exact multiplication: "       << time_diff1.count() << " s" << std::endl
              << "Time needed for approximate multiplication: " << time_diff2.count() << " s" << std::endl;
    std::cout << "enddeg+2" << std::endl;
}
