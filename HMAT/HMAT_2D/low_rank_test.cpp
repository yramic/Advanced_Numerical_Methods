#include "include/kernel.hpp"
#include "include/low_rank_app.hpp"
#include "include/global_interpolation_app.hpp"
#include "include/point.hpp"

#include <Eigen/Dense>
#include <cmath>
#include <ctime>
#include <iostream>
#include <chrono>

/*
 * Some of the definitions for the preprocessor must be deleted and also lots of comments that exist on the code.
 * This file will be cleaned from comments and not needed code in the final version because now everything is used for debugging the program.
 */


//#define circle
#define vector16
//#define random1


#define local
//#define global
//#define gsk
//#define gauss
int main() {

    // Input

    //std::cout << "Enter number of edges:" << std::endl;
    //unsigned n; std::cin >> n;

    // initializing n for testing
    unsigned n=20;
#ifdef vector16
    n = 16;
#endif
    //Eigen::VectorXd grid = Eigen::VectorXd::LinSpaced(n, 0., (n-1.)/n);
    std::vector<Point> PPoints; // initalizing Polygon Points properties
    PPoints.reserve(n);

    // x,y,v vextors used for testing(n=16)
    //std::vector<int> x = {14,77,11,38,81,54,7,48,17,30}, y = {62,45,46,87,37,20,72,21,30,90}, v = {59,7,23,71,34,50,48,54,30,60};
#ifdef vector16
    std::vector<int> x = {24,22,73,63,14,17,39,99,83,41,40,4,83,30,65,23}, y = {24,55,87,91,30,1,9,28,67,85,10,84,57,72,86,56}, v = {66,63,27,46,83,58,46,8,95,57,2,79,34,21,64,95};
#endif
    //std::vector<int> x = {18,89,39,39,26,53,4,44}, y = {63,49,66,27,46,4,54,59}, v = {1,2,3,4,5,6,7,8}; // vector with problem with setnearfar
    //std::vector<int> x = {80,27,17,80,23,86,99,89,70,43}, y = {65,74,12,88,44,56,77,31,41,88}, v = {1,2,3,4,5,6,7,8,9,10};
    //std::vector<int> x = {2,2,2,2,2,2,2,2,2,2}, y = {4,5,6,7,8,9,10,11,12,13}, v = {1,2,3,4,5,6,7,8,9,10};


    std::srand(std::time(0));   // initializing points properties randomly
    double tx,ty;
    int f=1;
    for (int i=0; i<n; i++){
        Point p;
        p.setId(i);
        p.setV(std::rand()%100);    // values 0-100

#ifdef random1
        tx = ((double)rand() / (double)(RAND_MAX));
        ty = ((double)rand() / (double)(RAND_MAX));
        double t1 = ((double)rand() / (double)(RAND_MAX));
        if(t1 < 0.5){
            tx = -tx;
        }
        double t2 = ((double)rand() / (double)(RAND_MAX));
        if(t2 < 0.5){
            ty = -ty;
        }
#endif

#ifdef circle
        double angle = M_PI*2/n;
        tx = std::cos(angle*i);
        ty = std::sin(angle*i);
#endif

        while (f) {             // checking if a point exists 2 times
            f=0;
            for (std::vector<Point>::iterator it=PPoints.begin(); it!=PPoints.end(); it++) {
                if (it->getX() == tx && it->getY() == ty) {
                    f = 1;
                }
            }
#ifdef random1
            tx = ((double)rand() / (double)(RAND_MAX));
            ty = ((double)rand() / (double)(RAND_MAX));
            t1 = ((double)rand() / (double)(RAND_MAX));
            if(t1 < 0.5){
                tx = -tx;
            }
            t2 = ((double)rand() / (double)(RAND_MAX));
            if(t2 < 0.5){
                ty = -ty;
            }
#endif

#ifdef circle
            double angle = M_PI*2/n;
            tx = std::cos(angle*i);
            ty = std::sin(angle*i);
#endif
        }
#if defined(random1) || defined(circle)
        p.setX(tx);
        p.setY(ty);
#endif


#ifdef vector16
        p.setX(x[i]);             // for testing
        p.setY(y[i]);             // for testing
        p.setV(v[i]);             // for testing
#endif
        PPoints.push_back(p);
    }
    for(std::vector<Point>::iterator iter = PPoints.begin(); iter< PPoints.end(); iter++) std::cout << (*iter).getX() << std::endl;
    std::cout << std::endl;
    for(std::vector<Point>::iterator iter = PPoints.begin(); iter< PPoints.end(); iter++) std::cout << (*iter).getY() << std::endl;
    std::cout << std::endl;

    Eigen::VectorXd    c(n);
    c= Eigen::VectorXd::Random(n);  // random initialization of vector c
    std::cout << c << std::endl << std::endl;
    // vector c initialization for testing
#ifdef vector16
    //Eigen::VectorXd c(n);
    c << 0.168182, -0.835255, 0.661581, -0.827464, 0.0385175, 0.559586, 0.745203, -0.624835, 0.273339, -0.327815, -0.528608, -0.320443, -0.731551, -0.600541, 0.474063, 0.98031;
#endif
    // printing the Polygon Points
    //for (std::vector<Point>::iterator it=PPoints.begin(); it!=PPoints.end(); it++)
    //    std::cout << "(" << it->getX() << "," << it->getY() << ") "; //<< it->getId() << ' ' << it->getV() << std::endl << std::flush;

    //std::cout << "Enter admissibility constant:" << std::endl;
    //double eta; std::cin >> eta;

    // initializing eta for testing
    double eta = 2;

    //std::cout << "Enter degree of interpolating polynomials:" << std::endl;
    //unsigned d; std::cin >> d;

    // initializing degree of interpolation for testing
    //unsigned d=1;

    //std::vector<unsigned> t = {1,2,3,5,7,10,20,40,80,100};
    std::vector<unsigned> t = {3};
    for(auto d : t){            // multiple runs of the program to check the behavior of the error for different degrees
    KernelGalerkin G;           // initialization of Galerkin kernel for 2d problem -1/(2*pi)*log||x-y||
    PolynomialKernel P;         // Polynomial Kernel initialization
    ConstantKernel CK(1);       // Constant Kernel initialization
    GlobalSmoothKernel gskernel;// Global Smooth Periodic Kernel initialization
    GaussKernel gkernel;        // Gauss Kernel initialization
    SingularKernel skernel;     // Singular Kernel initialization (log|x-y|)
    SingularKernelf skernelf;   // Singular Kernel initialization (log(1/|x-y|))

    // Compute exact matrix-vector product
#ifdef local
    auto start1 = std::chrono::system_clock::now();

    Eigen::MatrixXd M(n,n);
    for(int i=0; i<n; ++i)
        for(int j=0; j<n; ++j)
            M(i,j) = P(PPoints[i].getX(), PPoints[i].getY(), PPoints[j].getX(), PPoints[j].getY());
    Eigen::VectorXd f_exact = M * c;


    auto end1 = std::chrono::system_clock::now();
    std::chrono::duration<double> time_diff1 = end1 - start1;

    // Compute approximated matrix-vector product, given admissibility constant 'eta'

    auto start2 = std::chrono::system_clock::now();

    //LowRankApp lra(&P, PPoints, n);         // initialization of low rank approximation for BEM approx for matrix multiplication
    LowRankApp HMat(&P, PPoints, n, eta, d);
    Eigen::VectorXd f_approx = HMat.mvProd(c, eta, d);   // calculation of the low rank approximation

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

    std::cout << std::endl << "Number of points: " << n << std::endl;
    std::cout << "Degree: " << d << std::endl;

    // Compute approximation error
    std::cout << std::endl << "Local Interpolation" << std::endl;

    Eigen::VectorXd diff = f_exact - f_approx;

    // printing for testing
    std::cout << "f_exact" << std::endl;
    std::cout << f_exact << std::endl;

    std::cout << "f_approx" << std::endl;
    std::cout << f_approx << std::endl;

    std::cout << "diff" << std::endl;
    std::cout << diff << std::endl;
    std::cout << "Approximation error (l-inf norm): " << diff.lpNorm<Eigen::Infinity>() << std::endl
              << "Approximation error (l-2 norm): "   << diff.lpNorm<2>() << std::endl
              << "Relative Approximation error (l-2 norm): "   << diff.lpNorm<2>()/f_exact.lpNorm<2>() << std::endl
              << "Time needed for exact multiplication: "       << time_diff1.count() << " s" << std::endl
              << "Time needed for approximate multiplication: " << time_diff2.count() << " s" << std::endl;
#endif

#ifdef global
    // Global Interpolation
    std::cout << std::endl << "Global Interpolation" << std::endl;
    std::cout << "Number of points: "<< n << " Degree of interpolation: " << d << std::endl;
#ifdef gsk
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
              << "Relative Approximation error (l-2 norm): "   << diff_g.lpNorm<2>()/f_g_exact.lpNorm<2>() << std::endl
              << "Time needed for exact multiplication: "       << time_diff4.count() << " s" << std::endl
              << "Time needed for approximate multiplication: " << time_diff3.count() << " s" << std::endl;
#endif

#ifdef gauss
    std::cout << "Gauss" << std::endl;

    // Compute approximated matrix-vector product, given admissibility constant 'eta'
    auto start5 = std::chrono::system_clock::now();

    GlobalInterpolationApp gip_gkernel(&gkernel, PPoints, n);
    Eigen::VectorXd f_gg_approx = gip_gkernel.mvProd(c,d);

    auto end5 = std::chrono::system_clock::now();
    std::chrono::duration<double> time_diff5 = end5 - start5;
    // Compute exact matrix-vector product
    auto start6 = std::chrono::system_clock::now();

    Eigen::MatrixXd MGG(n,n);
    for(int i=0; i<n; ++i)
        for(int j=0; j<n; ++j)
            MGG(i,j) = gkernel(PPoints[i].getX(), PPoints[i].getY(), PPoints[j].getX(), PPoints[j].getY());
    Eigen::VectorXd f_gg_exact = MGG * c;

    auto end6 = std::chrono::system_clock::now();
    std::chrono::duration<double> time_diff6 = end6 - start6;
    // Compute approximation error
    Eigen::VectorXd diff_gg = f_gg_exact - f_gg_approx;

    // printing for testing
    std::cout << "f_gg_exact   f_gg_approx    diff_gg" << std::endl;
    for(int i=0; i<n; i++){
        std::cout << f_gg_exact(i) << "    " << f_gg_approx(i) << "   " << diff_gg(i) << std::endl;
    }
    std::cout << "Approximation error (l-inf norm): " << diff_gg.lpNorm<Eigen::Infinity>() << std::endl
              << "Approximation error (l-2 norm): "   << diff_gg.lpNorm<2>() << std::endl
              << "Relative Approximation error (l-2 norm): "   << diff_gg.lpNorm<2>()/f_gg_exact.lpNorm<2>() << std::endl
              << "Time needed for exact multiplication: "       << time_diff6.count() << " s" << std::endl
              << "Time needed for approximate multiplication: " << time_diff5.count() << " s" << std::endl;
#endif
#endif
    }
}
