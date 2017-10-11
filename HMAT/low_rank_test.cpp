#include "include/kernel.hpp"
#include "include/low_rank_app.hpp"
#include "include/point.hpp"

#include <Eigen/Dense>
#include <cmath>
#include <ctime>
#include <iostream>

int main() {

    // Input

    //std::cout << "Enter number of edges:" << std::endl;
    //unsigned n; std::cin >> n;
    //test
    unsigned n=16;
    //Eigen::VectorXd grid = Eigen::VectorXd::LinSpaced(n, 0., (n-1.)/n);
    std::vector<Point> PPoints; // initalizing Polygon Points properties
    PPoints.reserve(n);
    //std::vector<int> x = {14,77,11,38,81,54,7,48,17,30}, y = {62,45,46,87,37,20,72,21,30,90}, v = {59,7,23,71,34,50,48,54,30,60};
    std::vector<int> x = {24,22,73,63,14,17,39,99,83,41,40,4,83,30,65,23}, y = {24,55,87,91,30,1,9,28,67,85,10,84,57,72,86,56}, v = {66,63,27,46,83,58,46,8,95,57,2,79,34,21,64,95};

    std::srand(std::time(0));
    int tx,ty;
    int f=1;
    for (int i=0; i<n; i++){
        Point p;
        p.setId(i);
        /*p.setV(std::rand()%100);
        tx = std::rand()%100;
        ty = std::rand()%100;
        while (f) {
            f=0;
            for (std::vector<Point>::iterator it=PPoints.begin(); it!=PPoints.end(); it++) {
                if (it->getX() == tx && it->getY() == ty) {
                    f = 1;
                }
            }
            tx = std::rand()%100;
            ty = std::rand()%100;
        }
        //p.setX(std::rand()%n);
        //p.setY(std::rand()%n);
        p.setX(tx);
        p.setY(ty);*/
        p.setX(x[i]);
        p.setY(y[i]);
        p.setV(v[i]);
        PPoints.push_back(p);
    }
    Eigen::VectorXd    c = Eigen::VectorXd::Random(n);
    // test out
    for (std::vector<Point>::iterator it=PPoints.begin(); it!=PPoints.end(); it++)
        std::cout << ' ' << it->getX() << ' ' << it->getY() << ' ' << it->getId() << ' ' << it->getV() << std::endl << std::flush;

    //std::cout << "Enter admissibility constant:" << std::endl;
    //double eta; std::cin >> eta;

    // debug
    double eta = 2;

    //std::cout << "Enter degree of interpolating polynomials:" << std::endl;
    //unsigned d; std::cin >> d;
    // debug
    unsigned d=2;

    Kernel G(1.);

    // Compute exact matrix-vector product

    /*time_t start1; time(&start1);

    Eigen::MatrixXd M(n,n);
    for(int i=0; i<n; ++i)
        for(int j=0; j<n; ++j)
            M(i,j) = G(grid[i], grid[j]);
    Eigen::VectorXd f_exact = M * c;

    time_t end1; time(&end1);
    double time_diff1 = std::difftime(end1, start1);
    */
    // Compute approximated matrix-vector product, given admissibility constant 'eta'

    time_t start2; time(&start2);

    //LowRankApp lra(G, grid, grid);
    LowRankApp lra(G, PPoints);

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

    /*time_t end2; time(&end2);
    double time_diff2 = std::difftime(end2, start2);

    // Compute approximation error

    Eigen::VectorXd diff = f_exact - f_approx;

    std::cout << "Approximation error (l-inf norm): " << diff.lpNorm<Eigen::Infinity>() << std::endl
              << "Approximation error (l-2 norm): "   << diff.lpNorm<2>() << std::endl
              << "Time needed for exact multiplication: " << time_diff1 << std::endl
              << "Time needed for approximate multiplication: " << time_diff2 << std::endl;*/
}
