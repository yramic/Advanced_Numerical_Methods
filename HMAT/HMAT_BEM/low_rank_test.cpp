/***********************************************************************
 *                                                                     *
 * Code for Course "Advanced Numerical Methods for CSE"                *
 * (Prof. Dr. R. Hiptmair)                                             *
 * Author: Ioannis Magkanaris                                          *
 * Date: 11/2017                                                       *
 * (C) Seminar for Applied Mathematics, ETH Zurich                     *
 * This code can be freely used for non-commercial purposes as long    *
 * as this header is left intact.                                      *
 ***********************************************************************/
#include "include/kernel.hpp"
#include "include/low_rank_app.hpp"
#include "include/segment.hpp"

#include <Eigen/Dense>
#include <chrono>
#include <cmath>
#include <ctime>
#include <iostream>


int main() {

    // Input

    std::cout << "Enter gridsize:" << std::endl;
    unsigned n; std::cin >> n;
//    unsigned n = 1000;

    // initialization of segments
    std::vector<Segment> segments;
    segments.reserve(n);
    std::srand(std::time(0)); // initializing points properties randomly
    for(unsigned i=1; i<n; ++i) {
        Segment s;
        s.setId(i);
        s.setV(std::rand()%100); // values 0--n

        double angle = M_PI*2/n;
        Eigen::Vector2d a, b;
        a << std::cos(angle*(i-1.)), std::sin(angle*(i-1.));
        b << std::cos(angle* i),     std::sin(angle* i);

        s.setA(a);
        s.setB(b);
        segments.push_back(s);
    }

    std::cout << "Enter admissibility constant:" << std::endl;
    double eta; std::cin >> eta;
//    double eta = 0.5;

    std::cout << "Enter degree of interpolating polynomials:" << std::endl;
    unsigned q; std::cin >> q;
//    unsigned q = 2;

    KernelGalerkin G; // initialization of Galerkin kernel for 2d problem -1/(2*pi)*log||x-y||

    // Compute exact matrix

    Eigen::MatrixXd M(n,n);
    for(int i=0; i<n; ++i)
        for(int j=0; j<n; ++j)
            M(i,j) = G(segments[i].getA(), segments[i].getB(), segments[j].getA(), segments[j].getB());

    // Compute approximated matrix

    Eigen::MatrixXd Mtilde(n,n);
    for(int i=0; i<n; ++i) {
        LowRankApp HMat_tmp(&G, segments, eta, q); // initialization of low rank approximation for BEM approx for matrix multiplication
        Mtilde.col(i) = HMat_tmp.mvProd(Eigen::VectorXd::Unit(n,i));
    }

    // Compute approximation error

    Eigen::MatrixXd diff_M = M - Mtilde;

    std::cout << "Approximation error in Frobenius norm: " << diff_M.norm()/n << std::endl
              << "Approximation error in max norm: "       << diff_M.cwiseAbs().maxCoeff() << std::endl;
}
