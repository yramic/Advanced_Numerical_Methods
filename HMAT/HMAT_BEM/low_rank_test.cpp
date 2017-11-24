/***********************************************************************
 *                                                                     *
 * Code for Course "Advanced Numerical Methods for CSE"                *
 * (Prof. Dr. R. Hiptmair)                                             *
 * Author: Daniele Casati                                              *
 * Date: 11/2017                                                       *
 * (C) Seminar for Applied Mathematics, ETH Zurich                     *
 * This code can be freely used for non-commercial purposes as long    *
 * as this header is left intact.                                      *
 ***********************************************************************/
#include "include/kernel.hpp"
#include "include/low_rank_app.hpp"
#include "include/segment.hpp"
extern "C" {
#include "../BEM/CppHilbert/Library/source/gaussQuadrature.h"
}

#include <Eigen/Dense>
#include <chrono>
#include <cmath>
#include <ctime>
#include <iostream>


int main() {

    // Test

    std::vector<Segment> segments;
    {
        Segment s;
        s.setId(0);
        Eigen::Vector2d a, b;
        a << 0, 1;
        b << 0, 0;
        s.setA(a);
        s.setB(b);
        segments.push_back(s);
    }
    {
        Segment s;
        s.setId(1);
        Eigen::Vector2d a, b;
        a << 0, 0;
        b << 1, 0;
        s.setA(a);
        s.setB(b);
        segments.push_back(s);
    }
    {
        Segment s;
        s.setId(2);
        Eigen::Vector2d a, b;
        a << 1, 0;
        b << 1, 1;
        s.setA(a);
        s.setB(b);
        segments.push_back(s);
    }
    {
        Segment s;
        s.setId(3);
        Eigen::Vector2d a, b;
        a << 1, 1;
        b << 0, 1;
        s.setA(a);
        s.setB(b);
        segments.push_back(s);
    }
    unsigned n = segments.size();
    Node node(segments, 10);
    node.getRect();
    node.setV();

    KernelGalerkin G; // initialization of Galerkin kernel for 2d problem -1/(2*pi)*log||x-y||

    BlockCluster bc(&node, &node);
    bc.setMatrix(&G);
    std::cout << "VCVt:" << std::endl << bc.getVCV() << std::endl;

    Eigen::MatrixXd M(n,n);
    for(int i=0; i<n; ++i)
        for(int j=0; j<n; ++j)
            M(i,j) = G(segments[i].getA(), segments[i].getB(), segments[j].getA(), segments[j].getB());
    std::cout << "CppHilbert:" << std::endl << M << std::endl;

    // Input
/*
    std::cout << "Enter gridsize:" << std::endl;
    unsigned n; std::cin >> n;
//    unsigned n = 1000;

    // initialization of segments on a square
    unsigned n = 2; // segments per edge
    std::vector<Segment> segments;
    segments.reserve(n);
    std::srand(std::time(0)); // initializing points properties randomly
    double progress = 0.; double shift = 1./n;
    for(unsigned i=0; i<n; ++i) {
        {
            Segment s;
            s.setId(i*4);
            Eigen::Vector2d a, b;
            a << progress, 0.;
            b << progress+shift, 0.;
            s.setA(a);
            s.setB(b);
            segments.push_back(s);
        }
        {
            Segment s;
            s.setId(i*4+1);
            Eigen::Vector2d a, b;
            a << 0., progress;
            b << 0., progress+shift;
            s.setA(a);
            s.setB(b);
            segments.push_back(s);
        }
        {
            Segment s;
            s.setId(i*4+2);
            Eigen::Vector2d a, b;
            a << progress, 1.;
            b << progress+shift, 1.;
            s.setA(a);
            s.setB(b);
            segments.push_back(s);
        }
        {
            Segment s;
            s.setId(i*4+3);
            Eigen::Vector2d a, b;
            a << 1., progress;
            b << 1., progress+shift;
            s.setA(a);
            s.setB(b);
            segments.push_back(s);
        }
        progress += 1./n;
    }
    n = segments.size();
*/
/*
    // initialization of segments on a circle
    std::vector<Segment> segments;
    segments.reserve(n);
    std::srand(std::time(0)); // initializing points properties randomly
    for(unsigned i=0; i<n; ++i) {
        Segment s;
        s.setId(i);

        double angle = M_PI*2/n;
        Eigen::Vector2d a, b;
        a << std::cos(angle* i),     std::sin(angle* i);
        b << std::cos(angle*(i+1.)), std::sin(angle*(i+1.));

        s.setA(a);
        s.setB(b);
        segments.push_back(s);
    }

    Eigen::VectorXd c = Eigen::VectorXd::Random(n);

    std::cout << "Enter admissibility constant:" << std::endl;
    double eta; std::cin >> eta;
//    double eta = 0.5;

    std::cout << "Enter degree of interpolating polynomials:" << std::endl;
    unsigned q; std::cin >> q;
//    unsigned q = 2;

    KernelGalerkin G; // initialization of Galerkin kernel for 2d problem -1/(2*pi)*log||x-y||

    // Compute exact matrix-vector product

    auto start1 = std::chrono::system_clock::now();

    Eigen::MatrixXd M(n,n);
    for(int i=0; i<n; ++i)
        for(int j=0; j<n; ++j)
            M(i,j) = G(segments[i].getA(), segments[i].getB(), segments[j].getA(), segments[j].getB());
    Eigen::VectorXd f_exact = M * c;

    auto end1 = std::chrono::system_clock::now();
    std::chrono::duration<double> time_diff1 = end1 - start1;

    // Compute approximated matrix-vector product, given admissibility constant 'eta'

    auto start2 = std::chrono::system_clock::now();

    LowRankApp HMat(&G, segments, eta, q); // initialization of low rank approximation for BEM approx for matrix multiplication
    Eigen::VectorXd f_approx = HMat.mvProd(c); // calculation of the low rank approximation

    auto end2 = std::chrono::system_clock::now();
    std::chrono::duration<double> time_diff2 = end2 - start2;

    Eigen::VectorXd diff = f_exact - f_approx;

    // Compute approximation error

    std::cout << "Number of matrix operations performed for exact matrix: " << n*n << std::endl;

    std::cout << "Approximation error between f_exact and f_approx (l-inf norm of vector diff): " << diff.lpNorm<Eigen::Infinity>() << std::endl
              << "Approximation error between f_exact and f_approx (l-2 norm of vector diff): "   << diff.lpNorm<2>() << std::endl
              << "Relative Approximation error between f_exact and f_approx (l-2 norm of diff/l-2 norm of f_exact): " << diff.lpNorm<2>()/f_exact.lpNorm<2>() << std::endl
              << "Time needed for exact multiplication: "       << time_diff1.count() << " s" << std::endl
              << "Time needed for approximate multiplication: " << time_diff2.count() << " s" << std::endl;
*/
}
