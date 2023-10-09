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
#include <Eigen/Dense>
#include <chrono>
#include <cmath>
#include <ctime>
#include <iostream>

#include "include/global_interpolation_app.hpp"
#include "include/kernel.hpp"
#include "include/low_rank_app.hpp"
#include "include/point.hpp"

// #define random1
#define circle

#define local
// #define global

int main() {
  // Input

  std::cout << "Enter gridsize:" << std::endl;
  unsigned n;
  std::cin >> n;
  //    unsigned n = 1000;

  // initialization of points
  std::vector<Point> points;
  points.reserve(n);
  std::srand(std::time(0));  // initializing points properties randomly
  double tx, ty;
  for (unsigned i = 0; i < n; ++i) {
    Point p;
    p.setId(i);
    p.setV(std::rand() % 100);  // values 0--n

#ifdef random1
    tx = ((double)rand() / (double)(RAND_MAX));
    ty = ((double)rand() / (double)(RAND_MAX));
    double t1 = ((double)rand() / (double)(RAND_MAX));
    if (t1 < 0.5) {
      tx = -tx;
    }
    double t2 = ((double)rand() / (double)(RAND_MAX));
    if (t2 < 0.5) {
      ty = -ty;
    }
#endif

#ifdef circle
    double angle = M_PI * 2 / n;
    tx = std::cos(angle * i);
    ty = std::sin(angle * i);
#endif

    p.setX(tx);
    p.setY(ty);
    points.push_back(p);
  }

  Eigen::VectorXd c = Eigen::VectorXd::Random(n);

  std::cout << "Enter admissibility constant:" << std::endl;
  double eta;
  std::cin >> eta;
  //    double eta = 0.5;

  std::cout << "Enter degree of interpolating polynomials:" << std::endl;
  unsigned q;
  std::cin >> q;
  //    unsigned q = 2;

  KernelGalerkin G(1.);  // initialization of Galerkin kernel for 2d problem
                         // -1/(2*pi)*log||x-y||

  // Compute exact matrix-vector product

  auto start1 = std::chrono::system_clock::now();

  Eigen::MatrixXd M(n, n);
  for (int i = 0; i < n; ++i)
    for (int j = 0; j < n; ++j)
      M(i, j) = G(points[i].getX(), points[i].getY(), points[j].getX(),
                  points[j].getY());
  Eigen::VectorXd f_exact = M * c;

  auto end1 = std::chrono::system_clock::now();
  std::chrono::duration<double> time_diff1 = end1 - start1;

#ifdef local
  // Compute approximated matrix-vector product, given admissibility constant
  // 'eta'

  auto start2 = std::chrono::system_clock::now();

  LowRankApp HMat(&G, points, eta,
                  q);  // initialization of low rank approximation for BEM
                       // approx for matrix multiplication
  Eigen::VectorXd f_approx =
      HMat.mvProd(c);  // calculation of the low rank approximation

  auto end2 = std::chrono::system_clock::now();
  std::chrono::duration<double> time_diff2 = end2 - start2;

  Eigen::VectorXd diff = f_exact - f_approx;

  // Compute approximation error
  std::cout << std::endl << "Local Interpolation" << std::endl;

  std::cout << "Number of matrix operations performed for exact matrix: "
            << n * n << std::endl;

  std::cout << "Approximation error between f_exact and f_approx (l-inf norm "
               "of vector diff): "
            << diff.lpNorm<Eigen::Infinity>() << std::endl
            << "Approximation error between f_exact and f_approx (l-2 norm of "
               "vector diff): "
            << diff.lpNorm<2>() << std::endl
            << "Relative Approximation error between f_exact and f_approx (l-2 "
               "norm of diff/l-2 norm of f_exact): "
            << diff.lpNorm<2>() / f_exact.lpNorm<2>() << std::endl
            << "Time needed for exact multiplication: " << time_diff1.count()
            << " s" << std::endl
            << "Time needed for approximate multiplication: "
            << time_diff2.count() << " s" << std::endl;
#endif

#ifdef global
  // Compute approximated matrix-vector product, given admissibility constant
  // 'eta'

  auto start3 = std::chrono::system_clock::now();

  GlobalInterpolationApp gip_gskernel(&G, points);
  Eigen::VectorXd f_g_approx = gip_gskernel.mvProd(c, q);

  auto end3 = std::chrono::system_clock::now();
  std::chrono::duration<double> time_diff3 = end3 - start3;

  Eigen::VectorXd diff_g = f_exact - f_g_approx;

  // Compute approximation error
  std::cout << std::endl << "Global Interpolation" << std::endl;

  std::cout << "Number of matrix operations performed for exact matrix: "
            << n * n << std::endl;

  std::cout << "Approximation error between f_exact and f_approx (l-inf norm "
               "of vector diff): "
            << diff_g.lpNorm<Eigen::Infinity>() << std::endl
            << "Approximation error between f_exact and f_approx (l-2 norm of "
               "vector diff): "
            << diff_g.lpNorm<2>() << std::endl
            << "Relative Approximation error between f_exact and f_approx (l-2 "
               "norm of diff/l-2 norm of f_exact): "
            << diff_g.lpNorm<2>() / f_exact.lpNorm<2>() << std::endl
            << "Time needed for exact multiplication: " << time_diff1.count()
            << " s" << std::endl
            << "Time needed for approximate multiplication: "
            << time_diff3.count() << " s" << std::endl;
#endif
}
