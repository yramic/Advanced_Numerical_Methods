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
#include "include/point.hpp"
#include "include/uni-direct/block_cluster_Y.hpp"
#include "include/uni-direct/node_Y.hpp"
#include "include/naive_mul.hpp"

#include <Eigen/Dense>
#include <chrono>
#include <cmath>
#include <fstream>
#include <iomanip>
#include <iostream>


int main() {
  // Output filename
  std::string filename = "test_hmat_1d_error_uni.txt";
  std::ofstream myfile;
  myfile.open(filename);
  myfile << std::setw(15) << "#n" << std::setw(15) << "exact(us)"
         << std::setw(15) << "approx(us)" << std::endl;
  KernelCosine G(100.); // Kernel initialization
  for (unsigned N : {1,  2,  3,  4,  5,  6,  7,  8,  9,  10}) {
    // Number of collocation points
    unsigned n = std::pow(2,N) * 10;
    // Number of iterations for calculating the time
    unsigned iter = 100;
    double time_approx = 0., time_exact = 0.;
    // Collocation points
    Eigen::VectorXd grid = Eigen::VectorXd::LinSpaced(n, 0., 100.);
    // Vector to be multiplied
    Eigen::VectorXd c = Eigen::VectorXd::Random(n);
    std::vector<Point> GPoints; // initializing Grid Points properties
    GPoints.reserve(n);
    int k = 0;
    for (int i = 0; i < n; ++i) {
      Point p;
      p.setId(k);
      p.setX(grid[i]);
      k++;
      GPoints.push_back(p);
    }
    // Filling exact entries for matrix M
    Eigen::MatrixXd M(n, n);
    for (int i = 0; i < n; ++i)
      for (int j = 0; j < n; ++j)
        M(i, j) = G(GPoints[i].getX(), GPoints[j].getX());
    // Admissibility constant
    double eta = 0.1;
    // Degree of interpolation
    unsigned q = 2;
    // Constructor for Hierarchical matrix
    LowRankApp<BlockCluster_Y, Node_Y> HMat(&G, GPoints, eta, q/*, filename*/);
    // Iterations for timing
    for (unsigned i = 0; i < iter; ++i) {
      // Compute exact matrix-vector product using naive O(n^2) multiplication
      auto start1 = std::chrono::high_resolution_clock::now();
      Eigen::VectorXd f_exact = prod(M, c); // Calculating the product
      auto end1 = std::chrono::high_resolution_clock::now();
      auto time_diff1 =
          std::chrono::duration_cast<std::chrono::microseconds>(end1 - start1);
      time_exact += time_diff1.count();

      // Compute approximated matrix-vector product
      auto start2 = std::chrono::high_resolution_clock::now();
      Eigen::VectorXd f_approx = HMat.mvProd(c);
      auto end2 = std::chrono::high_resolution_clock::now();
      auto time_diff2 =
          std::chrono::duration_cast<std::chrono::microseconds>(end2 - start2);
      time_approx += time_diff2.count();;
      // std::cout << "Number of matrix operations performed for exact matrix: "
      // << n*n << std::endl;

      // Compute approximation error

      // Eigen::VectorXd diff = f_exact - f_approx;

      /*std::cout << "Approximation error (l-inf norm): " <<
         diff.lpNorm<Eigen::Infinity>() << std::endl
                << "Approximation error (l-2 norm): "   << diff.lpNorm<2>() <<
         std::endl
                << "Relative Approximation error (l-2 norm): "    <<
         diff.lpNorm<2>()/f_exact.lpNorm<2>() << std::endl
                << "Time needed for exact multiplication: "       <<
         time_diff1.count() << " s" << std::endl
                << "Time needed for approximate multiplication: " <<
         time_diff2.count() << " s" << std::endl;*/

      //    myfile << "time, " << n << ", " << std::setprecision(10) <<
      //    time_diff1.count() - time_diff2.count() << std::endl;

      // Alternative way to compute matrix error:
      // std::cout << "\nFor n = " << n << std::endl;
      /*Eigen::MatrixXd Mtilde(n,n);
      for(int i=0; i<n; ++i) {
          Mtilde.col(i) = HMat.mvProd(Eigen::VectorXd::Unit(n,i));
      }
      Eigen::MatrixXd diff_M = M - Mtilde;
      std::cout << "error_Frobenius, " << n << ", " << std::setprecision(10) <<
      diff_M.norm()/n << std::endl;
      std::cout << "error_max, "       << n << ", " << std::setprecision(10) <<
      diff_M.cwiseAbs().maxCoeff() << std::endl;*/
    }
    myfile << std::setw(15) << n << std::setw(15) << time_exact / iter
           << std::setw(15) << time_approx / iter << std::endl;
  }
}
