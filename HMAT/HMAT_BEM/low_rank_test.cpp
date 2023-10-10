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
  // Input

  //    std::cout << "Enter gridsize:" << std::endl;
  //    unsigned n; std::cin >> n;
  //    unsigned n = 1000;

  //    std::cout << "Enter admissibility constant:" << std::endl;
  //    double eta; std::cin >> eta;
  //    double eta = 0.5;

  std::cout << "Enter degree of interpolating polynomials:" << std::endl;
  unsigned q;
  std::cin >> q;
  //    unsigned q = 2;

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
  Node node(segments, q);
  node.getRect();
  node.setV();

  KernelGalerkin G;  // initialization of Galerkin kernel for 2d problem
                     // -1/(2*pi)*log||x-y||

  BlockCluster bc(&node, &node);
  bc.setMatrix(&G);
  Eigen::MatrixXd VCV = bc.getVCV();
  Eigen::MatrixXd VCVord(VCV.rows(), VCV.cols());
  for (unsigned i = 0; i < bc.getXNode()->getSegments().size(); ++i)
    for (unsigned j = 0; j < bc.getYNode()->getSegments().size(); ++j)
      VCVord(bc.getXNode()->getSegments()[i].getId(),
             bc.getYNode()->getSegments()[j].getId()) = VCV(i, j);
  std::cout << "VCVt:" << std::endl << VCVord << std::endl;

  Eigen::MatrixXd M(n, n);
  for (int i = 0; i < n; ++i)
    for (int j = 0; j < n; ++j)
      M(i, j) = G(segments[i].getA(), segments[i].getB(), segments[j].getA(),
                  segments[j].getB());
  std::cout << "CppHilbert:" << std::endl << M << std::endl;

  /*
      // initialization of segments on a square
      std::vector<Segment> segments;
      segments.reserve(4*n); // 'n' segments per edge
      std::srand(std::time(0)); // initializing segment properties randomly
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
  */
}
