/**
 * @file mfmg.h
 * @brief ADVNCSE homework MFMG code
 * @author Bob Schreiner, Jörg Nick
 * @date October 2023
 * @copyright Developed at SAM, ETH Zurich
 */

#ifndef MFMG_H_
#define MFMG_H_

#include <Eigen/Dense>
#include <iomanip>
#include <iostream>

namespace MFMG {
/** @brief Class implementing a multigrid method for a finite element
 * discretization of a scalar Dirichlet problem on the unit square
 */
/* SAM_LISTING_BEGIN_1 */
using GridFunction = Eigen::MatrixXd;
class DirichletBVPMultiGridSolver {
 public:
  DirichletBVPMultiGridSolver(double c, unsigned int L)
      : c_(c), L_(L) {}
  [[nodiscard]] GridFunction applyGridOperator(unsigned int level,
                                               const GridFunction &mu) const;
  [[nodiscard]] GridFunction directSolve(unsigned int level,
                                         const GridFunction &phi) const;
  void sweepGaussSeidel(unsigned int level, GridFunction &mu,
                        const GridFunction &phi) const;
  [[nodiscard]] GridFunction residual(unsigned int level,
                                      const GridFunction &mu,
                                      const GridFunction &phi) const;
  [[nodiscard]] GridFunction prolongate(unsigned int level,
                                        const GridFunction &gamma) const;
  [[nodiscard]] GridFunction restrict(unsigned int level,
                                      const GridFunction &rho) const;
  GridFunction multigridIteration(const GridFunction &mu,
                                  const GridFunction &phi,
                                  unsigned int L0 = 2) const;

 private:
  const double c_;        // Reaction coefficient
  const unsigned int L_;  // Level of the finest grid
};
/* SAM_LISTING_END_1 */

double estimateMGConvergenceRate(double c, unsigned int L, double tol = 1.0E-3);

void tabulateMGConvergenceRate();

}  // namespace MFMG

namespace MFMGLshape {
/** @brief Class implementing a multigrid method for a finite element
 * discretization of a scalar Dirichlet problem on the unit square
 */
/* SAM_LISTING_BEGIN_1 */
using GridFunction = Eigen::MatrixXd;
class DirichletBVPMultiGridSolver {
 public:
  DirichletBVPMultiGridSolver(double c, unsigned int L)
      : c_(c), L_(L){}
  [[nodiscard]] GridFunction applyGridOperator(unsigned int level,
                                               const GridFunction &mu) const;
  [[nodiscard]] GridFunction directSolve(unsigned int level,
                                         const GridFunction &phi) const;
  void sweepGaussSeidel(unsigned int level, GridFunction &mu,
                        const GridFunction &phi) const;
  [[nodiscard]] GridFunction residual(unsigned int level,
                                      const GridFunction &mu,
                                      const GridFunction &phi) const;
  [[nodiscard]] GridFunction prolongate(unsigned int level,
                                        const GridFunction &gamma) const;
  [[nodiscard]] GridFunction restrict(unsigned int level,
                                      const GridFunction &rho) const;
  GridFunction multigridIteration(const GridFunction &mu,
                                  const GridFunction &phi,
                                  unsigned int L0 = 2) const;

 private:
  const double c_;        // Reaction coefficient
  const unsigned int L_;  // Level of the finest grid
};
/* SAM_LISTING_END_1 */

double estimateMGConvergenceRate(double c, unsigned int L, double tol = 1.0E-3);

void tabulateMGConvergenceRate();

}  // namespace MFMGLshape

#endif  //MFMG_H_