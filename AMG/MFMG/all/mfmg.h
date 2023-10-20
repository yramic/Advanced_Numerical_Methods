/**
 * @file mfmg.h
 * @brief NPDE homework MFMG code
 * @author
 * @date October 2023
 * @copyright Developed at SAM, ETH Zurich
 */

#include <Eigen/Dense>
#include <iomanip>
#include <iostream>

namespace MFMG {
/** @brief Class implemnting a multigrid method for finite-difference
 * discretization of a scalar Dirichlet problem on the unit square
 */
/* SAM_LISTING_BEGIN_1 */
using GridFunction = Eigen::MatrixXd;
class DirichletBVPMultiGridSolver {
 public:
  DirichletBVPMultiGridSolver(double c, unsigned int L) : c_(c), L_(L) {}
  [[nodiscard]] GridFunction applyGridOperator(unsigned int level,
                                               const GridFunction &) const;
  [[nodiscard]] GridFunction directSolve(unsigned int level, const GridFunction &phi) const;
  void sweepGaussSeidel(unsigned int level, GridFunction &mu,
                        const GridFunction &phi);
  [[nodiscard]] GridFunction residual(unsigned int level,
                                      const GridFunction &mu,
                                      const GridFunction &phi) const;
  [[nodiscard]] GridFunction prolongate(unsigned int level,
                                        const GridFunction &gamma) const;
  [[nodiscard]] GridFunction restrict(const GridFunction &rho) const;
  void multigridIteration(GridFunction &mu, const GridFunction &phi,
                          unsigned int L0 = 2) const;

 private:
  const double c_;        // Reaction coefficient
  const unsigned int L_;  // Level of finest grid
};
/* SAM_LISTING_END_1 */

double estimateMGConvergenceRate(double c, unsigned int L, double tol = 1.0E-3);

void tabulateMGConvergenceRate();

  
}  // namespace MFMG
