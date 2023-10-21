/**
 * @file mfmgLshape.cpp
 * @brief NPDE homework MFMG code
 * @author
 * @date October 2023
 * @copyright Developed at SAM, ETH Zurich
 */

#include "mfmg.h"

namespace MFMG {

/* SAM_LISTING_BEGIN_1 */
GridFunction DirichletBVPMultiGridSolver::applyGridOperator(
    unsigned int level, const GridFunction &) const {
#if SOLUTION

#else
  // ************************************************************
  // Code for L-shaped case to be supplemented
  // ************************************************************
#endif
}
/* SAM_LISTING_END_1 */

/* SAM_LISTING_BEGIN_2 */
GridFunction DirichletBVPMultiGridSolver::directSolve(
    unsigned int level, const GridFunction &phi) const {
#if SOLUTION

#else
  // ************************************************************
  // Code for L-shaped case to be supplemented
  // ************************************************************
#endif
}
/* SAM_LISTING_END_2 */

/* SAM_LISTING_BEGIN_3 */
void DirichletBVPMultiGridSolver::sweepGaussSeidel(unsigned int level,
                                                   GridFunction &mu,
                                                   const GridFunction &phi) {
#if SOLUTION

#else
  // ************************************************************
  // Code for L-shaped case to be supplemented
  // ************************************************************
#endif
}
/* SAM_LISTING_END_3 */

/* SAM_LISTING_BEGIN_4 */
GridFunction DirichletBVPMultiGridSolver::residual(
    unsigned int level, const GridFunction &mu, const GridFunction &phi) const {
#if SOLUTION

#else
  // ************************************************************
  // Code for L-shaped case to be supplemented
  // ************************************************************
#endif
}
/* SAM_LISTING_END_4 */

/* SAM_LISTING_BEGIN_5 */
GridFunction DirichletBVPMultiGridSolver::prolongate(
    unsigned int level, const GridFunction &gamma) const {
#if SOLUTION

#else
  // ************************************************************
  // Code for L-shaped case to be supplemented
  // ************************************************************
#endif
}
/* SAM_LISTING_END_5 */

/* SAM_LISTING_BEGIN_6 */
GridFunction DirichletBVPMultiGridSolver::restrict(
    const GridFunction &rho) const {
#if SOLUTION

#else
  // ************************************************************
  // Code for L-shaped case to be supplemented
  // ************************************************************
#endif
}
/* SAM_LISTING_END_6 */

/* SAM_LISTING_BEGIN_7 */
void DirichletBVPMultiGridSolver::multigridIteration(GridFunction &mu,
                                                     const GridFunction &phi,
                                                     unsigned int L0) const {
#if SOLUTION

#else
  // ************************************************************
  // Code for L-shaped case to be supplemented
  // ************************************************************
#endif
}
/* SAM_LISTING_END_7 */

/* SAM_LISTING_BEGIN_8 */
double estimateMGConvergenceRate(double c, unsigned int L, double tol) {
#if SOLUTION

#else
  // ************************************************************
  // Code for L-shaped case to be supplemented
  // ************************************************************
#endif
}
/* SAM_LISTING_END_8 */

/* SAM_LISTING_BEGIN_9 */
void tabulateMGConvergenceRate() {
#if SOLUTION

#else
  // ************************************************************
  // Code for L-shaped case to be supplemented
  // ************************************************************
#endif
}
/* SAM_LISTING_END_9 */

}  // namespace MFMG
