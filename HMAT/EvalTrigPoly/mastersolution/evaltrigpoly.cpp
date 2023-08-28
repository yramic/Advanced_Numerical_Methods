/**
 * @file evaltrigpoly.cpp
 * @brief NPDE homework EvalTrigPoly code
 * @author Ralf Hiptmair
 * @date August 2023
 * @copyright Developed at SAM, ETH Zurich
 */

#include "evaltrigpoly.h"

namespace EvalTrigPoly {
/* SAM_LISTING_BEGIN_1 */
Eigen::VectorXcd evalTrigPolyEquid(const Eigen::VectorXcd &coeffs) {
  // **********************************************************************
  // TODO
  // **********************************************************************
}
/* SAM_LISTING_END_1 */

/* SAM_LISTING_BEGIN_2 */
Eigen::VectorXcd evalTrigPoly(const Eigen::VectorXcd &gamma,
                              const Eigen::VectorXd &x) {
  // **********************************************************************
  // TODO
  // **********************************************************************
}
/* SAM_LISTING_END_2 */

/* SAM_LISTING_BEGIN_3 */
std::vector<unsigned int> checkX(const Eigen::VectorXcd &gamma, ,
                                 const Eigen::VectorXd &x, double tau) {
  // **********************************************************************
  // TODO
  // **********************************************************************
}
/* SAM_LISTING_END_3 */

/* SAM_LISTING_BEGIN_4 */
Eigen::VectorXcd evalTrigPolyApprox(const Eigen::VectorXcd &gamma,
                                    const Eigen::VectorXd &x, unsigned int q) {
  // **********************************************************************
  // TODO
  // **********************************************************************
}
/* SAM_LISTING_END_4 */

}  // namespace EvalTrigPoly
