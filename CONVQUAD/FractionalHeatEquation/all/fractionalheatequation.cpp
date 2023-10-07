/**
 * @file fractionalheatequation.cpp
 * @brief NPDE homework FractionalHeatEquation code
 * @author Jörg Nick
 * @date October 2023
 * @copyright Developed at SAM, ETH Zurich
 */

#include "fractionalheatequation.h"

#include <memory>

namespace FractionalHeatEquation {

unsigned int SqrtsMplusA::solve_cnt{0};
unsigned int SqrtsMplusA::ludec_cnt{0};

SqrtsMplusA::SqrtsMplusA(unsigned int n, std::complex<double> s)
    : p_matrix_(std::make_unique<Eigen::SparseMatrix<std::complex<double>>>(
          n * n, n * n)) {
  (*p_matrix_).reserve(Eigen::RowVectorXi::Constant(n * n, 5));
  const double h = 1.0 / (1 + n);
  unsigned int rown = 0;
  for (int i = 0; i < n; ++i) {
    for (int j = 0; j < n; ++j, ++rown) {
      (*p_matrix_).insert(rown, rown) = 4.0 + h * h * std::sqrt(s);
      if (j > 0) (*p_matrix_).insert(rown, rown - 1) = -1.0;
      if (j < n - 1) (*p_matrix_).insert(rown, rown + 1) = -1.0;
      if (i > 0) (*p_matrix_).insert(rown, rown - n) = -1.0;
      if (i < n - 1) (*p_matrix_).insert(rown, rown + n) = -1.0;
    }
  }
}

/* SAM_LISTING_BEGIN_1 */
Eigen::VectorXd cqWeights(unsigned int M, double tau) {
  Eigen::VectorXd w(M + 1);

  return w;
}
/* SAM_LISTING_END_1 */

}  // namespace FractionalHeatEquation
