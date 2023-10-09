#ifndef CHEBYSHEVGAUSSQUAD
#define CHEBYSHEVGAUSSQUAD

#include <Eigen/Dense>
#include <cmath>
#include <utility>
#define USE_MATH
inline std::pair<Eigen::RowVectorXd, double> ChebyshevGaussQuad(
    unsigned int N) {
  Eigen::VectorXd nodes(N);
  double weight = M_PI / N;

  for (unsigned int i = 1; i < N + 1; ++i)
    nodes(i - 1) = cos((2 * i - 1) * M_PI / (2 * N));
  return std::make_pair(nodes, weight);
}

#endif
