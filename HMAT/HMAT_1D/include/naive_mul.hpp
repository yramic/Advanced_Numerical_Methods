#ifndef NAIVEMATVEC
#define NAIVEMATVEC

#include <Eigen/Dense>
#include <cassert>

inline Eigen::VectorXd prod(const Eigen::MatrixXd &M, const Eigen::VectorXd &c) {
  assert(M.cols()==c.rows());
  unsigned m = M.rows();
  unsigned n = M.cols();
  Eigen::VectorXd result(m);
  result.setZero();
  for (unsigned i = 0 ; i < m ; ++i) {
    for (unsigned j = 0 ; j < n ; ++j) {
      result(i) += M(i,j)*c(j);
    }
  }
  return result;
}

#endif
