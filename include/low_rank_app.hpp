#ifndef LOW_RANK_APP_HPP
#define LOW_RANK_APP_HPP

#include <Eigen/Dense>
#include <vector>
#include "BC.hpp"
#include "ctree.hpp"
#include "node.hpp"

// Approximate matrix-vector multiplication

Eigen::VectorXd mvProd(const std::vector<double>& x,
                       const std::vector<double>& y,
                       double admis_const, double kernel_const, int deg, const Eigen::VectorXd& c_);

#endif // LOW_RANK_APP_HPP
