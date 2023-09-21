/**
 * @file lowtriangtoeplitz.h
 * @brief ADVNCSE homework LowTriangToeplitz code
 * @author R. Hiptmair , Bob Schreiner
 * @date August 2023
 * @copyright Developed at SAM, ETH Zurich
 */
#ifndef LOWTRIANGTOEPLITZ_H_
#define LOWTRIANGTOEPLITZ_H_

#include <Eigen/Dense>
#include <cmath>
#include <iomanip>
#include <iostream>
#include <random>
#include <vector>

namespace LowTriangToeplitz {
Eigen::MatrixXcd toeplitz(const Eigen::VectorXcd& c, const Eigen::VectorXcd& r);
Eigen::VectorXcd pconvfft(const Eigen::VectorXcd& u, const Eigen::VectorXcd& x);
Eigen::VectorXcd toepMatVecMult(const Eigen::VectorXcd& c,
                                const Eigen::VectorXcd& r,
                                const Eigen::VectorXcd& x);
Eigen::VectorXcd ltpMult(const Eigen::VectorXcd& f, const Eigen::VectorXcd& g);
Eigen::VectorXcd ltpMultold(const Eigen::VectorXcd& f, const Eigen::VectorXcd& g);
Eigen::VectorXcd ltpSolve(const Eigen::VectorXcd& f, const Eigen::VectorXcd& y);
void test_accuracy_ltpMult();
void time_measure_ltpMult();
void test_accuracy_ltpSolve();
void time_measure_ltpSolve();

}  // namespace LowTriangToeplitz
#endif  // LOWTRIANGTOEPLITZ_H_