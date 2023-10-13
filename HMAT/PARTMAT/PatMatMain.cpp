/***********************************************************************
 *                                                                     *
 * Code for Course "Advanced Numerical Methods for CSE"                *
 * (Prof. Dr. R. Hiptmair)                                             *
 * Author: R.H.                                                        *
 * Date: Nov 7, 2017                                                   *
 * (C) Seminar for Applied Mathematics, ETH Zurich                     *
 * This code can be freely used for non-commercial purposes as long    *
 * as this header is left intact.                                      *
 ***********************************************************************/

// General includes
#include <Eigen/Dense>
#include <cmath>
#include <exception>
#include <iostream>
#include <vector>

namespace HMAT {

/* @brief Data structure for box of matrix partition
   Describes a matrix block through two index sets.
   The block itself is stored in rank-q factorized form.
*/

/* SAM_LISTING_BEGIN_0 */
// Rank-q matrix block in factorized form
template <int q>
struct FarFieldBlock {
  const std::vector<int> i_idx, j_idx;            // contained indices
  Eigen::Matrix<double, Eigen::Dynamic, q> U, V;  // low-rank factors
};

// Submatrix; no special structure assumed
struct NearFieldBlock {
  const std::vector<int> i_idx, j_idx;  // contained indices
  Eigen::MatrixXd Mloc;                 // matrix block
};
/* SAM_LISTING_END_0 */

/* SAM_LISTING_BEGIN_1 */
template <int q>
class PartMatrix {
 public:
  PartMatrix(size_t _n, size_t _m);
  // Matrix$\times$vector operation
  Eigen::VectorXd operator*(const Eigen::VectorXd &v) const;

 private:
  size_t m, n;  // dimensions of matrix
  std::vector<FarFieldBlock<q>> farField;
  std::vector<NearFieldBlock> nearField;
};
/* SAM_LISTING_END_1 */

template <int q>
PartMatrix<q>::PartMatrix(size_t _n, size_t _m) : n(_n), m(_m) {}

/* SAM_LISTING_BEGIN_2 */
// Partitioned $n\times m$-matrix split in near-field and
// far-field blocks, the latter of rank q
template <int q>
Eigen::VectorXd PartMatrix<q>::operator*(const Eigen::VectorXd &v) const {
  if (v.size() != m) throw(std::runtime_error("Size mismatch in *"));
  Eigen::VectorXd y(n);
  y.setZero();  // std::vector for returning result
  // Traverse far field boxes
  for (const FarFieldBlock<q> &B : farField) {
    // Get no. of x and y collocation points in box
    const size_t nB = B.i_idx.size();
    const size_t mB = B.j_idx.size();
    // Obtain values of argument std::vector corresponding to y-points
    Eigen::VectorXd tmp(mB);
    for (int j = 0; j < mB; j++) tmp(j) = v(B.j_idx[j]);
    // Multiply std::vector with low-rank matrix: Effort \cob{$\sharp I_k+\sharp J_k$}
    Eigen::VectorXd res(nB);
    res = B.U * (B.V.transpose() * tmp);
    // Accumlate result into components of result std::vector
    for (int i = 0; i < nB; i++) y(B.i_idx[i]) += res(i);
  }
  // Traverse near field boxes
  for (const NearFieldBlock &B : nearField) {
    // Get no. of x and y collocation points in box
    const size_t nB = B.i_idx.size();
    const size_t mB = B.j_idx.size();
    // Obtain values of argument std::vector corresponding to y-points
    Eigen::VectorXd tmp(mB);
    for (int j = 0; j < mB; j++) tmp(j) = v(B.j_idx[j]);
    // Multiply std::vector with local collocation matrix
    Eigen::VectorXd res(nB);
    res = B.Mloc * tmp;
    // Accumlate result into components of result std::vector
    for (int i = 0; i < nB; i++) y(B.i_idx[i]) += res(i);
  }
  return (y);  // (Move) return result vector
}
/* SAM_LISTING_END_2 */
}  // namespace HMAT

int main(int, char **) {
  HMAT::PartMatrix<7> PM(10, 10);
  Eigen::VectorXd v = Eigen::VectorXd::Random(10);
  Eigen::VectorXd y = PM * v;
  return (0);
}
