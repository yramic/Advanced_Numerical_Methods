/**
 * @file cfsplit.cpp
 * @brief NPDE homework CFSplit code
 * @author  R. Hiptmair
 * @date November 2023
 * @copyright Developed at SAM, ETH Zurich
 */

#include "cfsplit.h"

#include <filesystem>
#include <fstream>
#include <stdexcept>

namespace GalerkinConstruction {
// Copied from homework project "GalerkinConstruction"
Eigen::SparseMatrix<double> buildAH(
    const Eigen::SparseMatrix<double> &A,
    const Eigen::SparseMatrix<double, Eigen::RowMajor> &P) {
  Eigen::SparseMatrix<double> AH(P.cols(), P.cols());
  // define vector of triplets and reserve memory
  std::vector<Eigen::Triplet<double>> AH_trp{};  // For temporary COO format
  // For both A and P the inner dimensions are rows
  // loop over rows of A
  for (int k = 0; k < A.outerSize(); ++k) {
    for (Eigen::SparseMatrix<double>::InnerIterator l(A, k); l; ++l) {
      // loop over rows of P
      for (Eigen::SparseMatrix<double, Eigen::RowMajor>::InnerIterator j(
               P, l.row());
           j; ++j) {
        for (Eigen::SparseMatrix<double, Eigen::RowMajor>::InnerIterator i(P,
                                                                           k);
             i; ++i) {
          AH_trp.emplace_back(i.col(), j.col(),
                              j.value() * i.value() * l.value());
        }
      }
    }
  }
  // Initialize CCS matrix from COO format
  AH.setFromTriplets(AH_trp.begin(), AH_trp.end());
  AH.makeCompressed();
  return AH;
}
}  // namespace GalerkinConstruction

namespace CFSplit {

Eigen::SparseMatrix<double> readMMatrixFromFile(std::string filename) {
  // Default choices resides in the same directory as this source file
  const std::filesystem::path here = __FILE__;
  auto meshfilepath = here.parent_path() / filename;

  FILE *f = fopen(meshfilepath.string().c_str(), "r");
  if (f == NULL) {
    throw std::runtime_error("Cannot opem matrix file");
  }
  MM_typecode matcode;
  if (mm_read_banner(f, &matcode) != 0) {
    throw std::runtime_error("Could not process Matrix Market banner!");
  }
  std::cout << "Matrix market type = " << mm_typecode_to_str(matcode)
            << std::endl;
  if (!mm_is_sparse(matcode)) {
    throw std::runtime_error("Matrix must be sparse!");
  }
  if (mm_is_complex(matcode)) {
    throw std::runtime_error("Complex matrices not supported!");
  }
  int ret_code;
  int M, N, nz;
  if ((ret_code = mm_read_mtx_crd_size(f, &M, &N, &nz)) != 0) {
    throw std::runtime_error("Cannot read matrix size");
  }
  std::vector<Eigen::Triplet<double>> trp_vec;
  for (int i = 0; i < nz; i++) {
    int I, J;
    double val;
    fscanf(f, "%d %d %lg\n", &I, &J, &val);
    trp_vec.emplace_back(I, J, val);
  }
  Eigen::SparseMatrix<double> A(M, N);
  A.setFromTriplets(trp_vec.begin(), trp_vec.end());
  A.makeCompressed();
  return A;
}

/* SAM_LISTING_BEGIN_2 */
std::pair<std::vector<std::vector<RugeStuebenAMG::nodeidx>>,
          std::vector<std::vector<RugeStuebenAMG::nodeidx>>>
RugeStuebenAMG::getStrongConn(const Eigen::SparseMatrix<double> &A,
                              double tauC) const {
#if SOLUTION

#else
  // ************************************************************
  // Code to be supplemented
  // ************************************************************
#endif
}
/* SAM_LISTING_END_2 */

/* SAM_LISTING_BEGIN_3 */
void RugeStuebenAMG::setSmoothable(const Eigen::SparseMatrix<double> &A,
                                   double sigma,
                                   std::vector<NodeFlag> &nodeflags) const {
#if SOLUTION

#else
  // ************************************************************
  // Code to be supplemented
  // ************************************************************
#endif
}
/* SAM_LISTING_END_3 */

/* SAM_LISTING_BEGIN_4 */
std::vector<unsigned int> RugeStuebenAMG::connectivityMeasure(
    const std::vector<NodeFlag> &nodeflags,
    const std::vector<std::vector<nodeidx>> &Sj_star) const {
#if SOLUTION

#else
  // ************************************************************
  // Code to be supplemented
  // ************************************************************
#endif
}
/* SAM_LISTING_END_4 */

/* SAM_LISTING_BEGIN_5 */
std::vector<RugeStuebenAMG::NodeFlag> RugeStuebenAMG::CFsplit(
    const Eigen::SparseMatrix<double> &A,
    const std::vector<std::vector<nodeidx>> &Sj_star) const {
#if SOLUTION

#else
  // ************************************************************
  // Code to be supplemented
  // ************************************************************
#endif
}
/* SAM_LISTING_END_5 */

/* SAM_LISTING_BEGIN_6 */
Eigen::SparseMatrix<double, Eigen::RowMajor> RugeStuebenAMG::computeP(
    const Eigen::SparseMatrix<double> &A,
    const std::vector<std::vector<nodeidx>> &Si,
    const std::vector<NodeFlag> &nodeflags) const {
#if SOLUTION

#else
  // ************************************************************
  // Code to be supplemented
  // ************************************************************
#endif
}
/* SAM_LISTING_END_6 */

/* SAM_LISTING_BEGIN_7 */
void RugeStuebenAMG::iterate(const Eigen::VectorXd &phi,
                             Eigen::VectorXd &mu) const {
#if SOLUTION

#else
  // ************************************************************
  // Code to be supplemented
  // ************************************************************
#endif
}
/* SAM_LISTING_END_7 */

/* SAM_LISTING_BEGIN_9 */
double cvgRateAMG(const Eigen::SparseMatrix<double> &A, double tol) {
  double rate = 0.0;
#if SOLUTION

#else
  // ************************************************************
  // Code to be supplemented
  // ************************************************************
#endif
  return rate;
}
/* SAM_LISTING_END_9 */

}  // namespace CFSplit
