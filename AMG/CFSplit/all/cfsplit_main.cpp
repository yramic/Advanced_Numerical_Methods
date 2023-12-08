/**
 * @ file cfsplit_main.cpp
 * @ brief NPDE homework CFSplit MAIN FILE
 * @ author R. Hiptmair
 * @ date December 2023
 * @ copyright Developed at SAM, ETH Zurich
 */

#include <filesystem>

#include "cfsplit.h"

int main(int argc, char** argv) {
  std::cout << "Running code for ADVNCSE HW CFSplit" << std::endl;

  // Load file from data directory: default ../../data/galmat.mm
  // An alternative file in the data directory may be given as command line
  // argument
  std::string filename{"galmat.mm"};
  if (argc > 1) {
    filename = argv[1];
  }
  const std::filesystem::path here = __FILE__;
  auto mm_filepath = here.parent_path() / "../data" / filename;
  // Read matrix in Matrix Market Format
  std::cout << "Reading matrix from " << mm_filepath.string() << std::endl;
  Eigen::SparseMatrix<double> A =
      CFSplit::readMMatrixFromFile(mm_filepath.string());
  // Output matrix
  Eigen::VectorXd row_sums = A * Eigen::VectorXd::Constant(A.cols(), 1.0);
  int int_nd_cnt = 0;
  for (int k = 0; k < A.rows(); ++k) {
    if (std::abs(row_sums[k]) < 1.0E-6) {
      int_nd_cnt++;
    }
  }
  std::cout << "Info: A is " << A.rows() << " x " << A.cols()
            << "-matrix, nnz = " << A.nonZeros() << ", " << int_nd_cnt
            << " fully interior nodes" << std::endl;
  int L = 0;                            // Level count
  std::vector<unsigned long int> Adim;  // Number of d.o.f. on coarser levels
  std::vector<unsigned long int> nnz;  // Number of nnz of coarse level matrices
  CFSplit::RugeStuebenAMG AMG(
      A, 10, 0.25, 1.5,
      [&L, &Adim, &nnz](
          const Eigen::SparseMatrix<double>& Al,
          const Eigen::SparseMatrix<double>& /*P*/,
          const std::vector<CFSplit::RugeStuebenAMG::NodeFlag>& /*nodeflag*/)
          -> void {
        L++;
        Adim.push_back(Al.rows());
        nnz.push_back(Al.nonZeros());
      });
  std::cout << "AMG hierarchy:" << std::endl;
  std::cout << "Level " << L << ": A is " << A.rows() << " x " << A.cols()
            << "-matrix, nnz = " << A.nonZeros() << std::endl;
  for (int l = 0; l < L; ++l) {
    std::cout << "Level " << L - l - 1 << ": A is " << Adim[l] << " x "
              << Adim[l] << "-matrix, nnz = " << nnz[l] << std::endl;
  }

  std::cout << "Convergence rate = " << CFSplit::cvgRateAMG(A) << std::endl;
  return 0;
}
