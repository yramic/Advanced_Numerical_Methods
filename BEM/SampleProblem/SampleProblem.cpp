#include <iostream>
#include <iomanip>

#include <Eigen/Dense>
#include <Eigen/Sparse>

using namespace Eigen;

/* @brief Compute the CCS format of matrix $A$
 * \param[in] A An $n \times n$ matrix without columns only made by zeros
 * \param[out] val Vector of nonzero values of $A$
 * \param[out] row\_ind Row indices of each element of 'val'
 * \param[out] col\_ptr Indices of the elements in 'val' which start a column of $A$
 */
/* SAM_LISTING_BEGIN_0 */
void CCS(const MatrixXd & A, VectorXd & val, VectorXd & row_ind, VectorXd & col_ptr)
{
	// Number of rows and columns
	int m = A.rows();
	int n = A.cols();

	// Number of nonzero entries
	int nnz = 0;
	for(int i=0; i<m; ++i) { // Row iterator
		for(int j=0; j<n; ++j) { // Col iterator
			if(A(i,j) != 0) {
				++nnz;
			}
		}
	}

	// Initialization
	val.resize(nnz);
	row_ind.resize(nnz);
	col_ptr.resize(n);

#if SOLUTION
	// Store $A$ in CCS format
	int index = 0;
	for(int j=0; j<n; ++j) { // Col iterator

		// Update 'col\_ptr'
		col_ptr(j) = index;

		for(int i=0; i<m; ++i) { // Row iterator
			if(A(i,j) != 0) {
				// Record the value to 'val'
				val(index) = A(i,j);
				// Record the row index to 'row\_ind'
				row_ind(index) = i;
				++index;
			}
		}
	}
#else // TEMPLATE
    // TODO: compute the CCS format of matrix $A$
#endif // TEMPLATE
}
/* SAM_LISTING_END_0 */

int main() {
    // Initialization
    unsigned int n = 6;
    MatrixXd A(n,n);
	// Poisson matrix
    A <<  4, -1,  0, -1,  0,  0,
         -1,  4, -1,  0, -1,  0,
          0, -1,  4,  0,  0, -1,
         -1,  0,  0,  4, -1,  0,
          0, -1,  0, -1,  4, -1,
		  0,  0, -1,  0, -1,  4;

    // Test 'CCS'
	VectorXd val_1, row_ind_1, col_ptr_1;
      CCS(A, val_1, row_ind_1, col_ptr_1);

/* @brief Compute the CCS format of matrix $A$ using Eigen methods
 */
/* SAM_LISTING_BEGIN_1 */
	double *  val_2;
	int * row_ind_2;
	int * col_ptr_2;

	// Caste $A$ to 'SparseMatrix' and compress it to store the matrix in CCS format
	SparseMatrix<double> As = A.sparseView();
	As.makeCompressed();

	val_2 = As.valuePtr(); // Pointer to values
	row_ind_2 = As.innerIndexPtr(); // Pointer to indices
	col_ptr_2 = As.outerIndexPtr(); // Pointer to first indices of each inner vector
/* SAM_LISTING_END_1 */

	// Verify that the solutions are the same
	// Compute l2-norm of the differences between the CCS vectors
	double diff_val = 0, diff_row_ind = 0, diff_col_ptr = 0;
	for(int i=0; i<val_1.size(); ++i) {
		diff_val += pow(val_1(i) - *(val_2+i), 2);
		diff_row_ind += pow(row_ind_1(i) - *(row_ind_2+i), 2);
	}
	for(int i=0; i<col_ptr_1.size(); ++i) {
		diff_col_ptr += pow(col_ptr_1(i) - *(col_ptr_2+i), 2);
	}
	std::cout << "l2-norm of the difference between val = " << sqrt(diff_val) << std::endl;
	std::cout << "l2-norm of the difference between row_ind = " << sqrt(diff_row_ind) << std::endl;
	std::cout << "l2-norm of the difference between col_ptr = " << sqrt(diff_col_ptr) << std::endl;
}
