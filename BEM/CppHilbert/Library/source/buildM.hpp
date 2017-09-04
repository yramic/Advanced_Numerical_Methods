///////////////////////////////////////////////////////////////////////////////
/// \file buildM.hpp
/// \brief This file provides the function to calculate the mass matrices for
///        low order Galerkin discretization
///
///  This file is part of the HILBERT program package for the numerical
///  solution of the Laplace equation with mixed boundary conditions by use of
///  BEM in 2D.
///
///  C++ adaptation for ANCSE17 of HILBERT V3.1 TUWien 2009-2013
///////////////////////////////////////////////////////////////////////////////


#include <Eigen/Sparse>

/**
 *  Assembles a mass-type matrix M for P0 x S1.
 *
 *  The entries of the (nE x nC)-matrix M for P0 x S1 read \f$ M_{ij} = \int_{Ei}
 *  \phi_j ds, \f$, where \f$ \phi_j \f$ is the S1-hat function associated with
 *  the node zj. The output M is a sparse matrix.
 *
 * @param[out] M
 * @param[in] coordinates  (nC x 2) matrix of the mesh nodes {z1,...,zn}.
 * @param[in] elements  (nE x 2) Matrix representing a partition {E1,...,Em} of a
 *                  boundary \f$\Gamma\f$.
 */
void computeM01(Eigen::SparseMatrix<double>& M, const Eigen::MatrixXd& coordinates,
            Eigen::MatrixXi& elements)
{

    int nE = elements.rows();
    // Define triplets vector
    typedef Eigen::Triplet<double> triplet_t;
    std::vector<triplet_t> triplets;
    triplets.reserve(2*nE);

    // traverse elements
    for(int i=0; i<nE; i++){
        // identify element's vertices
        int aidx = elements(i,0);
        int bidx = elements(i,1);
        const Eigen::Vector2d& a = coordinates.row(aidx);
        const Eigen::Vector2d& b = coordinates.row(bidx);

        // Fill triplets with the contribution corresponding to their associated
        // basis functions
        double h = (b-a).norm();
        triplets.push_back(triplet_t(i, aidx, h/2.));
        triplets.push_back(triplet_t(i, bidx, h/2.));
    }

    M.setFromTriplets(triplets.begin(), triplets.end());
}

/**
 *  Assembles a mass-type matrix M for S1 x S1.
 *
 *  The entries of the (nC x nC)-matrix M for S1 x S1 read \f$ M_{ij} = \int_{
 *  supp \phi_i} \int_{supp \phi_j} \phi_i \phi_j ds, \f$, where \f$ \phi_i \f$
 *  is the S1-hat function associated with the node zi. The output M is a sparse matrix.
 *
 * @param[out] M
 * @param[in] coordinates  (nC x 2) matrix of the mesh nodes {z1,...,zn}.
 * @param[in] elements  (nE x 2) Matrix representing a partition {E1,...,Em} of a
 *                  boundary \f$\Gamma\f$.
 */
void computeM11(Eigen::SparseMatrix<double>& M, const Eigen::MatrixXd& coordinates,
            Eigen::MatrixXi& elements)
{

    int nE = elements.rows();
    int nC = coordinates.rows();
    // Define triplets vector
    typedef Eigen::Triplet<double> triplet_t;
    std::vector<triplet_t> triplets;
    triplets.reserve(3*nE);

    // traverse elements
    for(int i=0; i<nE; i++){
        // identify element's vertices
        int aidx = elements(i,0);
        int bidx = elements(i,1);
        const Eigen::Vector2d& a = coordinates.row(aidx);
        const Eigen::Vector2d& b = coordinates.row(bidx);

        // Fill triplets with the contribution corresponding to their associated
        // basis functions
        double h = (b-a).norm();
        triplets.push_back(triplet_t(aidx, aidx, h/3.));
        triplets.push_back(triplet_t(aidx, bidx, h/6.));
        triplets.push_back(triplet_t(bidx, aidx, h/6.));
        triplets.push_back(triplet_t(bidx, bidx, h/3.));
    }

    M.setFromTriplets(triplets.begin(), triplets.end());
}

