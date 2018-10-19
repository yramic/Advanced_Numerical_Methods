//------------------------------------------------------------------------------
// <preamble>
// </preamble>
//------------------------------------------------------------------------------

//! @file   matrix.hpp
//! @author Lars Kielhorn
//! @date   2013

#ifndef ETH_LINALG_MATRIX_HPP
#define ETH_LINALG_MATRIX_HPP

//------------------------------------------------------------------------------
#ifdef EIGEN3
#include "matrix_eigen3.hpp"
#elif ETH_MATRIX
#include "matrix_eth.hpp"
#endif

#endif // ETH_LINALG_FIXED_SIZE_MATRIX_HPP
