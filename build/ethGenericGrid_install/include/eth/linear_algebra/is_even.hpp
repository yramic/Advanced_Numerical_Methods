//------------------------------------------------------------------------------
// <preamble>
// </preamble>
//------------------------------------------------------------------------------

//! @file   matrix_eth_operations.hpp
//! @author Lars Kielhorn
//! @date   2013

// syste includes --------------------------------------------------------------
#include <cstddef>

#ifndef ETH_LINALG_IS_EVEN_HPP
#define ETH_LINALG_IS_EVEN_HPP


//------------------------------------------------------------------------------
namespace eth {
  namespace linalg {
    namespace fixed {

      namespace /* anonymous */ {

        template< bool COND >
        struct static_if {
          static const bool value = true;
        };
        
        template< >
        struct static_if< false > {
          static const bool value = false;
        };
        
      } // end namespace anonymous


      template< std::size_t N >
      struct isEven {
        static const bool value = static_if< N == 2*(N/2) >::value;
      };

    } // end namespace fixed
  } // end namespace linalg
} // end namespace eth


#endif // ETH_LINALG_IS_EVEN_HPP
