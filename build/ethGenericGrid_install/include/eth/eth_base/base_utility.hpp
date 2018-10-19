/********************************************************************
  created:	2013/07/23
  created:	23:7:2013   10:28
  filename: 	base_utility.hpp
  author:		Raffael Casagrande
  
  purpose:	Provide some very general utility functions...
*********************************************************************/


#ifndef ETH_BASE_UTILITY_HPP
#define ETH_BASE_UTILITY_HPP

namespace eth {
  namespace base {
    
    //--------------------------------------------------------------------------
    /** \function makeConst 
     * \brief Turn a mutable into a constant expression. E.g., this may be 
     * useful when one wants to establish a const_iterator on a mutable 
     * container by using the auto-keyword
     */
    //--------------------------------------------------------------------------
    template< typename T >
    const T& makeConst( T& ref ) { return static_cast< const T& >( ref ); }

  }
}

#endif // ETH_BASE_UTILITY_HPP
