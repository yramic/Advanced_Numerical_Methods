/********************************************************************
  created:	2013/07/23
  created:	23:7:2013   15:04
  filename: enum_ref_el_types.hpp
  author:		Lars Kielhorn
  
  purpose:	Define all supported Reference Element types
*********************************************************************/

#ifndef ETH_GRID_ENUM_REF_EL_TYPES_HPP
#define ETH_GRID_ENUM_REF_EL_TYPES_HPP

// system includes -------------------------------------------------------------
#include <string>

//------------------------------------------------------------------------------
namespace eth {
  namespace base {

    enum class RefElType : char;
    /// Return the name of the given RefElType...
    std::string getRefElName(const RefElType& ef);
    /// Output to a stream...
    std::ostream& operator<<( std::ostream& out, const RefElType& ef );

  } // end namespace base
} // end namespace eth


//------------------------------------------------------------------------------
/** \enum RefElType
 * \brief	List all possible reference elements which are in principal
 * 			supported by this interface. The class ReferenceElement
 * 			provides more information about each type.
 * \remark This list must be extended if other types of reference 
 * 		   elements are supported by an implementation.
 */
enum class eth::base::RefElType : char {
  POINT,		//!< One node
    SEGMENT,	//!< LineSegment [0,1]
    TRIA,		//!< triangular element type
    QUAD,		//!< quadrilateral element type
    TETRA,		//!< tetrahedral element type
    HEXA,		//!< hexahedral element type
    PRISM,		//!< prism element type
    PYRAMID		//!< Walk like an egyptian: the pyramid type
    }; // end enum RefElType


#endif // ETH_GRID_ENUM_REF_EL_TYPES_HPP
