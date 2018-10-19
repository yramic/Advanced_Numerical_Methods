/********************************************************************
  created:	2013/06/25
  created:	25:6:2013   15:04
  filename: 	intersection.hpp
  author:		Raffael Casagrande
  
  purpose:	Provide abstract interface for an intersection of one
        entity with another.
*********************************************************************/

#ifndef intersection_HPP__
#define intersection_HPP__

namespace eth{
  namespace grid{
    template<class VIEW_TRAITS>
    class Intersection;
  }
}

// Implementation
//////////////////////////////////////////////////////////////////////////
#include "eth_base.hpp"
#include "geometry.hpp"
#include "entity_pointer.hpp"

/**
 * \brief	Represents a codim=1 geometric object which lies at the intersection
 * 			between two elements (codim=0) or an element and the domain boundary
 * \tparam	VIEW_TRAITS A class which models \ref ViewTraitsDoc. This determines
 * 						the underlying implementation class which can be 
 * 						different for leaf and level views.
 * 
 * On matching grids, an Intersection always belongs to exactly one entity of
 * codim=1 (all interfaces are conforming). However in the case of non-matching
 * grids this is not the case and an Intersection can be just a subset of a face
 * of an element.
 * 
 * Using IntersectionIterator 's one can iterate through the boundaries of a 
 * given element. The union of all such intersections is the boundary of the
 * element.
 * 
 * The method boundary() returns true if and only if the intersection involves
 * the domain's boundary. In this case one can use the Grid::boundaryIndexSet() to
 * obtain a unique index for the corresponding boundary segment in the macro 
 * grid. See its documentation for further information.
 * 
 * \note Intersection objects are immutable and cannot be copied/assigned to. 
 * 		 Use IntersectionIterator 's to pass around references.
 * \note The lifetime of an Intersection generally ends when the user has no
 * 		 IntersectionIterators pointing to it anymore. I.e. Intersections can
 * 		 be created on the fly.		 
 * 
 * 
 * <H3> Note to the Interface-implementator </H3>
 * This interface uses the CRTP in much the same way as the 
 * eth::grid::Entity class does. The implementation must therefore contain 
 * exactly the same methods. See the documentation of the Entity class for a
 * discussion CRTP vs. Engine.
 *
 */
template<class GRID_TRAITS>
class eth::grid::Intersection {
public:
  //////////////////////////////////////////////////////////////////////////
  // typedefs and static constants
  //////////////////////////////////////////////////////////////////////////
  
  typedef GRID_TRAITS gridTraits_t;
  typedef typename gridTraits_t::intersection_t impl_t;
  typedef typename gridTraits_t::size_type size_type;

  // Some geometric stuff:
  /// The dimension of the space in which the mesh is embedded.
  static const int dimWorld = gridTraits_t::dimWorld;
  /// dimension of the reference element.
  static const int dimRefEl = gridTraits_t::dimMesh-1;

  /// type of global coordinates.
  typedef typename gridTraits_t::template fixedSizeMatrix_t<dimWorld,1> globalCoord_t;
  /// type of reference element local coordinates.
  typedef typename gridTraits_t::template fixedSizeMatrix_t<dimRefEl,1> localCoord_t;


  //////////////////////////////////////////////////////////////////////////
  // Methods
  //////////////////////////////////////////////////////////////////////////

  // Forbid Assignment operator:
  Intersection& operator=(const Intersection& rhs) = delete;

  /**
   * \brief	Return the type of the underlying reference element of this
   * 			intersection.
   * \remark	this is only a shortcut for geometry().refElType().
   */
  eth::base::RefElType refElType() const {
    CHECK_CRTP_RETURN(asImp().refElType());
  }

  /**
  * \brief	Return the geometric realization of this intersection.
  * \remark	The returned geometry object is valid for at least as long as this
  * 			intersection remains valid. (i.e. as long as there is an
  * 			IntersectionIterator pointing to it).
  */
  Geometry<gridTraits_t, typename gridTraits_t::template geometry_t<1>>
    geometry() const {
      CHECK_CRTP_RETURN(asImp().geometry());
  }



  /// Is this an intersection with the boundary of the domain?
  bool boundary() const {
    CHECK_CRTP_RETURN(asImp().boundary());
  }

  /// is this intersection shared with another element? (we have always: boundary() != neighbor)
  bool neighbor() const {
    CHECK_CRTP_RETURN(asImp().neighbor());
  }

  /**
  * \brief	Returns an entity pointer to the inside of this intersection.
  */
  EntityPointer<gridTraits_t,typename gridTraits_t::template entityPointer_t<0>> 
    inside() const {
    CHECK_CRTP_RETURN(asImp().inside());
  }

  /**
  * \brief	Returns an entity pointer to the other element with which
  * 			this intersection is shared.
  * \warning If neighbor() is false, respectively boundary() is true,
  * 		   this is an intersection with the boundary and the result of this
  * 		   function is undefined.
  */
  EntityPointer<gridTraits_t,typename gridTraits_t::template entityPointer_t<0>>
    outside() const {
      CHECK_CRTP_RETURN(asImp().outside());
  }

  /**
   * \brief	Is this intersection conforming?
   * 
   * The result of this method is:
   * \code
   * (inside()->subEntity<1>(indexInInside()) == outside()->subEntity<1>(indexInOutside()) || boundary()
   * \endcode
   */
  bool isConforming() const {
    CHECK_CRTP_RETURN(asImp().isConforming());
  }

  /**
   * \brief	This intersection is geometrically part of a face of the inside
   * 			Element. This method returns the (local) index of this face as 
   * 			specified by the generic ReferenceElement class. 
   */
  size_type indexInInside() const {
    CHECK_CRTP_RETURN(asImp().indexInInside());
  }

  /**
   * \brief	This intersection is geometrically part of a face of the outside
   * 			element. This method returns the (local) index of this face
   * 			as specified by the generic Reference element class.
   * \warning The behaviour of this method is unspecified if boundary() 
   * 			is true.
   */
  size_type indexInOutside() const {
    CHECK_CRTP_RETURN(asImp().indexInOutside());
  }

  /**
   * \brief	A mapping which maps from the intersections local coordinates
   * 			to the inside() 's element local coordinates.
   * \note	The returned geometry object is valid for at least as long as this
   * 			intersection remains valid. (i.e. as long as there is an
   * 			IntersectionIterator pointing to it).
   * \note the returned geometry object is not const, otherwise it would not 
   *       be possible to MOVE it around!
   */
  Geometry<gridTraits_t, typename gridTraits_t::template localGeometry_t<gridTraits_t::dimMesh-1,gridTraits_t::dimMesh>>
    geometryInInside() const {
    CHECK_CRTP_RETURN(asImp().geometryInInside());
  }

  /**
   * \brief	A mapping which maps from the intersections local coordinates
   * 			to the outside() 's element local coordinates.
   * \note	The returned geometry object is valid for at least as long as this
   * 			intersection remains valid. (i.e. as long as there is an
   * 			IntersectionIterator pointing to it).
   * \warning	The behaviour of the method is undefined if boundary() == true.
   * \note the returned geometry object is not const, otherwise it would not 
   *       be possible to MOVE it around!
   */
  Geometry<gridTraits_t,typename gridTraits_t::template localGeometry_t<gridTraits_t::dimMesh-1,gridTraits_t::dimMesh>>
    geometryInOutside() const {
    CHECK_CRTP_RETURN(asImp().geometryInOutside());
  }


  /**
  *  \brief Return a normal vector (length not necessarily 1) at the given
  * 		  local coordinates that points from the inside to outside.
  */
  globalCoord_t outerNormal(const localCoord_t& local) const {
    CHECK_CRTP_RETURN(asImp().outerNormal(local));
  }

  /**
  *  \brief Return a normal vector (length == 1) at the given
  * 		  local coordinates that points from the inside to outside.
  */
  globalCoord_t unitOuterNormal(const localCoord_t& local) const {
    CHECK_CRTP_RETURN(asImp().unitOuterNormal(local));
  }

  /**
  *  \brief Return a normal vector (length not necessarily = 1) with length
  *  		  geometry().integrationElement(local). The normal points from
  *  		  inside to the outside.
  */
  globalCoord_t integrationOuterNormal(const localCoord_t& local) const {
    CHECK_CRTP_RETURN(asImp().integrationOuterNormal(local));
  }


  /**
  *  \brief Return a normal vector (length == 1) at the geometry() 's 
  *  		  center (Geometry::center() method) pointing from inside to outside
  */
  globalCoord_t centerUnitOuterNormal() const {
    CHECK_CRTP_RETURN(asImp().centerUnitOuterNormal());
  }

protected:
  /// Forbid copy construction (at least at this level)
  Intersection(const Intersection& /* other */) {  }
  /// Default Constructor:
  Intersection() {}

private:
  /// CRTP trick
  impl_t& asImp() { return static_cast<impl_t&> (*this);}
  const impl_t& asImp()const {return static_cast<const impl_t&>(*this);}
};

#endif // intersection_HPP__
