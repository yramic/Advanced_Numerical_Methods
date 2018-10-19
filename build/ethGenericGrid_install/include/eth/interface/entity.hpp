/********************************************************************
  created:	2013/06/18
  created:	18:6:2013   15:20
  filename: 	entity.h
  author:		Raffael Casagrande
  
  purpose:	Provide abstract interface for an entity
*********************************************************************/
/*! \file */

#ifndef _HPP_ECDC6483_88C4_4D97_A430_89DDF4712264
#define _HPP_ECDC6483_88C4_4D97_A430_89DDF4712264

namespace eth{
  namespace grid{
    template<class GRID_TRAITS, int CODIM>
    class Entity;

    // Forward declarations:
    template<class,class> class EntityIterator;
  }
}

// Implementation
//////////////////////////////////////////////////////////////////////////

#include "entity_seed.hpp"
#include "eth_base/CRTP_check.hpp"
#include "geometry.hpp"
#include "entity_collection.hpp"
#include "eth_base/ref_el_types.hpp"
#include <memory>

/**
 * \brief	Abstract interface for a n-dimensional entity.  Entities provide 
 * 			"algebraic" information about the grid, i.e. how the different
 * 			elements are connected with each other. Furthermore they provide
 * 			access to geometry objects which provide the actual element 
 * 			transformation.
 * 			
 * \tparam GRID_TRAITS a class which models \ref GridTraitsDoc
 * \tparam CODIM The codimension of the entity (w.r.t. to dimMesh!)
 * 
 * \note	Entities are always immutable and cannot be copied or assigned to.
 * 			If you need a pointer, use EntityPointer s.
 *
 * \note	The lifetime of an entity generally ends when the user has no 
 * 			EntityPointers that reference the entity anymore. I.e. it's
 * 			possible that the entities are just created on the fly. The 
 * 			same thing holds for the Geometry objects: As soon as the "owning"
 * 			entity is destructed, the geometry object becomes invalid.
 * 
 * <H3> Note for Grid Implementators </H3>
 * This interface class uses the CRTP pattern and you should implement all the
 * methods exactly as they are implemented in this interface. 
 * 
 * Note that this interface class has a protected copy constructor and therefore
 * the actual implementation class may still be copyable (and thus be stored in
 * a STL Container).
 * 
 * In principal this interface could have been realized with the engine pattern
 * as well but there are several advantages when using the CRTP Pattern:
 * - The implementation class can obtain a reference to the interface. This is  
 *	 particularly handy in the implementation of IndexSet / IdSet . 
 * - The interface class has automatically the same lifetime as the underlying   
 *   implementation (i.e. the wrapper class must not be destroyed manually). 
 *   This is especially handy for the implementation of the EntityPointers which
 *   return a REFERENCE to an entity. If the engine pattern was used, the 
 *   EntityPointer class would have needed to obtain the engine wrapper somehow
 *   (this wrapper must be valid for as long as the entity is valid!). Naturally
 *   the wrapper would have been stored with the underlying Entity 
 *   Implementation but this would have complicated programming much more.
 *   CRTP is more natural.
 *   
 * \note If entities were copyable/assignable it would be dangerous to use the 
 * 		 CRTP pattern
 * 		 because of the slicing problem: \code
 * 		 Entity<gridTraits_t,CODIM> e = EntityImplementation(); \endcode
 * 		 would cause the actual implementation to be stripped off e!
 */
template<class GRID_TRAITS, int CODIM>
class eth::grid::Entity {
public:
  typedef GRID_TRAITS gridTraits_t;
  typedef typename gridTraits_t::template entity_t<CODIM> impl_t;
  static const int codim = CODIM;
  typedef typename gridTraits_t::size_type size_type;

  // Forbid Assignment operator:
  Entity& operator=(const Entity&) = delete;

  /**
   * \brief	Return the geometric realization of the entity. 
   * \remark The returned object is valid as long as this entity remains valid.
   */
  Geometry<gridTraits_t,typename gridTraits_t::template geometry_t<codim>>
    geometry() const {
    CHECK_CRTP_RETURN(asImp().geometry());
  }

  /**
   * \brief	Return a memory efficient representation of this entity.
   * 			The returned seed can be converted into an Entity Pointer by
   * 			eth::grid::Grid::entityPointer function.
   */
  const EntitySeed<gridTraits_t,codim> seed() const {
    CHECK_CRTP_RETURN(asImp().seed());
  }

  /**
   * \brief	Return type of underlying reference element.
   * \remark	This is only a shortcut for geometry().type().
   */
  eth::base::RefElType refElType() const {
    CHECK_CRTP_RETURN(asImp().refElType());
  }

  /// Return the zero-based level of this entity.
  size_type level() const {
    CHECK_CRTP_RETURN(asImp().level());
  }

protected:
  /// Forbid copy (at least of interface class).
  Entity(const Entity&)  {};
  /// Default Constructor...
  Entity() {};

private:
  /// CRTP trick
  impl_t& asImp() { return static_cast<impl_t&> (*this);}
  const impl_t& asImp()const {return static_cast<const impl_t&>(*this);}
};

////////////////////////////////////////////////////////////////////////////


namespace eth {
  namespace grid {
/**
 * \brief	Entity of codim 0 (has additional functionality!)
 */
template<class GRID_TRAITS>
class Entity<GRID_TRAITS,0> {
public:
  typedef GRID_TRAITS gridTraits_t;
  typedef typename gridTraits_t::template entity_t<0> impl_t;
  static const int codim = 0;
  typedef typename gridTraits_t::size_type size_type;

  // Delete Assignment operator:
  Entity& operator=(const Entity&) = delete;
  


  /**
    * \brief	Return the geometric realization of the entity. 
    * \remark The returned object is valid as long as this entity remains valid.
    */
    Geometry<gridTraits_t,
    typename gridTraits_t::template geometry_t<codim>> geometry() const 
  {
    CHECK_CRTP_RETURN(asImp().geometry());
  }

  /**
   * \brief	Return a memory efficient representation of this entity.
   * 			The returned seed can be converted into an Entity Pointer by
   * 			eth::grid::Grid::entityPointer function.
   */
  const EntitySeed<gridTraits_t,codim> seed() const {
    CHECK_CRTP_RETURN(asImp().seed());
  }

  /**
   * \brief	Return type of underlying reference element.
   * \remark	This is only a shortcut for geometry().type().
   */
  eth::base::RefElType refElType() const {
    CHECK_CRTP_RETURN(asImp().refElType());
  }

  /// Refinement level of this entity (zero-based)
  size_type level() const {
    CHECK_CRTP_RETURN(asImp().level());
  }

  //////////////////////////////////////////////////////////////////////////
  // Additional functions for codim=0
  //////////////////////////////////////////////////////////////////////////
  

  /// Return number of sub entities of given codimension (is here for convenience,
  /// is a duplicate of ReferenceElements::numSubEntities(refElType(),CODIM)
  template<int CODIM>
  size_type countSubEntities() const {
    static_assert(CODIM>= 0 && CODIM<=gridTraits_t::dimMesh,
      "CODIM is negative or CODIM is greater than mesh dimension.");
    CHECK_CRTP_RETURN(asImp().template countSubEntities<CODIM>());
  }


  /**
   * \brief	Return entity pointer to i-th sub entity of codimension CODIM.
   * \param	subIndex	Zero-based index of the sub entity. They are ordered
   * 						according to the rules specified in the 
   * 						ReferenceElement class.
   */
  template<int CODIM>
  EntityPointer<gridTraits_t,
                typename gridTraits_t::template entityPointer_t<CODIM>> 
    subEntity(size_type subIndex) const 
  {
    static_assert(CODIM>= 0 && CODIM<=gridTraits_t::dimMesh,
      "CODIM is negative or CODIM is greater than mesh dimension.");
    CHECK_CRTP_RETURN(asImp().template subEntity<CODIM>(subIndex));
  }


  /**
   * \brief Return the orientation (true indicates positive, 
   * false indicates negative orientation) for the i-th sub entity of codimension CODIM.
   * \param subIndex Zero-based index of the sub-entity. They are ordered according
   *                 to the rules specified in the ReferenceElement class.
   */
  template< int CODIM >
  inline bool orientationSign( size_type subIndex ) const
  {
    static_assert( CODIM > 0 && CODIM < gridTraits_t::dimMesh,
                   "CODIM is zero or CODIM equals the mesh dimension." );
    CHECK_CRTP_RETURN(asImp().template orientationSign<CODIM>(subIndex));
  }


  /**
   * \brief Return the shift of orientation (E.g. 2 of a triangle with indices {2,1,0} w.r.t. {0,1,2})
   * \param subIndex Zero-based index of the sub-entity. They are ordered according
   *                 to the rules specified in the ReferenceElement class.
   */
  template< int CODIM >
  inline int orientationShift( size_type subIndex ) const
  {
    static_assert( CODIM > 0 && CODIM < gridTraits_t::dimMesh,
                   "CODIM is zero or CODIM equals the mesh dimension." );
    CHECK_CRTP_RETURN(asImp().template orientationShift<CODIM>(subIndex));
  }

  

  /**
   * \brief	Returns the father entity on the next coarser grid. If the entity
   * 			is part of the macro grid (coarses level) the result is undefined.
   */
  EntityPointer<gridTraits_t,
                typename gridTraits_t::template entityPointer_t<0>> father()
    const 
  {
    CHECK_CRTP_RETURN(asImp().father());
  }

  /// Return true if entity has a father which can be accessed through father() method.
  inline bool hasFather() const {
    CHECK_CRTP_RETURN(asImp().hasFather());
  }

  /// Return true if this is a leaf entity.
  inline bool isLeaf() const {
    CHECK_CRTP_RETURN(asImp().isLeaf());
  }

  /**
  * \brief	Get Information about how this entity is embedded geometrically
  * 			in the fathers entity.
  * 
  * This method gives us a Geometry object which maps from the reference
  * element of this entity to the reference element of the fathers entity.
  * This can be used e.g. to interpolate to the finer mesh level...
  */
  Geometry<gridTraits_t, typename gridTraits_t::template 
    localGeometry_t<gridTraits_t::dimMesh,gridTraits_t::dimMesh>> geometryInFather() const 
  {
      CHECK_CRTP_RETURN(asImp().geometryInFather());
  }

  /**
   * \brief	Access to the refinement children of this entity.
   * \param	maxLevel	Include only entities with level <= maxLevel.
   * \return	An entity collection which contains all the refinement children
   * 			of this entity up to the specified level).
   * \note	The returned collection is easily traversed with range-based
   * 			for loops.
   * 
   * The returned entities contain all children entities and possibly children
   * of the children etc. E.g. if a TRIG is subdivided into four TRIG children
   * and the first child is again refined into four children, the entity
   * collection will contain 8 elements (codim==0).
   */
  const 
  EntityCollection<EntityIterator<gridTraits_t,
                                  typename gridTraits_t::hierarchicIterator_t>>
    children(size_type maxLevel) const 
  {
    CHECK_CRTP_RETURN(asImp().children(maxLevel));
  }


  /**
   * \brief	True if entity has been created during last call to Grid::adapt().
   * 			This flag is reset when Grid::postAdapt() is called.
   */
  bool isNew() const {
    CHECK_CRTP_RETURN(asImp().isNew());
  }

  /**
   * \brief	true if this entity might vanish in next adapt call (because
   * 			of coarsening). If it returns false, the entity will still exist
   * 			after the call to Grid::adapt().
   * \note	Grid.adapt() will reset the mightVanish flags.
   */
  bool mightVanish() const {
    CHECK_CRTP_RETURN(asImp().mightVanish());
  }

protected:
  /// Forbid copy (at least of interface class).
  Entity(const Entity&)  {};
  /// Default Constructor...
  Entity() {};
  /// Default destructor (must be protected in CRTP!!!)
  ~Entity() {};


private:
  /// CRTP trick
  impl_t& asImp() { return static_cast<impl_t&> (*this);}
  const impl_t& asImp()const {return static_cast<const impl_t&>(*this);}
};
  } // end namespace grid
} // end namespace eth

#endif // _HPP_ECDC6483_88C4_4D97_A430_89DDF4712264
