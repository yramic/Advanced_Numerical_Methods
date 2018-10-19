/********************************************************************
  created:	2013/06/19
  created:	19:6:2013   16:15
  filename: 	IdSet.hpp
  author:		Raffael Casagrande
  
  purpose:	Provide abstract interface for IdSet
*********************************************************************/

#ifndef IdSet_HPP__
#define IdSet_HPP__


namespace eth{
  namespace grid{
    template<class GRID_TRAITS>
    class IdSet;
  }
}


// Implementation
//////////////////////////////////////////////////////////////////////////

#include "eth_base/CRTP_check.hpp"

/**
 * \brief	Map Entities and boundary interfaces to unique, \e persistent id's.
 *         Possibly in constant 
 * 			time (e.g. use Unordered containers or something else...)
 * \tparam GRID_TRAITS A model of \ref GridTraitsDoc . The GridTraits determine
 * 					   the underlying implementation as well as the idType_t.
 * 
 * Every Grid implementation provides exactly one idSet implementation through
 * its gridTraits_t::idSet_t type. It provides essentially one mapping:
 * \f[ m : EI \mapsto \mathbb{I} \f]
 * where \f$ EI \f$ is the set of \e all entities (any codimension, level) joined
 * with the set of all \e boundary intersections and
 * \f$ \mathbb{I} = \{ i \in \mathbb{N}^+ | i>=0 \wedge i<\text{maxSize()-1} \} \f$
 * 	   represents the set of all idTypes_t (See below for more
 * info concerning \f$ \mathbb{I} \f$).
 * 
 * This mapping has the following important properties:
 * - The mapped Id is persistent, i.e. the id of an entity/boundary intersection  
 *   remains the same,   
 *	 even if the grid is changing (main difference to IndexSet 's)
 * - All copies of an entity share the same id.  
 *  (see  \cite bastian2008generic , Definition 16). That is if an entity
 *	 	 is an exact copy of its father, the father and child have the same id.
 * - The Id of a Boundary Intersection equals the id of the corresponding  
 *     boundary intersection on the coarse grid. 
 * - In contrast to IndexSet 's, id's must not be consecutive and zero based,  
 *   they are however assumed to lie in the range \f$ [0,\text{maxSize}()-1] \f$
 * - The idType must be an integral, unsigned type and is specified by   
 *   gridTraits::idType_t (see \ref GridTraitsDoc )
 * - The mapping is unique, that is two elements \f$ a,b \in EI \f$ have never the same id:   
 *   \f$ (m(a) \neq m(b) \f$ unless one is an exact copy of its father (see above).
 *
 *  \note An IdSet is obtained from Grid::idSet() and remains valid as along as 
 *  	  the grid remains valid.
 *
 *
 * 
 * <H3> Note for Grid-implementor </H3>
 * The IdSet interface class uses the CRTP pattern much in the same way as the
 * Entity, IndexSet classes do.
 * 
 * The size which is returned by the maxSize() member function should be as 
 * small as possible because this is the size for which the user must acquire 
 * an Array to attach additional information to entities. 
 * An id set is also allowed to reuse ids of previously deleted entities/intersections
 * and to assign them to newly appearing entities/intersections.
 * 
 */
template<class GRID_TRAITS>
class eth::grid::IdSet {
public:
  typedef GRID_TRAITS gridTraits_t;
  typedef typename gridTraits_t::idSet_t impl_t;
  /// The type to which entities are mapped.
  typedef typename gridTraits_t::idType_t idType_t;
  typedef typename gridTraits_t::size_type size_type;

  static_assert(std::is_integral<idType_t>::value, 
    "idType_t must be an integral data type.");
  static_assert(std::is_unsigned<idType_t>::value,
    "idType_t must be an unsigned data type.");

  /// Forbid assignment
  IdSet& operator=(const IdSet& rhs) = delete;

  /// Map an Entity with codimension=CODIM to its id.
  template<int CODIM>
  idType_t id(const Entity<gridTraits_t,CODIM>& entity) const {
    CHECK_CRTP_RETURN(asImp().id(entity));
  }
  
  /**
   * \brief Map a boundary Intersection to its unique, persistent id.
   * \param i The Boundary Intersection which should be mapped.
   * \return  Zero-based id (see description of this class)
   * 
   * \note This id has nothing to do with the insertion index of a boundary
   *       segment in the GridFactory class.
   */
  idType_t id(const Intersection<gridTraits_t>& i) const {
    ETH_ASSERT_MSG(i.boundary(),"Only boundary intersections are mapped.");
    CHECK_CRTP_RETURN(asImp().id(i));
  }

  /**
   * \brief	Map a subentity of given element to its id.
   * \tparam	CODIM   	The codim of the subentity.
   * \param	entity  	The element from which the subEntity is from.
   * \param	subIndex	Zero-based index of the sub Entity.
   * \remark	This method is the same as calling \code
   * 			id(*entity->subEntity<CODIM>(subIndex)) \endcode
   * 			(not exactly legal code). This method is only here
   * 			for convenience and to support possibly faster implementations.
   */
  template<int CODIM>
  idType_t subId(const Entity<gridTraits_t,0>& entity, size_type subIndex) const {
    CHECK_CRTP_RETURN(asImp().template subId<CODIM>(entity,subIndex));
  }

  /**
   * \brief	maxSize()-1 is the maximum id to which an entity or intersection
   *        can be mapped.
   * \warning maxSize() may change when the grid has changed.
   */
  idType_t maxSize() const {
    CHECK_CRTP_RETURN(asImp().maxSize());
  }



protected:
  /// Forbid copy on Interface level:
  IdSet(const IdSet& /* other */) {
  }
  /// Default constructor:
  IdSet() {}
  /// Virtual destructor (CRTP!)
  virtual ~IdSet() {}

private:
  /// CRTP trick
  impl_t& asImp() { return static_cast<impl_t&> (*this);}
  const impl_t& asImp()const {return static_cast<const impl_t&>(*this);}
};

#endif // IdSet_HPP__
