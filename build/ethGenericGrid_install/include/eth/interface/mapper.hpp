/********************************************************************
  created:	2013/07/03
  created:	3:7:2013   15:14
  filename: 	Mapper.hpp
  author:		Raffael Casagrande
  
  purpose:	Define Abstract interface for a mapper.
*********************************************************************/
#ifndef Mapper_HPP__
#define Mapper_HPP__


namespace eth{
  namespace grid{
    template<class GRID_TRAITS, class MAPPER_IMPL>
    class Mapper;
  }
}

// Implementation
//////////////////////////////////////////////////////////////////////////
#include "entity.hpp"
#include "eth_base/CRTP_check.hpp"


/**
 * \brief	Enforces the abstract interface for mappers (through CRTP) 
 * 			but is in practice almost never used by the user. Instead the user
 * 			should instantiate the concrete implementations directly 
 * 			(e.g. GenericMapper )
 * 
 * Mappers are always built on top of IndexSets and provide the functionality
 * to map only a small subset of all entities of a mesh to a consecutive and
 * zero-based indices (IndexSets contain the basic ingredients to support such
 * a functionality, but do not provide it directly). In particular index Mappers
 * allow the user to define a subset \f$ E' \subseteq E \f$ of all entities
 * \f$ E \f$ in a Leaf and/or Level View. The mapper then maps all 
 * 	   \f$ e \in E' \f$ to an index:
 * \f[ m : E' \mapsto \left[0,|E'|-1\right] \f]
 * 	   
 * The mapped index \f$ m(e) \f$ can then be used to attach user defined data
 * to the grid.
 * 
 * \warning A mapper becomes typically invalid when the grid changes.
 */
template<class GRID_TRAITS, class MAPPER_IMPL>
class eth::grid::Mapper {
public:
  typedef GRID_TRAITS gridTraits_t;
  typedef MAPPER_IMPL impl_t;
  /// This is the type to which we map the entities
  typedef typename gridTraits_t::size_type size_type;

  /**
  * \brief	Map entity to array index. 
  * \param	e	Reference to entity with codim=CODIM
  * \return	An index in the range 0... Max number of entities in Set \f$E\f$
  * \warning	if \f$ e \not \in E \f$ the behavior is undefined and may not
  * 			even be noticed if NDEBUG is set (in this case a wrong array 
  * 			index may be returned...)
  */
  template<int CODIM>
  size_type map(const Entity<gridTraits_t,CODIM>& e) const {
    CHECK_CRTP_RETURN(asImp().map(e));
  }

  /**
   * \brief	Map subentity of an element to its array index.
   * \param	e	 	An element of which the subEntity should be considered
   * \param	i	 	Zero-based index of the subEntity (as specified by 
   * 					ReferenceElement)
   * \param	codim	The codimension of the subEntity.
   * \return	An index in the range 0... Max number of entities in Set \f$E\f$
   * \warning	if the subEntity \f$\not \in E \f$ the behavior is undefined and
   * 			may not	even be noticed if NDEBUG is set (in this case a wrong
   * 			array index may be returned...)
   */
  size_type subMap(const Entity<gridTraits_t,0>& e, int codim, size_type i) const {
    CHECK_CRTP_RETURN(asImp().template subMap(e,codim,i));
  }

  template<int CODIM>
  size_type subMap(const Entity<gridTraits_t,0>& e, size_type i ) const {
    CHECK_CRTP_RETURN(asImp().template subMap<CODIM>(e,i));
  }


  /// Return \f$ |E| \f$.
  size_type size() const {
    CHECK_CRTP_RETURN(asImp().size());
  }

  /// Is the given Entity contained in \f$ E \f$ ?
  template<int CODIM>
  bool contains(const Entity<gridTraits_t,CODIM>& e) const {
    CHECK_CRTP_RETURN(asImp().contains(e)); 
  }

  /**
   * \brief	Is the given subEntity contained in \f$ E \f$ ?
   * \param	e	 	The element of which we take the subEntity.
   * \param	i	 	Zero-based index of the subEntity
   * \param	codim	The codim of the subEntity
   */
  bool subContains(const Entity<gridTraits_t,0>& e, int /* codim */,size_type /* i */) const {
    CHECK_CRTP_RETURN(asImp().contains(e));
  }

  /// Virtual Destructor (Needed for CRTP!!!)
  virtual ~Mapper() {}

protected:
  /// Default Constructor
  Mapper() {};
  /// Forbid Copy Construction (CRTP has slicing problem!)
  Mapper(const Mapper& /* other */) {
  }
  /// Forbid Assignment (CRTP has slicing problem!).
  Mapper& operator=(const Mapper& /* rhs */) {
  }

private:
  /// CRTP trick
  impl_t& asImp() { return static_cast<impl_t&> (*this);}
  const impl_t& asImp()const {return static_cast<const impl_t&>(*this);}
};


#endif // Mapper_HPP__
