/********************************************************************
  created:	2013/07/03
  created:	3:7:2013   18:27
  filename: 	single_type_mapper.h
  author:		Raffael Casagrande
  
  purpose:	Provide the implementation for a Mapper mapping from
        entities with a given ReferenceElementType(s) to
        consecutive, zero based indices.
*********************************************************************/

#ifndef ETH_GRID_GENERIC_MAPPER_HPP
#define ETH_GRID_GENERIC_MAPPER_HPP

#include "eth_base/ref_el_types.hpp"
#include "eth_base/apply_numeric.hpp"

namespace eth{
  namespace grid{
    template<class GRID_TRAITS, eth::base::RefElType ... TYPES>
    class GenericMapper;
  }
}

// Implementation
//////////////////////////////////////////////////////////////////////////
#include "boost/assert.hpp"
#include "entity.hpp"
#include "mapper.hpp"



/**
 * \brief	This is a generic mapper which allows one to map a set of entities
 * 			of given ReferenceElementTypes to a zero-based and consecutive index.
 * 			(see Mapper for general info about mappers.)
 * \tparam	IS_IMPL	The type of index set on which the mapper is build.
 * \tparam	TYPES	A lis of all ReferenceElement Types which should be mapped
 * 					by this mapper (variadic template!)
 * 
 * Let \f$ E \f$ denote the set of all entities in the given index set and
 * denote by \f$ E^t \subset E \f$ the set of all entities with Reference
 * Element type \f$ t \f$. This mapper assigns then a consecutive
 * 	   and zero based index to the following subset:
 * \f[ \left\{ e \in E^t | t \in \text{TYPES} \right\} \f].
 */
template<class IS_IMPL, eth::base::RefElType... TYPES>
class eth::grid::GenericMapper : 
  public Mapper<typename IS_IMPL::gridTraits_t,
  eth::grid::GenericMapper<IS_IMPL,TYPES...>>  {
public:
  typedef Mapper<typename IS_IMPL::gridTraits_t,
    eth::grid::GenericMapper<IS_IMPL,TYPES...>> base_t;
  typedef IS_IMPL indexSetImpl_t;
  using typename base_t::gridTraits_t;
  using typename base_t::size_type;

  explicit GenericMapper(const IndexSet<gridTraits_t,indexSetImpl_t>& is) 
    : indexSet_(is), size_(0) {
      for(auto t : is.refElTypes())
        size_ += is.size(t);
  }

  template<int CODIM>
  size_type map(const Entity<gridTraits_t,CODIM>& e) const {
    return offset_<TYPES...>(e.refElType()) + indexSet_.index(e);
  }


  template<int CODIM>
  size_type subMap(const Entity<gridTraits_t,0>& e, size_type i) const {
    return offset_<TYPES...>
      (eth::base::ReferenceElements::getSubEntityType(e.refElType(),CODIM,i)) +
      indexSet_.template subIndex<CODIM>(e,i);
  }

  size_type subMap(const Entity<gridTraits_t,0>& e, int codim,size_type i) const {
    return eth::base::applyNumeric<subMapFunctor,0,base_t::gridTraits_t::dimMesh,
      const Entity<gridTraits_t,0>&,size_type,const GenericMapper&>(codim,e,i,*this);
  }

  size_type size() const {
    return size_;
  }

  template<int CODIM>
  bool contains(const Entity<gridTraits_t,CODIM>& e) const {
    return contains_<TYPES...>(e.refElType());
  }

  bool subContains(const Entity<gridTraits_t,0>& e, int codim, size_type i) const {
    return contains_<TYPES...>
      (eth::base::ReferenceElements::getSubEntityType(e.refElType(),codim,i));
  } 

  virtual ~GenericMapper() {}

private:
  const IndexSet<gridTraits_t,indexSetImpl_t>& indexSet_;
  /// Total number of indices.
  size_type size_;

private:
  // Utility Functions
  //////////////////////////////////////////////////////////////////////////
  /// function which assigns an int to the RefElType of the entity.
  template<eth::base::RefElType T_, eth::base::RefElType T2_, eth::base::RefElType... TYPES_>
  int offset_(eth::base::RefElType t) const {
    if(t==T_) return 0;
    else return indexSet_.size(T_)+ offset_<T2_,TYPES_...>(t);
  }
  /// If there are no more types left, this mapper is called.
  template<eth::base::RefElType T_>
  int offset_(eth::base::RefElType t) const {
    if(t==T_) return 0;
    else 
    ETH_VERIFY_MSG(false, "The given RefElType is not mapped by this Mapper");
    return 0; // only here to suppress warnings...
  }

  template<eth::base::RefElType T_, eth::base::RefElType T2_, eth::base::RefElType... TYPES_>
  bool contains_(eth::base::RefElType t) const {
    if(t==T_) return true;
    else return contains_<T2_,TYPES_...>(t);
  }
  template<eth::base::RefElType T_>
  bool contains_(eth::base::RefElType t) const {
    if(t==T_) return true;
    else return false;
  }

  template<int CODIM>
  struct subMapFunctor {
    typedef size_type returnType_t;
    static returnType_t invoke(const Entity<gridTraits_t,0>& e,size_type subIndex,
      const GenericMapper& mapper) {
      return mapper.subMap<CODIM>(e,subIndex);
    } 
  };
};



#endif // ETH_GRID_GENERIC_MAPPER_HPP
