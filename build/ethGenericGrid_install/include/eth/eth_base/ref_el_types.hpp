/********************************************************************
  created:	2013/07/23
  created:	23:7:2013   15:04
  filename: 	ref_el_types.hpp
  author:		Raffael Casagrande
  
  purpose:	Define all Reference Elements which are considered in this
            project.
*********************************************************************/

#ifndef ETH_GRID_REF_EL_TYPE_HPP
#define ETH_GRID_REF_EL_TYPE_HPP


#include <iostream>
#include <string>
#include "numeric_types.hpp"
#include "integer_list.hpp"
#include <vector>
#include <array>
#include "ETH_ASSERT.hpp"
#include <boost/lexical_cast.hpp>
#include "type_traits_ref_el.hpp"
#include "enum_ref_el_types.hpp"

//------------------------------------------------------------------------------
namespace eth {
  namespace base {

    // forward declaration.
    template<RefElType TYPE>
    class ReferenceElement;

  } // end namespace base
} // end namespace eth



namespace eth {
  namespace base {

    template<RefElType... TYPES>
    struct RefElList{
    private:
      // forward declarations:
      template<int DIM,RefElType... TYPES_>
      struct onlyDimHelper;

    public:
      /// how many types are stored in this RefElList?
      static const unsigned_t numTypes = sizeof...(TYPES);
      /// return type i (zero based)
      static constexpr RefElType getRefElType(unsigned_t i) {
        return getRefElType_<TYPES...>(i);
      }

      /// create a new RefElList with the given type added at front.
      template<RefElType TYPE>
      using pushFront_t = RefElList<TYPE,TYPES...>;

      /// create a new RefElList from this one that only contains reference elements with the given dimension.
      template<int DIM>
      using onlyDim_t = typename onlyDimHelper<DIM,TYPES...>::truncatedList_t;

      /// with this, a user can obtain the ReferenceElements as a paramter pack.
      template<template<RefElType...> class TEMPLATE>
      using getParamPack = TEMPLATE<TYPES...>;

    private:
      template<RefElType TYPE_1, RefElType TYPE_2, RefElType... TYPES_>
      static constexpr RefElType getRefElType_(unsigned_t i) {
        return i==0 ? TYPE_1 : getRefElType_<TYPE_2,TYPES_...>(i-1);
      }
      template<RefElType TYPE_1>
      static constexpr RefElType getRefElType_(unsigned_t i) {
        return i==0 ? TYPE_1 : throw std::logic_error("This reference element is not contained in this list");
      }

      // Now comes some sophisticated template magic which is used to provide a
      // sublist of all entities of a given dimension:
      template<int DIM, bool isDim, RefElType... TYPES_>
      struct onlyDimHelperBase;

      template<int DIM, RefElType TYPE_1, RefElType TYPE_2, RefElType... TYPES_>
      struct onlyDimHelperBase<DIM,false,TYPE_1,TYPE_2,TYPES_...> : 
        public onlyDimHelperBase<DIM,eth::base::ReferenceElement<TYPE_2>::dimension==DIM,TYPE_2, TYPES_...> {
      };

      template<int DIM, RefElType TYPE_1, RefElType TYPE_2, RefElType... TYPES_>
      struct onlyDimHelperBase<DIM,true,TYPE_1,TYPE_2,TYPES_...> {
        //typedef RefElList<TYPE_1,TYPE_2, TYPES_...> truncatedList_t;
        typedef typename onlyDimHelperBase<DIM,eth::base::ReferenceElement<TYPE_2>::dimension==DIM,TYPE_2,TYPES_...>::truncatedList_t::template pushFront_t<TYPE_1> truncatedList_t;
      };

      template<int DIM, RefElType TYPE_>
      struct onlyDimHelperBase<DIM,false,TYPE_> {
        typedef RefElList<> truncatedList_t;
      };

      template<int DIM, RefElType TYPE_>
      struct onlyDimHelperBase<DIM,true,TYPE_> {
        typedef RefElList<TYPE_> truncatedList_t;
      };

      template<int DIM, RefElType TYPE_1, RefElType TYPE_2, RefElType... TYPES_>
      struct onlyDimHelper<DIM,TYPE_1,TYPE_2,TYPES_...> : 
        public onlyDimHelperBase<DIM,eth::base::ReferenceElement<TYPE_1>::dimension==DIM,TYPE_1,TYPE_2, TYPES_...> {
          typedef onlyDimHelperBase<DIM,eth::base::ReferenceElement<TYPE_1>::dimension==DIM,TYPE_1,TYPE_2, TYPES_...> base_t;
          static_assert(base_t::truncatedList_t::numTypes>0, "This list of RefElTypes is empty.");
      };

      template<int DIM, RefElType TYPE_>
      struct onlyDimHelper<DIM,TYPE_> : public onlyDimHelperBase<DIM,eth::base::ReferenceElement<TYPE_>::dimension==DIM,TYPE_> {
        typedef onlyDimHelperBase<DIM,eth::base::ReferenceElement<TYPE_>::dimension==DIM,TYPE_> base_t;
        static_assert(base_t::truncatedList_t::numTypes>0, "This list of RefElTypes is empty.");
      };

      // end of template magic.
    };

    /// A list of all reference elements which exist.
    typedef RefElList<RefElType::POINT,
      RefElType::SEGMENT,RefElType::TRIA,RefElType::QUAD,RefElType::TETRA,
      RefElType::HEXA,RefElType::PRISM,RefElType::PYRAMID> AllRefElTypesList_t;


    /**
   * \brief	This is a utility class which provides runtime access to an
   * 			arbitrary ReferenceElement. This class has the same
   * 			methods as a ReferenceElement<TYPE> template class but
   * 			instead of having to specify the Type at compile time, we can
   * 			specify the type at runtime. In essence this class
   * 			is just a wrapper and forwards all the methods to the rigth
   * 			template class.
   */
  class ReferenceElements {
  public:
    /// This is a static class, forbid construction
    ReferenceElements() = delete;

    /// Type used to count nodes and other things
    typedef eth::base::unsigned_t						size_type;
    /// Type used to group node/edge/face numbers together
    typedef std::vector<size_type>				indexVec_t;
    /// Type to group RefElTypes together...
    typedef std::vector<eth::base::RefElType>	typeVec_t;

    
    
    
    /// Return the type of the given reference element.
    static int getDimension(RefElType type) {
      return apply<dimFunctor_>(type);
    }

    /**
    * \brief	Number of SubEntities of the given codimension?
    * \param type	The Reference Element type.	
    * \param codim	The Codimension of the SubEntities with respect
    * 					to dimension of this Reference Element.
    */
    static int numSubEntities(RefElType type, int codim) {
      return apply<numSubEntitiesFunctor_,int>(type,codim);
    }

    /**
    * \brief	Get the type of i-th subentity of given codimension
    * \param	type	The Reference Element type.
    * \param	codim	The codim with respect to dimension of TRIA (2)
    * \param	i	 	Zero-based (sub)index of the subEntity.
    */
    static RefElType getSubEntityType(RefElType type, int codim, size_type i) {
      return apply<getSubEntityTypeFunctor_,int,size_type>(type,codim,i);
    }

    /**
     * \brief	Get zero-based node indices for the given subEntity.
     * \param   type	The type of the ReferenceElement.
     * \param	codim	The codim of the subEntity (with respect to
     * 					dimension of this element)
     * \param	i	 	Zero-based index of the subEntity of the given codim
     * \return	The sub entity nodes.
     */
    static inline const indexVec_t getSubEntityNodes
      (RefElType type,int codim, size_type i) {
      return apply<getSubEntityNodesFunctor_,int,size_type>(type,codim,i);
    }

    /**
     * \brief	Get (local) coordinates of the given node.
     * \param	type Reference Element type.
     * \param	i	Zero-based index of the node.
     * \return	A vector containing coordinates of the node.
     */
    static inline const Eigen::Matrix<coord_t, Eigen::Dynamic,1> getNodeCoord
      (RefElType type, size_type i) {
      return apply<getNodeCoordFunctor_,size_type>(type,i);
    }


  public:
    // Template magic: This function takes a TEMPLATE(!) which in turn
    // accepts a RefElType as template parameters and then calls
    // invoke() on it with the correct number of arguments.
    // Unfortunately the return type must still be specified in the template
    // :-( Hopefully intel fixes the bug soon.
    template<template<RefElType> class FUNCTOR, class... ARGS>
    static auto apply(RefElType type, ARGS... args) -> 
      typename FUNCTOR<RefElType::POINT>::returnType_t{
        // The following would actually be much nicer, but doesn't seem to 
        // work with intel :-(
        //decltype(FUNCTOR<RefElType::POINT>::invoke(temp))  {
        switch (type){
        case RefElType::POINT:
          return FUNCTOR<RefElType::POINT>::invoke(args...);
        case RefElType::SEGMENT:
          return FUNCTOR<RefElType::SEGMENT>::invoke(args...);
        case RefElType::TRIA:
          return FUNCTOR<RefElType::TRIA>::invoke(args...);
        case RefElType::QUAD:
          return FUNCTOR<RefElType::QUAD>::invoke(args...);
        case RefElType::TETRA:
          return FUNCTOR<RefElType::TETRA>::invoke(args...);
        case RefElType::HEXA:
          return FUNCTOR<RefElType::HEXA>::invoke(args...);
        case RefElType::PRISM:
          return FUNCTOR<RefElType::PRISM>::invoke(args...);
        case RefElType::PYRAMID:
          return FUNCTOR<RefElType::PYRAMID>::invoke(args...);
        default:
          ETH_ASSERT_MSG(false,"Unknown element type");
        }
        // This return statement is only here to suppress warnings in 
        // compiler...
        return FUNCTOR<RefElType::POINT>::invoke(args...);
    }

    // The following is a specialization of the above method which
    // has switch statements only for the selected dimension:
    template<int DIM, template<RefElType> class FUNCTOR, class... ARGS>
    static auto applyDim(RefElType type, ARGS... args) -> 
      typename FUNCTOR<RefElType::SEGMENT>::returnType_t {
        static_assert(DIM>= 0 && DIM<4, "Only dimension 0-3 are supported.");
        // We have to use a very dirty trick because C++ doesn't allow us
        // to partially specialize template member functions of a template
        // class, so we use overload resolution instead:
        return applyDimHelper<FUNCTOR, ARGS...>(Dummy<DIM>(),type,args...);
    }

    /**
     * \brief The most generical apply: You can specify which ReferenceElements
     *        should be considered and provide your own functor...
     * \tparam FUNCTOR  A template which accepts a RefElType template parameter.
     *                  with a \e static invoke(ARGS...) method (which will
     *                  be called by apply)
     * \tparam LIST_OF_REF_EL_TYPES Must be a RefElList class which has as
     *              describes the RefElTypes which are possibly input through
     *              the type parameter at runtime.
     * \tparam ARGS  The types of the additional arguments which are passed on
     *               to the static invoke method of the FUNCTOR.
     * \param type  The runtime RefElType
     * \param args  Variable arguments providing the arguments.
     * \return  .
     */
    template<template<RefElType> class FUNCTOR, class LIST_OF_REF_EL_TYPES, class... ARGS>
    static auto applyUniversal(RefElType type, ARGS... args) ->
    typename FUNCTOR<LIST_OF_REF_EL_TYPES::getRefElType(0)>::returnType_t {
      static_assert(LIST_OF_REF_EL_TYPES::numTypes>0, "The List of Reference elements must contain at least one element.");
      static const int numTypes = LIST_OF_REF_EL_TYPES::numTypes;

      return applyUniversalHelper<FUNCTOR,LIST_OF_REF_EL_TYPES,numTypes-1, ARGS...>::invoke(type,args...);
    }

    

  private:
    // Functors for all the above methods, forward the call to the actual
    // Reference Element.
    template<RefElType TYPE> struct dimFunctor_ {
      typedef int returnType_t;
      static int invoke(){
        return ReferenceElement<TYPE>::dimension;
    }};

    template<RefElType TYPE> struct numSubEntitiesFunctor_ {
      typedef size_type returnType_t;
      static size_type invoke(int codim) {
        return ReferenceElement<TYPE>::numSubEntities(codim);
      }
    };

    template<RefElType TYPE> struct getSubEntityTypeFunctor_{
      typedef RefElType returnType_t;
      static returnType_t invoke(int codim, size_type i) {
        return ReferenceElement<TYPE>::getSubEntityType(codim,i);
      }
    };

    template<RefElType TYPE> struct getSubEntityNodesFunctor_{
      typedef indexVec_t returnType_t;
      static returnType_t invoke(int codim, size_type i) {
        return ReferenceElement<TYPE>::getSubEntityNodes(codim,i);
      }
    };

    template<RefElType TYPE> struct getNodeCoordFunctor_ {
      typedef Eigen::Matrix<coord_t,Eigen::Dynamic,1> returnType_t;
      static returnType_t invoke(size_type i) {
        return ReferenceElement<TYPE>::getNodeCoord(i);
      }
    };


    // Very dirty trick to make the apply routine which is
    // templated on the dim actually work: (x is actually not used except
    // to distiguish between different dimensions.)
    // For dim = 0:
    template<int DUMMY> struct Dummy {};

    template<template<RefElType> class FUNCTOR, class... ARGS>
    static auto applyDimHelper(Dummy<0>, RefElType type, ARGS... args) -> 
      typename FUNCTOR<RefElType::POINT>::returnType_t {
        ETH_ASSERT_MSG(type == RefElType::POINT, (std::string("RefElType")+
          getRefElName(type)+ std::string("doesn't have codim 0.")).c_str());
        static_cast<void>( type );
        return FUNCTOR<RefElType::POINT>::invoke(args...);
    }
    // For dim = 1
    template<template<RefElType> class FUNCTOR, class... ARGS>
    static auto applyDimHelper(Dummy<1>, RefElType type, ARGS... args) -> 
      typename FUNCTOR<RefElType::SEGMENT>::returnType_t {
        switch(type) {
        case RefElType::SEGMENT:
          return FUNCTOR<RefElType::SEGMENT>::invoke(args...);
        default:
          ETH_VERIFY_MSG(false,"Unknown reference element type or reference element type is not of dim 1");
        }
        // This return statement is only here to suppress warnings in 
        // compiler...
        return FUNCTOR<RefElType::SEGMENT>::invoke(args...);
    }
    template<template<RefElType> class FUNCTOR, class... ARGS>
    static auto applyDimHelper(Dummy<2>, RefElType type, ARGS... args) -> 
      typename FUNCTOR<RefElType::TRIA>::returnType_t {
        switch(type) {
        case RefElType::TRIA:
          return FUNCTOR<RefElType::TRIA>::invoke(args...);
        case RefElType::QUAD:
          return FUNCTOR<RefElType::QUAD>::invoke(args...);
        default:
          ETH_VERIFY_MSG(false,"Unknown reference element type or reference element type is not of dim 2");
        }
        // This return statement is only here to suppress warnings in 
        // compiler...
        return FUNCTOR<RefElType::TRIA>::invoke(args...);
    }
    template<template<RefElType> class FUNCTOR, class... ARGS>
    static auto applyDimHelper(Dummy<3>,RefElType type, ARGS... args) -> 
      typename FUNCTOR<RefElType::TETRA>::returnType_t {
        switch(type) {
        case RefElType::TETRA:
          return FUNCTOR<RefElType::TETRA>::invoke(args...);
        case RefElType::PYRAMID:
          return FUNCTOR<RefElType::PYRAMID>::invoke(args...);
        case RefElType::PRISM:
          return FUNCTOR<RefElType::PRISM>::invoke(args...);
        case RefElType::HEXA:
          return FUNCTOR<RefElType::HEXA>::invoke(args...);
        default:
          ETH_VERIFY_MSG(false,"Unknown reference element type or reference element type is not of dim 3");
        }
        // This return statement is only here to suppress warnings in 
        // compiler...
        return FUNCTOR<RefElType::TETRA>::invoke(args...);
    }

    /// Universal helpers which allows us to make a switch statement for only
    /// a selected range of Reference Elements.
    template<template<RefElType> class FUNCTOR, class LIST_OF_REF_EL_TYPES, int I, class... ARGS>
    struct applyUniversalHelper {
      static auto invoke(RefElType type, ARGS... args)
        -> typename FUNCTOR<LIST_OF_REF_EL_TYPES::getRefElType(0)>::returnType_t {
          return LIST_OF_REF_EL_TYPES::getRefElType(I) == type ? 
            FUNCTOR<LIST_OF_REF_EL_TYPES::getRefElType(I)>::invoke(args...) :
            applyUniversalHelper<FUNCTOR,LIST_OF_REF_EL_TYPES,I-1,ARGS...>::invoke(type,args...);
      }
    };
    
    template<template<RefElType> class FUNCTOR, class LIST_OF_REF_EL_TYPES, class... ARGS>
    struct applyUniversalHelper<FUNCTOR,LIST_OF_REF_EL_TYPES,0,ARGS...> {
      static auto invoke(RefElType type, ARGS... args)
        -> typename FUNCTOR<LIST_OF_REF_EL_TYPES::getRefElType(0)>::returnType_t {
          ETH_ASSERT_MSG(LIST_OF_REF_EL_TYPES::getRefElType(0) == type,
            (std::string("The reference element type ") + getRefElName(type) +
             std::string(" is not contained in the list which is passed to applyUniversal.")).c_str())
            static_cast<void>( type );
            return FUNCTOR<LIST_OF_REF_EL_TYPES::getRefElType(0)>::invoke(args...);
      }
    };
  };



  template<RefElType TYPE>
  class ReferenceElement {
  private:
    /// the underlying type traits class
    typedef detail::typeTraits_<TYPE> typeTraits_t;

  public:
    /// The dimension of this reference element.
    static constexpr int dimension = typeTraits_t::dim;

    /// An integer list which decodes the building blocks of the reference elemtnt
    typedef typename typeTraits_t::entity_list entity_list;
    /// Type used to count nodes and other things
    typedef eth::base::unsigned_t					size_type;
    /// Type used to group node/edge/face numbers together
    typedef std::vector<size_type>				indexVec_t;
    /// Type to group RefElTypes together...
    typedef std::vector<eth::base::RefElType>	typeVec_t;
    /// Type used to specify coordinates of nodes...
    typedef fixedMatrix_t<dimension,1>	coordinate_t;
    /// Type to specify all corner coordinates
    typedef fixedMatrix_t<dimension,
                          typeTraits_t::template numSubEntities<dimension>() > cornerCoordinates_t;

    /// Forbid instantiation of this class, this class is static!
    ReferenceElement() = delete;

  private:
    // The underlying static constants which really define the behavior
    // of this class. (Are implemented differently for different RefElTypes)
    
    /// subEntityTypes_[i][j] = Type of subentity with codim i, index j
    static const std::vector<typeVec_t>					subEntityTypes_;
    /// subEntityNodes_[i][j][k] = Index of k-th node of subentity with codim i, index j.
    static const std::vector<std::vector<indexVec_t>>	subEntityCorners_;
    /// nodeCoord_[i][k] = k-th coordinate of i-th node.
    static const fixedMatrix_t<dimension,
                               typeTraits_t::template numSubEntities<dimension>()>	cornerCoord_;

  public:

    /**
    * \brief	Number of SubEntities of the given codimension (compile-time)
    * \param  codim The Codimension of the SubEntities with respect
    * 				      to dimension of this Reference Element.
    */
    template< int CODIM >
    static constexpr size_type numSubEntities( ) {
      //ETH_ASSERT_MSG(CODIM <= dimension,
      //  (eth::base::getRefElName(TYPE) + 
      //  std::string(" has no subentity of codim ") +
      //  boost::lexical_cast<std::string>(CODIM)).c_str());
      return typeTraits_t::template numSubEntities<CODIM>();
    }

    /**
    * \brief	Number of SubEntities of the given codimension (run-time)
    * \param  codim The Codimension of the SubEntities with respect
    * 				      to dimension of this Reference Element.
    */
    static inline size_type numSubEntities( int codim ) {
      ETH_ASSERT_MSG(codim <= dimension,
        (eth::base::getRefElName(TYPE) + 
        std::string(" has no subentity of codim ") +
        boost::lexical_cast<std::string>(codim)).c_str());
      return typeTraits_t::numSubEntities(codim);
    }


    /**
    * \brief	Get the type of i-th subentity of given codimension
    * \param	codim	The codim with respect to dimension of TRIA (2)
    * \param	i	 	Zero-based (sub)index of the subEntity.
    */
    static inline RefElType getSubEntityType(int codim, size_type i) {
      ETH_ASSERT_MSG(codim <= dimension,
        (eth::base::getRefElName(TYPE) + 
        std::string(" has no subentity of codim ") +
        boost::lexical_cast<std::string>(codim)).c_str());
      ETH_ASSERT_MSG(i<ReferenceElement<TYPE>::numSubEntities(codim),
        (eth::base::getRefElName(TYPE) + 
        std::string(" has only ") + 
        boost::lexical_cast<std::string>(numSubEntities(codim)) +
        std::string(" subentities of codim ") + 
        boost::lexical_cast<std::string>(codim) + 
        std::string(" but you requested type of ") +
        boost::lexical_cast<std::string>(i) +
        std::string("-th subentity.")).c_str());
      return subEntityTypes_[codim][i];
    }

    /**
     * \brief	Get zero-based node indices for the given subEntity.
     * \param	codim	The codim of the subEntity (with respect to
     * 					dimension of this element)
     * \param	i	 	Zero-based index of the subEntity of the given codim
     * \return	The sub entity nodes.
     */
    static inline const indexVec_t getSubEntityNodes(int codim, size_type i) {
      ETH_ASSERT_MSG(codim <= dimension,
        (eth::base::getRefElName(TYPE) + 
        std::string(" has no subentity of codim ") +
        boost::lexical_cast<std::string>(codim)).c_str());
      ETH_ASSERT_MSG(i < numSubEntities(codim),
        (eth::base::getRefElName(TYPE) + 
        std::string(" has only ") + 
        boost::lexical_cast<std::string>(numSubEntities(codim)) +
        std::string(" subentities of codim ") + 
        boost::lexical_cast<std::string>(codim) + 
        std::string(" but you requested nodes of ") +
        boost::lexical_cast<std::string>(i) +
        std::string("-th subentity.")).c_str());
      return subEntityCorners_[codim][i];
    }


    /**
     * \brief	Get (local) coordinates of the given node.
     * \param	i	Zero-based index of the node.
     * \return	A vector containing the values of the nodes.
     */
    static inline const coordinate_t getNodeCoord(size_type i) {
      ETH_ASSERT_MSG(i < numSubEntities(dimension),
        (eth::base::getRefElName(TYPE) + 
        std::string(" has only ") + 
        boost::lexical_cast<std::string>(numSubEntities(dimension)) +
        std::string(" nodes, but you requested coordinates of node") + 
        boost::lexical_cast<std::string>(i)).c_str());
      return cornerCoord_.col(i);
    }


    /**
     * \brief	Get (local) coordinates of all corner nodes.
     * \return	A matrix containing the corner nodes columnwise.
     */
    static inline const cornerCoordinates_t& getNodeCoords() {
      return cornerCoord_;
    }


  };
  }
}

// include instantiations
// #include "ref_el_types.tpp"


#endif // ETH_GRID_REF_EL_TYPE_HPP
