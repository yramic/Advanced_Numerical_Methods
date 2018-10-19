/********************************************************************
  created:	2013/07/15
  created:	15:7:2013   16:24
  filename: 	entity_types.hpp
  author:		Raffael Casagrande
  
  purpose:	Enumerate all possible element types and mappings
*********************************************************************/



#ifndef entity_types_hybrid_HPP__
#define entity_types_hybrid_HPP__

#include <string>
#include <iostream>
#include "ref_el_types.hpp"

namespace eth {
  namespace base {

      enum class GeometryType : char {
        POINT,
        SEGEMENT2,
        TRIA3,
        QUAD4,
        TET4,
        HEX8,
        PRISM6,
        PYRAMID5
      };

      /// Return the name of the Geometry Type as a string...
      std::string getGeometryTypeName(GeometryType type);

      /// Print a Geometry type to a stream...
      std::ostream& operator<<(std::ostream& os, GeometryType type);


      //!\brief  Additional info about Geometry Types, values are saved in cpp
      //!        file
      template<GeometryType TYPE>
      struct GeometryTypeInfo;

      /// Provide runtime info about the Geometry Types...
      struct GeometryTypeInfos {
        typedef eth::base::unsigned_t size_type;

        static int numNodes(GeometryType type) {
          return apply<numNodesFunctor_>(type);
        }

        static int dimension(GeometryType type) {
          return apply<dimensionFunctor_>(type);
        }

        static eth::base::RefElType refElType(GeometryType type) {
          return apply<refElTypeFunctor_>(type);
        }

        static size_type numSubEntites(GeometryType type, int codim) {
          return apply<numSubEntitiesFunctor_>(type,codim);
        }

        static std::vector<size_type> getSubEntityNodes
          (GeometryType type, int codim, size_type i) {
          return apply<getSubEntityNodesFunctor_>(type,codim,i);
        }

        static GeometryType getSubEntityType(GeometryType type, int codim, size_type i) {
          return apply<getSubEntityTypeFunctor_>(type,codim,i);
        }




      public:
      
        // Template magic: This function takes a TEMPLATE(!) which in turn
        // accepts a GeometryType as template parameter and then calls
        // invoke() on it with the correct number of arguments.
        // Unfortunately the return type must still be specified in the template
        // :-( Hopefully intel fixes the bug soon.
        template<template<GeometryType> class FUNCTOR, class... ARGS>
        static auto apply(GeometryType type, ARGS... args) -> 
          typename FUNCTOR<GeometryType::SEGEMENT2>::returnType_t{
            // The following would actually be much nicer, but doesn't seem to 
            // work with intel :-(
            //typename decltype(FUNCTOR<GeometryType::SEGEMENT2>::invoke(args...))  {
            switch(type) {
            case GeometryType::SEGEMENT2:
              return FUNCTOR<GeometryType::SEGEMENT2>::invoke(args...);
            case GeometryType::TRIA3:
              return FUNCTOR<GeometryType::TRIA3>::invoke(args...);
            case GeometryType::QUAD4:
              return FUNCTOR<GeometryType::QUAD4>::invoke(args...);
            case GeometryType::TET4:
              return FUNCTOR<GeometryType::TET4>::invoke(args...);
            case GeometryType::PYRAMID5:
              return FUNCTOR<GeometryType::PYRAMID5>::invoke(args...);
            case GeometryType::PRISM6:
              return FUNCTOR<GeometryType::PRISM6>::invoke(args...);
            case GeometryType::HEX8:
              return FUNCTOR<GeometryType::HEX8>::invoke(args...);
            default:
              ETH_VERIFY_MSG(false,"Unknown geometry type");
            }
            // This return statement is only here to suppress warnings in 
            // compiler...
            return FUNCTOR<GeometryType::SEGEMENT2>::invoke(args...);
        }

        // The following is a specialization of the above method which
        // has switch statements only for the selected dimension:
        template<int DIM, template<GeometryType> class FUNCTOR, class... ARGS>
        static auto applyDim(GeometryType type, ARGS... args) -> 
          typename FUNCTOR<GeometryType::SEGEMENT2>::returnType_t {
            static_assert(DIM>= 0 && DIM<4, "Only dimension 0-3 are supported.");
            // We have to use a very dirty trick because C++ doesn't allow us
            // to partially specialize template member functions of a template
            // class, so we use overload resolution instead:
            return applyDimHelper<FUNCTOR, ARGS...>(Dummy<DIM>(),type,args...);
        }



      private:
        // Functors for all the methods which forward the call to the actual
        // reference element.
        template<GeometryType TYPE> struct numNodesFunctor_ {
          typedef int returnType_t;
          static int invoke() {
            return GeometryTypeInfo<TYPE>::numNodes;
        }};
        
        template<GeometryType TYPE> struct dimensionFunctor_ {
          typedef int returnType_t;
          static int invoke() {
            return GeometryTypeInfo<TYPE>::dimension;
        }};

        template<GeometryType TYPE> struct refElTypeFunctor_ {
          typedef eth::base::RefElType returnType_t;
          static eth::base::RefElType invoke() {
            return GeometryTypeInfo<TYPE>::refElType;
        }};

        template<GeometryType TYPE> struct numSubEntitiesFunctor_ {
          typedef size_type returnType_t;
          static returnType_t invoke(int codim) {
            return GeometryTypeInfo<TYPE>::numSubEntities(codim);
          }
        };

        template<GeometryType TYPE> struct getSubEntityNodesFunctor_ {
          typedef std::vector<size_type> returnType_t;
          static returnType_t invoke(int codim, size_type i) {
            return GeometryTypeInfo<TYPE>::getSubEntityNodes(codim,i);
          }
        };

        template<GeometryType TYPE> struct getSubEntityTypeFunctor_ {
          typedef GeometryType returnType_t;
          static returnType_t invoke(int codim, size_type i) {
            return GeometryTypeInfo<TYPE>::getSubEntityType(codim,i);
          }
        };

      private:
        // Very dirty trick to make the apply routine which is
        // templated on the dim actually work: (x is actually not used except
        // to distiguish between different dimensions.)
        // For dim = 0:
        template<int DUMMY> struct Dummy {};

        template<template<GeometryType> class FUNCTOR, class... ARGS>
        static auto applyDimHelper(Dummy<0>, GeometryType type, ARGS... args) -> 
          typename FUNCTOR<GeometryType::POINT>::returnType_t {
            ETH_ASSERT_MSG(type == GeometryType::POINT, (std::string("Geometrytype")+
              getGeometryTypeName(type)+ std::string("doesn't have codim 0.")).c_str());
            static_cast<void>( type );
            return FUNCTOR<GeometryType::POINT>::invoke(args...);
        }
        // For dim = 1
        template<template<GeometryType> class FUNCTOR, class... ARGS>
        static auto applyDimHelper(Dummy<1>, GeometryType type, ARGS... args) -> 
          typename FUNCTOR<GeometryType::SEGEMENT2>::returnType_t {
            switch(type) {
            case GeometryType::SEGEMENT2:
              return FUNCTOR<GeometryType::SEGEMENT2>::invoke(args...);
            default:
              ETH_VERIFY_MSG(false,"Unknown geometry type or geometry type is not of dim 1");
            }
            // This return statement is only here to suppress warnings in 
            // compiler...
            return FUNCTOR<GeometryType::SEGEMENT2>::invoke(args...);
        }
        template<template<GeometryType> class FUNCTOR, class... ARGS>
        static auto applyDimHelper(Dummy<2>, GeometryType type, ARGS... args) -> 
          typename FUNCTOR<GeometryType::TRIA3>::returnType_t {
            switch(type) {
            case GeometryType::TRIA3:
              return FUNCTOR<GeometryType::TRIA3>::invoke(args...);
            case GeometryType::QUAD4:
              return FUNCTOR<GeometryType::QUAD4>::invoke(args...);
            default:
              ETH_VERIFY_MSG(false,"Unknown geometry type or geometry type is not of codim 2");
            }
            // This return statement is only here to suppress warnings in 
            // compiler...
            return FUNCTOR<GeometryType::TRIA3>::invoke(args...);
        }
        template<template<GeometryType> class FUNCTOR, class... ARGS>
        static auto applyDimHelper(Dummy<3>,eth::base::GeometryType type, ARGS... args) -> 
          typename FUNCTOR<GeometryType::TET4>::returnType_t {
            switch(type) {
            case GeometryType::TET4:
              return FUNCTOR<GeometryType::TET4>::invoke(args...);
            case GeometryType::PYRAMID5:
              return FUNCTOR<GeometryType::PYRAMID5>::invoke(args...);
            case GeometryType::PRISM6:
              return FUNCTOR<GeometryType::PRISM6>::invoke(args...);
            case GeometryType::HEX8:
              return FUNCTOR<GeometryType::HEX8>::invoke(args...);
            default:
              ETH_VERIFY_MSG(false,"Unknown geometry type or geometry type is not of codim 0");
            }
            // This return statement is only here to suppress warnings in 
            // compiler...
            return FUNCTOR<GeometryType::TET4>::invoke(args...);
        }

        
      };



      // Implementation for concrete types:
      //////////////////////////////////////////////////////////////////////////
      
      //! \brief Base class for all flat GeometryTypeInfo objects which 
      //! correspond directly to a eth::base::ReferenceElement
      template<eth::base::RefElType TYPE>
      class FlatGeometryTypeInfoT {
      public:
        typedef eth::base::unsigned_t size_type;
        // just take values from ReferenceElement
        static inline const std::vector<size_type> getSubEntityNodes
          (int codim, size_type i) {
            return eth::base::ReferenceElement<TYPE>::getSubEntityNodes(codim,i);
        }

        static inline size_type numSubEntities(int codim) {
          return eth::base::ReferenceElement<TYPE>::numSubEntities(codim);
        }

        static inline GeometryType getSubEntityType(int codim, size_type i) {
          eth::base::RefElType retType = eth::base::ReferenceElement<TYPE>::
            getSubEntityType(codim,i);
          switch (retType)
          {
          case eth::base::RefElType::POINT: return GeometryType::POINT;
          case eth::base::RefElType::SEGMENT: return GeometryType::SEGEMENT2;
          case eth::base::RefElType::TRIA: return GeometryType::TRIA3;
          case eth::base::RefElType::QUAD: return GeometryType::QUAD4;
          case eth::base::RefElType::TETRA: return GeometryType::TET4;
          case eth::base::RefElType::PYRAMID: return GeometryType::PYRAMID5;
          case eth::base::RefElType::PRISM: return GeometryType::PRISM6;
          case eth::base::RefElType::HEXA: return GeometryType::HEX8;
          default:
            ETH_VERIFY_MSG(false,"Unknown Element type");
            break;
          }
          // Make compiler happy:
          return GeometryType::POINT;
        }
      };
       
       
      // SEGMENT2:
      template<> struct GeometryTypeInfo<GeometryType::SEGEMENT2> 
      : public FlatGeometryTypeInfoT<eth::base::RefElType::SEGMENT>  {
        typedef eth::base::unsigned_t size_type;
        static const  int dimension = 1;
        static const  int numNodes = 2;
        static const eth::base::RefElType refElType = eth::base::RefElType::SEGMENT;

      };

      // TRIA3:
      template<> struct GeometryTypeInfo<GeometryType::TRIA3>
      : public FlatGeometryTypeInfoT<eth::base::RefElType::TRIA> {
        typedef eth::base::unsigned_t size_type;
        static const  int dimension = 2;
        static const  int numNodes = 3;
        static const eth::base::RefElType refElType = eth::base::RefElType::TRIA;
      };

      // QUAD4:
      template<> struct GeometryTypeInfo<GeometryType::QUAD4>
      : public FlatGeometryTypeInfoT<eth::base::RefElType::QUAD> {
        typedef eth::base::unsigned_t size_type;
        static const  int dimension = 2;
        static const  int numNodes = 4;
        static const eth::base::RefElType refElType = eth::base::RefElType::QUAD;
      };

      // TET4:
      template<> struct GeometryTypeInfo<GeometryType::TET4> 
      : public FlatGeometryTypeInfoT<eth::base::RefElType::TETRA> {
        typedef eth::base::unsigned_t size_type;
        static const  int dimension = 3;
        static const  int numNodes = 4;
        static const eth::base::RefElType refElType = eth::base::RefElType::TETRA;

        // just take values from ReferenceElement
        static inline const std::vector<size_type> getSubEntityNodes
          (int codim, size_type i) {
            return eth::base::ReferenceElement<eth::base::RefElType::TETRA>::
              getSubEntityNodes(codim,i);
        }
      };

      // PYRAMID5:
      template<> struct GeometryTypeInfo<GeometryType::PYRAMID5> 
      : public FlatGeometryTypeInfoT<eth::base::RefElType::PYRAMID> {
        typedef eth::base::unsigned_t size_type;
        static const  int dimension = 3;
        static const  int numNodes = 5;
        static const eth::base::RefElType refElType = eth::base::RefElType::PYRAMID;

        // just take values from ReferenceElement
        static inline const std::vector<size_type> getSubEntityNodes
          (int codim, size_type i) {
            return eth::base::ReferenceElement<eth::base::RefElType::PYRAMID>::
              getSubEntityNodes(codim,i);
        }
      };

      // PRISM6:
      template<> struct GeometryTypeInfo<GeometryType::PRISM6> 
      : public FlatGeometryTypeInfoT<eth::base::RefElType::PRISM> {
        typedef eth::base::unsigned_t size_type;
        static const  int dimension = 3;
        static const  int numNodes = 6;
        static const eth::base::RefElType refElType = eth::base::RefElType::PRISM;

        // just take values from ReferenceElement
        static inline const std::vector<size_type> getSubEntityNodes
          (int codim, size_type i) {
            return eth::base::ReferenceElement<eth::base::RefElType::PRISM>::
              getSubEntityNodes(codim,i);
        }
      };

      // HEX8:
      template<> struct GeometryTypeInfo<GeometryType::HEX8>
      : public FlatGeometryTypeInfoT<eth::base::RefElType::HEXA> {
        typedef eth::base::unsigned_t size_type;
        static const  int dimension = 3;
        static const  int numNodes = 8;
        static const eth::base::RefElType refElType = eth::base::RefElType::HEXA;

        // just take values from ReferenceElement
        static inline const std::vector<size_type> getSubEntityNodes
          (int codim, size_type i) {
            return eth::base::ReferenceElement<eth::base::RefElType::HEXA>::
              getSubEntityNodes(codim,i);
        }
      };
  } // namespace base
} // namespace eth


#endif // entity_types_HPP__
