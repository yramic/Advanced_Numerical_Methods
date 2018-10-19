/********************************************************************
	created:	2013/06/27
	created:	27:6:2013   14:34
	filename: 	RefElTypes.h
	author:		Raffael Casagrande
	
	purpose:	List all ReferenceElement types which are supported
				by the interface.
*********************************************************************/

/*! \file refElTypes.hpp */

#ifndef RefElTypes_HPP__
#define RefElTypes_HPP__

#include <iostream>
#include "../eth_base/ref_el_types.hpp"

namespace eth{
	namespace grid{
		/// Use RefElTypes also in interface namespace...
		typedef eth::base::RefElType RefElType;


		template<RefElType TYPE>
		class ReferenceElement;

		class ReferenceElements;
	}
}

// Implementation
//////////////////////////////////////////////////////////////////////////


#include <vector>
#include <boost/assert.hpp>
#include <boost/lexical_cast.hpp>
#include <string>
#include <cassert>
#include <array>
#include <cstdarg>
#include "../eth_base/numeric_types.hpp"


namespace eth {
	namespace grid{

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
		 * \return	A vector containing the values of the nodes.
		 */
		static inline const std::vector<double> getNodeCoord
			(RefElType type, size_type i) {
			return apply<getNodeCoordFunctor_,size_type>(type,i);
		}


	protected:
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
					BOOST_ASSERT_MSG(false,"Unknown element type");
				}
				// This return statement is only here to suppress warnings in 
				// compiler...
				return FUNCTOR<RefElType::POINT>::invoke(args...);
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
			static const returnType_t invoke(int codim, size_type i) {
				return ReferenceElement<TYPE>::getSubEntityNodes(codim,i);
			}
		};

		template<RefElType TYPE> struct getNodeCoordFunctor_ {
			typedef std::vector<double> returnType_t;
			static const returnType_t invoke(size_type i) {
				return std::vector<double>
					(ReferenceElement<TYPE>::getNodeCoord(i).begin(), 
					ReferenceElement<TYPE>::getNodeCoord(i).end());
			}
		};
	};



	template<RefElType TYPE>
	class ReferenceElement {
	private:
		// Helper traits class to determine dimension of Reference Elements...
    // intel compiler support explicit specialization but g++, clang 
    // deal only with partial specialization of member classes. Hence, we add
    // the artificial template parameter DUMMY. Not nice! :-(
		template<RefElType TYPE2,typename DUMMY = int>
		struct typeTraits_;

		template<typename DUMMY> struct typeTraits_<RefElType::POINT,DUMMY>
		{static const int dim = 0;};
		template<typename DUMMY> struct typeTraits_<RefElType::SEGMENT,DUMMY>
		{static const int dim = 1;};
		template<typename DUMMY> struct typeTraits_<RefElType::TRIA,DUMMY>
		{ static const int dim = 2; };
		template<typename DUMMY> struct typeTraits_<RefElType::QUAD,DUMMY>
		{ static const int dim = 2; };
		template<typename DUMMY> struct typeTraits_<RefElType::TETRA,DUMMY>
		{ static const int dim = 3; };
		template<typename DUMMY> struct typeTraits_<RefElType::HEXA,DUMMY>
		{ static const int dim = 3; };
		template<typename DUMMY> struct typeTraits_<RefElType::PYRAMID,DUMMY>
		{ static const int dim = 3; };
		template<typename DUMMY> struct typeTraits_<RefElType::PRISM,DUMMY>
		{ static const int dim = 3; };

	public:
		/// The dimension of this reference element.
		static const int dimension = typeTraits_<TYPE>::dim;

		/// Type used to count nodes and other things
		typedef eth::base::unsigned_t					size_type;
		/// Type used to group node/edge/face numbers together
		typedef std::vector<size_type>				indexVec_t;
		/// Type to group RefElTypes together...
		typedef std::vector<eth::base::RefElType>	typeVec_t;
		/// Type used to specify coordinates of nodes...
		typedef std::array<double, dimension>		coordinate_t;

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
		static const std::vector<coordinate_t>				cornerCoord_;

	public:

		/**
		* \brief	Number of SubEntities of the given codimension?
		* \param CODIM	The Codimension of the SubEntities with respect
		* 					to dimension of this Reference Element.
		*/
		static inline size_type numSubEntities(int codim) {
			BOOST_ASSERT_MSG(codim <= dimension,
				(eth::base::getRefElName(TYPE) + 
				std::string(" has no subentity of codim ") +
				boost::lexical_cast<std::string>(codim)).c_str());
			return subEntityTypes_[codim].size();
		}

		/**
		* \brief	Get the type of i-th subentity of given codimension
		* \param	codim	The codim with respect to dimension of TRIA (2)
		* \param	i	 	Zero-based (sub)index of the subEntity.
		*/
		static inline RefElType getSubEntityType(int codim, size_type i) {
			BOOST_ASSERT_MSG(codim <= dimension,
				(eth::base::getRefElName(TYPE) + 
				std::string(" has no subentity of codim ") +
				boost::lexical_cast<std::string>(codim)).c_str());
			BOOST_ASSERT_MSG(i<ReferenceElement<TYPE>::numSubEntities(codim),
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
			BOOST_ASSERT_MSG(codim <= dimension,
				(eth::base::getRefElName(TYPE) + 
				std::string(" has no subentity of codim ") +
				boost::lexical_cast<std::string>(codim)).c_str());
			BOOST_ASSERT_MSG(i < numSubEntities(codim),
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
			BOOST_ASSERT_MSG(i < numSubEntities(dimension),
				(eth::base::getRefElName(TYPE) + 
				std::string(" has only ") + 
				boost::lexical_cast<std::string>(numSubEntities(dimension)) +
				std::string(" nodes, but you requested coordinates of node") + 
				boost::lexical_cast<std::string>(i)).c_str());
			return cornerCoord_[i];
		}


	};
	}
}

#endif // RefElTypes_HPP__
