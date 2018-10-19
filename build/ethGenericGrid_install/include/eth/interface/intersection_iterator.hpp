/********************************************************************
	created:	2013/06/21
	created:	21:6:2013   14:43
	filename: 	intersectionIterator.hpp
	author:		Raffael Casagrande
	
	purpose:	Provide abstract interface for IntersectionIterators
*********************************************************************/

#ifndef intersectionIterator_HPP__
#define intersectionIterator_HPP__


namespace eth{
	namespace grid{
		template<class VIEW_TRAITS>
		class IntersectionIterator;
	}
}

// Implementation
//////////////////////////////////////////////////////////////////////////
#include "intersection.hpp"

/**
 * \brief	Intersection iterator which can iterate over all Intersections of an
 * 			entity with other entities or the boundary.
 * \tparam VIEW_TRAITS the view traits class (see \ref ViewTraitsDoc) which
 * 					   determines the underlying implementation class.
 * 
 * The user can obtain an Intersection Iterator from a GridView object. An
 * IntersectionIterator can iterate through all intersections of a given
 * element (codim=0) with other elements or the boundary.
 * 
 * \remark Intersections are only defined for elements (entities with codim=0)
 * \remark Intersection Iterators are always symmetric: Let \b I1 be an intersection iterator
 * 		   that points to an intersection with inner element \b e1 and outer
 * 		   element \b e2. Then there exists another intersection iterator \b I2 
 * 		   of same type (Leaf/Level) as \b I1 that points to an intersection with inner element
 * 		   \b e2 and outer element \b e1 (cf. the discussion below on Leaf/Level
 * 		   intersection iterators)
 * 
 * <H3> Leaf and Level Intersection Iterators </H3>
 * Let \b e be an Entity of codimension zero. Then there exist two ways two 
 * access the intersections of \b e with other elements resp. the boundary:
 * - Level intersection iterators result by calling GridView::intersections on a 
 *  Level GridView and traverse only entities on the same level. 
 * - Leaf intersection iterators are obtained by calling GridView::intersections 
 * on a Leaf GridView and traverse only leaf entities. If \b e is not a leaf,
 * the GridView::intersections returns an empty Collection (i.e. begin and
 * end iterator are equal).
 * 
 * TODO: Add a picture visualizing this...
 * 
 * <H3> Implementation of specific functions </H3>
 * For clarity we list here all the functions which must be implemented by the 
 * underlying engine class:
 * 
 * Implementation function                                      | Corresponding Interface function(s)
 * -------------------------------------------------------------|-------------------------------------------------
 * II_IMPL(const II_IMPL& other) (Copy constructor)             | IntersectionIterator(const IntersectionIterator& other)
 * void assign(const II_IMPL& rhs)	(assignment)			        	| IntersectionIterator& operator=(const IntersectionIterator& rhs)
 * void increment()											                      	| IntersectionIterator& operator++()
 * const gridTraits_t::intersection_t& dereference() const   		| operator* and operator ->
 * bool operator==(const II_IMPL& other) const			        		| operator== and operator!=
 * 
 * Where II_IMPL stands for the underlying engine class 
 * (viewTraits_t::intersectionIterator_t).
 */
template<class VIEW_TRAITS>
class eth::grid::IntersectionIterator {
public:
	/// The viewTraits (see \ref ViewTraitsDoc)
	typedef VIEW_TRAITS viewTraits_t;
	/// The underlying implementation of this intersection iterator
	typedef typename viewTraits_t::intersectionIterator_t impl_t;
	/// the traits class which describes the generic Grid Implementation
	typedef typename viewTraits_t::gridTraits_t gridTraits_t;
	typedef typename gridTraits_t::size_type size_type;


	/// Copy constructor 
	// (this is default copy constructor, but it is kept for clarity in documentation)
	IntersectionIterator(const IntersectionIterator& other) 
		: impl_(other.impl_) {}

	/// Assignment operator
	// (this is the default assignment operator, it is kept for clarity in documentation)
	IntersectionIterator& operator=(const IntersectionIterator& rhs) {
		impl_.assign(rhs.impl_);
		return *this;
	}

	/// Implementation constructor
	explicit IntersectionIterator(impl_t impl) : impl_(std::move(impl)) {}


	/**
	 * \brief	Increment operator; Moves the intersection iterator to the next
	 * 			intersection.
	 * \return	A viewTraits::Intersection object which contains information
	 * 			about the interseciton of the given entity with other entities.
	 * \remark This method must not be implemented in the engine class, instead
	 * 		   implement void increment();
	 */
	IntersectionIterator& operator++() {
		impl_.increment();
		return *this;
	}

	/**
	 * \brief	Indirection operator. Returns a reference to the Intersection to
	 * 			which this pointer points.
	 * \remark Do not implement this method in the engine class, instead
	 * 		   implement const viewTraits_t::intersection & dereference();
	 */
	const Intersection<gridTraits_t>& operator*() const {
		static_assert(std::is_reference<decltype(impl_.dereference())>::value, 
			"intersectionIterator::dereference() in implementation must return a reference!");
		return impl_.dereference();
	}

	/**
	 * \brief	Member dereference operator.
	 * \remark Do not implement this method in the engine class, instead
	 * 		   implement const viewTraits_t::intersection& dereference()
	 */
	const Intersection<gridTraits_t>* operator->() const {
		static_assert(std::is_reference<decltype(impl_.dereference())>::value, 
			"intersectionIterator::dereference() in implementation must return a reference!");
		return &impl_.dereference();
	}

	/**
	 * \brief	Equality operator, checks if two IntersectionIterators point to
	 * 			the same Intersection.
	 *
	 * \return	true if both intersection iterators point to the same
	 * 			intersection.
	 * \remark implent this method in the engine class (but operator!= is not 
	 * 		   required)
	 */
	bool operator==(const IntersectionIterator& rhs) const {
		return (impl_ == rhs.impl_);
	}

	/**
	 * \brief	Inequality operator.
	 * \remark You do not need to implement this method in the engine class,
	 * 		   operator== is just negated.
	 */
	bool operator!=(const IntersectionIterator& rhs) const {
		return !(impl_ == rhs.impl_);
	}

private:
	/// underlying engine class.
	impl_t impl_;

};

#endif // intersectionIterator_HPP__
