/********************************************************************
	created:	2013/06/18
	created:	18:6:2013   15:01
	filename: 	entitySeed.hpp
	author:		Raffael Casagrande
	
	purpose:	Provide an abstract interface for Entity seed...
*********************************************************************/

#ifndef entitySeed_HPP__
#define entitySeed_HPP__

namespace eth{
	namespace grid{
		template<class GRID_TRAITS, int CODIM>
		class EntitySeed;
	}
}

// Implementation
//////////////////////////////////////////////////////////////////////////


/**
 * \brief	Store a reference to an entity with a minimal memory footprint.
 * 
 * An EntitySeed can be obtained from an Entity object by calling Entity::seed().
 * The seed can then be turned into an EntityPointer by calling
 * Grid::entityPointer(). 
 *
 * In contrast to the related entity object, the EntitySeed consumes minimal
 * memory (e.g. just an index) and can be copied and assigned to as needed. This
 * is perfect for storing entities e.g. in a std::vector (which needs a copy
 * constructor)
 * In comparison to an entity pointer a EntitySeed should consume even less 
 * memory. On the other hand an EntitySeed cannot be used directly to access the
 * underlying entity, for this an EntityPointer is much more suitable.
 * 
 * <H3> Notes for Grid-Implementor </H3>
 * - The EntitySeed class uses the engine pattern and stores an 
 * implementation object internally (not a reference!). Whenever the EntitySeed
 * is copied or assigned to, the underlying implementation is copied/assigned
 *  as well.
 * - Use the EntitySeed::impl() method in the Grid implementation to create  
 *	 an entityPointer.
 */
template<class GRID_TRAITS, int CODIM>
class eth::grid::EntitySeed {
public:
	/// A model of GridTraits
	typedef GRID_TRAITS gridTraits_t;
	/// Underlying implementation
	typedef typename GRID_TRAITS::template entitySeed_t<CODIM> impl_t;
	/// Codimension of the entity
	static const int codim = CODIM;

	/// Standard Copy constructor (defined for documentation)
	EntitySeed(const EntitySeed& other) : impl_(other.impl_) {}

	/// Standard Assignment operator (defined for documentation)
	EntitySeed& operator=(const EntitySeed& rhs) {
		impl_=rhs.impl_;
		return *this;
	}

	/**
	 * \brief	Return the underlying implementation, this should not be used by
	 * 			the user of the interface. It is only here for the interface
	 * 			implementator.
	 */
	const impl_t& impl() const {
		return impl_;
	}


	/// Implementation constructor (to create wrapper)
	EntitySeed(const impl_t& implementation) :
		impl_(implementation) {
	}

private:
	impl_t impl_;
};


#endif // entitySeed_HPP__
