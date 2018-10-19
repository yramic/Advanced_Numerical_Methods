/********************************************************************
  created:	2013/06/17
  created:	17:6:2013   17:45
  filename: 	entityCollection.h
  author:		Raffael Casagrande
  
  purpose:	A minimal interface for traversing a collection of entities.
*********************************************************************/

#ifndef entityCollection_HPP__
#define entityCollection_HPP__




namespace eth{
  namespace grid{
    template<class ITERATOR>
    class EntityCollection;
  }
}

#include <boost/utility.hpp>

/**
*	\brief A minimal interface for traversing a collection of entities
*		   
*	\remarks This class is most of all used to group a begin and end iterator
*			 together such that the new c++ range based for loop applies.
*/
template<class ITERATOR>
class eth::grid::EntityCollection {
public:
  /// the type of iterators returned by begin() and end().
  typedef ITERATOR iterator_t;

  /**
   * \brief	Create an entity collection with the 
   *
   * \param	begin	The begin.
   * \param	end  	The end.
   */
  EntityCollection(iterator_t begin, iterator_t end) 
    : begin_(std::move(begin)), end_(std::move(end)){

  }

  /// Returns a iterator to the first element of the collection
  iterator_t begin() const {
    return begin_;
  }

  /// Returns an iterator one past the last element.
  iterator_t end() const {
    return end_;
  }

private:
  iterator_t begin_;
  iterator_t end_;
};


#endif // entityCollection_HPP__