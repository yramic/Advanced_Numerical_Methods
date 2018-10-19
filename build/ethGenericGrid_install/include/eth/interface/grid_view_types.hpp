/********************************************************************
  created:	2013/06/21
  created:	21:6:2013   11:47
  filename: 	GridViewTypes.hpp
  author:		Raffael Casagrande
  
  purpose:	Define the possible GridView Types through enum...
*********************************************************************/
#ifndef GridViewTypes_HPP__
#define GridViewTypes_HPP__

namespace eth {
  namespace grid{
    /// Defines the possible View types.
    enum class GridViewTypes : int {
      LevelView = 1,
      LeafView = 2
    };
  }
}


#endif // GridViewTypes_HPP__