/********************************************************************
  created:	2013/07/12
  created:	12:7:2013   10:49
  filename: 	storage_types.hpp
  author:		Raffael Casagrande
  
  purpose:	Define possible ways in which an Engine class can
        store the underlying implementation...
*********************************************************************/


namespace eth{
  namespace grid {
    //! \brief Possible types of storage (only important for implementor)
    //			(This is so far only used by Geometry Wrapper).
    enum class StorageType {
      UniquePtr,
      Reference
    };

    //! \brief A Helper class which provides types and functions for 
    // different storage strategies.
    template<StorageType TYPE, class IMPL>
    struct StorageHelper;

    template<class IMPL>
    struct StorageHelper<eth::grid::StorageType::Reference,IMPL>;

    template<class IMPL>
    struct StorageHelper<eth::grid::StorageType::UniquePtr,IMPL>;

  }
}

// Implementation
//////////////////////////////////////////////////////////////////////////

#include <memory>

template<class IMPL>
struct eth::grid::StorageHelper<eth::grid::StorageType::Reference,IMPL>
{
  using storageType_t = const IMPL&;
  
  static inline const IMPL& ref(const IMPL& object)
  {
    return object;
  }

};

template<class IMPL>
struct eth::grid::StorageHelper<eth::grid::StorageType::UniquePtr,IMPL>
{  
  using storageType_t = std::unique_ptr<const IMPL>;
  
  static inline const IMPL& ref(const std::unique_ptr<const IMPL>& object) {
    return *object;
  }

};
