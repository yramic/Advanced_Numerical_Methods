/********************************************************************
  created:	2013/06/24
  created:	24:6:2013   8:59
  filename: 	~vsCB32.h
  author:		Raffael Casagrande
  
  purpose:	Provide a macro to check that a function has been
        implemented in the derived class (in the context
        of CRTP pattern)
*********************************************************************/
#ifndef CRTP_check_HPP__
#define CRTP_check_HPP__


#include <boost/assert.hpp>
#include <cassert>
#include <string>
#include <boost/preprocessor.hpp>
#include <type_traits>



namespace eth {
  namespace base {
    namespace detail {
      // Helper to produce better error messages.
      inline std::string CRTPError(char const* method) {
        return std::string("The CRTP method ") + std::string(method) + 
          std::string(" is not implemented in derived class");
      }
      // Helper to determine the correct type to use to temporarily store a value
      template<class T>
      struct typeHelper___;
      template<class T2>
      struct typeHelper___<T2&> {
        static constexpr bool isRValue = false;
        typedef T2& typeToDeclare_t;
      };
      template<class T2>
      struct typeHelper___<T2&&> {
        static constexpr bool isRValue = true;
        typedef T2&& typeToDeclare_t;
      };
      template<class T2>
      struct typeHelper___ {
        static constexpr bool isRValue = false;
        typedef T2 typeToDeclare_t;
      };
    }
  }
}


#ifdef CHECK_CRTP
#undef CHECK_CRTP
#endif

// This macro checks if the CRTP method has been implemented in the subclass. If
// not an error is thrown. (Note that this macro effectively calls the method
// in the subclass!)
// Usage: e.g. CHECK_CRTP(asImp().size(5)) --> Just put a regular function call
// in the brackets. If you also want to return the value, use
// CHECK_CRTP_RETURN(asImp().size(5)) --> calls return as well.
#ifdef NDEBUG
#define CHECK_CRTP(__interface_method__) __interface_method__;
#else
#define CHECK_CRTP(__interface_method__) \
  {\
    static bool UNIQUE_LABEL ## __COUNTER__ = false; \
    if(UNIQUE_LABEL ## __COUNTER__ == true) \
      BOOST_ASSERT_MSG(false, eth::base::detail::CRTPError(#__interface_method__).c_str()); \
    UNIQUE_LABEL ## __COUNTER__ = true; \
    __interface_method__; \
    UNIQUE_LABEL ## __COUNTER__ = false; \
  }
#endif

#ifdef CHECK_CRTP_RETURN
#undef CHECK_CRTP_RETURN
#endif

#ifdef NDEBUG
#define CHECK_CRTP_RETURN(__interface_method__)  return __interface_method__;
#else
#define CHECK_CRTP_RETURN(__interface_method__) \
{\
  static bool UNIQUE_LABEL ## __COUNTER__ = false; \
  if(UNIQUE_LABEL ## __COUNTER__ == true) \
  BOOST_ASSERT_MSG(false, eth::base::detail::CRTPError(#__interface_method__).c_str()); \
  UNIQUE_LABEL ## __COUNTER__ = true; \
  typedef eth::base::detail::typeHelper___<decltype(__interface_method__)> helper_t; \
    typename helper_t::typeToDeclare_t result = __interface_method__; \
    UNIQUE_LABEL ## __COUNTER__ = false; \
    return result; \
}
#endif



#endif // CRTP_check_HPP__
