/********************************************************************
  created:	2013/07/17
  created:	17:7:2013   16:00
  filename: 	ETH_ASSERT.hpp
  author:		Raffael Casagrande
  
  purpose:	Define ASSERTION Macros that can be called in the code..
*********************************************************************/
#ifndef ETH_BASE_ASSERT_HPP
#define ETH_BASE_ASSERT_HPP

#include <boost/assert.hpp>
#include <iostream>
// Some boost specific functions:
//  IDE's like Visual Studio perform better if output goes to std::cout or
//  some other stream, so allow user to configure output stream:
 #ifndef BOOST_ASSERT_MSG_OSTREAM
 # define BOOST_ASSERT_MSG_OSTREAM std::cerr
 #endif

 namespace eth {
  namespace base {
    inline void assertion_failed_msg(char const * expr, char const * msg, char const * function,
      char const * file, long line)
    {
      BOOST_ASSERT_MSG_OSTREAM
        << "***** Internal Program Error - assertion (" << expr << ") failed in "
        << function << ":\n"
        << file << '(' << line << "): " << msg << std::endl;
      std::abort();
    }
  }
 }

// ETH_ASSERT : Assert if NDEBUG is undefined:
//////////////////////////////////////////////////////////////////////////

#ifdef ETH_ASSERT
#undef ETH_ASSERT
#endif

#define ETH_ASSERT(__condition__) \
  BOOST_ASSERT(__condition__)


// ETH_ASSERT_MSH : Assert if NDEBUG is undefined and throw message
//////////////////////////////////////////////////////////////////////////

#ifdef ETH_ASSERT_MSG
#undef ETH_ASSERT_MSG
#endif

#define ETH_ASSERT_MSG(__condition__, __message__) \
  BOOST_ASSERT_MSG(__condition__, __message__);

// ETH_VERIFY : Always check, if condition is not satisfied program stops
//////////////////////////////////////////////////////////////////////////

#ifdef ETH_VERIFY
#undef ETH_VERIFY
#endif

#define ETH_VERIFY(__condition__)  ((__condition__) \
  ? ((void)0) \
  : ::eth::base::assertion_failed_msg(#__condition__,"", \
  BOOST_CURRENT_FUNCTION, __FILE__, __LINE__))
  

// ETH_VERIFY_MESSAGE : Always check and show message
//////////////////////////////////////////////////////////////////////////

#ifdef ETH_VERIFY_MSG
#undef ETH_VERIFY_MSG
#endif

#define ETH_VERIFY_MSG(__condition__, __message__)  ((__condition__) \
  ? ((void)0) \
  : ::eth::base::assertion_failed_msg(#__condition__, __message__, \
  BOOST_CURRENT_FUNCTION, __FILE__, __LINE__)) 

  //BOOST_CURRENT_FUNCTION, __FILE__, __LINE__))

#endif // ETH_BASE_ASSERT_HPP
