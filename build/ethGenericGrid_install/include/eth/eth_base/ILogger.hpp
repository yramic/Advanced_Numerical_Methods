/********************************************************************
  created:	2013/09/24
  created:	24:9:2013   10:28
  filename: 	ILogger.hpp
  author:		Raffael Casagrande
  
  purpose:	Provide an abstract (minimalistic) interface for a
            Logger which can be passed around to record debug information.
*********************************************************************/
//          Copyright Raffael Casagrande 2013 - 2015.
// Distributed under the Boost Software License, Version 1.0.
//    (See accompanying file LICENSE_1_0.txt or copy at
//          http://www.boost.org/LICENSE_1_0.txt)


#ifndef _HPP_CECE2435_4DD7_43EE_8B64_013DC456FFBE
#define _HPP_CECE2435_4DD7_43EE_8B64_013DC456FFBE

// declaration:
//////////////////////////////////////////////////////////////////////////
namespace eth {
  namespace base {
    class ILogger;
  }
}

// Definition
//////////////////////////////////////////////////////////////////////////

#include <string>
#include "numeric_types.hpp"

namespace eth {
  namespace base {
    
    /**
     * \brief Abstract interface for any logger which can be passed around
     *        to write log entries. This is used e.g. to make performance tests
     *        or to write output to std::cout during a simulation.
     * \note  You should use the c++ macro LOGGER_ENTRY to write a log entry in
     *        general because this allows us to automatically record the 
     *        line number and the file from which the log-entry comes.
     */
    class ILogger {
    public:
      /**
       * \brief Write an entry to the logger.
       * \param message     The message which should appear.
       * \param importance  How important is this message? (1=most important,
       *                    5 = least important).  
       * \note You should generally prefer the LOGGER_ENTRY macro over this method
       *       because it will automatically record the line number and the file
       *       from which the log-entry is written.
       */
      virtual void entry(const std::string& message, unsigned_t importance)=0;

      /**
       * \brief Write an entry to the logger.
       * \param message     The message which should appear
       * \param importance  How important is this message? (1=most important,
       *                    5 = least important).
       * \param fileName    Filename of the c++ source file from which this 
       *                    log entry is written (usually set through c++ macro)
       * \param lineNumber  The line number in the c++ source file from which
       *                    log entry is written (usually set through c++ macro)
       */
      virtual void entry(const std::string& message, unsigned_t importance,
        const std::string& fileName, unsigned_t lineNumber)=0;

      /// virtual destructor.
      inline virtual ~ILogger(){}
    };
  }
}


// Define LOGGER_ENTRY
//////////////////////////////////////////////////////////////////////////
// This macro takes the following form:
// LOGGER_ENTRY(logger,message,importance)
// logger     : an eth::base::ILogger* object (a pointer)
// message    : The message which should be recorded (must be convertible to std::string)
// importance : the importance of this message (1 = very important, 5 = doesn't matter)
//
// If logger==nullptr, nothing will happen.
// 
// Through the preprocessor macro LOGGER_LEVEL one can determine which log entries are
// actually written out:
// LOGGER_LEVEL = 0 : No entries are written.
// LOGGER_LEVEL = 3 : all entries up to and including importance 3 are written.
//
// If LOGGER_LEVEL is not set, we set it here to 5:

#ifndef LOGGER_LEVEL
#define LOGGER_LEVEL 5
#endif

#ifdef LOGGER_ENTRY
#undef LOGGER_ENTRY
#endif

#define LOGGER_ENTRY(__logger__,__message__,__importance__) \
{\
  if(__logger__!=nullptr && __importance__ <= LOGGER_LEVEL) { \
    __logger__->entry(__message__,__importance__,__FILE__,__LINE__);\
  }\
}


#endif // _HPP_CECE2435_4DD7_43EE_8B64_013DC456FFBE
