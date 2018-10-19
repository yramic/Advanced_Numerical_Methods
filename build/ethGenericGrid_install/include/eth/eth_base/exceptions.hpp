/********************************************************************
  created:	2013/07/16
  created:	16:7:2013   18:10
  filename: 	exceptions.hpp
  author:		Raffael Casagrande
  
  purpose:	Provide a set of exception classes which can be thrown at
            runtime...
*********************************************************************/

#ifndef eth_exceptions_HPP__
#define eth_exceptions_HPP__


#include <boost/exception/all.hpp>
#include <string>

namespace eth {
  namespace base {

    struct EthException : virtual boost::exception, virtual std::exception {
      
      EthException(const std::string& what) : errorMessage_(what) {

      }

      std::string getErrorMessage() const {
        return errorMessage_;
      }

      virtual const char* what() const noexcept {
        return errorMessage_.c_str();
      }

    private:
      std::string errorMessage_;  
    };

    struct GridException : EthException {
      GridException(const std::string& what) : EthException(what) {

      }
    };
  }
}


#endif // exceptions_HPP__
