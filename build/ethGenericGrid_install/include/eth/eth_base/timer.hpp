/********************************************************************
  created:	2013/07/23
  created:	23:7:2013   15:03
  filename: 	timer.hpp
  author:		Raffael Casagrande
  
  purpose:	Provide a basic interface for timers using the new c++11 features.
*********************************************************************/

#ifndef ETH_TIMER_HPP
#define ETH_TIMER_HPP

// system includes -------------------------------------------------------------
#include <cstddef>
#include <chrono>

// own includes ----------------------------------------------------------------
#include "numeric_types.hpp"

//------------------------------------------------------------------------------
namespace eth {
  namespace base {

    enum class TimeUnit : eth::base::unsigned_t;

    template< enum TimeUnit >
    class Timer;
    
    template< enum TimeUnit tu >
    std::ostream& operator<<( std::ostream& out, const Timer< tu >& timer );

    std::ostream& operator<<( std::ostream& out, TimeUnit tu );
  }
}


//------------------------------------------------------------------------------
/** \enum TimeUnit
 * \brief TimeUnit defines the most common time units which may be used for
 * time measurement operations. Timer-objects are intended to be defined via
 * this enumerator class.
 */
//------------------------------------------------------------------------------
enum class eth::base::TimeUnit : eth::base::unsigned_t 
{ hour, min, sec, milli, micro, nano };

//------------------------------------------------------------------------------
namespace /* anonymous */ {
  template< enum eth::base::TimeUnit TU > struct TimerTraits { };

  //----------------------------------------------------------------------------
  template< >
  struct TimerTraits< eth::base::TimeUnit::hour > {
    typedef std::ratio<3600,1> unit_t;
    static const std::string unit;
  };
  //----------------------------------------------------------------------------
  template< >
  struct TimerTraits< eth::base::TimeUnit::min > {
    typedef std::ratio<60,1> unit_t;
    static const std::string unit;
  };
  //----------------------------------------------------------------------------
  template< >
  struct TimerTraits< eth::base::TimeUnit::sec > {
    typedef std::ratio<1,1> unit_t;
    static const std::string unit;
  };
  //----------------------------------------------------------------------------
  template< >
  struct TimerTraits< eth::base::TimeUnit::milli > {
    typedef std::ratio<1,1000> unit_t;
    static const std::string unit;
  };
  //----------------------------------------------------------------------------
  template< >
  struct TimerTraits< eth::base::TimeUnit::micro > {
    typedef std::ratio<1,1000000> unit_t;
    static const std::string unit;
  };
  //----------------------------------------------------------------------------
  template< >
  struct TimerTraits< eth::base::TimeUnit::nano > {
    typedef std::ratio<1,1000000000> unit_t;
    static const std::string unit;
  };

  const std::string TimerTraits<eth::base::TimeUnit::hour >::unit = "h";
  const std::string TimerTraits<eth::base::TimeUnit::min  >::unit = "min";
  const std::string TimerTraits<eth::base::TimeUnit::sec  >::unit = "sec";
  const std::string TimerTraits<eth::base::TimeUnit::milli>::unit = "milli-sec";
  const std::string TimerTraits<eth::base::TimeUnit::micro>::unit = "micro-sec";
  const std::string TimerTraits<eth::base::TimeUnit::nano >::unit = "nano-sec";
} // end namespace anonymous


//------------------------------------------------------------------------------
/** \class Timer
 * \brief Simple, leightweight, and easy-to-use timer.
 *
 * \tparam TU Time unit as it is defined in the scoped enumerator TimeUnit
 */
//------------------------------------------------------------------------------
template< enum eth::base::TimeUnit TU = eth::base::TimeUnit::sec >
class eth::base::Timer {
private:
  /** @name typedefs concerning the used clock and time data type */
  //@{
#ifdef NO_INITIALIZER_LIST_SUPPORT
  typedef std::chrono::system_clock         clock_t;
#else
  typedef std::chrono::steady_clock          clock_t;
#endif
  typedef std::chrono::time_point< clock_t > time_point_t;
  typedef TimerTraits< TU >                  timer_traits;
  //@}

  /** @name Timer stores begin-, and end-time as well as an status flag */
  //@{
  //! instantiation time, either set by constructor or by restart()-member
  time_point_t clock_begin_;
  //! on stop() the final time point is stored herein
  time_point_t clock_end_;
  //! a flag indicating whether the timer has stopped
  bool         is_stopped_;
  //@}

  //! access to time data is granted only via a stream
  template< enum eth::base::TimeUnit tu >
  friend 
  std::ostream& operator<<( std::ostream& out, const Timer< tu >& timer );

public:
  /** @name constructors, assignments, etc... */
  //@{
  //! only default constructor is provided -> instantiates private data
  Timer( ) { this -> restart( ); }
  //! no copies
  Timer( const Timer& ) = delete;
  //! no assignments
  const Timer& operator=( const Timer& ) = delete;
  //@}

  //! stop the timer to fix time duration until output is performed
  void stop( ) { clock_end_ = clock_t::now( ); is_stopped_ = true; }

  //! restart the timer
  void restart( ) { clock_begin_ = clock_t::now( ); is_stopped_ = false; }

  //! get the time unit
  static constexpr eth::base::TimeUnit unit( ) { return TU; }
  
  //! get time-measurement as a pair < double, time-unit >
  std::pair< double, std::string > get( ) const
  {
    const auto diff = ( is_stopped_ ? clock_end_ - clock_begin_ 
                        : clock_t::now( ) - clock_begin_ );
    typedef std::chrono::duration< double,
                                   typename timer_traits::unit_t > duration_t;
    const duration_t duration( diff );

    return std::make_pair( static_cast<double>( duration.count() ), 
                           timer_traits::unit );
  }

private:
  //! compute the time duration and output to stream
  std::ostream& write_( std::ostream& out ) const
  {
    out << this -> get().first << " " << timer_traits::unit;
    return out;
  }

}; // end class Timer


//------------------------------------------------------------------------------
template< enum eth::base::TimeUnit tu >
std::ostream& eth::base::operator<<( std::ostream& out, 
                                     const eth::base::Timer<tu>& timer )
{
  return timer.write_( out );
}

//------------------------------------------------------------------------------
// some timer typedefs
namespace eth {
  namespace base {
    //! time measurements in hours
    typedef Timer< TimeUnit::hour  > HourTimer;
    //! time measurements in minutes
    typedef Timer< TimeUnit::min   > MinTimer;
    //! time measurements in seconds
    typedef Timer< TimeUnit::sec   > SecTimer;
    //! time measurements in milliseconds
    typedef Timer< TimeUnit::milli > MilliTimer;
    //! time measurements in microseconds
    typedef Timer< TimeUnit::micro > MicroTimer;
    //! time measurements in nanoseconds
    typedef Timer< TimeUnit::nano  > NanoTimer;
  }
}
#endif // ETH_TIMER_HPP
