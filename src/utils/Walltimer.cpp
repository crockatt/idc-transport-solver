//------------------------------------------------------------------------------------------------------------
//! \file   utils/Walltimer.cpp
//! \brief  Implementation of Walltimer class.
//!
//! \author Michael Crockatt
//! \date   May 2017
//------------------------------------------------------------------------------------------------------------

# include <cmath>
# include <cstring>
# include <ctime>

# include "utils/CLog.hpp"
# include "utils/global.hpp"
# include "utils/Walltimer.hpp"


//------------------------------------------------------------------------------------------------------------
//! \brief  Obtains the current walltime as a double-precision floating point value of seconds.
//!
//! Obtains the current wall-clock time as given by \c CLOCK_MONOTONIC. The value is scaled to a
//! double-precision floating point value given in units of seconds.
//!
//! \param[in,out]  wcTime      Pointer to double value to place current time in.
//------------------------------------------------------------------------------------------------------------
static void getWalltime_d (

    double * const wcTime
) {

    struct timespec clockTime;

    clock_gettime( CLOCK_MONOTONIC, &clockTime );

    *wcTime = clockTime.tv_sec + 1.0e-9 * clockTime.tv_nsec;
}


//============================================================================================================
//=== CONSTRUCTORS AND DESTRUCTOR ============================================================================
//============================================================================================================


//------------------------------------------------------------------------------------------------------------
//! \brief  Constructs a zero-initialized Walltimer object.
//------------------------------------------------------------------------------------------------------------
Walltimer::Walltimer( void ) :

    start_time( std::nan("0") ),
    stop_time( std::nan("0") ),
    elapsed_time( 0.0 ),
    is_running( false )
{}


//============================================================================================================
//=== INTERFACE ROUTINES =====================================================================================
//============================================================================================================


//------------------------------------------------------------------------------------------------------------
//! \brief  Clears the Walltimer object.
//------------------------------------------------------------------------------------------------------------
Walltimer & Walltimer::Clear( void ) {

    PRINT_STATUS( "Clearing timer %p.\n", this )

    this->start_time    = std::nan( "0" );
    this->stop_time     = std::nan( "0" );
    this->elapsed_time  = 0.0;
    this->is_running    = false;

    return *this;
}


//------------------------------------------------------------------------------------------------------------
//! \brief  Starts timing.
//------------------------------------------------------------------------------------------------------------
Walltimer & Walltimer::Start( void ) {

    PRINT_STATUS( "Starting timer %p.\n", this )

    if ( this->is_running ) {

        PRINT_WARNING( "Walltimer already running. Ignoring start.\n" )
        return *this;
    }

    getWalltime_d( &this->start_time );
    this->is_running = true;

    return *this;
}


//------------------------------------------------------------------------------------------------------------
//! \brief  Stops the given Walltimer object, adding the time elapsed since \pp{this->start_time} to
//!         \pp{this->elapsed_time}.
//------------------------------------------------------------------------------------------------------------
Walltimer & Walltimer::Stop( void ) {

    PRINT_STATUS( "Stopping timer %p.\n", this )

    if ( !this->is_running ) {

        PRINT_WARNING( "Timer not running. Ignoring stop.\n" )
        return *this;
    }

    getWalltime_d( &this->stop_time );

    this->elapsed_time += this->stop_time - this->start_time;

    this->is_running    = false;
    this->start_time    = std::nan( "0" );

    return *this;
}


//------------------------------------------------------------------------------------------------------------
//! \brief  Returns the time elapsed in seconds.
//------------------------------------------------------------------------------------------------------------
const double & Walltimer::GetElapsedTime( void ) const {

    return this->elapsed_time;
}


//------------------------------------------------------------------------------------------------------------
//! \brief  Sets the elapsed time stored using the given value.
//!
//! The given value is assumed to be in units of seconds.
//------------------------------------------------------------------------------------------------------------
void Walltimer::SetElapsedTime(

    const double set_time
) {

    this->elapsed_time = set_time;
}
