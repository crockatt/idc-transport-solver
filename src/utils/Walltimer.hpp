//------------------------------------------------------------------------------------------------------------
//! \file   utils/Walltimer.hpp
//! \brief  Header file for #Walltimer class implementation.
//!
//! \author Michael Crockatt
//! \date   January 2018
//------------------------------------------------------------------------------------------------------------

# ifndef __WALLTIMER_HPP__
# define __WALLTIMER_HPP__


class Walltimer {

private:

    double start_time;
    double stop_time;
    double elapsed_time;
    bool is_running;

public:

    // --- Constructors and destructor. --------------------------------------------------------------- //

    Walltimer( void );

    Walltimer( const Walltimer & ) = delete;
    Walltimer( const Walltimer && ) = delete;

    ~Walltimer( void ) = default;

    // --- Operator overloads. ------------------------------------------------------------------------ //

    Walltimer & operator=( const Walltimer & ) = delete;
    Walltimer & operator=( const Walltimer && ) = delete;

    // --- Interface routines. ------------------------------------------------------------------------ //

    //!
    //! \brief  Clears all accumulated time from the object.
    //!
    Walltimer & Clear( void );

    //!
    //! \brief  Starts the timer.
    //!
    Walltimer & Start( void );

    //!
    //! \brief  Stops the timer and adds the elapsed time since last start to the total time elapsed.
    //!
    Walltimer & Stop( void );

    //!
    //! \brief  Returns true if the timer is running and false otherwise.
    //!
    bool IsRunning( void ) const {  return this->is_running;  }

    //!
    //! \brief  Returns the total time elapsed in seconds.
    //!
    const double & GetElapsedTime( void ) const;

    //!
    //! \brief  Sets the elapsed time stored using the given value.
    //!
    void SetElapsedTime( const double );
};


# endif // ifndef __WALLTIMER_HPP__
