//------------------------------------------------------------------------------------------------------------
//! \file   utils/CLog.hpp
//! \brief  Header for basic C logging interface.
//!
//! \author Michael Crockatt
//! \date   July 2017
//------------------------------------------------------------------------------------------------------------

# ifndef __CLOG_HPP__
# define __CLOG_HPP__


# include <cmath>
# include <cstdint>
# include <cstdio>


//
// Define macros for ANSI color strings.
//
# if defined (COLOR_TERM) || defined (DOXYCOMPILE)

    //!
    //! \brief  ANSI string to turn terminal text red.
    //!
    //! \see    \c #COLOR_TERM
    //!
    # define __RED "\x1B[31m"

    //!
    //! \brief  ANSI string to turn terminal text yellow.
    //!
    //! \see    \c #COLOR_TERM
    //!
    # define __YEL "\x1B[33m"

    //!
    //! \brief  ANSI string to turn terminal text green.
    //!
    //! \see    \c #COLOR_TERM
    //!
    # define __GRN "\x1B[32m"

    //!
    //! \brief  ANSI string to turn terminal text cyan.
    //!
    //! \see    \c #COLOR_TERM
    //!
    # define __CYN "\x1B[36m"

    //!
    //! \brief  ANSI string to turn terminal text blue.
    //!
    //! \see    \c #COLOR_TERM
    //!
    # define __BLU "\x1B[34m"

    //!
    //! \brief  ANSI string to turn terminal text magenta.
    //!
    //! \see    \c #COLOR_TERM
    //!
    # define __MAG "\x1B[35m"

    //!
    //! \brief  ANSI string to reset text color to normal.
    //!
    //! \see    \c #COLOR_TERM
    //!
    # define __RESET "\x1B[0m"

# else // if defined (COLOR_TERM)

    # define __RED
    # define __YEL
    # define __GRN
    # define __CYN
    # define __BLU
    # define __MAG
    # define __RESET

# endif // if defined (COLOR_TERM)


//============================================================================================================
//=== MACROS AND FUNCTION DEFINITIONS FOR GLOBAL LOGGING INTERFACE ===========================================
//============================================================================================================


//------------------------------------------------------------------------------------------------------------
//! \brief  Wrapper macro for calling printError().
//!
//! Passes the message to print through to PrintError() adding the macros \c \_\_FILE__ and \c \_\_LINE__ to
//! the function call in the appropriate locations.
//!
//! \param[in]      ...         Arguments for message to print, formatted as for \c printf().
//!
//! \see    \c #PRINT_WARNING
//! \see    \c #PRINT_LOG
//! \see    \c #PRINT_STATUS
//! \see    \c #PRINT_NOTE
//!
//! \see    PrintError()
//------------------------------------------------------------------------------------------------------------
# define PRINT_ERROR(...) \
    PrintError( __FILE__, __LINE__, __VA_ARGS__ );


//------------------------------------------------------------------------------------------------------------
//! \brief  Prints an error message to stderr and #Global::logfp.
//------------------------------------------------------------------------------------------------------------
void PrintError( const char * file, const int line, const char * message, ... );


//------------------------------------------------------------------------------------------------------------
//! \brief  Wrapper macro for calling PrintWarning().
//!
//! Passes the message to print through to PrintWarning() adding the macros \c \_\_FILE__ and \c \_\_LINE__ to
//! the function call in the appropriate locations.
//!
//! \param[in]      ...         Arguments for message to print, formatted as for \c printf().
//!
//! \see    \c #PRINT_ERROR
//! \see    \c #PRINT_LOG
//! \see    \c #PRINT_STATUS
//! \see    \c #PRINT_NOTE
//!
//! \see    PrintWarning()
//------------------------------------------------------------------------------------------------------------
# define PRINT_WARNING(...) \
    PrintWarning( __FILE__, __LINE__, __VA_ARGS__ );


//------------------------------------------------------------------------------------------------------------
//! \brief  Prints an error message to stderr and #Global::logfp.
//------------------------------------------------------------------------------------------------------------
void PrintWarning( const char * file, const int line, const char * message, ... );


//------------------------------------------------------------------------------------------------------------
//! \brief  Wrapper macro for calling PrintLog().
//!
//! Passes the message to print through to PrintLog() adding the macros \c \_\_FILE__ and \c \_\_LINE__ to
//! the function call in the appropriate locations.
//!
//! \param[in]      ...         Arguments for message to print, formatted as for \c printf().
//!
//! \see    \c #PRINT_ERROR
//! \see    \c #PRINT_WARNING
//! \see    \c #PRINT_STATUS
//! \see    \c #PRINT_NOTE
//!
//! \see    PrintLog()
//------------------------------------------------------------------------------------------------------------
# define PRINT_LOG(...) \
    PrintLog( __FILE__, __LINE__, __VA_ARGS__ );


//------------------------------------------------------------------------------------------------------------
//! \brief  Prints an error message to stderr and #Global::logfp.
//------------------------------------------------------------------------------------------------------------
void PrintLog( const char * file, const int line, const char * message, ... );


//------------------------------------------------------------------------------------------------------------
//! \brief  Wrapper macro for calling PrintStatus().
//!
//! Passes the message to print through to PrintStatus() adding the macros \c \_\_FILE__ and \c \_\_LINE__ to
//! the function call in the appropriate locations.
//!
//! \attention  If \c #LOGLEVEL < 3, then this macro has no effect.
//!
//! \param[in]      ...         Arguments for message to print, formatted as for \c printf().
//!
//! \see    \c #PRINT_ERROR
//! \see    \c #PRINT_WARNING
//! \see    \c #PRINT_LOG
//! \see    \c #PRINT_NOTE
//!
//! \see    PrintStatus()
//------------------------------------------------------------------------------------------------------------
# if LOGLEVEL >= 3
    # define PRINT_STATUS(...) \
        PrintStatus( __FILE__, __LINE__, __VA_ARGS__ );
# else
    # define PRINT_STATUS(...)
# endif


//------------------------------------------------------------------------------------------------------------
//! \brief  Prints an error message to stderr and #Global::logfp.
//------------------------------------------------------------------------------------------------------------
void PrintStatus( const char * file, const int line, const char * message, ... );


//------------------------------------------------------------------------------------------------------------
//! \brief  Wrapper macro for calling PrintNote().
//!
//! Passes the message to print through to PrintNote() adding the macros \c \_\_FILE__ and \c \_\_LINE__ to
//! the function call in the appropriate locations.
//!
//! \attention  If \c #LOGLEVEL < 2, then this macro has no effect.
//!
//! \param[in]      ...         Arguments for message to print, formatted as for \c printf().
//!
//! \see    \c #PRINT_ERROR
//! \see    \c #PRINT_WARNING
//! \see    \c #PRINT_LOG
//! \see    \c #PRINT_STATUS
//!
//! \see    PrintNote()
//------------------------------------------------------------------------------------------------------------
# if LOGLEVEL >= 2
    # define PRINT_NOTE(...) \
        PrintNote( __FILE__, __LINE__, __VA_ARGS__ );
# else
    # define PRINT_NOTE(...)
# endif


//------------------------------------------------------------------------------------------------------------
//! \brief  Prints an error message to stderr and #Global::logfp.
//------------------------------------------------------------------------------------------------------------
void PrintNote( const char * file, const int line, const char * message, ... );


//============================================================================================================
//=== GLOBAL VARIABLES =======================================================================================
//============================================================================================================


//!
//! \brief  A hack so global variables only need to be declared in this file.
//!
# ifdef NO_EXTERN_CLOG
    # define EXTERN_CLOG
# else
    # define EXTERN_CLOG extern
# endif


//!
//! \brief  Contains global variables.
//!
namespace Global {


//!
//! \brief  File pointer for global log file.
//!
EXTERN_CLOG FILE * logfp;


} // namespace Global



# endif // ifndef __CLOG_HPP__
