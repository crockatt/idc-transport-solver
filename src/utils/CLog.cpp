//------------------------------------------------------------------------------------------------------------
//! \file   utils/CLog.cpp
//! \brief  Implementation of basic C logging interface.
//!
//! \author Michael Crockatt
//! \date   July 2017
//------------------------------------------------------------------------------------------------------------

# include <cstdarg>
# include <cstdint>
# include <cstdlib>
# include <cstring>

//!
//! \brief  Trick to declare global variables only in this file and \c extern elsewhere.
//!
# define NO_EXTERN_CLOG

# include "utils/CLog.hpp"

# if defined (ENABLE_MPI)
    # include "utils/global.hpp"
# endif


using namespace Global;


//============================================================================================================
//=== FORMATTED TEXT OUTPUT ROUTINES =========================================================================
//============================================================================================================


//------------------------------------------------------------------------------------------------------------
//! \brief  Prints an error message to stderr and #Global::logfp.
//!
//! Error messages are labeled in red (if the macro \c #COLOR_TERM is defined) and are always accompanied by
//! the filename and line number of the location where the error occurred.
//!
//! \attention  The macro \c #PRINT_ERROR is provided as a wrapper for calling this function.
//!
//! \param[in]      file        Name of file where function is called from (usually just \c \_\_FILE__).
//! \param[in]      line        Line number where function is called from (usually just \c \_\_LINE__).
//! \param[in]      message     C-string containing message to print, possibly with format specifiers.
//! \param[in]      ...         Arguments specifying data to print.
//!
//! \see    \c #PRINT_ERROR
//------------------------------------------------------------------------------------------------------------
void PrintError (

    const char * file,
    const int line,
    const char * message,
    ...
) {

    va_list args;
    va_start( args, message );

    char * buff = new char[ std::strlen(message) + 128 ];

    std::sprintf( buff, __RED "ERROR: " __RESET "[%s:%d] "
        # if defined (ENABLE_MPI)
            "(%d/%d)  ",
        # else
            " ",
        # endif
            file, line
        # if defined (ENABLE_MPI)
            , MPI_rank, MPI_num_ranks
        # endif
            );

    std::strcat( buff, message );

    if ( logfp != nullptr ) {

        va_list args2;
        va_copy( args2, args );
        std::vfprintf( logfp, buff, args2 );
        va_end( args2 );
    }

    vfprintf( stderr, buff, args );
    va_end( args );
    delete [] buff;
}


//------------------------------------------------------------------------------------------------------------
//! \brief  Prints a warning message to stderr and #Global::logfp.
//!
//! Warning messages are labeled in magenta (if the macro \c #COLOR_TERM is defined) and are always
//! accompanied by the filename and line number of the location where the warning was generated.
//!
//! \attention  The macro \c #PRINT_WARNING is provided as a wrapper for calling this function.
//!
//! \param[in]      file        Name of file where function is called from (usually just \c \_\_FILE__).
//! \param[in]      line        Line number where function is called from (usually just \c \_\_LINE__).
//! \param[in]      message     C-string containing message to print, possibly with format specifiers.
//! \param[in]      ...         Arguments specifying data to print
//!
//! \see    \c #PRINT_WARNING
//------------------------------------------------------------------------------------------------------------
void PrintWarning (

    const char * file,
    const int line,
    const char * message,
    ...
) {

    va_list args;
    va_start( args, message );

    char * buff = new char[ std::strlen(message) + 128 ];

    std::sprintf( buff, __MAG "WARNING: " __RESET "[%s:%d] "
        # if defined (ENABLE_MPI)
            "(%d/%d)  ",
        # else
            " ",
        # endif
            file, line
        # if defined (ENABLE_MPI)
            , MPI_rank, MPI_num_ranks
        # endif
            );

    std::strcat( buff, message );

    if ( logfp != nullptr ) {

        va_list args2;
        va_copy( args2, args );
        std::vfprintf( logfp, buff, args2 );
        va_end( args2 );
    }

    std::vfprintf( stderr, buff, args );
    va_end( args );
    delete [] buff;
}


//------------------------------------------------------------------------------------------------------------
//! \brief  Prints a log message to stdout and #Global::logfp.
//!
//! If \c #LOGLEVEL is defined to be at least 2, then the log message is labeled in blue (if the macro
//! \c #COLOR_TERM is defined) and is accompanied by the filename and line number of the location where the
//! log message was generated.
//!
//! If \c #LOGLEVEL is less than 2, then the message is printed as given.
//!
//! \attention  The macro \c #PRINT_LOG is provided as a wrapper for calling this function.
//!
//! \param[in]      file        Name of file where function is called from (usually just \c \_\_FILE__).
//! \param[in]      line        Line number where function is called from (usually just \c \_\_LINE__).
//! \param[in]      message     C-string containing message to print, possibly with format specifiers.
//! \param[in]      ...         Arguments specifying data to print
//!
//! \see    \c #PRINT_LOG
//------------------------------------------------------------------------------------------------------------
void PrintLog (

    const char * file,
    const int line,
    const char * message,
    ...
) {

# if defined (ENABLE_MPI) && LOGLEVEL <= 2
    if ( MPI_rank != 0 ) {  return;  }
# endif

    va_list args;
    va_start( args, message );

    char * buff = new char[ std::strlen(message) + 128 ];
    buff[0] = '\0';

# if LOGLEVEL >= 2

    std::sprintf( buff, __BLU "LOG: " __RESET "[%s:%d] "
        # if defined (ENABLE_MPI)
            "(%d/%d)  ",
        # else
            " ",
        # endif
            file, line
        # if defined (ENABLE_MPI)
            , MPI_rank, MPI_num_ranks
        # endif
            );

# else // if LOGLEVEL >= 2

    (void) file;
    (void) line;

# endif // if LOGLEVEL >= 2

    std::strcat( buff, message );

    if ( logfp != nullptr ) {

        va_list args2;
        va_copy( args2, args );
        std::vfprintf( logfp, buff, args2 );
        va_end( args2 );
    }

    std::vfprintf( stdout, buff, args );
    va_end( args );
    delete [] buff;
}


//------------------------------------------------------------------------------------------------------------
//! \brief  Prints a status message to stderr and #Global::logfp.
//!
//! Status messages are labeled in green (if the macro \c #COLOR_TERM is defined) and are always accompanied
//! by the filename and line number of the location where the status message was generated.
//!
//! \attention  The macro \c #PRINT_STATUS is provided as a wrapper for calling this function.
//!
//! \param[in]      file        Name of file where function is called from (usually just \c \_\_FILE__).
//! \param[in]      line        Line number where function is called from (usually just \c \_\_LINE__).
//! \param[in]      message     C-string containing message to print, possibly with format specifiers.
//! \param[in]      ...         Arguments specifying data to print
//!
//! \see    \c #PRINT_STATUS
//------------------------------------------------------------------------------------------------------------
void PrintStatus (

    const char * file,
    const int line,
    const char * message,
    ...
) {

    va_list args;
    va_start( args, message );

    char * buff = new char[ std::strlen(message) + 128 ];

    std::sprintf( buff, __GRN "STATUS: " __RESET "[%s:%d] "
        # if defined (ENABLE_MPI)
            "(%d/%d)  ",
        # else
            " ",
        # endif
            file, line
        # if defined (ENABLE_MPI)
            , MPI_rank, MPI_num_ranks
        # endif
            );

    std::strcat( buff, message );

    if ( logfp != nullptr ) {

        va_list args2;
        va_copy( args2, args );
        std::vfprintf( logfp, buff, args2 );
        va_end( args2 );
    }

    std::vfprintf( stderr, buff, args );
    va_end( args );
    delete [] buff;
}


//------------------------------------------------------------------------------------------------------------
//! \brief  Prints a message indicating a potential performance impact to stderr and #Global::logfp.
//!
//! Performance notes are labeled in yellow (if the macro \c #COLOR_TERM is defined) and are always
//! accompanied by the filename and line number of the location where the performance note was generated.
//!
//! \attention  The macro \c #PRINT_NOTE is provided as a wrapper for calling this function.
//!
//! \param[in]      file        Name of file where function is called from (usually just \c \_\_FILE__).
//! \param[in]      line        Line number where function is called from (usually just \c \_\_LINE__).
//! \param[in]      message     C-string containing message to print, possibly with format specifiers.
//! \param[in]      ...         Arguments specifying data to print
//!
//! \see    \c #PRINT_NOTE
//------------------------------------------------------------------------------------------------------------
void PrintNote (

    const char * file,
    const int line,
    const char * message,
    ...
) {

    va_list args;
    va_start( args, message );

    char * buff = new char[ std::strlen(message) + 128 ];

    std::sprintf( buff, __YEL "NOTE: " __RESET "[%s:%d] "
        # if defined (ENABLE_MPI)
            "(%d/%d)  ",
        # else
            " ",
        # endif
            file, line
        # if defined (ENABLE_MPI)
            , MPI_rank, MPI_num_ranks
        # endif
            );

    std::strcat( buff, message );

    if ( logfp != nullptr ) {

        va_list args2;
        va_copy( args2, args );
        std::vfprintf( logfp, buff, args2 );
        va_end( args2 );
    }

    std::vfprintf( stderr, buff, args );
    va_end( args );
    delete [] buff;
}
