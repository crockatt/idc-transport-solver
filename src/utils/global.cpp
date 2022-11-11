//------------------------------------------------------------------------------------------------------------
//! \file   utils/global.cpp
//! \brief  Implementation of global functions. Mostly just a couple basic routines.
//!
//! \author Michael Crockatt, Kris Garrett
//! \date   January 2018
//------------------------------------------------------------------------------------------------------------

# include <algorithm>
# include <cmath>
# include <cstdarg>
# include <cstdlib>
# include <cstring>

# if defined (ENABLE_MPI)
    # include <mpi.h>
# endif

//!
//! \brief  Trick to declare global variables only in this file and \c extern elsewhere.
//!
# define NO_EXTERN_GLOBAL

# include "utils/global.hpp"
# include "utils/Quadrule/Quadrule.hpp"


using namespace Global;
using namespace Quadrule;


namespace Global {

    //!
    //! \brief  Used to map from a (dimension, edge) index pair to the associated member of the BoundaryEdge
    //!         enum class.
    //!
    //! Used to determine whether certain edges should be included in operations using a dimension-independent
    //! implementation.
    //!
    extern const BoundaryEdge DimEdgeIdx_to_BoundaryEdge [SPACE_DIMS][2] =
        {
            { BoundaryEdge::X_Min, BoundaryEdge::X_Max }
        # if SPACE_DIMS >= 2
        ,   { BoundaryEdge::Y_Min, BoundaryEdge::Y_Max }
        # endif
        # if SPACE_DIMS == 3
        ,   { BoundaryEdge::Z_Min, BoundaryEdge::Z_Max }
        # endif
        };

    extern const int col_width = 20;
}


//------------------------------------------------------------------------------------------------------------
//! \brief  Allocates memory for and computes the full tensor of triple product integrals of Legendre
//!         polynomials.
//!
//! Triple product integrals are evaluated by Gaussian quadrature. Any values below a specified error
//! tolerance chosen to correspond with machine precision (currently 5e-16) are set to floating point zero.
//!
//! \attention  Values are stored at the global pointer #Global::TPI.
//!
//! \see    #Global::TPI
//------------------------------------------------------------------------------------------------------------
void ComputeTPI ( void ) {

    const double tol = 5e-16;


    // --- Allocate memory for storing tensor of product integrals. ----------------------------------- //

    const size_t sizeof_tensor = (DG_degree + 1)*(DG_degree + 1)*(DG_degree + 1) * sizeof(double);

    delete [] TPI;
    TPI = new double[ sizeof_tensor ];

    std::memset( TPI, 0, sizeof_tensor );


    // --- Compute a Gauss-Legendre quadrature rule which computes all integrals exactly. ------------- //

    const int64_t numGLnodes = std::max( 2 * DG_degree, int64_t(2) );

    double * const GLnodes   = new double[ numGLnodes ];
    double * const GLweights = new double[ numGLnodes ];

    ComputeQuadrature( numGLnodes, GLnodes, GLweights, NodesType::GaussLegendre );


    // --- Compute the coefficients of the tensor. ---------------------------------------------------- //

    for ( int64_t a = 0; a <= DG_degree; ++a ) {
    for ( int64_t b = 0; b <= DG_degree; ++b ) {
    for ( int64_t c = 0; c <= DG_degree; ++c ) {

        // Compute integral using quadrature rule.
        for ( int64_t i = 0; i < numGLnodes; ++i )
            TPI[ ITPI(a,b,c) ] += GLweights[i] * Legendre( a, GLnodes[i] )
                                               * Legendre( b, GLnodes[i] )
                                               * Legendre( c, GLnodes[i] );

        // Clean elements of tensor which are numerically zero.
        if ( std::abs( TPI[ ITPI(a,b,c) ] ) < tol ) {  TPI[ ITPI(a,b,c) ] = 0.0;  }

    }}}


    // --- Clean up. ---------------------------------------------------------------------------------- //

    delete [] GLnodes;
    delete [] GLweights;
}


//============================================================================================================
//=== ADDITIONAL SIMPLE ROUTINES =============================================================================
//============================================================================================================


//------------------------------------------------------------------------------------------------------------
//! \brief  Reads a null-terminated C-style string from a binary file at the given file pointer into an
//!         std::string object.
//!
//! \param[out]     str     String to store result in.
//! \param[in,out]  fp      File pointer from which to read the string.
//------------------------------------------------------------------------------------------------------------
void ReadStringFromFile (

    std::string & str,
# if defined (ENABLE_MPI)
    MPI_File & fp
# else
    std::FILE * & fp
# endif
) {

    str.clear();
    char c = '0';

    while (
            # if defined (ENABLE_MPI)
                MPI_SUCCESS == MPI_File_read( fp, &c, 1, MPI_CHAR, MPI_STATUS_IGNORE )
            # else
                std::fread( &c, sizeof(char), 1, fp )
            # endif
            &&  c != '\0'
    ) {
        str.push_back( c );
    }
}


//------------------------------------------------------------------------------------------------------------
//! \brief  Writes the contents of an std::string object to a binary file at the given file pointer as a
//!         null-terminated C-style string.
//!
//! \param[in]      str     String to write to disk.
//! \param[in,out]  fp      File pointer at which to write the string.
//------------------------------------------------------------------------------------------------------------
void WriteStringToFile (

    const std::string & str,
# if defined (ENABLE_MPI)
    MPI_File & fp
# else
    std::FILE * & fp
# endif
) {

    const char null_byte = '\0';

# if defined (ENABLE_MPI)

    for ( const char & c : str )
        MPI_File_write( fp, &c, 1, MPI_CHAR, MPI_STATUS_IGNORE );

    MPI_File_write( fp, &null_byte, 1, MPI_CHAR, MPI_STATUS_IGNORE );

# else // if defined (ENABLE_MPI)

    for ( const char & c : str )
        std::fwrite( &c, sizeof(char), 1, fp );

    std::fwrite( &null_byte, sizeof(char), 1, fp );

# endif // if defined (ENABLE_MPI)
}


//------------------------------------------------------------------------------------------------------------
//! \brief  Returns the degree of the DG approximation for the spatial variables based on the contents of the
//!         given ParameterList.
//!
//! \param[in]  plist       List of parameters to decide degree of approximation from.
//------------------------------------------------------------------------------------------------------------
int64_t GetDGDegreeX (

    const ParameterList & plist
) {

    int64_t DG_degree_x = -1;

    try {  DG_degree_x = plist.GetValue<int64_t>( "dg_degree"   );  } catch (...) {/* empty */}
    try {  DG_degree_x = plist.GetValue<int64_t>( "dg_degree_x" );  } catch (...) {/* empty */}

    if ( DG_degree_x < 0 ) {

        std::string error_message =   std::string()
                                    + "Failed to determine degree of DG spatial discretization from "
                                    + "ParameterList in "
                                    + std::string(__func__)
                                    + ".\n";

        throw std::runtime_error( error_message );
    }

    return DG_degree_x;
}


//------------------------------------------------------------------------------------------------------------
//! \brief  Returns the degree of the DG approximation in time based on the contents of the given
//!         ParameterList.
//!
//! \param[in]  plist       List of parameters to decide degree of approximation from.
//------------------------------------------------------------------------------------------------------------
int64_t GetDGDegreeT (

    const ParameterList & plist
) {

    int64_t DG_degree_t = -1;

    try {  DG_degree_t = plist.GetValue<int64_t>( "dg_degree"   );  } catch (...) {/* empty */}
    try {  DG_degree_t = plist.GetValue<int64_t>( "dg_degree_t" );  } catch (...) {/* empty */}

    if ( DG_degree_t < 0 ) {

        std::string error_message =   std::string()
                                    + "Failed to determine degree of DG temporal discretization from "
                                    + "ParameterList in "
                                    + std::string(__func__)
                                    + ".\n";

        throw std::runtime_error( error_message );
    }

    return DG_degree_t;
}
