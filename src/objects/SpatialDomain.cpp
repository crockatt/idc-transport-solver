//------------------------------------------------------------------------------------------------------------
//! \file   objects/SpatialDomain.cpp
//! \brief  Implementation of SpatialDomain class.
//!
//! \author Michael M. Crockatt
//! \date   November 2017
//------------------------------------------------------------------------------------------------------------


# include <cinttypes>
# include <initializer_list>
# include <iomanip>
# include <sstream>

# if defined (ENABLE_MPI)
    # include <mpi.h>
# endif

# include "objects/SpatialDomain.hpp"
# include "utils/CLog.hpp"


//============================================================================================================
//=== CONSTRUCTORS AND DESTRUCTOR ============================================================================
//============================================================================================================


//------------------------------------------------------------------------------------------------------------
//! \brief  Constructs a SpatialDomain with default (and generally physically invalid) values.
//------------------------------------------------------------------------------------------------------------
SpatialDomain::SpatialDomain ( void ) :

    nx_var {},
    ax_var {},
    bx_var {}
{
    for ( int dim = 0; dim < SPACE_DIMS; ++dim )
        this->dx_var[dim] = (this->bx_var[dim] - this->ax_var[dim]) / this->nx_var[dim];
}


//------------------------------------------------------------------------------------------------------------
//! \brief  Creates a SpatialDomain object using the provided values.
//------------------------------------------------------------------------------------------------------------
SpatialDomain::SpatialDomain (

    const int64_t (& nx) [SPACE_DIMS],
    const double  (& ax) [SPACE_DIMS],
    const double  (& bx) [SPACE_DIMS]
) {

    for ( int64_t dim = 0; dim < SPACE_DIMS; ++dim ) {

        this->nx_var[dim] = nx[dim];
        this->ax_var[dim] = ax[dim];
        this->bx_var[dim] = bx[dim];

        this->dx_var[dim] = (this->bx_var[dim] - this->ax_var[dim]) / this->nx_var[dim];
    }
}


//------------------------------------------------------------------------------------------------------------
//! \brief  Creates a SpatialDomain object using the provided values.
//------------------------------------------------------------------------------------------------------------
SpatialDomain::SpatialDomain (

    std::initializer_list<int64_t> nx,
    std::initializer_list<double>  ax,
    std::initializer_list<double>  bx
) {

# if defined (STRICT_CHECK)

    if ( nx.size() != SPACE_DIMS ) {

        std::string error_message = "std::initializer_list<int64_t> nx has wrong size in '"
                                    + std::string(__func__)
                                    + "': expected '"
                                    + std::to_string(SPACE_DIMS)
                                    + "' got '"
                                    + std::to_string( nx.size() )
                                    + "'.";

        PRINT_ERROR( error_message.c_str() )
        throw std::out_of_range( error_message );
    }

    if ( ax.size() != SPACE_DIMS ) {

        std::string error_message = "std::initializer_list<int64_t> ax has wrong size in '"
                                    + std::string(__func__)
                                    + "': expected '"
                                    + std::to_string(SPACE_DIMS)
                                    + "' got '"
                                    + std::to_string( nx.size() )
                                    + "'.";

        PRINT_ERROR( error_message.c_str() )
        throw std::out_of_range( error_message );
    }

    if ( bx.size() != SPACE_DIMS ) {

        std::string error_message = "std::initializer_list<int64_t> bx has wrong size in '"
                                    + std::string(__func__)
                                    + "': expected '"
                                    + std::to_string(SPACE_DIMS)
                                    + "' got '"
                                    + std::to_string( nx.size() )
                                    + "'.";

        PRINT_ERROR( error_message.c_str() )
        throw std::out_of_range( error_message );
    }

# endif // if defined (STRICT_CHECK)

    int64_t dim = 0;

    for ( auto & val : nx ) {  this->nx_var[dim] = val;  ++dim;  }      dim = 0;
    for ( auto & val : ax ) {  this->ax_var[dim] = val;  ++dim;  }      dim = 0;
    for ( auto & val : bx ) {  this->bx_var[dim] = val;  ++dim;  }

    for ( dim = 0; dim < SPACE_DIMS; ++dim )
        this->dx_var[dim] = (this->bx_var[dim] - this->ax_var[dim]) / this->nx_var[dim];
}


//------------------------------------------------------------------------------------------------------------
//! \brief  Attempts to construct a SpatialDomain object using values read from a ParameterList.
//!
//! An exception is thrown if the required values cannot be found in the provided ParameterList.
//------------------------------------------------------------------------------------------------------------
SpatialDomain::SpatialDomain (

    const ParameterList & input_list
) {

    // Read mesh information for X dimension.
    input_list.GetValue( "nx", this->nx_var[0] );
    input_list.GetValue( "ax", this->ax_var[0] );
    input_list.GetValue( "bx", this->bx_var[0] );

# if SPACE_DIMS >= 2

    // Read mesh information for Y dimension.
    input_list.GetValue( "ny", this->nx_var[1] );
    input_list.GetValue( "ay", this->ax_var[1] );
    input_list.GetValue( "by", this->bx_var[1] );

# endif // if SPACE_DIMS >= 2
# if SPACE_DIMS == 3

    // Read mesh information for Z dimension.
    input_list.GetValue( "nz", this->nx_var[2] );
    input_list.GetValue( "az", this->ax_var[2] );
    input_list.GetValue( "bz", this->bx_var[2] );

# endif // if SPACE_DIMS == 3

    // Compute values for mesh sizes.
    for ( int64_t dim = 0; dim < SPACE_DIMS; ++dim )
        this->dx_var[dim] = (this->bx_var[dim] - this->ax_var[dim]) / this->nx_var[dim];
}


//============================================================================================================
//=== STATIC MEMBER FUNCTIONS ================================================================================
//============================================================================================================


//------------------------------------------------------------------------------------------------------------
//! \brief  Two SpatialDomain objects "are matching" if all of their parameters are equal.
//!
//! \param[in]      first       The first of the two SpatialDomain objects to compare.
//! \param[in]      second      The second of the two SpatialDomain objects to compare.
//!
//! \return     Returns true if the SpatialDomain objects "are matching" and false otherwise.
//------------------------------------------------------------------------------------------------------------
bool SpatialDomain::AreMatching (

    const SpatialDomain & first,
    const SpatialDomain & second
) {

    bool result = true;

    for ( int i = 0; i < SPACE_DIMS; ++i ) {

        result &= ( first.nx(i) == second.nx(i) );
        result &= ( first.ax(i) == second.ax(i) );
        result &= ( first.bx(i) == second.bx(i) );
        result &= ( first.dx(i) == second.dx(i) );
    }

    return result;
}


//============================================================================================================
//=== ADDITIONAL MEMBER FUNCTIONS ============================================================================
//============================================================================================================


//------------------------------------------------------------------------------------------------------------
//! \brief  Returns a string describing the values stored in the object for printing.
//!
//! \return Returns a string describing the values stored in the object for printing.
//------------------------------------------------------------------------------------------------------------
std::string SpatialDomain::PrintString ( void ) const {

    std::stringstream print_stream;

    for ( int64_t dim = 0; dim < SPACE_DIMS; ++dim ) {

        print_stream << "Dim " << dim << ": "
                     << std::setw(7) << std::right
                     << nx(dim)

                     << std::scientific << std::setprecision(4) << std::right
                     << " ["
                     << std::setw(11) << ax(dim)
                     << ","
                     << std::setw(11) << bx(dim)
                     << "] "

                     << std::setw(11) << dx(dim)

                     << "  \t";
    }

    return print_stream.str();
}


//------------------------------------------------------------------------------------------------------------
//! \brief  Reads parameters for a SpatialDomain object from the given file pointer.
//!
//! \param[in]  fp          File pointer at which to begin reading object from.
//! \param[in]  version     Integer for file format to use. This option in here only to enable conversion of
//!                         old output files to the current format. Generally the default value should be
//!                         used.
//!
//! \see    SpatialDomain::WriteToDisk()
//------------------------------------------------------------------------------------------------------------
void SpatialDomain::ReadFromDisk (

# if defined (ENABLE_MPI)
    MPI_File & fp,
# else
    std::FILE * & fp,
# endif
    const int64_t version // = 0
) {

# if defined (ENABLE_MPI)

    //!
    //! \brief   Macro to check return values from \c MPI_File_read for errors.
    //!
    # define CHECK_READ \
        if ( ret != MPI_SUCCESS ) {  goto read_error;  }

    int ret;

# else // if defined (ENABLE_MPI)

    //!
    //! \brief  Macro to check return values from \c std::fread for errors.
    //!
    # define CHECK_READ \
        if ( !ret ) {  goto read_error;  }

    //!
    //! \brief  Macro to check return values from \c std::fseek for errors.
    //!
    # define CHECK_SEEK \
        if ( ret ) {  goto read_error;  }

    size_t ret;

# endif // if defined (ENABLE_MPI)

    for ( int64_t dim = 0; dim < SPACE_DIMS; ++dim ) {

        switch ( version ) {

        default:
        /*  Current file format version.
         *
         *  For each dimension, read in order:
         *
         *      Variable        Type            Count
         *     -------------   -------------   -------
         *      nx_var          int64_t         1
         *      ax_var          double          1
         *      bx_var          double          1
         */
        {
        # if defined (ENABLE_MPI)

            ret = MPI_File_read( fp, &this->nx_var[dim], 1, MPI_INT64_T, MPI_STATUS_IGNORE );  CHECK_READ
            ret = MPI_File_read( fp, &this->ax_var[dim], 1, MPI_DOUBLE,  MPI_STATUS_IGNORE );  CHECK_READ
            ret = MPI_File_read( fp, &this->bx_var[dim], 1, MPI_DOUBLE,  MPI_STATUS_IGNORE );  CHECK_READ

        # else // if defined (ENABLE_MPI)

            ret = std::fread( &this->nx_var[dim], sizeof(int64_t), 1, fp );  CHECK_READ
            ret = std::fread( &this->ax_var[dim], sizeof(double),  1, fp );  CHECK_READ
            ret = std::fread( &this->bx_var[dim], sizeof(double),  1, fp );  CHECK_READ

        # endif // if defined (ENABLE_MPI)
        } break;

        case 1:
        /*  Reads old version compatible with commit da940f016bcf76717b4a058fe354e77d0b937a8c (and similar).
         *
         *  This format predates changes to the DomainDecomposition implementation that introduced the
         *  SpatialDomain class. Extra information that is no longer stored in the new file format must be
         *  skipped over.
         *
         *  For each dimension, read in order:
         *
         *      Variable        Type            Count
         *     -------------   -------------   -------
         *      nx_var          int64_t         1
         *      ax_var          double          1
         *      bx_var          double          1
         *
         *  Then skip in order:
         *
         *      Variable                    Type            Count
         *     -------------------------   -------------   -------
         *      dx_var                      double          1
         *      MPI_num_blocks_var          int             1
         *      MPI_block_coords_var        int             1
         *      MPI_ax_var                  double          1
         *      MPI_bx_var                  double          1
         *
         *  Total bytes skipped per dimension:  32 Bytes
         */
        {
        # if defined (ENABLE_MPI)

            ret = MPI_File_read( fp, &this->nx_var[dim], 1, MPI_INT64_T, MPI_STATUS_IGNORE );  CHECK_READ
            ret = MPI_File_read( fp, &this->ax_var[dim], 1, MPI_DOUBLE,  MPI_STATUS_IGNORE );  CHECK_READ
            ret = MPI_File_read( fp, &this->bx_var[dim], 1, MPI_DOUBLE,  MPI_STATUS_IGNORE );  CHECK_READ

            ret = MPI_File_seek( fp, 32, MPI_SEEK_CUR );    CHECK_READ

        # else // if defined (ENABLE_MPI)

            ret = std::fread( &this->nx_var[dim], sizeof(int64_t), 1, fp );  CHECK_READ
            ret = std::fread( &this->ax_var[dim], sizeof(double),  1, fp );  CHECK_READ
            ret = std::fread( &this->bx_var[dim], sizeof(double),  1, fp );  CHECK_READ

            ret = std::fseek( fp, 32, SEEK_CUR );   CHECK_SEEK

        # endif // if defined (ENABLE_MPI)
        } break;

        } // end switch.

        this->dx_var[dim] = (this->bx_var[dim] - this->ax_var[dim]) / this->nx_var[dim];
    }

    return;


# undef CHECK_READ

read_error:
    throw std::ios_base::failure( "Failed to read SpatialDomain object from file pointer.\n" );
}


//------------------------------------------------------------------------------------------------------------
//! \brief  Writes the parameters of a SpatialDomain object to the given file pointer.
//!
//! \param[in]      fp      File pointer at which to write object to.
//!
//! \see    SpatialDomain::ReadFromDisk()
//------------------------------------------------------------------------------------------------------------
void SpatialDomain::WriteToDisk (

# if defined (ENABLE_MPI)
    MPI_File & fp
# else
    std::FILE * & fp
# endif

) const {

    for ( int64_t dim = 0; dim < SPACE_DIMS; ++dim ) {

    # if defined (ENABLE_MPI)

        MPI_File_write( fp, &this->nx_var[dim], 1, MPI_INT64_T, MPI_STATUS_IGNORE );
        MPI_File_write( fp, &this->ax_var[dim], 1, MPI_DOUBLE,  MPI_STATUS_IGNORE );
        MPI_File_write( fp, &this->bx_var[dim], 1, MPI_DOUBLE,  MPI_STATUS_IGNORE );

    # else // if defined (ENABLE_MPI)

        std::fwrite( &this->nx_var[dim], sizeof(int64_t), 1, fp );
        std::fwrite( &this->ax_var[dim], sizeof(double),  1, fp );
        std::fwrite( &this->bx_var[dim], sizeof(double),  1, fp );

    # endif // if defined (ENABLE_MPI)
    }
}
