//------------------------------------------------------------------------------------------------------------
//! \file   objects/DomainDecomposition.cpp
//! \brief  Implementation of DomainDecomposition class.
//!
//! \author Michael M. Crockatt
//! \date   April 2018
//------------------------------------------------------------------------------------------------------------

# include <algorithm>
# include <cinttypes>
# include <ios>
# include <vector>

# if defined (ENABLE_MPI)
    # include <mpi.h>
# endif

# include "objects/SpatialDomain.hpp"
# include "objects/DomainDecomposition.hpp"
# include "utils/CLog.hpp"


//============================================================================================================
//=== CONSTRUCTORS, DESTRUCTORS, AND ASSOCIATED HELPER ROUTINES ==============================================
//============================================================================================================


//------------------------------------------------------------------------------------------------------------
//! \brief  Constructs an empty DomainDecomposition object with zero-initialized parameters.
//------------------------------------------------------------------------------------------------------------
DomainDecomposition::DomainDecomposition ( void ) :

    global_domain {},
    MPI_num_blocks_arr {},
    MPI_block_coords_arr {},
    local_subdomains { nullptr }
{}


//------------------------------------------------------------------------------------------------------------
//! \brief  Constructs a DomainDecomposition object using the given values.
//!
//! \param[in]  global_domain_in    Contains the global mesh information.
//! \param[in]  decompose_mpi_ranks (optional) <br>
//!                                 Specifies whether or not a domain decomposition should be applied to the
//!                                 global mesh using the topology information stored in
//!                                 Global::MPI_cart_comm. Default value is true. The value false is used by
//!                                 MPI implementations of binary read/write routines to index into output
//!                                 files (which are structured without domain decomposition).
//------------------------------------------------------------------------------------------------------------
DomainDecomposition::DomainDecomposition (

    const SpatialDomain & global_domain_in,
    const bool decompose_mpi_ranks // = true
) :
    global_domain{ global_domain_in },
    MPI_num_blocks_arr {},
    MPI_block_coords_arr {},
    local_subdomains{ nullptr }
{
    // Set domain decomposition topology.
    SetTopology( decompose_mpi_ranks );

    // Compute information for subdomains of domain decomposition.
    ComputeSubdomains();
}


//------------------------------------------------------------------------------------------------------------
//! \brief  Constructs a DomainDecomposition object using values read from a ParameterList for the global
//!         mesh information.
//!
//! \param[in]  input_list          ParameterList to read global mesh information from.
//! \param[in]  decompose_mpi_ranks (optional) <br>
//!                                 Specifies whether or not a domain decomposition should be applied to the
//!                                 global mesh using the topology information stored in
//!                                 Global::MPI_cart_comm. Default value is true. The value false is used by
//!                                 MPI implementations of binary read/write routines to index into output
//!                                 files (which are structured without domain decomposition).
//------------------------------------------------------------------------------------------------------------
DomainDecomposition::DomainDecomposition (

    const ParameterList & input_list,
    const bool decompose_mpi_ranks // = true
) :
    // Construct SpatialDomain object for global mesh information and delegate construction.
    DomainDecomposition( SpatialDomain( input_list ), decompose_mpi_ranks )
{}


//------------------------------------------------------------------------------------------------------------
//! \brief  Sets values for DomainDecomposition::MPI_num_blocks_arr and
//!         DomainDecomposition::MPI_block_coords_arr from the Cartesian communicator Global::MPI_cart_comm
//!         if the input argument is true.
//!
//! \param[in]  decompose_mpi_ranks Specifies whether or not a domain decomposition should be applied to the
//!                                 global mesh using the topology information stored in
//!                                 Global::MPI_cart_comm.
//------------------------------------------------------------------------------------------------------------
void DomainDecomposition::SetTopology (

    const bool decompose_mpi_ranks
) {

    // Initialize default values for sequential (non-MPI) implementation.
    for ( int64_t dim = 0; dim < SPACE_DIMS; ++dim ) {

        this->MPI_num_blocks_arr[dim] = 1;
        this->MPI_block_coords_arr[dim] = 0;
    }

# if defined (ENABLE_MPI)

    if ( decompose_mpi_ranks ) {

        int periods [SPACE_DIMS];

        // Retrieve topology information from MPI communicator for MPI implementations.
        int MPI_err = MPI_Cart_get( Global::MPI_cart_comm, SPACE_DIMS,
                                    this->MPI_num_blocks_arr, periods, this->MPI_block_coords_arr );

        if ( MPI_err ) {

            char err_str [MPI_MAX_ERROR_STRING];
            int err_len;

            MPI_Error_string( MPI_err, err_str, &err_len );
            std::string error_string = std::string( err_str );
            std::replace( error_string.begin(), error_string.end(), '\n', ' ' );

            PRINT_ERROR( "MPI_Cart_get returned error '%s'.\nAborting run.\n", error_string.c_str() )
            MPI_Abort( MPI_COMM_WORLD, EXIT_FAILURE );
        }
    }

# else // if defined (ENABLE_MPI)

    // Avoid compiler warning if MPI is not enabled.
    (void) decompose_mpi_ranks;

# endif // if defined (ENABLE_MPI)

}


//------------------------------------------------------------------------------------------------------------
//! \brief  Copy constructor for DomainDecomposition class.
//!
//! \param[in]  that    Object containing global mesh information to use for construction.
//------------------------------------------------------------------------------------------------------------
DomainDecomposition::DomainDecomposition (

    const DomainDecomposition & that
) :
    DomainDecomposition( that.global_domain )
{}


//------------------------------------------------------------------------------------------------------------
//! \brief  Allocates memory and computes subdomains based on global parameters.
//------------------------------------------------------------------------------------------------------------
void DomainDecomposition::ComputeSubdomains ( void ) {

    int64_t num_subdomains = 1;

    for ( int64_t dim = 0; dim < SPACE_DIMS; ++dim )
        num_subdomains *= MPI_num_blocks_arr[dim];

    local_subdomains = new SpatialDomain[ num_subdomains ];

    // --- Compute bounds and number of cells across each dimension. ---------------------------------- //

    const int64_t nx_per_block [] =
        {       global_domain.nx(0) / MPI_num_blocks(0)
        # if SPACE_DIMS >= 2
            ,   global_domain.nx(1) / MPI_num_blocks(1)
        # endif
        # if SPACE_DIMS == 3
            ,   global_domain.nx(2) / MPI_num_blocks(2)
        # endif
        };

    const int64_t nx_leftover [] =
        {       global_domain.nx(0) % MPI_num_blocks(0)
        # if SPACE_DIMS >= 2
            ,   global_domain.nx(1) % MPI_num_blocks(1)
        # endif
        # if SPACE_DIMS == 3
            ,   global_domain.nx(2) % MPI_num_blocks(2)
        # endif
        };

    std::vector<int64_t> nx [SPACE_DIMS];
    std::vector<double>  ax [SPACE_DIMS];
    std::vector<double>  bx [SPACE_DIMS];

    for ( int64_t dim = 0; dim < SPACE_DIMS; ++dim ) {

        for ( int64_t blk_idx = 0; blk_idx < MPI_num_blocks(dim); ++blk_idx ) {

            nx[dim].push_back( nx_per_block[dim] + ( nx_leftover[dim] > blk_idx ? 1 : 0 ) );

            if ( blk_idx == 0 )
                ax[dim].push_back( global_domain.ax(dim) );
            else
                ax[dim].push_back( bx[dim].back() );

            bx[dim].push_back( ax[dim].back() + nx[dim].back() * global_domain.dx(dim) );
        }
    }

    // --- Initialize subdomain objects using tensor product. ----------------------------------------- //

# if SPACE_DIMS == 1

    for ( int64_t i = 0; i < MPI_num_blocks(0); ++i ) {

        Subdomain_NonConst(i) = SpatialDomain( { nx[0].at(i) }, { ax[0].at(i) }, { bx[0].at(i) } );
    }

# elif SPACE_DIMS == 2

    for ( int64_t i = 0; i < MPI_num_blocks(0); ++i ) {
    for ( int64_t j = 0; j < MPI_num_blocks(1); ++j ) {

        Subdomain_NonConst(i,j) = SpatialDomain( { nx[0].at(i), nx[1].at(j) },
                                                 { ax[0].at(i), ax[1].at(j) },
                                                 { bx[0].at(i), bx[1].at(j) } );
    }}

# elif SPACE_DIMS == 3

    for ( int64_t i = 0; i < MPI_num_blocks(0); ++i ) {
    for ( int64_t j = 0; j < MPI_num_blocks(1); ++j ) {
    for ( int64_t k = 0; k < MPI_num_blocks(2); ++k ) {

        Subdomain_NonConst(i,j,k) = SpatialDomain( { nx[0].at(i), nx[1].at(j), nx[2].at(k) },
                                                   { ax[0].at(i), ax[1].at(j), ax[2].at(k) },
                                                   { bx[0].at(i), bx[1].at(j), bx[2].at(k) } );
    }}}

# endif // if SPACE_DIMS == ?
}


//------------------------------------------------------------------------------------------------------------
//! \brief  Destructor for DomainDecomposition class.
//------------------------------------------------------------------------------------------------------------
DomainDecomposition::~DomainDecomposition ( void ) {

    delete [] local_subdomains;
    local_subdomains = nullptr;
}


//------------------------------------------------------------------------------------------------------------
//! \brief  Determines offsets for handling of boundary and halo cells in operations.
//!
//! Offsets indicate whether the ghost cells along the specified local domain edge are included (1) or not
//! (0) in the operation.
//!
//! The first index specifies spatial dimension.
//! The second index specifies either:
//!     - index 0:  offset for lower bound of loop
//!     - index 1:  offset for upper bound of loop
//! in the given spatial dimension.
//!
//! \param[in]  op_domain   Domain over which operation is to be applied.
//! \param[out] offset      Array containing loop offsets.
//------------------------------------------------------------------------------------------------------------
void DomainDecomposition::DetermineOffsets (

    const OpDomain op_domain,
    int64_t (& offset) [SPACE_DIMS][2]

) const {

    /*
     *  dim = spatial dimension: { { 0, X }, { 1, Y }, { 2, Z } }.
     *  bound = boundary edge: { { 0, lower }, { 1, upper } }.
     */
    for ( int64_t dim = 0; dim < SPACE_DIMS; ++dim ) {
    for ( int64_t bound = 0; bound <= 1; ++bound ) {

        if (/*  Include edge in operation if one of the following conditions is satisfied:
             *
             *  1.  (a) Boundaries ARE periodic; AND
             *      (b) interior halo cells ARE included in the operation.
             */
                (    Global::periodic
                  && BitmaskHasAll( op_domain, OpDomain::Halo )
                )
            /*
             *  2.  (a) The given boundary IS reflecting; AND
             *      (b) interior halo cells ARE included in the operation.
             */
             || (    BitmaskHasAll( Global::reflecting_boundaries,
                                    Global::DimEdgeIdx_to_BoundaryEdge[dim][bound] )
                  && BitmaskHasAll( op_domain, OpDomain::Halo )
                )
            /*
             *  3.  (a) Boundaries ARE NOT periodic; AND
             *      (b) the given boundary IS NOT reflecting; AND
             *      (c) the halo cells of the given edge for the current MPI rank DO lie on the physical
             *          domain boundary (NOTE: this is always true if MPI is not enabled); AND
             *      (d) boundary halo cells ARE included in the operation.
             */
             || (    !Global::periodic
                  && !BitmaskHasAny( Global::reflecting_boundaries,
                                     Global::DimEdgeIdx_to_BoundaryEdge[dim][bound] )
            # if defined (ENABLE_MPI)
                  && this->MPI_block_coords(dim) ==
                        (   bound == 0
                          ? 0                               // Lower bound.
                          : this->MPI_num_blocks(dim) - 1   // Upper bound.
                        )
            # endif
                  && BitmaskHasAll( op_domain, OpDomain::Boundary )
                )
            # if defined (ENABLE_MPI)
            /*
             *  4.  (a) Boundaries ARE NOT periodic; AND
             *      (b) the given boundary IS NOT reflecting; AND
             *      (c) the halo cells of the given edge for the current MPI rank DO NOT lie on the physical
             *          domain boundary; AND
             *      (d) interior halo cells ARE included in the operation.
             */
             || (    !Global::periodic
                  && !BitmaskHasAny( Global::reflecting_boundaries,
                                     Global::DimEdgeIdx_to_BoundaryEdge[dim][bound] )
                  && this->MPI_block_coords(dim) !=
                        (   bound == 0
                          ? 0                               // Lower bound.
                          : this->MPI_num_blocks(dim) - 1   // Upper bound.
                        )
                  && BitmaskHasAll( op_domain, OpDomain::Halo )
                )
        # endif // if defined (ENABLE_MPI)
        )    {  offset[dim][bound] = 1;  }
        else {  offset[dim][bound] = 0;  }
    }}
}


//============================================================================================================
//=== OPERATOR OVERLOADS =====================================================================================
//============================================================================================================


//------------------------------------------------------------------------------------------------------------
//! \brief  Copy assignment operator for DomainDecomposition class.
//------------------------------------------------------------------------------------------------------------
DomainDecomposition & DomainDecomposition::operator= (

    const DomainDecomposition & that
) {

    if ( AreMatching( that, DomainDecomposition{} ) ) {

        this->Zero();

    } else {

        delete [] this->local_subdomains;
        this->local_subdomains = nullptr;

        this->global_domain = that.global_domain;

        for ( int64_t dim = 0; dim < SPACE_DIMS; ++dim ) {

            this->MPI_num_blocks_arr[dim] = that.MPI_num_blocks(dim);
            this->MPI_block_coords_arr[dim] = that.MPI_block_coords(dim);
        }

        this->ComputeSubdomains();
    }

    return *this;
}


//============================================================================================================
//=== STATIC MEMBER FUNCTIONS ================================================================================
//============================================================================================================


//------------------------------------------------------------------------------------------------------------
//! \brief  Two DomainDecomposition objects "are matching" if all of their parameters are equal.
//!
//! \param[in]      first       The first of the two DomainDecomposition objects to compare.
//! \param[in]      second      The second of the two DomainDecomposition objects to compare.
//!
//! \return     Returns true if the DomainDecomposition objects "are matching" and false otherwise.
//------------------------------------------------------------------------------------------------------------
bool DomainDecomposition::AreMatching (

    const DomainDecomposition & first,
    const DomainDecomposition & second
) {

    bool result = SpatialDomain::AreMatching( first.global_domain, second.global_domain );

    for ( int i = 0; i < SPACE_DIMS; ++i ) {

        result &= ( first.MPI_num_blocks(i) == second.MPI_num_blocks(i) );
        result &= ( first.MPI_block_coords(i) == second.MPI_block_coords(i) );
    }

    return result;
}


//============================================================================================================
//=== ADDITIONAL MEMBER FUNCTIONS ============================================================================
//============================================================================================================


//------------------------------------------------------------------------------------------------------------
//! \brief  Zeros all parameters of a DomainDecomposition object.
//------------------------------------------------------------------------------------------------------------
DomainDecomposition & DomainDecomposition::Zero ( void ) {

    delete [] this->local_subdomains;
    this->local_subdomains = nullptr;

    this->global_domain = SpatialDomain{};

    for ( int64_t dim = 0; dim < SPACE_DIMS; ++dim ) {

        this->MPI_num_blocks_arr[dim] = 0;
        this->MPI_block_coords_arr[dim] = 0;
    }

    return *this;
}


//------------------------------------------------------------------------------------------------------------
//! \brief  Prints a summary of the object configuration to the logging interface.
//!
//! \param[in]  prefix  String prepended to each line of output, intended to be used for controlling
//!                     indentation of output.
//------------------------------------------------------------------------------------------------------------
void DomainDecomposition::Print (

    const std::string & prefix // = "  "

) const {

    PRINT_STATUS( "Executing DomainDecomposition::%s.\n", __func__ )

    PRINT_LOG( "\n" )

# if defined (ENABLE_MPI)

    if ( Global::MPI_rank == 0 ) {

    # if SPACE_DIMS == 1

        PRINT_LOG( "%s%-*s  [ %d ]\n", prefix.c_str(), Global::col_width, "MPI subdomains:", MPI_num_blocks(0) )

    # elif SPACE_DIMS == 2

        PRINT_LOG( "%s%-*s  [ %d, %d ]\n", prefix.c_str(), Global::col_width, "MPI subdomains:",
                   MPI_num_blocks(0), MPI_num_blocks(1) )

    # elif SPACE_DIMS == 3

        PRINT_LOG( "%s%-*s  [ %d, %d, %d ]\n", prefix.c_str(), Global::col_width, "MPI subdomains:",
                   MPI_num_blocks(0), MPI_num_blocks(1), MPI_num_blocks(2) )

    # endif // if SPACE_DIMS == ?

        PRINT_LOG( "\n%sGlobal Domain:\n", prefix.c_str() )

# endif // if defined (ENABLE_MPI)

        PRINT_LOG( "%s%s\n\n", prefix.c_str(), global_domain.PrintString().c_str() )

# if defined (ENABLE_MPI)

    }

    MPI_Barrier( Global::MPI_cart_comm );

# if SPACE_DIMS == 1

    PRINT_LOG( "%s(%d/%d)  %s  [ %d ]\n", prefix.c_str(),
               Global::MPI_rank, Global::MPI_num_ranks,
               "MPI Coords:",
               MPI_block_coords(0) )

# elif SPACE_DIMS == 2

    PRINT_LOG( "%s(%d/%d)  %s  [ %d, %d ]\n", prefix.c_str(),
               Global::MPI_rank, Global::MPI_num_ranks,
               "MPI Coords:",
               MPI_block_coords(0), MPI_block_coords(1) )

# elif SPACE_DIMS == 3

    PRINT_LOG( "%s(%d/%d)  %s  [ %d, %d, %d ]\n", prefix.c_str(),
               Global::MPI_rank, Global::MPI_num_ranks,
               "MPI Coords:",
               MPI_block_coords(0), MPI_block_coords(1), MPI_block_coords(2) )

# endif // if SPACE_DIMS == ?

    MPI_Barrier( Global::MPI_cart_comm );

    if ( Global::MPI_rank == 0 ) {

        PRINT_LOG( "\n" )
        PRINT_LOG( "%sLocal domains:\n", prefix.c_str() )

    # if SPACE_DIMS == 1

        for ( int64_t i = 0; i < MPI_num_blocks(0); ++i ) {

            PRINT_LOG( "%s%s\n", prefix.c_str(), Subdomain(i).PrintString().c_str() )
        }

        PRINT_LOG( "\n" )

    # elif SPACE_DIMS == 2

        for ( int64_t i = 0; i < MPI_num_blocks(0); ++i ) {
        for ( int64_t j = 0; j < MPI_num_blocks(1); ++j ) {

            PRINT_LOG( "%s%s\n", prefix.c_str(), Subdomain(i,j).PrintString().c_str() )
        }
            PRINT_LOG( "\n" )
        }

    # elif SPACE_DIMS == 3
        # warning "DomainDecomposition::Print not implemented for 3 spatial dimensions."
    # endif
    }

# endif // if defined (ENABLE_MPI)
}


//------------------------------------------------------------------------------------------------------------
//! \brief  Computes the indices (block and local) associated with the spatial cell in which the given point
//!         lies.
//!
//! \param[in]  point           Point to locate cell for.
//! \param[out] block_coords    The indices of the subdomain in which the point lies.
//! \param[out] local_coords    The indices of the spatial cell local to the subdomain in which the point
//!                             lies.
//------------------------------------------------------------------------------------------------------------
void DomainDecomposition::PointToCell (

    const double (& point) [SPACE_DIMS],
    int (& block_coords) [SPACE_DIMS],
    int (& local_coords) [SPACE_DIMS]

) const {

    for ( int64_t dim = 0; dim < SPACE_DIMS; ++dim ) {  block_coords[dim] = 0;  }

    for ( int64_t dim = 0; dim < SPACE_DIMS; ++dim ) {

        for ( block_coords[dim] = 0;
              point[dim] > this->Subdomain( block_coords ).bx(dim);
              ++block_coords[dim]
        ) {/* empty */}

        local_coords[dim] = (point[dim] - this->Subdomain(block_coords).ax(dim)) / this->dx(dim);
    }
}


//------------------------------------------------------------------------------------------------------------
//! \brief  Computes the indices (block and local) associated with the spatial cell in which the given point
//!         lies.
//!
//! \param[in]  point           Point to locate cell for.
//! \param[out] block_coords    The indices of the subdomain in which the point lies.
//! \param[out] local_coords    The indices of the spatial cell local to the subdomain in which the point
//!                             lies.
//------------------------------------------------------------------------------------------------------------
void DomainDecomposition::PointToCell (

    std::initializer_list<double> point,
    int (& block_coords) [SPACE_DIMS],
    int (& local_coords) [SPACE_DIMS]

) const {

# if defined (STRICT_CHECK)

    if ( point.size() != SPACE_DIMS ) {

        std::string error_message = "std::initializer_list<int> nx has wrong size in '"
                                    + std::string(__func__)
                                    + "': expected '"
                                    + std::to_string(SPACE_DIMS)
                                    + "' got '"
                                    + std::to_string( point.size() )
                                    + "'.";

        PRINT_ERROR( error_message.c_str() )
        throw std::out_of_range( error_message );
    }

# endif // if defined (STRICT_CHECK)

    double point_arr [SPACE_DIMS] = {};
    int64_t dim = 0;

    for ( auto & val : point ) {  point_arr[dim] = val;  ++dim;  }

    PointToCell( point_arr, block_coords, local_coords );
}


//------------------------------------------------------------------------------------------------------------
//! \brief  Constructs a DomainDecomposition object form the values stored at the given file pointer using the
//!         provided subdomain decomposition parameters.
//!
//! \param[in]  fp                  File pointer at which to begin reading object from.
//! \param[in]  decompose_mpi_ranks (optional) <br>
//!                                 Specifies whether or not a domain decomposition should be applied to the
//!                                 global mesh using the topology information stored in
//!                                 Global::MPI_cart_comm. Default value is true. The value false is used by
//!                                 MPI implementations of binary read/write routines to index into output
//!                                 files (which are structured without domain decomposition).
//! \param[in]  version             Integer for file format to use. This option in here only to enable
//!                                 conversion of old output files to the current format. Generally the
//!                                 default value should be used.
//!
//! \see    DomainDecomposition::WriteToDisk()
//------------------------------------------------------------------------------------------------------------
void DomainDecomposition::ReadFromDisk (

# if defined (ENABLE_MPI)
    MPI_File & fp,
# else
    std::FILE * & fp,
# endif

    const bool decompose_mpi_ranks, // = true
    const int64_t version // = 0
) {

    delete [] this->local_subdomains;
    this->local_subdomains = nullptr;

    // Read global mesh information from disk.
    this->global_domain.ReadFromDisk( fp, version );

    // Set domain decomposition topology information.
    SetTopology( decompose_mpi_ranks );

    // Compute information for subdomains of domain decomposition.
    ComputeSubdomains();
}


//------------------------------------------------------------------------------------------------------------
//! \brief  Writes the parameters in a DomainDecomposition object to the given file pointer.
//!
//! \param[in]      fp      File pointer at which to write object to.
//!
//! \see    DomainDecomposition::ReadFromDisk()
//------------------------------------------------------------------------------------------------------------
void DomainDecomposition::WriteToDisk (

# if defined (ENABLE_MPI)
    MPI_File & fp
# else
    std::FILE * & fp
# endif

) const {

    global_domain.WriteToDisk( fp );
}
