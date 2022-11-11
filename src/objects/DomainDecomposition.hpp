//------------------------------------------------------------------------------------------------------------
//! \file   objects/DomainDecomposition.hpp
//! \brief  Header for DomainDecomposition class.
//!
//! \author Michael M. Crockatt
//! \date   February 2018
//------------------------------------------------------------------------------------------------------------

# ifndef __DOMAIN_DECOMPOSITION_HPP__
# define __DOMAIN_DECOMPOSITION_HPP__


# include <cstddef>
# include <cstdint>
# include <cstdio>
# include <initializer_list>

# include "objects/SpatialDomain.hpp"
# include "utils/bitmask_operators.hpp"
# include "utils/CLog.hpp"
# include "utils/global.hpp"


// Forward declarations for friend init.
namespace RKDG {
    class OrdinateFlux;
    class CrossSection;
}


//------------------------------------------------------------------------------------------------------------
//! \brief  Used to specify where operations are to be performed on objects using DomainDecomposition mesh
//!         descriptions.
//------------------------------------------------------------------------------------------------------------
enum class OpDomain {

    //!
    //! Specifies that the operation should apply to all cells on the interior of the local subdomain.
    //!
    Interior = 0x01,

    //!
    //! Specifies that the operation should apply to all cells on the boundary of the physical domain where,
    //! e.g., an inflow condition is defined. This does _not_ include halo cells from other subdomains or
    //! halo cells that lie along periodic boundary edges.
    //!
    Boundary = 0x02,

    //!
    //! Specifies that the operation should apply to all halo cells from other subdomains. This does _not_
    //! include halo cells on the boundary of the physical domain (i.e., inflow boundaries), but _does_
    //! include halo cells that lie along periodic boundary edges.
    //!
    Halo = 0x04,

    //!
    //! Specifies that the operation should be applied to all cells, including interior, halo, and boundary
    //! cells.
    //!
    All = 0x07

};

ENABLE_BITMASK_OPERATORS( OpDomain )


//------------------------------------------------------------------------------------------------------------
//! \brief  Class for storing and managing spatial parameters such as bounds, mesh sizes, and MPI domain
//!         decomposition configuration.
//!
//------------------------------------------------------------------------------------------------------------
class DomainDecomposition {

public:

    //========================================================================================================
    //=== CONSTRUCTORS, DESTRUCTOR, AND RECONFIGURATION ROUTINES =============================================
    //========================================================================================================

    //--------------------------------------------------------------------------------------------------------
    //! \brief  Constructs an empty DomainDecomposition object with zero-initialized parameters.
    //--------------------------------------------------------------------------------------------------------
    DomainDecomposition( void );

    //--------------------------------------------------------------------------------------------------------
    //! \brief  Constructs a DomainDecomposition object using the given values.
    //--------------------------------------------------------------------------------------------------------
    DomainDecomposition(
        const SpatialDomain &,
        const bool = true
    );

    //--------------------------------------------------------------------------------------------------------
    //! \brief  Constructs a DomainDecomposition object using values read from a ParameterList for the global
    //!         mesh information.
    //--------------------------------------------------------------------------------------------------------
    DomainDecomposition(
        const ParameterList &,
        const bool = true
    );

    //--------------------------------------------------------------------------------------------------------
    //! \brief  Copy constructor for DomainDecomposition class.
    //--------------------------------------------------------------------------------------------------------
    DomainDecomposition( const DomainDecomposition & );

    //--------------------------------------------------------------------------------------------------------
    //! \brief  Destructor for DomainDecomposition class.
    //--------------------------------------------------------------------------------------------------------
    virtual ~DomainDecomposition( void );


    //========================================================================================================
    //=== OPERATOR OVERLOADS =================================================================================
    //========================================================================================================

    //--------------------------------------------------------------------------------------------------------
    //! \brief  Copy assignment operator for DomainDecomposition class.
    //--------------------------------------------------------------------------------------------------------
    DomainDecomposition & operator=( const DomainDecomposition & );


    //========================================================================================================
    //=== STATIC MEMBER FUNCTIONS ============================================================================
    //========================================================================================================

    //--------------------------------------------------------------------------------------------------------
    //! \brief  Two DomainDecomposition objects "are matching" if all of their parameters are equal.
    //--------------------------------------------------------------------------------------------------------
    static bool AreMatching(
        const DomainDecomposition &,
        const DomainDecomposition &
    );


    //========================================================================================================
    //=== ADDITIONAL MEMBER FUNCTIONS ========================================================================
    //========================================================================================================

    //--------------------------------------------------------------------------------------------------------
    //! \brief  Returns a const reference to the number of spatial cells in the specified dimension of the
    //!         local subdomain.
    //!
    //! \param[in]  dim     Index of spatial dimension to return value for.
    //--------------------------------------------------------------------------------------------------------
    inline const int64_t & nx (

        const int64_t dim

    ) const {

    # if defined (STRICT_CHECK)
        CHECK_DIM_VALUE(dim)
    # endif

        return this->Subdomain( this->MPI_block_coords_arr ).nx(dim);
    }

    //--------------------------------------------------------------------------------------------------------
    //! \brief  Returns a const reference to an array containing the number spatial cells in each dimension of
    //!         the local subdomain.
    //--------------------------------------------------------------------------------------------------------
    inline auto nx( void ) const -> const int64_t (&) [SPACE_DIMS] {

        return this->Subdomain( this->MPI_block_coords_arr ).nx();
    }

    //--------------------------------------------------------------------------------------------------------
    //! \brief  Returns a const reference to the lower limit of the local spatial domain in the specified
    //!         dimension.
    //!
    //! \param[in]  dim     Index of spatial dimension to return value for.
    //--------------------------------------------------------------------------------------------------------
    inline const double & ax (

        const int64_t dim

    ) const {

    # if defined (STRICT_CHECK)
        CHECK_DIM_VALUE(dim)
    # endif

        return this->Subdomain( this->MPI_block_coords_arr ).ax(dim);
    }

    //--------------------------------------------------------------------------------------------------------
    //! \brief  Returns a const reference to an array containing the lower limit of the local subdomain in
    //!         each dimension.
    //--------------------------------------------------------------------------------------------------------
    inline auto ax( void ) const -> const double (&) [SPACE_DIMS] {

        return this->Subdomain( this->MPI_block_coords_arr ).ax();
    }

    //--------------------------------------------------------------------------------------------------------
    //! \brief  Returns a const reference to the upper limit of the local spatial domain in the specified
    //!         dimension.
    //!
    //! \param[in]  dim     Index of spatial dimension to return value for.
    //--------------------------------------------------------------------------------------------------------
    inline const double & bx (

        const int64_t dim

    ) const {

    # if defined (STRICT_CHECK)
        CHECK_DIM_VALUE(dim)
    # endif

        return this->Subdomain( this->MPI_block_coords_arr ).bx(dim);
    }

    //--------------------------------------------------------------------------------------------------------
    //! \brief  Returns a const reference to an array containing the upper limit of the local subdomain in
    //!         each dimension.
    //--------------------------------------------------------------------------------------------------------
    inline auto bx( void ) const -> const double (&) [SPACE_DIMS] {

        return this->Subdomain( this->MPI_block_coords_arr ).bx();
    }

    //--------------------------------------------------------------------------------------------------------
    //! \brief  Returns a const reference to the spatial mesh size in the specified dimension.
    //!
    //! \param[in]  dim     Index of spatial dimension to return value for.
    //--------------------------------------------------------------------------------------------------------
    inline const double & dx (

        const int64_t dim

    ) const {

    # if defined (STRICT_CHECK)
        CHECK_DIM_VALUE(dim)
    # endif

        return this->global_domain.dx(dim);
    }

    //--------------------------------------------------------------------------------------------------------
    //! \brief  Returns a const reference to an array containing the spatial mesh size in each dimension.
    //--------------------------------------------------------------------------------------------------------
    inline auto dx( void ) const -> const double (&) [SPACE_DIMS] {  return this->global_domain.dx();  }

    //--------------------------------------------------------------------------------------------------------
    //! \brief  Returns a const reference to the number of subdomains in the specified dimension.
    //!
    //! \param[in]  dim     Index of spatial dimension to return value for.
    //--------------------------------------------------------------------------------------------------------
    inline const int & MPI_num_blocks (

        const int64_t dim

    ) const {

    # if defined (STRICT_CHECK)
        CHECK_DIM_VALUE(dim)
    # endif

        return this->MPI_num_blocks_arr[dim];
    }

    //--------------------------------------------------------------------------------------------------------
    //! \brief  Returns a const reference to an array containing the number of subdomains in each dimension.
    //--------------------------------------------------------------------------------------------------------
    inline auto MPI_num_blocks( void ) const -> const int (&) [SPACE_DIMS] {

        return this->MPI_num_blocks_arr;
    }

    //--------------------------------------------------------------------------------------------------------
    //! \brief  Returns a const reference to the coordinate of the local subdomain in the specified dimension.
    //!
    //! \param[in]  dim     Index of spatial dimension to return value for.
    //--------------------------------------------------------------------------------------------------------
    inline const int & MPI_block_coords (

        const int64_t dim

    ) const {

    # if defined (STRICT_CHECK)
        CHECK_DIM_VALUE(dim)
    # endif

        return this->MPI_block_coords_arr[dim];
    }

    //--------------------------------------------------------------------------------------------------------
    //! \brief  Returns a const reference to an array containing the coordinates of the local subdomain.
    //--------------------------------------------------------------------------------------------------------
    inline auto MPI_block_coords( void ) const -> const int (&) [SPACE_DIMS] {

        return this->MPI_block_coords_arr;
    }

    //--------------------------------------------------------------------------------------------------------
    //! \brief  Returns a const reference to the total number of spatial cells in the global spatial mesh in
    //!         the specified dimension.
    //!
    //! \param[in]  dim     Index of spatial dimension to return value for.
    //--------------------------------------------------------------------------------------------------------
    inline const int64_t & global_nx (

        const int64_t dim

    ) const {

    # if defined (STRICT_CHECK)
        CHECK_DIM_VALUE(dim)
    # endif

        return this->global_domain.nx(dim);
    }

    //--------------------------------------------------------------------------------------------------------
    //! \brief  Returns a const reference to an array containing the number of spatial cells in each dimension
    //!         of the global spatial domain.
    //--------------------------------------------------------------------------------------------------------
    inline auto global_nx( void ) const -> const int64_t (&) [SPACE_DIMS] {

        return this->global_domain.nx();
    }

    //--------------------------------------------------------------------------------------------------------
    //! \brief  Returns a const reference to the lower limit of the global spatial domain in the specified
    //!         dimension.
    //!
    //! \param[in]  dim     Index of spatial dimension to return value for.
    //--------------------------------------------------------------------------------------------------------
    inline const double & global_ax (

        const int64_t dim

    ) const {

    # if defined (STRICT_CHECK)
        CHECK_DIM_VALUE(dim)
    # endif

        return this->global_domain.ax(dim);
    }

    //--------------------------------------------------------------------------------------------------------
    //! \brief  Returns a const reference to an array containing the lower limit of the global spatial domain
    //!         in each dimension.
    //--------------------------------------------------------------------------------------------------------
    inline auto global_ax( void ) const -> const double (&) [SPACE_DIMS] {

        return this->global_domain.ax();
    }

    //--------------------------------------------------------------------------------------------------------
    //! \brief  Returns a const reference to the upper limit of the global spatial domain in the specified
    //!         dimension.
    //!
    //! \param[in]  dim     Index of spatial dimension to return value for.
    //--------------------------------------------------------------------------------------------------------
    inline const double & global_bx (

        const int64_t dim

    ) const {

    # if defined (STRICT_CHECK)
        CHECK_DIM_VALUE(dim)
    # endif

        return this->global_domain.bx(dim);
    }

    //--------------------------------------------------------------------------------------------------------
    //! \brief  Returns a const reference to an array containing the upper limit of the global spatial domain
    //!         in each dimension.
    //--------------------------------------------------------------------------------------------------------
    inline auto global_bx( void ) const -> const double (&) [SPACE_DIMS] {

        return this->global_domain.bx();
    }

    //--------------------------------------------------------------------------------------------------------
    //! \brief  Returns a const reference to the SpatialDomain object for the subdomain with given
    //!         coordinates.
    //!
    //! \param[in]  coords  Array containing coordinates of subdomain to retrieve object for.
    //--------------------------------------------------------------------------------------------------------
    inline const SpatialDomain & Subdomain (

        const int (& coords) [SPACE_DIMS]

    ) const {

        return Subdomain
                    # if SPACE_DIMS == 1
                        ( coords[0] );
                    # elif SPACE_DIMS == 2
                        ( coords[0], coords[1] );
                    # elif SPACE_DIMS == 3
                        ( coords[0], coords[1], coords[2] );
                    # endif
    }

    //--------------------------------------------------------------------------------------------------------
    //! \brief  Returns a const reference to the SpatialDomain object for the subdomain with given
    //!         coordinates.
    //!
    //! \param[in]  i   First coordinate of subdomain to retrieve object for.
    //! \param[in]  j   Second coordinate of subdomain to retrieve object for.
    //--------------------------------------------------------------------------------------------------------
    inline const SpatialDomain & Subdomain (

            const int64_t i
    # if SPACE_DIMS >= 2 || defined (DOXYCOMPILE)
        ,   const int64_t j
    # endif
    # if SPACE_DIMS == 3
        ,   const int64_t k
    # endif

    ) const {

        return const_cast<DomainDecomposition *>( this )->Subdomain_NonConst
                                                                        # if SPACE_DIMS == 1
                                                                            (i);
                                                                        # elif SPACE_DIMS == 2
                                                                            (i,j);
                                                                        # elif SPACE_DIMS == 3
                                                                            (i,j,k);
                                                                        # endif
    }

    //--------------------------------------------------------------------------------------------------------
    //! \brief  Zeros all parameters of a DomainDecomposition object.
    //--------------------------------------------------------------------------------------------------------
    virtual DomainDecomposition & Zero( void );

    //--------------------------------------------------------------------------------------------------------
    //! \brief  Prints a summary of the object configuration to the logging interface.
    //--------------------------------------------------------------------------------------------------------
    virtual void Print( const std::string & = "  " ) const;

    //--------------------------------------------------------------------------------------------------------
    //! \brief  Computes the indices (block and local) associated with the spatial cell in which the given
    //!         point lies.
    //--------------------------------------------------------------------------------------------------------
    void PointToCell(
        const double (& point) [SPACE_DIMS],
        int (& block_coords) [SPACE_DIMS],
        int (& local_coords) [SPACE_DIMS]
    ) const;

    //--------------------------------------------------------------------------------------------------------
    //! \brief  Computes the indices (block and local) associated with the spatial cell in which the given
    //!         point lies.
    //--------------------------------------------------------------------------------------------------------
    void PointToCell(
        std::initializer_list<double> point,
        int (& block_coords) [SPACE_DIMS],
        int (& local_coords) [SPACE_DIMS]
    ) const;

    //--------------------------------------------------------------------------------------------------------
    //! \brief  Determines offsets for handling of boundary and halo cells in operations.
    //--------------------------------------------------------------------------------------------------------
    void DetermineOffsets(
        const OpDomain op_domain,
        int64_t (& offset) [SPACE_DIMS][2]
    ) const;

    //--------------------------------------------------------------------------------------------------------
    //! \brief  Performs initialization of cross sections, initial condition, and sources.
    //--------------------------------------------------------------------------------------------------------
    friend void init(
        RKDG::OrdinateFlux &,
        RKDG::OrdinateFlux &,
        RKDG::CrossSection &,
        RKDG::CrossSection &,
        DomainDecomposition &
    );


protected:

    //========================================================================================================
    //=== INTERNAL VARIABLES =================================================================================
    //========================================================================================================

    SpatialDomain global_domain;            //!< Contains information about the global spatial domain.

    int MPI_num_blocks_arr [SPACE_DIMS];    //!< Number of blocks in each dimension of domain decomposition.
    int MPI_block_coords_arr [SPACE_DIMS];  //!< Coordinates of local block in global domain decomposition.

    SpatialDomain * local_subdomains;       //!< Contains information about each subdomain in the spatial domain decomposition.


    //========================================================================================================
    //=== PROTECTED MEMBER FUNCTIONS =========================================================================
    //========================================================================================================

    //--------------------------------------------------------------------------------------------------------
    //! \brief  Sets values for DomainDecomposition::MPI_num_blocks_arr and
    //!         DomainDecomposition::MPI_block_coords_arr from the Cartesian communicator
    //!         Global::MPI_cart_comm if the input argument is true.
    //--------------------------------------------------------------------------------------------------------
    void SetTopology( const bool );

    //--------------------------------------------------------------------------------------------------------
    //! \brief  Returns a non-const reference to the SpatialDomain object for the subdomain with given
    //!         coordinates.
    //!
    //! \param[in]  i   First coordinate of subdomain to retrieve object for.
    //! \param[in]  j   Second coordinate of subdomain to retrieve object for.
    //--------------------------------------------------------------------------------------------------------
    inline SpatialDomain & Subdomain_NonConst (

            const int64_t i
    # if SPACE_DIMS >= 2 || defined (DOXYCOMPILE)
        ,   const int64_t j
    # endif
    # if SPACE_DIMS == 3
        ,   const int64_t k
    # endif
    ) {

    # if defined (STRICT_CHECK)

        if ( this->local_subdomains == nullptr ) {

            std::string error_message = "Pointer this->local_subdomains is NULL in '" + std::string(__func__)
                                        + "'.\n";

            PRINT_ERROR( error_message.c_str() )
            throw std::runtime_error( error_message );
        }

        if ( i < 0 || i >= MPI_num_blocks(0) ) {

            std::string error_message = "Domain block index '"
                                        + std::to_string(i)
                                        + "' not in valid range [0,"
                                        + std::to_string( MPI_num_blocks(0) )
                                        + ") for dimension 0 in "
                                        + std::string(__func__)
                                        + ".\n";

            PRINT_ERROR( error_message.c_str() )
            throw std::out_of_range( error_message );
        }

    # if SPACE_DIMS >= 2

        if ( j < 0 || j >= MPI_num_blocks(1) ) {

            std::string error_message = "Domain block index '"
                                        + std::to_string(j)
                                        + "' not in valid range [0,"
                                        + std::to_string( MPI_num_blocks(1) )
                                        + ") for dimension 0 in "
                                        + std::string(__func__)
                                        + ".\n";

            PRINT_ERROR( error_message.c_str() )
            throw std::out_of_range( error_message );
        }

    # endif // if SPACE_DIMS >= 2
    # if SPACE_DIMS == 3

        if ( k < 0 || k >= MPI_num_blocks(2) ) {

            std::string error_message = "Domain block index '"
                                        + std::to_string(k)
                                        + "' not in valid range [0,"
                                        + std::to_string( MPI_num_blocks(2) )
                                        + ") for dimension 0 in "
                                        + std::string(__func__)
                                        + ".\n";

            PRINT_ERROR( error_message.c_str() )
            throw std::out_of_range( error_message );
        }

    # endif // if SPACE_DIMS == 3
    # endif // if defined (STRICT_CHECK)

    # if SPACE_DIMS == 1

        return this->local_subdomains[ i ];

    # elif SPACE_DIMS == 2

        return this->local_subdomains[ j + MPI_num_blocks(1) * i ];

    # elif SPACE_DIMS == 3

        return this->local_subdomains[ k + MPI_num_blocks(2) * ( j + MPI_num_blocks(1) * i ) ];

    # endif // if SPACE_DIMS == ?
    }

    //--------------------------------------------------------------------------------------------------------
    //! \brief  Allocates memory and computes subdomains based on global parameters.
    //--------------------------------------------------------------------------------------------------------
    void ComputeSubdomains( void );

    //--------------------------------------------------------------------------------------------------------
    //! \brief  Constructs a DomainDecomposition object form the values stored at the given file pointer using
    //!         the provided subdomain decomposition parameters.
    //--------------------------------------------------------------------------------------------------------
    void ReadFromDisk(
    # if defined (ENABLE_MPI)
        MPI_File &,
    # else
        std::FILE * &,
    # endif
        const bool decompose_mpi_ranks = true,
        const int64_t version = 0
    );

    //--------------------------------------------------------------------------------------------------------
    //! \brief  Writes the parameters in a DomainDecomposition object to the given file pointer.
    //--------------------------------------------------------------------------------------------------------
    void WriteToDisk(
    # if defined (ENABLE_MPI)
        MPI_File &
    # else
        std::FILE * &
    # endif
    ) const;

};


# endif // ifndef __DOMAIN_DECOMPOSITION_HPP__
