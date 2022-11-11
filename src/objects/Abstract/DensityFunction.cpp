//------------------------------------------------------------------------------------------------------------
//! \file   objects/Abstract/DensityFunction.cpp
//! \brief  Implementation of Abstract::DensityFunction class.
//!
//! \author Michael Crockatt
//! \date   February 2018
//------------------------------------------------------------------------------------------------------------

# include <cinttypes>

# if defined (ENABLE_BELOS_TPETRA) || defined (DOXYCOMPILE)
    # include <Tpetra_MultiVector.hpp>

    # include "utils/BelosTpetraInterface.hpp"
# endif

# include "objects/Abstract/DensityFunction.hpp"
# include "utils/CLog.hpp"


using namespace Abstract;


//============================================================================================================
//=== CONSTRUCTORS, DESTRUCTOR, AND ASSOCIATED HELPER ROUTINES ===============================================
//============================================================================================================


//------------------------------------------------------------------------------------------------------------
//! \brief  Constructs an empty Abstract::DensityFunction object with zero-initialized parameters.
//------------------------------------------------------------------------------------------------------------
DensityFunction::DensityFunction ( void ) :

    DomainDecomposition(),
    halo_cells_dirty { false },
    density { nullptr },
    sizeof_density {},
    dimof_density {},
    DOF_per_cell {}
{
    PRINT_STATUS( "Executing Abstract::DensityFunction::%s.\n", __func__ )
}


//------------------------------------------------------------------------------------------------------------
//! \brief  Constructs an Abstract::DensityFunction object with the given set of parameters.
//!
//! \param[in]      spatial_params  Contains the parameters of the spatial discretization.
//------------------------------------------------------------------------------------------------------------
DensityFunction::DensityFunction (

    const DomainDecomposition & spatial_params
) :
    DomainDecomposition( spatial_params )
{
    PRINT_STATUS( "Executing Abstract::DensityFunction::%s.\n", __func__ )
}


//------------------------------------------------------------------------------------------------------------
//! \brief  Destructor for class Abstract::DensityFunction.
//!
//! Deallocates all internal memory before destroying the object.
//!
//! \see    Abstract::DensityFunction::Deallocate()
//------------------------------------------------------------------------------------------------------------
DensityFunction::~DensityFunction ( void ) {

    PRINT_STATUS( "Executing Abstract::DensityFunction::%s.\n", __func__ )

    Deallocate();
}


//------------------------------------------------------------------------------------------------------------
//! \brief  Sets parameters and allocates memory at internal pointers.
//------------------------------------------------------------------------------------------------------------
void DensityFunction::Allocate ( void ) {

    PRINT_STATUS( "Executing Abstract::DensityFunction::%s.\n", __func__ )

    SetDensityDimensions();

    // Allocate memory for density array.
    this->density = (double *)
        # if defined (ENABLE_HWLOC) && defined (USE_HWLOC_ALLOC)

            hwloc_alloc_membind( Global::machine_topology, this->sizeof_density, Global::active_core_mask,
                                 HWLOC_MEMBIND_FIRSTTOUCH, HWLOC_MEMBIND_PROCESS );

        # elif defined (USE_ALIGNED_ALLOC)
            aligned_alloc( ALIGNED_ALLOC_ALIGNMENT, this->sizeof_density );
        # else
            std::malloc( this->sizeof_density );
        # endif

    ZeroDensity();
}


//------------------------------------------------------------------------------------------------------------
//! \brief  Deallocates memory at all internal pointers.
//!
//! \see    DomainDecomposition::Zero()
//------------------------------------------------------------------------------------------------------------
void DensityFunction::Deallocate ( void ) {

    PRINT_STATUS( "Executing Abstract::DensityFunction::%s.\n", __func__ )

    if ( density != nullptr ) {

        # if defined (ENABLE_HWLOC) && defined (USE_HWLOC_ALLOC)
            hwloc_free( Global::machine_topology, density, sizeof_density );
        # else
            std::free( density );
        # endif

        density = nullptr;
    }

    this->DOF_per_cell = 0;
    this->sizeof_density = 0;
    this->dimof_density = 0;
}


//------------------------------------------------------------------------------------------------------------
//! \brief  Reconfigures an Abstract::DensityFunction object with a different set of parameters.
//!
//! \param[in]      spatial_params  Contains the parameters of the spatial discretization.
//!
//! \see    Abstract::DensityFunction::Allocate()
//! \see    Abstract::DensityFunction::Deallocate()
//------------------------------------------------------------------------------------------------------------
DensityFunction & DensityFunction::Reconfigure (

    const DomainDecomposition & spatial_params
) {

    PRINT_STATUS( "Executing Abstract::DensityFunction::%s.\n", __func__ )

    Deallocate();

    DomainDecomposition::operator=( spatial_params );

    Allocate();

    return *this;
}


//------------------------------------------------------------------------------------------------------------
//! \brief  Sets appropriate values for object parameters.
//!
//! \attention  Assumes that DensityFunction::DOF_per_cell has been properly initialized by derived class.
//!
//------------------------------------------------------------------------------------------------------------
void DensityFunction::SetDensityDimensions ( void ) {

    PRINT_STATUS( "Executing Abstract::DensityFunction::%s.\n", __func__ )

    this->dimof_density = 1;

    for ( int64_t dim = 0; dim < SPACE_DIMS; ++dim ) {

        // Force dimof_density to be zero if nx(dim) is for any dimension.
        this->dimof_density *= ( this->nx(dim) > 0 ? 1 : 0 );

        this->dimof_density *= this->nx(dim) + 2;
    }

    this->dimof_density *= this->DOF_per_cell;
    this->sizeof_density = this->dimof_density * sizeof(double);
}


//============================================================================================================
//=== STATIC MEMBER FUNCTIONS ================================================================================
//============================================================================================================


//------------------------------------------------------------------------------------------------------------
//! \brief  Returns true if the objects have matching DomainDecomposition parameters and the same number
//!         of degrees of freedom per cell.
//!
//! \param[in]  first   Objects to compare.
//! \param[in]  second  Objects to compare.
//!
//! \return     Returns true if:
//!                 1.  \code DomainDecomposition::AreMatching( first, second ) \endcode returns
//!                     true, and
//!                 2.  \c first.DOF_per_cell is equal to \c second.DOF_per_cell.
//------------------------------------------------------------------------------------------------------------
bool DensityFunction::AreCompatible (

    const DensityFunction & first,
    const DensityFunction & second
) {

    PRINT_STATUS( "Executing Abstract::DensityFunction::%s.\n", __func__ )

    return      DomainDecomposition::AreMatching( first, second )
            &&  first.DOF_per_cell == second.DOF_per_cell;
}


//------------------------------------------------------------------------------------------------------------
//! \brief  Computes the \f$ \ell_1 \f$ distance between the vectors stored at DensityFunction::density,
//!         excluding ghost cells.
//!
//! \param[in]  first   Objects to compute distance between.
//! \param[in]  second  Objects to compute distance between.
//!
//! \return     Returns the \f$ \ell_1 \f$ distance between the two vectors.
//------------------------------------------------------------------------------------------------------------
double DensityFunction::l1Dist (

    const DensityFunction & first,
    const DensityFunction & second
) {

    PRINT_STATUS( "Executing Abstract::DensityFunction::%s.\n", __func__ )

# if defined (STRICT_CHECK)

    if ( !AreCompatible( first, second )) {

        std::string error_message =   "Objects passed to Abstract::DensityFunction::"
                                    + std::string(__func__)
                                    + " are not compatible: operation cannot be applied.\n";

        PRINT_ERROR( error_message.c_str() )
        throw std::invalid_argument( error_message );
    }

# endif // if defined (STRICT_CHECK)

    double dist = 0.0;

    const int64_t (& nx) [SPACE_DIMS] = first.nx();

# if SPACE_DIMS == 1

    # pragma omp parallel for schedule(static) reduction(+:dist)
    for ( int64_t i = 1; i <= nx[0];              ++i ) {
    for ( int64_t l = 0; l <  first.DOF_per_cell; ++l ) {

        dist += std::abs( first(i,l) - second(i,l) );
    }}

# elif SPACE_DIMS == 2

    # pragma omp parallel for collapse(2) schedule(static) reduction(+:dist)
    for ( int64_t i = 1; i <= nx[0];              ++i ) {
    for ( int64_t j = 1; j <= nx[1];              ++j ) {
    for ( int64_t l = 0; l <  first.DOF_per_cell; ++l ) {

        dist += std::abs( first(i,j,l) - second(i,j,l) );
    }}}

# elif SPACE_DIMS == 3

    # pragma omp parallel for collapse(3) schedule(static) reduction(+:dist)
    for ( int64_t i = 1; i <= nx[0];              ++i ) {
    for ( int64_t j = 1; j <= nx[1];              ++j ) {
    for ( int64_t k = 1; k <= nx[2];              ++k ) {
    for ( int64_t l = 0; l <  first.DOF_per_cell; ++l ) {

        dist += std::abs( first(i,j,k,l) - second(i,j,k,l) );
    }}}}

# endif // if SPACE_DIMS == ?

# if defined (ENABLE_MPI)
    MPI_Allreduce( MPI_IN_PLACE, &dist, 1, MPI_DOUBLE, MPI_SUM, Global::MPI_cart_comm );
# endif

    return dist;
}


//------------------------------------------------------------------------------------------------------------
//! \brief  Computes the \f$ \ell_2 \f$ distance between the vectors stored at DensityFunction::density,
//!         excluding ghost cells.
//!
//! \param[in]  first   Objects to compute distance between.
//! \param[in]  second  Objects to compute distance between.
//!
//! \return     Returns the \f$ \ell_2 \f$ distance between the two vectors.
//------------------------------------------------------------------------------------------------------------
double DensityFunction::l2Dist (

    const DensityFunction & first,
    const DensityFunction & second
) {

    PRINT_STATUS( "Executing Abstract::DensityFunction::%s.\n", __func__ )

# if defined (STRICT_CHECK)

    if ( !AreCompatible( first, second )) {

        std::string error_message =   "Objects passed to Abstract::DensityFunction::"
                                    + std::string(__func__)
                                    + " are not compatible: operation cannot be applied.\n";

        PRINT_ERROR( error_message.c_str() )
        throw std::invalid_argument( error_message );
    }

# endif // if defined (STRICT_CHECK)

    double dist = 0.0;

    const int64_t (& nx) [SPACE_DIMS] = first.nx();

# if SPACE_DIMS == 1

    # pragma omp parallel for schedule(static) reduction(+:dist)
    for ( int64_t i = 1; i <= nx[0];              ++i ) {
    for ( int64_t l = 0; l <  first.DOF_per_cell; ++l ) {

        const double temp = first(i,l) - second(i,l);
        dist += temp * temp;
    }}

# elif SPACE_DIMS == 2

    # pragma omp parallel for collapse(2) schedule(static) reduction(+:dist)
    for ( int64_t i = 1; i <= nx[0];              ++i ) {
    for ( int64_t j = 1; j <= nx[1];              ++j ) {
    for ( int64_t l = 0; l <  first.DOF_per_cell; ++l ) {

        const double temp = first(i,j,l) - second(i,j,l);
        dist += temp * temp;
    }}}

# elif SPACE_DIMS == 3

    # pragma omp parallel for collapse(3) schedule(static) reduction(+:dist)
    for ( int64_t i = 1; i <= nx[0];              ++i ) {
    for ( int64_t j = 1; j <= nx[1];              ++j ) {
    for ( int64_t k = 1; k <= nx[2];              ++k ) {
    for ( int64_t l = 0; l <  first.DOF_per_cell; ++l ) {

        const double temp = first(i,j,k,l) - second(i,j,k,l);
        dist += temp * temp;
    }}}}

# endif // if SPACE_DIMS == ?

# if defined (ENABLE_MPI)
    MPI_Allreduce( MPI_IN_PLACE, &dist, 1, MPI_DOUBLE, MPI_SUM, Global::MPI_cart_comm );
# endif

    return std::sqrt( dist );
}


//------------------------------------------------------------------------------------------------------------
//! \brief  Computes the \f$ \ell_{\infty} \f$ distance between the vectors stored at
//!         DensityFunction::density, excluding ghost cells.
//!
//! \param[in]  first   Objects to compute distance between.
//! \param[in]  second  Objects to compute distance between.
//!
//! \return     Returns the \f$ \ell_{\infty} \f$ distance between the two vectors.
//------------------------------------------------------------------------------------------------------------
double DensityFunction::linfDist (

    const DensityFunction & first,
    const DensityFunction & second
) {

    PRINT_STATUS( "Executing Abstract::DensityFunction::%s.\n", __func__ )

# if defined (STRICT_CHECK)

    if ( !AreCompatible( first, second )) {

        std::string error_message =   "Objects passed to Abstract::DensityFunction::"
                                    + std::string(__func__)
                                    + " are not compatible: operation cannot be applied.\n";

        PRINT_ERROR( error_message.c_str() )
        throw std::invalid_argument( error_message );
    }

# endif // if defined (STRICT_CHECK)

    double dist = 0.0;

    const int64_t (& nx) [SPACE_DIMS] = first.nx();

# if SPACE_DIMS == 1

    # if _OPENMP >= 201107
        # pragma omp parallel for schedule(static) reduction(max:dist)
    # endif
    for ( int64_t i = 1; i <= nx[0];              ++i ) {
    for ( int64_t l = 0; l <  first.DOF_per_cell; ++l ) {

        dist = std::max( dist, std::abs( first(i,l) - second(i,l) ) );
    }}

# elif SPACE_DIMS == 2

    # if _OPENMP >= 201107
        # pragma omp parallel for collapse(2) schedule(static) reduction(max:dist)
    # endif
    for ( int64_t i = 1; i <= nx[0];              ++i ) {
    for ( int64_t j = 1; j <= nx[1];              ++j ) {
    for ( int64_t l = 0; l <  first.DOF_per_cell; ++l ) {

        dist = std::max( dist, std::abs( first(i,j,l) - second(i,j,l) ) );
    }}}

# elif SPACE_DIMS == 3

    # if _OPENMP >= 201107
        # pragma omp parallel for collapse(3) schedule(static) reduction(max:dist)
    # endif
    for ( int64_t i = 1; i <= nx[0];              ++i ) {
    for ( int64_t j = 1; j <= nx[1];              ++j ) {
    for ( int64_t k = 1; k <= nx[2];              ++k ) {
    for ( int64_t l = 0; l <  first.DOF_per_cell; ++l ) {

        dist = std::max( dist, std::abs( first(i,j,k,l) - second(i,j,k,l) ) );
    }}}}

# endif // if SPACE_DIMS == ?

# if defined (ENABLE_MPI)
    MPI_Allreduce( MPI_IN_PLACE, &dist, 1, MPI_DOUBLE, MPI_MAX, Global::MPI_cart_comm );
# endif

    return dist;
}


//------------------------------------------------------------------------------------------------------------
//! \brief  Performs a deep copy of one object to another.
//!
//! \param[out]     dest        Destination for copy.
//! \param[in]      source      Source for copy.
//! \param[in]      op_domain   (optional) <br>
//!                             Domain over which to apply copy operation. Defaults to OpDomain::All.
//------------------------------------------------------------------------------------------------------------
void DensityFunction::Copy (

    DensityFunction & dest,
    const DensityFunction & source,
    const OpDomain op_domain        // = OpDomain::All
) {
    dest.Copy( source, op_domain );
}


//------------------------------------------------------------------------------------------------------------
//! \brief  \f$ \vec{y} \gets \beta \vec{y} + \alpha \vec{x} \f$.
//!
//! Performs an operation on two DensityFunction objects corresponding to the BLAS routines \c ?axpy.
//!
//! \param[in]      alpha       The scalar \f$ \alpha \f$.
//! \param[in]      beta        The scalar \f$ \beta \f$.
//! \param[in]      x           The vector \f$ \vec{x} \f$.
//! \param[in,out]  y           The vector \f$ \vec{y} \f$.
//------------------------------------------------------------------------------------------------------------
void DensityFunction::AXPY (

    const double alpha,
    const double beta,
    const DensityFunction & x,
    DensityFunction & y
) {

    PRINT_STATUS( "Executing Abstract::DensityFunction::%s.\n", __func__ )

# if defined (STRICT_CHECK)

    if ( !AreCompatible( x, y )) {

        std::string error_message =   "Objects passed to Abstract::DensityFunction::"
                                    + std::string(__func__)
                                    + " are not compatible: operation cannot be applied.\n";

        PRINT_ERROR( error_message.c_str() )
        throw std::invalid_argument( error_message );
    }

# endif // if defined (STRICT_CHECK)

    Global::TMR_SD_axpy.Start();

    const int64_t (& nx) [SPACE_DIMS] = y.nx();

# if SPACE_DIMS == 1

    # pragma omp parallel for schedule(static)
    for ( int64_t i = 0; i <= nx[0] + 1;      ++i ) {
    for ( int64_t l = 0; l <  y.DOF_per_cell; ++l ) {

        y(i,l) = beta * y(i,l) + alpha * x(i,l);
    }}

# elif SPACE_DIMS == 2

    # pragma omp parallel for collapse(2) schedule(static)
    for ( int64_t i = 0; i <= nx[0] + 1;      ++i ) {
    for ( int64_t j = 0; j <= nx[1] + 1;      ++j ) {
    for ( int64_t l = 0; l <  y.DOF_per_cell; ++l ) {

        y(i,j,l) = beta * y(i,j,l) + alpha * x(i,j,l);
    }}}

# elif SPACE_DIMS == 3

    # pragma omp parallel for collapse(3) schedule(static)
    for ( int64_t i = 0; i <= nx[0] + 1;      ++i ) {
    for ( int64_t j = 0; j <= nx[1] + 1;      ++j ) {
    for ( int64_t k = 0; k <= nx[2] + 1;      ++k ) {
    for ( int64_t l = 0; l <  y.DOF_per_cell; ++l ) {

        y(i,j,k,l) = beta * y(i,j,k,l) + alpha * x(i,j,k,l);
    }}}}

# endif // if SPACE_DIMS == ?

    y.halo_cells_dirty |= x.halo_cells_dirty;

    Global::TMR_SD_axpy.Stop();
}


//============================================================================================================
//=== ADDITIONAL MEMBER FUNCTIONS ============================================================================
//============================================================================================================


//------------------------------------------------------------------------------------------------------------
//! \brief  Prints a summary of the object configuration to the logging interface.
//!
//! \param[in]  prefix  String prepended to each line of output, intended to be used for controlling
//!                     indentation of output.
//------------------------------------------------------------------------------------------------------------
void DensityFunction::Print (

    const std::string & prefix // = "  "

) const {

    PRINT_STATUS( "Executing Abstract::DensityFunction::%s.\n", __func__ )

    this->DomainDecomposition::Print( prefix );

    PRINT_LOG( "%s%-*s % " PRId64 "\n", prefix.c_str(),
               Global::col_width, "Dim of density:", this->dimof_density )

    PRINT_LOG( "%s%-*s  %zu\n", prefix.c_str(),
               Global::col_width, "Size (Bytes):", this->sizeof_density )

    PRINT_LOG( "%s%-*s  %s\n", prefix.c_str(),
               Global::col_width, "Halo dirty:", ( this->HalosAreDirty() ? "true" : "false" ) )

    PRINT_LOG( "\n" )

    PRINT_LOG( "%s%-*s % .2e\n", prefix.c_str(),
               Global::col_width, "l1 norm:", this->l1Norm() )

    PRINT_LOG( "%s%-*s % .2e\n", prefix.c_str(),
               Global::col_width, "l2 norm:", this->l2Norm() )

    PRINT_LOG( "%s%-*s % .2e\n", prefix.c_str(),
               Global::col_width, "linf norm:", this->linfNorm() )
}


//------------------------------------------------------------------------------------------------------------
//! \brief  Deallocates memory at all internal pointers and zeros all parameters of an
//!         Abstract::DensityFunction object.
//------------------------------------------------------------------------------------------------------------
DensityFunction & DensityFunction::Zero ( void ) {

    PRINT_STATUS( "Executing Abstract::DensityFunction::%s.\n", __func__ )

    Deallocate();
    DomainDecomposition::Zero();

    return *this;
}


//------------------------------------------------------------------------------------------------------------
//! \brief  Zeros coefficients at the pointer RKDG::DensityFunction::density.
//!
//! \param[in]  op_domain   (optional) <br>
//!                         Domain over which to apply operation. Defaults to OpDomain::All.
//------------------------------------------------------------------------------------------------------------
DensityFunction & DensityFunction::ZeroDensity (

    const OpDomain op_domain    // = OpDomain::All
) {

    PRINT_STATUS( "Executing Abstract::DensityFunction::%s.\n", __func__ )

    Global::TMR_SD_zero.Start();

    /* Determine offsets for handling boundary and halo cells.
     *
     * Offsets indicate whether the ghost cells along the specified local domain edge are included (1) or not
     * (0) in the operation.
     *
     * First index specifies spatial dimension.
     * Second index specifies either:
     *      index 0:    offset for lower bound of loop
     *      index 1:    offset for upper bound of loop
     * in the given spatial dimension.
     */
    int64_t offset[][2] = {     { 0, 0 }
                        # if SPACE_DIMS >= 2
                            ,   { 0, 0 }
                        # endif
                        # if SPACE_DIMS == 3
                            ,   { 0, 0 }
                        # endif
                          };

    this->DetermineOffsets( op_domain, offset );

    /* Update status of halo cells in object.
     *
     * If we zero both the interior and halo components, then we know the synchronization is up to date. If
     * we do one, but not both, then we know that the ranks are unsynchronized. Otherwise, we leave the
     * current state as it is.
     */
    if ( BitmaskHasAll( op_domain, OpDomain::Interior | OpDomain::Halo ) )
        this->MarkAllHalosClean();
    else if ( BitmaskHasAny( op_domain, OpDomain::Interior | OpDomain::Halo ) )
        this->MarkAllHalosDirty();

# if SPACE_DIMS == 1

    // --- Operation includes interior of domain. ----------------------------------------------------- //

    if ( BitmaskHasAll( op_domain, OpDomain::Interior ) ) {

        # pragma omp parallel for schedule(static)
        for ( int64_t i = 1 - offset[0][0]; i <= this->nx(0) + offset[0][1]; ++i ) {
        for ( int64_t l = 0; l < this->DOF_per_cell; ++l ) {

            (*this)(i,l) = 0.0;
        }}

    // --- Handle boundary/halo only as a special case. ----------------------------------------------- //

    } else {

        if ( offset[0][0] ) {

            for ( int64_t l = 0; l < this->DOF_per_cell; ++l )
                (*this)(0,l) = 0.0;
        }

        if ( offset[0][1] ) {

            for ( int64_t l = 0; l < this->DOF_per_cell; ++l )
                (*this)( this->nx(0) +1, l) = 0.0;
        }
    }

# elif SPACE_DIMS == 2

    // --- Operation includes interior of domain. ----------------------------------------------------- //

    if ( BitmaskHasAll( op_domain, OpDomain::Interior ) ) {

        const int64_t (& nx) [SPACE_DIMS] = this->nx();

        # pragma omp parallel for collapse(2) schedule(static)
        for ( int64_t i = 1 - offset[0][0]; i <= nx[0] + offset[0][1]; ++i ) {
        for ( int64_t j = 1 - offset[1][0]; j <= nx[1] + offset[1][1]; ++j ) {
        for ( int64_t l = 0; l < this->DOF_per_cell; ++l ) {

            (*this)(i,j,l) = 0.0;
        }}}

    // --- Handle boundary/halo only as a special case. ----------------------------------------------- //

    } else
    # pragma omp parallel
    {
        const int64_t (& nx) [SPACE_DIMS] = this->nx();

        if ( offset[0][0] ) {

            # pragma omp for collapse(2) schedule(static)
            for ( int64_t j = 0; j <= nx[1] + 1;          ++j ) {
            for ( int64_t l = 0; l <  this->DOF_per_cell; ++l ) {

                (*this)(0,j,l) = 0.0;
            }}
        }

        if ( offset[0][1] ) {

            # pragma omp for collapse(2) schedule(static)
            for ( int64_t j = 0; j <= nx[1] + 1;          ++j ) {
            for ( int64_t l = 0; l <  this->DOF_per_cell; ++l ) {

                (*this)( this->nx(0) +1, j,l) = 0.0;
            }}
        }

        if ( offset[1][0] ) {

            # pragma omp for collapse(2) schedule(static)
            for ( int64_t i = 0; i <= nx[0] + 1;          ++i ) {
            for ( int64_t l = 0; l <  this->DOF_per_cell; ++l ) {

                (*this)(i,0,l) = 0.0;
            }}
        }

        if ( offset[1][1] ) {

            # pragma omp for collapse(2) schedule(static)
            for ( int64_t i = 0; i <= nx[0] + 1;          ++i ) {
            for ( int64_t l = 0; l <  this->DOF_per_cell; ++l ) {

                (*this)(i, this->nx(1) +1, l) = 0.0;
            }}
        }
    }

# elif SPACE_DIMS == 3

    // --- Operation includes interior of domain. ----------------------------------------------------- //

    if ( BitmaskHasAll( op_domain, OpDomain::Interior ) ) {

        const int64_t (& nx) [SPACE_DIMS] = this->nx();

        # pragma omp parallel for collapse(3) schedule(static)
        for ( int64_t i = 1 - offset[0][0]; i <= nx[0] + offset[0][1]; ++i ) {
        for ( int64_t j = 1 - offset[1][0]; j <= nx[1] + offset[1][1]; ++j ) {
        for ( int64_t k = 1 - offset[2][0]; k <= nx[2] + offset[2][1]; ++k ) {
        for ( int64_t l = 0; l < this->DOF_per_cell; ++l ) {

            (*this)(i,j,k,l) = 0.0;
        }}}}

    // --- Handle boundary/halo only as a special case. ----------------------------------------------- //

    } else
    # pragma omp parallel
    {
        const int64_t (& nx) [SPACE_DIMS] = this->nx();

        if ( offset[0][0] ) {

            # pragma omp for collapse(3) schedule(static)
            for ( int64_t j = 0; j <= nx[1] + 1;          ++j ) {
            for ( int64_t k = 0; k <= nx[2] + 1;          ++k ) {
            for ( int64_t l = 0; l <  this->DOF_per_cell; ++l ) {

                (*this)(0,j,k,l) = 0.0;
            }}}
        }

        if ( offset[0][1] ) {

            # pragma omp for collapse(3) schedule(static)
            for ( int64_t j = 0; j <= nx[1] + 1;          ++j ) {
            for ( int64_t k = 0; k <= nx[2] + 1;          ++k ) {
            for ( int64_t l = 0; l <  this->DOF_per_cell; ++l ) {

                (*this)( this->nx(0) +1, j,k,l) = 0.0;
            }}}
        }

        if ( offset[1][0] ) {

            # pragma omp for collapse(3) schedule(static)
            for ( int64_t i = 0; i <= nx[0] + 1;          ++i ) {
            for ( int64_t k = 0; k <= nx[2] + 1;          ++k ) {
            for ( int64_t l = 0; l <  this->DOF_per_cell; ++l ) {

                (*this)(i,0,k,l) = 0.0;
            }}}
        }

        if ( offset[1][1] ) {

            # pragma omp for collapse(3) schedule(static)
            for ( int64_t i = 0; i <= nx[0] + 1;          ++i ) {
            for ( int64_t k = 0; k <= nx[2] + 1;          ++k ) {
            for ( int64_t l = 0; l <  this->DOF_per_cell; ++l ) {

                (*this)(i, this->nx(1) +1, k,l) = 0.0;
            }}}
        }

        if ( offset[2][0] ) {

            # pragma omp for collapse(3) schedule(static)
            for ( int64_t i = 0; i <= nx[0] + 1;          ++i ) {
            for ( int64_t j = 0; j <= nx[1] + 1;          ++j ) {
            for ( int64_t l = 0; l <  this->DOF_per_cell; ++l ) {

                (*this)(i,j,0,l) = 0.0;
            }}}
        }

        if ( offset[2][1] ) {

            # pragma omp for collapse(3) schedule(static)
            for ( int64_t i = 0; i <= nx[0] + 1;          ++i ) {
            for ( int64_t j = 0; j <= nx[1] + 1;          ++j ) {
            for ( int64_t l = 0; l <  this->DOF_per_cell; ++l ) {

                (*this)(i,j, this->nx(2) +1, l) = 0.0;
            }}}
        }
    }

# endif // if SPACE_DIMS == ?

    Global::TMR_SD_zero.Stop();

    return *this;
}


//------------------------------------------------------------------------------------------------------------
//! \brief  Synchronizes the values of all halo cells between MPI ranks and across periodic boundaries.
//!
//! \todo   Implement 3D exchange.
//------------------------------------------------------------------------------------------------------------
void DensityFunction::SynchronizeHalos ( void ) {

    PRINT_STATUS( "Executing Abstract::DensityFunction::%s.\n", __func__ )

    Global::TMR_SynchronizeHalos.Start();

# if defined (ENABLE_MPI)
# if SPACE_DIMS == 1

    int source_rank, dest_rank;
    int block_length = this->DOFPerCell();

    // X+ shift.
    MPI_Cart_shift( Global::MPI_cart_comm, 0, 1, &source_rank, &dest_rank );

    MPI_Sendrecv( this->PointerAtCell( this->nx(0) ), block_length, MPI_DOUBLE, dest_rank,   0,
                  this->PointerAtCell( 0           ), block_length, MPI_DOUBLE, source_rank, 0,
                  Global::MPI_cart_comm, MPI_STATUS_IGNORE );

    // X- shift.
    MPI_Cart_shift( Global::MPI_cart_comm, 0, -1, &source_rank, &dest_rank );

    MPI_Sendrecv( this->PointerAtCell( 1               ), block_length, MPI_DOUBLE, dest_rank,   1,
                  this->PointerAtCell( this->nx(0) + 1 ), block_length, MPI_DOUBLE, source_rank, 1,
                  Global::MPI_cart_comm, MPI_STATUS_IGNORE );

# elif SPACE_DIMS == 2

    MPI_Datatype Y_row = MPI_DATATYPE_NULL;
    MPI_Datatype X_row = MPI_DATATYPE_NULL;

    int source_rank, dest_rank;
    int block_length = this->DOFPerCell();

    // --- Transfers along X dimension. --------------------------------------------------------------- //

    int stride_Y = this->CellStride(1);

    MPI_Type_vector( this->nx(1), block_length, stride_Y, MPI_DOUBLE, &Y_row );
    MPI_Type_commit( &Y_row );

    // X+ shift.
    MPI_Cart_shift( Global::MPI_cart_comm, 0, 1, &source_rank, &dest_rank );

    MPI_Sendrecv( this->PointerAtCell( this->nx(0), 1 ), 1, Y_row, dest_rank,   0,
                  this->PointerAtCell( 0,           1 ), 1, Y_row, source_rank, 0,
                  Global::MPI_cart_comm, MPI_STATUS_IGNORE );

    // X- shift.
    MPI_Cart_shift( Global::MPI_cart_comm, 0, -1, &source_rank, &dest_rank );

    MPI_Sendrecv( this->PointerAtCell( 1,               1 ), 1, Y_row, dest_rank,   1,
                  this->PointerAtCell( this->nx(0) + 1, 1 ), 1, Y_row, source_rank, 1,
                  Global::MPI_cart_comm, MPI_STATUS_IGNORE );

    // --- Transfers along Y dimension. --------------------------------------------------------------- //

    int stride_X = this->CellStride(0);

    MPI_Type_vector( this->nx(0), block_length, stride_X, MPI_DOUBLE, &X_row );
    MPI_Type_commit( &X_row );

    // Y+ shift.
    MPI_Cart_shift( Global::MPI_cart_comm, 1, 1, &source_rank, &dest_rank );

    MPI_Sendrecv( this->PointerAtCell( 1, this->nx(1) ), 1, X_row, dest_rank,   2,
                  this->PointerAtCell( 1, 0           ), 1, X_row, source_rank, 2,
                  Global::MPI_cart_comm, MPI_STATUS_IGNORE );

    // Y- shift.
    MPI_Cart_shift( Global::MPI_cart_comm, 1, -1, &source_rank, &dest_rank );

    MPI_Sendrecv( this->PointerAtCell( 1, 1               ), 1, X_row, dest_rank,   3,
                  this->PointerAtCell( 1, this->nx(1) + 1 ), 1, X_row, source_rank, 3,
                  Global::MPI_cart_comm, MPI_STATUS_IGNORE );

    // --- Clean up. ---------------------------------------------------------------------------------- //

    MPI_Type_free( &X_row );
    MPI_Type_free( &Y_row );

# elif SPACE_DIMS == 3

    # warning "Implementation of Abstract::DensityFunction::SynchronizeHalos incomplete."

# endif // if SPACE_DIMS == ?
# else // if defined (ENABLE_MPI)

    if ( Global::periodic ) {
    # if SPACE_DIMS == 1

        for ( int64_t l = 0; l < this->DOF_per_cell; ++l ) {

            (*this)( 0,               l) = (*this)( this->nx(0), l);
            (*this)( this->nx(0) + 1, l) = (*this)( 1,           l);
        }

    # elif SPACE_DIMS == 2

        const int64_t (& nx) [SPACE_DIMS] = this->nx();

        # pragma omp parallel for collapse(2) schedule(static)
        for ( int64_t i = 1; i <= nx[0];              ++i ) {
        for ( int64_t l = 0; l <  this->DOF_per_cell; ++l ) {

            (*this)(i, 0,               l) = (*this)(i, this->nx(1), l);
            (*this)(i, this->nx(1) + 1, l) = (*this)(i, 1,           l);
        }}

        # pragma omp parallel for collapse(2) schedule(static)
        for ( int64_t j = 1; j <= nx[1];              ++j ) {
        for ( int64_t l = 0; l <  this->DOF_per_cell; ++l ) {

            (*this)( 0,               j,l) = (*this)( this->nx(0), j,l);
            (*this)( this->nx(0) + 1, j,l) = (*this)( 1,           j,l);
        }}

    # elif SPACE_DIMS == 3

        # warning "Implementation of Abstract::DensityFunction::SynchronizeHalos incomplete."

    # endif // if SPACE_DIMS == ?
    }

# endif // if defined (ENABLE_MPI)

    this->MarkAllHalosClean();

    Global::TMR_SynchronizeHalos.Stop();
}


//------------------------------------------------------------------------------------------------------------
//! \brief  Returns the number of degrees of freedom per spatial mesh element. This includes spatial,
//!         temporal, and angular degrees of freedom (when applicable).
//------------------------------------------------------------------------------------------------------------
int64_t DensityFunction::DOFPerCell ( void ) const {

    PRINT_STATUS( "Executing Abstract::DensityFunction::%s.\n", __func__ )

    return this->DOF_per_cell;
}


//------------------------------------------------------------------------------------------------------------
//! \brief  Returns the total number of degrees of freedom used in the underlying local representation of the
//!         object, either with or without including halo elements.
//------------------------------------------------------------------------------------------------------------
int64_t DensityFunction::GetLocalArrayDim (

    const bool with_halo // = false

) const {

    PRINT_STATUS( "Executing Abstract::DensityFunction::%s.\n", __func__ )

    int64_t array_dim = 1;

    for ( int64_t dim = 0; dim < SPACE_DIMS; ++dim )
        array_dim *= this->nx(dim) + ( with_halo ? 2 : 0 );

    array_dim *= this->DOF_per_cell;

    return array_dim;
}


//------------------------------------------------------------------------------------------------------------
//! \brief  Returns the total number of degrees of freedom used in the underlying global representation of the
//!         object, either with or without including halo elements.
//------------------------------------------------------------------------------------------------------------
int64_t DensityFunction::GetGlobalArrayDim (

    const bool with_halo // = false

) const {

    PRINT_STATUS( "Executing Abstract::DensityFunction::%s.\n", __func__ )

    int64_t array_dim = 1;

    for ( int64_t dim = 0; dim < SPACE_DIMS; ++dim )
        array_dim *= this->global_nx(dim) + ( with_halo ? 2 : 0 );

    array_dim *= this->DOF_per_cell;

    return array_dim;
}


//------------------------------------------------------------------------------------------------------------
//! \brief  Returns the state of \c Abstract::DensityFunction::halo_cells_dirty.
//------------------------------------------------------------------------------------------------------------
bool DensityFunction::HalosAreDirty ( void ) const {

    return this->halo_cells_dirty;
}


//------------------------------------------------------------------------------------------------------------
//! \brief  Marks all halo state variables as dirty.
//------------------------------------------------------------------------------------------------------------
void DensityFunction::MarkAllHalosDirty ( void ) {

    this->halo_cells_dirty = true;
}


//------------------------------------------------------------------------------------------------------------
//! \brief  Marks all halo state variables as clean.
//------------------------------------------------------------------------------------------------------------
void DensityFunction::MarkAllHalosClean ( void ) {

    this->halo_cells_dirty = false;
}


//------------------------------------------------------------------------------------------------------------
//! \brief  Scales the coefficients of the object by the given scalar.
//!
//! \param[in]  scal    Scalar to scale coefficients by.
//------------------------------------------------------------------------------------------------------------
void DensityFunction::Scale (

    const double scal

) const {

    PRINT_STATUS( "Scaling coefficients of RKDG::OrdinateFlux object.\n" )

    Global::TMR_AF_scal.Start();

    const int64_t (& nx) [SPACE_DIMS] = this->nx();

# if SPACE_DIMS == 1

    # pragma omp parallel for schedule(static)
    for ( int64_t i = 0; i <= nx[0] + 1;          ++i ) {
    for ( int64_t l = 0; l <  this->DOF_per_cell; ++l ) {

        (*this)(i,l) *= scal;
    }}

# elif SPACE_DIMS == 2

    # pragma omp parallel for collapse(2) schedule(static)
    for ( int64_t i = 0; i <= nx[0] + 1;          ++i ) {
    for ( int64_t j = 0; j <= nx[1] + 1;          ++j ) {
    for ( int64_t l = 0; l <  this->DOF_per_cell; ++l ) {

        (*this)(i,j,l) *= scal;
    }}}

# elif SPACE_DIMS == 3

    # pragma omp parallel for collapse(3) schedule(static)
    for ( int64_t i = 0; i <= nx[0] + 1;          ++i ) {
    for ( int64_t j = 0; j <= nx[1] + 1;          ++j ) {
    for ( int64_t k = 0; k <= nx[2] + 1;          ++k ) {
    for ( int64_t l = 0; l <  this->DOF_per_cell; ++l ) {

        (*this)(i,j,k,l) *= scal;
    }}}}

# endif // if SPACE_DIMS == ?

    Global::TMR_AF_scal.Stop();
}


//------------------------------------------------------------------------------------------------------------
//! \brief  Computes the \f$ \ell_1 \f$ norm of the vector stored at DensityFunction::density, excluding ghost
//!         cells.
//!
//! \return     Returns the \f$ \ell_1 \f$ norm of the coefficient vector.
//------------------------------------------------------------------------------------------------------------
double DensityFunction::l1Norm ( void ) const {

    PRINT_STATUS( "Executing Abstract::DensityFunction::%s.\n", __func__ )

    double norm = 0.0;

    const int64_t (& nx) [SPACE_DIMS] = this->nx();

# if SPACE_DIMS == 1

    # pragma omp parallel for schedule(static) reduction(+:norm)
    for ( int64_t i = 1; i <= nx[0];              ++i ) {
    for ( int64_t l = 0; l <  this->DOF_per_cell; ++l ) {

        norm += std::abs( (*this)(i,l) );
    }}

# elif SPACE_DIMS == 2

    # pragma omp parallel for collapse(2) schedule(static) reduction(+:norm)
    for ( int64_t i = 1; i <= nx[0];              ++i ) {
    for ( int64_t j = 1; j <= nx[1];              ++j ) {
    for ( int64_t l = 0; l <  this->DOF_per_cell; ++l ) {

        norm += std::abs( (*this)(i,j,l) );
    }}}

# elif SPACE_DIMS == 3

    # pragma omp parallel for collapse(3) schedule(static) reduction(+:norm)
    for ( int64_t i = 1; i <= nx[0];              ++i ) {
    for ( int64_t j = 1; j <= nx[1];              ++j ) {
    for ( int64_t k = 1; k <= nx[2];              ++k ) {
    for ( int64_t l = 0; l <  this->DOF_per_cell; ++l ) {

        norm += std::abs( (*this)(i,j,k,l) );
    }}}}

# endif // if SPACE_DIMS == ?

# if defined (ENABLE_MPI)
    MPI_Allreduce( MPI_IN_PLACE, &norm, 1, MPI_DOUBLE, MPI_SUM, Global::MPI_cart_comm );
# endif

    return norm;
}


//------------------------------------------------------------------------------------------------------------
//! \brief  Computes the \f$ \ell_2 \f$ norm of the vector stored at DensityFunction::density, excluding ghost
//!         cells.
//!
//! \return     Returns the \f$ \ell_2 \f$ norm of the coefficient vector.
//------------------------------------------------------------------------------------------------------------
double DensityFunction::l2Norm ( void ) const {

    PRINT_STATUS( "Executing Abstract::DensityFunction::%s.\n", __func__ )

    double norm = 0.0;

    const int64_t (& nx) [SPACE_DIMS] = this->nx();

# if SPACE_DIMS == 1

    # pragma omp parallel for schedule(static) reduction(+:norm)
    for ( int64_t i = 1; i <= nx[0];              ++i ) {
    for ( int64_t l = 0; l <  this->DOF_per_cell; ++l ) {

        const double tempVal = (*this)(i,l);
        norm += tempVal * tempVal;
    }}

# elif SPACE_DIMS == 2

    # pragma omp parallel for collapse(2) schedule(static) reduction(+:norm)
    for ( int64_t i = 1; i <= nx[0];              ++i ) {
    for ( int64_t j = 1; j <= nx[1];              ++j ) {
    for ( int64_t l = 0; l <  this->DOF_per_cell; ++l ) {

        const double tempVal = (*this)(i,j,l);
        norm += tempVal * tempVal;
    }}}

# elif SPACE_DIMS == 3

    # pragma omp parallel for collapse(3) schedule(static) reduction(+:norm)
    for ( int64_t i = 1; i <= nx[0];              ++i ) {
    for ( int64_t j = 1; j <= nx[1];              ++j ) {
    for ( int64_t k = 1; k <= nx[2];              ++k ) {
    for ( int64_t l = 0; l <  this->DOF_per_cell; ++l ) {

        const double tempVal = (*this)(i,j,k,l);
        norm += tempVal * tempVal;
    }}}}

# endif // if SPACE_DIMS == ?

# if defined (ENABLE_MPI)
    MPI_Allreduce( MPI_IN_PLACE, &norm, 1, MPI_DOUBLE, MPI_SUM, Global::MPI_cart_comm );
# endif

    return std::sqrt( norm );
}


//------------------------------------------------------------------------------------------------------------
//! \brief  Computes the \f$ \ell_{\infty} \f$ norm of the vector stored at DensityFunction::density,
//!         excluding ghost cells.
//!
//! \return     Returns the \f$ \ell_1 \f$ norm of the coefficient vector.
//------------------------------------------------------------------------------------------------------------
double DensityFunction::linfNorm ( void ) const {

    PRINT_STATUS( "Executing Abstract::DensityFunction::%s.\n", __func__ )

    double norm = 0.0;

    const int64_t (& nx) [SPACE_DIMS] = this->nx();

# if SPACE_DIMS == 1

    # if _OPENMP >= 201107
        # pragma omp parallel for schedule(static) reduction(max:norm)
    # endif
    for ( int64_t i = 1; i <= nx[0];              ++i ) {
    for ( int64_t l = 0; l <  this->DOF_per_cell; ++l ) {

        norm = std::max( norm, std::abs( (*this)(i,l) ) );
    }}

# elif SPACE_DIMS == 2

    # if _OPENMP >= 201107
        # pragma omp parallel for collapse(2) schedule(static) reduction(max:norm)
    # endif
    for ( int64_t i = 1; i <= nx[0];              ++i ) {
    for ( int64_t j = 1; j <= nx[1];              ++j ) {
    for ( int64_t l = 0; l <  this->DOF_per_cell; ++l ) {

        norm = std::max( norm, std::abs( (*this)(i,j,l) ) );
    }}}

# elif SPACE_DIMS == 3

    # if _OPENMP >= 201107
        # pragma omp parallel for collapse(3) schedule(static) reduction(max:norm)
    # endif
    for ( int64_t i = 1; i <= nx[0];              ++i ) {
    for ( int64_t j = 1; j <= nx[1];              ++j ) {
    for ( int64_t k = 1; k <= nx[2];              ++k ) {
    for ( int64_t l = 0; l <  this->DOF_per_cell; ++l ) {

        norm = std::max( norm, std::abs( (*this)(i,j,k,l) ) );
    }}}}

# endif // if SPACE_DIMS == ?

# if defined (ENABLE_MPI)
    MPI_Allreduce( MPI_IN_PLACE, &norm, 1, MPI_DOUBLE, MPI_MAX, Global::MPI_cart_comm );
# endif

    return norm;
}


//------------------------------------------------------------------------------------------------------------
//! \brief  Performs a deep copy of one object to another.
//!
//! \param[in]  source      Object to copy to \c this.
//! \param[in]  op_domain   (optional) <br>
//!                         Domain over which to apply copy operation. Defaults to OpDomain::All.
//------------------------------------------------------------------------------------------------------------
void DensityFunction::Copy (

    const DensityFunction & source,
    const OpDomain op_domain            // = OpDomain::All
) {

    PRINT_STATUS( "Executing Abstract::DensityFunction::%s.\n", __func__ )

# if defined (STRICT_CHECK)

    if ( !AreCompatible( *this, source )) {

        std::string error_message =   "Objects passed to Abstract::DensityFunction::"
                                    + std::string(__func__)
                                    + " are not compatible: operation cannot be applied.\n";

        PRINT_ERROR( error_message.c_str() )
        throw std::invalid_argument( error_message );
    }

# endif // if defined (STRICT_CHECK)

    Global::TMR_AF_copy.Start();

    /* Determine offsets for handling boundary and halo cells.
     *
     * Offsets indicate whether the ghost cells along the specified local domain edge are included (1) or not
     * (0) in the operation.
     *
     * First index specifies spatial dimension.
     * Second index specifies either:
     *      index 0:    offset for lower bound of loop
     *      index 1:    offset for upper bound of loop
     * in the given spatial dimension.
     */
    int64_t offset[][2] = {     { 0, 0 }
                        # if SPACE_DIMS >= 2
                            ,   { 0, 0 }
                        # endif
                        # if SPACE_DIMS == 3
                            ,   { 0, 0 }
                        # endif
                          };

    this->DetermineOffsets( op_domain, offset );

    // Update status of halo cells in object.
    if (/*
         *  Halos are marked dirty unless BOTH or NEITHER of interior and halo cells are included in the
         *  operation.
         */
             BitmaskHasAny( op_domain, OpDomain::Interior | OpDomain::Halo )
         && !BitmaskHasAll( op_domain, OpDomain::Interior | OpDomain::Halo )
    ) {
        this->MarkAllHalosDirty();

    } else {
    /*  In here: op_domain has BOTH interior and halo or NEITHER.
     *
     *  Only update state of current object if copying interior values; i.e., do not update
     *  object's state if only copying boundary conditions.
     */

        if ( BitmaskHasAll( op_domain, OpDomain::Interior | OpDomain::Halo ) )
            this->halo_cells_dirty = source.halo_cells_dirty;
    }

# if SPACE_DIMS == 1

    // --- Operation includes interior of domain. ----------------------------------------------------- //

    if ( BitmaskHasAll( op_domain, OpDomain::Interior ) ) {

        # pragma omp parallel for schedule(static)
        for ( int64_t i = 1 - offset[0][0]; i <= this->nx(0) + offset[0][1]; ++i ) {
        for ( int64_t l = 0; l < this->DOF_per_cell; ++l ) {

            (*this)(i,l) = source(i,l);
        }}

    // --- Handle boundary/halo only as a special case. ----------------------------------------------- //

    } else {

        if ( offset[0][0] ) {

            for ( int64_t l = 0; l < this->DOF_per_cell; ++l )
                (*this)(0,l) = source(0,l);
        }

        if ( offset[0][1] ) {

            for ( int64_t l = 0; l < this->DOF_per_cell; ++l )
                (*this)( this->nx(0) +1, l) = source( this->nx(0) +1, l);
        }
    }

# elif SPACE_DIMS == 2

    // --- Operation includes interior of domain. ----------------------------------------------------- //

    if ( BitmaskHasAll( op_domain, OpDomain::Interior ) ) {

        const int64_t (& nx) [SPACE_DIMS] = this->nx();

        # pragma omp parallel for collapse(2) schedule(static)
        for ( int64_t i = 1 - offset[0][0]; i <= nx[0] + offset[0][1]; ++i ) {
        for ( int64_t j = 1 - offset[1][0]; j <= nx[1] + offset[1][1]; ++j ) {
        for ( int64_t l = 0; l < this->DOF_per_cell; ++l ) {

            (*this)(i,j,l) = source(i,j,l);
        }}}

    // --- Handle boundary/halo only as a special case. ----------------------------------------------- //

    } else
    # pragma omp parallel
    {
        const int64_t (& nx) [SPACE_DIMS] = this->nx();

        if ( offset[0][0] ) {

            # pragma omp for collapse(2) schedule(static)
            for ( int64_t j = 0; j <= nx[1] + 1;          ++j ) {
            for ( int64_t l = 0; l <  this->DOF_per_cell; ++l ) {

                (*this)(0,j,l) = source(0,j,l);
            }}
        }

        if ( offset[0][1] ) {

            # pragma omp for collapse(2) schedule(static)
            for ( int64_t j = 0; j <= nx[1] + 1;          ++j ) {
            for ( int64_t l = 0; l <  this->DOF_per_cell; ++l ) {

                (*this)( this->nx(0) +1, j,l) = source( this->nx(0) +1, j,l);
            }}
        }

        if ( offset[1][0] ) {

            # pragma omp for collapse(2) schedule(static)
            for ( int64_t i = 0; i <= nx[0] + 1;          ++i ) {
            for ( int64_t l = 0; l <  this->DOF_per_cell; ++l ) {

                (*this)(i,0,l) = source(i,0,l);
            }}
        }

        if ( offset[1][1] ) {

            # pragma omp for collapse(2) schedule(static)
            for ( int64_t i = 0; i <= nx[0] + 1;          ++i ) {
            for ( int64_t l = 0; l <  this->DOF_per_cell; ++l ) {

                (*this)(i, this->nx(1) +1, l) = source(i, this->nx(1) +1, l);
            }}
        }
    }

# elif SPACE_DIMS == 3

    // --- Operation includes interior of domain. ----------------------------------------------------- //

    if ( BitmaskHasAll( op_domain, OpDomain::Interior ) ) {

        const int64_t (& nx) [SPACE_DIMS] = this->nx();

        # pragma omp parallel for collapse(3) schedule(static)
        for ( int64_t i = 1 - offset[0][0]; i <= nx[0] + offset[0][1]; ++i ) {
        for ( int64_t j = 1 - offset[1][0]; j <= nx[1] + offset[1][1]; ++j ) {
        for ( int64_t k = 1 - offset[2][0]; k <= nx[2] + offset[2][1]; ++k ) {
        for ( int64_t l = 0; l < this->DOF_per_cell; ++l ) {

            (*this)(q,i,j,k,d,e,f,s) = source(q,i,j,k,d,e,f,s);
        }}}}

    // --- Handle boundary/halo only as a special case. ----------------------------------------------- //

    } else
    # pragma omp parallel
    {
        const int64_t (& nx) [SPACE_DIMS] = this->nx();

        if ( offset[0][0] ) {

            # pragma omp for collapse(2) schedule(static)
            for ( int64_t j = 0; j <= nx[1] + 1;          ++j ) {
            for ( int64_t k = 0; k <= nx[2] + 1;          ++k ) {
            for ( int64_t l = 0; l <  this->DOF_per_cell; ++l ) {

                (*this)(0,j,k,l) = source(0,j,k,l);
            }}}
        }

        if ( offset[0][1] ) {

            # pragma omp for collapse(2) schedule(static)
            for ( int64_t j = 0; j <= nx[1] + 1;          ++j ) {
            for ( int64_t k = 0; k <= nx[2] + 1;          ++k ) {
            for ( int64_t l = 0; l <  this->DOF_per_cell; ++l ) {

                (*this)( this->nx(0) +1, j,k,l) = source( this->nx(0) +1, j,k,l);
            }}}
        }

        if ( offset[1][0] ) {

            # pragma omp for collapse(2) schedule(static)
            for ( int64_t i = 0; i <= nx[0] + 1;          ++i ) {
            for ( int64_t k = 0; k <= nx[2] + 1;          ++k ) {
            for ( int64_t l = 0; l <  this->DOF_per_cell; ++l ) {

                (*this)(i,0,k,l) = source(i,0,k,l);
            }}}
        }

        if ( offset[1][1] ) {

            # pragma omp for collapse(2) schedule(static)
            for ( int64_t i = 0; i <= nx[0] + 1;          ++i ) {
            for ( int64_t k = 0; k <= nx[2] + 1;          ++k ) {
            for ( int64_t l = 0; l <  this->DOF_per_cell; ++l ) {

                (*this)(i, this->nx(1) +1, k,l) = source(i, this->nx(1) +1, k,l);
            }}}
        }

        if ( offset[2][0] ) {

            # pragma omp for collapse(2) schedule(static)
            for ( int64_t i = 0; i <= nx[0] + 1;          ++i ) {
            for ( int64_t j = 0; j <= nx[1] + 1;          ++j ) {
            for ( int64_t l = 0; l <  this->DOF_per_cell; ++l ) {

                (*this)(i,j,0,l) = source(i,j,0,l);
            }}}
        }

        if ( offset[2][1] ) {

            # pragma omp for collapse(2) schedule(static)
            for ( int64_t i = 0; i <= nx[0] + 1;          ++i ) {
            for ( int64_t j = 0; j <= nx[1] + 1;          ++j ) {
            for ( int64_t l = 0; l <  this->DOF_per_cell; ++l ) {

                (*this)(i,j, this->nx(2) +1, l) = source(i,j, this->nx(2) +1, l);
            }}}
        }
    }

# endif // if SPACE_DIMS == ?

    Global::TMR_AF_copy.Stop();
}


//------------------------------------------------------------------------------------------------------------
//! \brief  Updates the values in an external object with an operation of the form
//!         \f$ \vec{y} \gets \beta \vec{y} + \alpha \vec{x} \f$ where \f$ \vec{x} \f$ is \c this.
//!
//! \param[out] dest    Object to pack values into.
//! \param[in]  alpha   The scalar \f$ \alpha \f$.
//! \param[in]  beta    The scalar \f$ \beta \f$.
//!
//! \see    DensityFunction::UnpackPointer()
//------------------------------------------------------------------------------------------------------------
template< typename ViewType >
void DensityFunction::PackPointer (

    ViewType & dest,
    const double alpha,
    const double beta

) const {

    PRINT_STATUS( "Executing Abstract::DensityFunction::%s.\n", __func__ )

    Global::TMR_SD_pack.Start();

    const int64_t (& nx) [SPACE_DIMS] = this->nx();

# if SPACE_DIMS == 1

    # pragma omp parallel for schedule(static)
    for ( int64_t i = 1; i <= nx[0];              ++i ) {
    for ( int64_t l = 0; l <  this->DOF_per_cell; ++l ) {

        dest( this->IndexNoGC(i,l) ) = beta * dest( this->IndexNoGC(i,l) ) + alpha * (*this)(i,l);
    }}

# elif SPACE_DIMS == 2

    # pragma omp parallel for collapse(2) schedule(static)
    for ( int64_t i = 1; i <= nx[0];              ++i ) {
    for ( int64_t j = 1; j <= nx[1];              ++j ) {
    for ( int64_t l = 0; l <  this->DOF_per_cell; ++l ) {

        dest( this->IndexNoGC(i,j,l) ) = beta * dest( this->IndexNoGC(i,j,l) ) + alpha * (*this)(i,j,l);
    }}}

# elif SPACE_DIMS == 3

    # pragma omp parallel for collapse(3) schedule(static)
    for ( int64_t i = 1; i <= nx[0];              ++i ) {
    for ( int64_t j = 1; j <= nx[1];              ++j ) {
    for ( int64_t k = 1; k <= nx[2];              ++k ) {
    for ( int64_t l = 0; l <  this->DOF_per_cell; ++l ) {

        dest( this->IndexNoGC(i,j,k,l) ) = beta * dest( this->IndexNoGC(i,j,k,l) ) + alpha * (*this)(i,j,k,l);
    }}}}

# endif // if SPACE_DIMS == ?

    Global::TMR_SD_pack.Stop();
}


//------------------------------------------------------------------------------------------------------------
//! \brief  Unpacks the contents of the provided external data structure into the present object.
//!
//! Ghost cells (which are not present in the external vector) are set to zero.
//!
//! \param[in]  src     Object to unpack values from.
//!
//! \see    DensityFunction::PackPointer()
//------------------------------------------------------------------------------------------------------------
template< typename ViewType >
void DensityFunction::UnpackPointer (

    const ViewType & src
) {

    PRINT_STATUS( "Executing Abstract::DensityFunction::%s.\n", __func__ )

    this->ZeroDensity();

    Global::TMR_SD_pack.Start();

    const int64_t (& nx) [SPACE_DIMS] = this->nx();

# if SPACE_DIMS == 1

    # pragma omp parallel for schedule(static)
    for ( int64_t i = 1; i <= nx[0];              ++i ) {
    for ( int64_t l = 0; l <  this->DOF_per_cell; ++l ) {

        (*this)(i,l) = src( this->IndexNoGC(i,l) );
    }}

# elif SPACE_DIMS == 2

    # pragma omp parallel for collapse(2) schedule(static)
    for ( int64_t i = 1; i <= nx[0];              ++i ) {
    for ( int64_t j = 1; j <= nx[1];              ++j ) {
    for ( int64_t l = 0; l <  this->DOF_per_cell; ++l ) {

        (*this)(i,j,l) = src( this->IndexNoGC(i,j,l) );
    }}}

# elif SPACE_DIMS == 3

    # pragma omp parallel for collapse(3) schedule(static)
    for ( int64_t i = 1; i <= nx[0];              ++i ) {
    for ( int64_t j = 1; j <= nx[1];              ++j ) {
    for ( int64_t k = 1; k <= nx[2];              ++k ) {
    for ( int64_t l = 0; l <  this->DOF_per_cell; ++l ) {

        (*this)(i,j,k,l) = src( this->IndexNoGC(i,j,k,l) );
    }}}}

# endif // if SPACE_DIMS == ?

    if (    Global::periodic
    # if defined (ENABLE_MPI)
         || Global::MPI_num_ranks > 1
    # endif
    ) {  this->MarkAllHalosDirty();  }

    Global::TMR_SD_pack.Stop();
}


//============================================================================================================
//=== PETSC INTERFACE ROUTINES ===============================================================================
//============================================================================================================

# if defined (ENABLE_PETSC) || defined (DOXYCOMPILE)


//------------------------------------------------------------------------------------------------------------
//! \brief  Updates the values in a PETScVec object with an operation of the form
//!         \f$ \vec{y} \gets \beta \vec{y} + \alpha \vec{x} \f$ where \f$ \vec{x} \f$ is \c this.
//!
//! \param[out]     dest    PETScVec to pack.
//! \param[in]      alpha   The scalar \f$ \alpha \f$.
//! \param[in]      beta    The scalar \f$ \beta \f$.
//!
//! \see    DensityFunction::UnpackExternalVector()
//------------------------------------------------------------------------------------------------------------
template<>
void DensityFunction::PackExternalVector< PETScVec, double >(

    PETScVec & dest,
    const double alpha,
    const double beta,
    const size_t

) const {

    PRINT_STATUS( "Executing Abstract::DensityFunction::%s.\n", __func__ )

    double * dest_array;

    Global::TMR_PETSc.Start();
    VecGetArray( dest, &dest_array );
    Global::TMR_PETSc.Stop();

    MVWrapper<double *> dest_view( dest_array );

    this->PackPointer( dest_view, alpha, beta );

    Global::TMR_PETSc.Start();
    VecRestoreArray( dest, &dest_array );
    Global::TMR_PETSc.Stop();
}


//------------------------------------------------------------------------------------------------------------
//! \brief  Unpacks the coefficients stored in a PETScVec into the given object.
//!
//! Ghost cells (which are not present in the PETScVec) are set to zero.
//!
//! \param[in]      source      PETScVec to unpack.
//!
//! \see    DensityFunction::PackExternalVector()
//------------------------------------------------------------------------------------------------------------
template<>
void DensityFunction::UnpackExternalVector<PETScVec>(

    const PETScVec & source,
    const size_t
) {

    PRINT_STATUS( "Executing Abstract::DensityFunction::%s.\n", __func__ )

    double * source_array;

    Global::TMR_PETSc.Start();
    VecGetArray( source, &source_array );
    Global::TMR_PETSc.Stop();

    this->ZeroDensity();

    MVWrapper<double *> source_view( source_array );

    this->UnpackPointer( source_view );

    Global::TMR_PETSc.Start();
    VecRestoreArray( source, &source_array );
    Global::TMR_PETSc.Stop();
}


# endif // if defined (ENABLE_PETSC)


//============================================================================================================
//=== EPETRA INTERFACE ROUTINES ==============================================================================
//============================================================================================================

# if defined (ENABLE_AZTEC_EPETRA) || defined (DOXYCOMPILE)


//------------------------------------------------------------------------------------------------------------
//! \brief  Updates the values in an Epetra_MultiVector object with an operation of the form
//!         \f$ \vec{y} \gets \beta \vec{y} + \alpha \vec{x} \f$ where \f$ \vec{x} \f$ is \c this.
//!
//! \param[out]     dest    Epetra_MultiVector to pack.
//! \param[in]      alpha   The scalar \f$ \alpha \f$.
//! \param[in]      beta    The scalar \f$ \beta \f$.
//! \param[in]      column  (optional) <br>
//!                         Column of Epetra_MultiVector to update. Defaults to 0 (first column).
//!
//! \see    DensityFunction::UnpackExternalVector()
//------------------------------------------------------------------------------------------------------------
template<>
void DensityFunction::PackExternalVector< Epetra_MultiVector, double >(

    Epetra_MultiVector & dest,
    const double alpha,
    const double beta,
    const size_t column         // = 0

) const {

    PRINT_STATUS( "Executing Abstract::DensityFunction::%s.\n", __func__ )

    double * dest_array = nullptr;
    int LD_dest = -1;

    dest.ExtractView( &dest_array, &LD_dest );

    MVWrapper<double *> dest_view( dest_array, LD_dest, column );

    this->PackPointer( dest_view, alpha, beta );
}


//------------------------------------------------------------------------------------------------------------
//! \brief  Unpacks the coefficients stored in an Epetra_MultiVector into the given object.
//!
//! Ghost cells (which are not present in the Epetra_MultiVector \pp{source}) are set to zero.
//!
//! \param[in]      source  Epetra_MultiVector object to unpack coefficients from.
//! \param[in]      column  (optional) <br>
//!                         Column of Epetra_MultiVector to unpack. Defaults to 0 (first column).
//!
//! \see    DensityFunction::PackExternalVector()
//------------------------------------------------------------------------------------------------------------
template<>
void DensityFunction::UnpackExternalVector<Epetra_MultiVector>(

    const Epetra_MultiVector & source,
    const size_t column         // = 0
) {

    PRINT_STATUS( "Executing Abstract::DensityFunction::%s.\n", __func__ )

    this->ZeroDensity();

    double * source_array = nullptr;
    int LD_source = -1;

    source.ExtractView( &source_array, &LD_source );

    MVWrapper<double *> source_view( source_array, LD_source, column );

    this->UnpackPointer( source_view );
}


# endif // if defined (ENABLE_AZTEC_EPETRA)


//============================================================================================================
//=== TPETRA INTERFACE ROUTINES ==============================================================================
//============================================================================================================

# if defined (ENABLE_BELOS_TPETRA) || defined (DOXYCOMPILE)


//------------------------------------------------------------------------------------------------------------
//! \brief  Updates the values in a TpetraVec object with an operation of the form
//!         \f$ \vec{y} \gets \beta \vec{y} + \alpha \vec{x} \f$. Where \f$ \vec{x} \f$ is \c this.
//!
//! \param[out]     y       TpetraVec object to update.
//! \param[in]      alpha   The scalar \f$ \alpha \f$.
//! \param[in]      beta    The scalar \f$ \beta \f$.
//! \param[in]      column  (optional) <br>
//!                         Column of TpetraVec (MultiVector) to update. Defaults to 0 (first column).
//!
//! \see    DensityFunction::UnpackExternalVector()
//------------------------------------------------------------------------------------------------------------
template<>
void DensityFunction::PackExternalVector< TpetraVec, TpetraScalar >(

    TpetraVec & y,
    const TpetraScalar alpha,
    const TpetraScalar beta,
    const size_t column         // = 0

) const {

    PRINT_STATUS( "Executing Abstract::DensityFunction::%s.\n", __func__ )

    TpetraVec::dual_view_type::t_host y_view = y.getLocalView<Kokkos::HostSpace>();

    MVWrapper<TpetraVec::dual_view_type::t_host> y_view_wrapper( y_view, 0, column );

    this->PackPointer( y_view_wrapper, alpha, beta );
}


//------------------------------------------------------------------------------------------------------------
//! \brief  Unpacks the coefficients stored in a TpetraVec into the given DensityFunction object.
//!
//! Ghost cells (which are not present in the TpetraVec \pp{source}) are set to zero.
//!
//! \param[in]      source  TpetraVec object to unpack coefficients from.
//! \param[in]      column  (optional) <br>
//!                         Column of TpetraVec (MultiVector) to unpack. Defaults to 0 (first column).
//!
//! \see    DensityFunction::PackExternalVector()
//------------------------------------------------------------------------------------------------------------
template<>
void DensityFunction::UnpackExternalVector<TpetraVec>(

    const TpetraVec & source,
    const size_t column         // = 0
) {

    PRINT_STATUS( "Executing Abstract::DensityFunction::%s.\n", __func__ )

    this->ZeroDensity();

    TpetraVec::dual_view_type::t_host source_view = source.getLocalView<Kokkos::HostSpace>();

    MVWrapper<TpetraVec::dual_view_type::t_host> source_view_wrapper( source_view, 0, column );

    this->UnpackPointer( source_view_wrapper );
}


# endif // if defined (ENABLE_BELOS_TPETRA)
