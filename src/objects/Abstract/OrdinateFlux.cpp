//------------------------------------------------------------------------------------------------------------
//! \file   objects/Abstract/OrdinateFlux.cpp
//! \brief  Implementation of Abstract::OrdinateFlux class.
//!
//! \author Michael Crockatt
//! \date   May 2017
//------------------------------------------------------------------------------------------------------------


# include "objects/Abstract/OrdinateFlux.hpp"
# include "utils/CLog.hpp"


using namespace Abstract;
using namespace Quadrule;


//============================================================================================================
//=== CONSTRUCTORS, DESTRUCTOR, AND ASSOCIATED HELPER ROUTINES ===============================================
//============================================================================================================


//------------------------------------------------------------------------------------------------------------
//! \brief  Constructs an empty Abstract::OrdinateFlux object with zero-initialized parameters.
//------------------------------------------------------------------------------------------------------------
OrdinateFlux::OrdinateFlux ( void ) :

    DomainDecomposition(),
    DensityFunction(),
    upwind_halo_cells_dirty { false }
{
    PRINT_STATUS( "Executing Abstract::OrdinateFlux::%s.\n", __func__ )
}


//------------------------------------------------------------------------------------------------------------
//! \brief  Constructs an Abstract::OrdinateFlux with the given parameters.
//!
//! \param[in]  ang_order           Order of discrete ordinates quadrature of specified type.
//! \param[in]  symmetric_reduce    (optional) <br/>
//!                                 Specifies whether or not the number of ordinates should be reduced by
//!                                 symmetry. Only used for 2D problems. Defaults to false.
//! \param[in]  ordinate_type       (optional) <br/>
//!                                 Type of discrete ordinate quadrature to use.
//!                                 Defaults to Gauss-Legendre quadrature in 1D and Chebyshev-Legendre
//!                                 quadrature in multi-D.
//------------------------------------------------------------------------------------------------------------
OrdinateFlux::OrdinateFlux (

    const int64_t ang_order,
    const bool symmetric_reduce,        // = false
    const OrdinateType ordinate_type    // Defaulted.
) :
    DomainDecomposition(),
    DensityFunction(),
    OrdinateSet( ang_order, symmetric_reduce, ordinate_type )
{
    PRINT_STATUS( "Executing Abstract::OrdinateFlux::%s.\n", __func__ )
}


//------------------------------------------------------------------------------------------------------------
//! \brief  Reconfigures an Abstract::OrdinateFlux object with a different set of parameters.
//!
//! \param[in]  ang_order           Order of discrete ordinates quadrature of specified type.
//! \param[in]  symmetric_reduce    (optional) <br/>
//!                                 Specifies whether or not the number of ordinates should be reduced by
//!                                 symmetry. Only used for 2D problems. Defaults to false.
//! \param[in]  ordinate_type       (optional) <br/>
//!                                 Type of discrete ordinate quadrature to use.
//!                                 Defaults to Gauss-Legendre quadrature in 1D and Chebyshev-Legendre
//!                                 quadrature in multi-D.
//------------------------------------------------------------------------------------------------------------
OrdinateFlux & OrdinateFlux::Reconfigure (

    const int64_t ang_order,
    const bool symmetric_reduce,        // = false
    const OrdinateType ordinate_type    // Defaulted.
) {

    PRINT_STATUS( "Executing Abstract::OrdinateFlux::%s.\n", __func__ )

    Deallocate();

    this->OrdinateSet::Reconfigure( ang_order, symmetric_reduce, ordinate_type );

    Allocate();

    return *this;
}


//------------------------------------------------------------------------------------------------------------
//! \brief  Reconfigures an Abstract::OrdinateFlux object with a different set of parameters.
//!
//! \param[in]      ordinate_set        OrdinateSet object containing parameters to reconfigure object with.
//------------------------------------------------------------------------------------------------------------
OrdinateFlux & OrdinateFlux::Reconfigure (

    const OrdinateSet & ordinate_set
) {

    PRINT_STATUS( "Executing Abstract::OrdinateFlux::%s.\n", __func__ )

    Deallocate();

    this->OrdinateSet::Reconfigure( ordinate_set );

    Allocate();

    return *this;
}


//============================================================================================================
//=== STATIC MEMBER FUNCTIONS ================================================================================
//============================================================================================================


//------------------------------------------------------------------------------------------------------------
//! \brief  Performs a deep copy of one object to another.
//!
//! \param[out]     dest        Destination for copy.
//! \param[in]      source      Source for copy.
//! \param[in]      op_domain   (optional) <br>
//!                             Domain over which to apply copy operation. Defaults to OpDomain::All.
//------------------------------------------------------------------------------------------------------------
void OrdinateFlux::Copy (

    OrdinateFlux & dest,
    const OrdinateFlux & source,
    const OpDomain op_domain        // = OpDomain::All
) {

    PRINT_STATUS( "Executing Abstract::OrdinateFlux::%s.\n", __func__ )

    dest.Copy( source, op_domain );
}


//------------------------------------------------------------------------------------------------------------
//! \brief  \f$ \vec{y} \gets \beta \vec{y} + \alpha \vec{x} \f$.
//!
//! Performs an operation on two OrdinateFlux objects corresponding to the BLAS routines \c ?axpy.
//!
//! \param[in]      alpha       The scalar \f$ \alpha \f$.
//! \param[in]      beta        The scalar \f$ \beta \f$.
//! \param[in]      x           The vector \f$ \vec{x} \f$.
//! \param[in,out]  y           The vector \f$ \vec{y} \f$.
//------------------------------------------------------------------------------------------------------------
void OrdinateFlux::AXPY (

    const double alpha,
    const double beta,
    const OrdinateFlux & x,
    OrdinateFlux & y
) {

    PRINT_STATUS( "Executing Abstract::OrdinateFlux::%s.\n", __func__ )

    DensityFunction::AXPY( alpha, beta, x, y );

    y.upwind_halo_cells_dirty |= x.upwind_halo_cells_dirty;
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
void OrdinateFlux::Print (

        const std::string & prefix // = "  "

) const {

    PRINT_STATUS( "Executing Abstract::OrdinateFlux::%s.\n", __func__ )

    this->OrdinateSet::Print( prefix );
    this->Abstract::DensityFunction::Print( prefix );

    PRINT_LOG( "\n" )

    PRINT_LOG( "%s%-*s  %s\n", prefix.c_str(),
               Global::col_width, "Upwind dirty:", this->upwind_halo_cells_dirty ? "true" : "false" )
}


//------------------------------------------------------------------------------------------------------------
//! \brief  Synchronizes the values of all halo cells between MPI ranks and across periodic boundaries.
//------------------------------------------------------------------------------------------------------------
void OrdinateFlux::SynchronizeHalos ( void ) {

    PRINT_STATUS( "Executing Abstract::OrdinateFlux::%s.\n", __func__ )

    this->ReflectBoundaries();
    this->DensityFunction::SynchronizeHalos();
}


//------------------------------------------------------------------------------------------------------------
//! \brief  Returns the number of spatial degrees of freedom per mesh element. This includes spatial and
//!         temporal, but not angular degrees of freedom; i.e., this is the number of degrees of freedom for
//!         each mesh cell and ordinate.
//------------------------------------------------------------------------------------------------------------
int64_t OrdinateFlux::SpatialDOFPerCell ( void ) const {

    PRINT_STATUS( "Executing Abstract::OrdinateFlux::%s.\n", __func__ )

# if defined (STRICT_CHECK)

    if ( this->DOFPerCell() % this->nq() ) {

        std::string error_message = "Number of ordinates (nq) does not evenly divide DOF count "
                                    "(DOFPerCell): something must have gone wrong.";

        throw std::invalid_argument( error_message );
    }

# endif // if defined (STRICT_CHECK)

    return this->DOFPerCell() / this->nq();
}


//------------------------------------------------------------------------------------------------------------
//! \brief  Returns the state of \c Abstract::OrdinateFlux::upwind_halo_cells_dirty.
//------------------------------------------------------------------------------------------------------------
bool OrdinateFlux::UpwindHalosAreDirty ( void ) const {

    PRINT_STATUS( "Executing Abstract::OrdinateFlux::%s.\n", __func__ )

    return this->upwind_halo_cells_dirty;
}


//------------------------------------------------------------------------------------------------------------
//! \brief  Marks all halo state variables as dirty.
//------------------------------------------------------------------------------------------------------------
void OrdinateFlux::MarkAllHalosDirty ( void ) {

    PRINT_STATUS( "Executing Abstract::OrdinateFlux::%s.\n", __func__ )

    this->upwind_halo_cells_dirty = true;

    // Call up class hierarchy.
    this->Abstract::DensityFunction::MarkAllHalosDirty();
}


//------------------------------------------------------------------------------------------------------------
//! \brief  Marks all halo state variables as clean.
//------------------------------------------------------------------------------------------------------------
void OrdinateFlux::MarkAllHalosClean ( void ) {

    PRINT_STATUS( "Executing Abstract::OrdinateFlux::%s.\n", __func__ )

    this->upwind_halo_cells_dirty = false;

    // Call up class hierarchy.
    this->Abstract::DensityFunction::MarkAllHalosClean();
}


//------------------------------------------------------------------------------------------------------------
//! \brief  Performs a deep copy of one object to another.
//------------------------------------------------------------------------------------------------------------
void OrdinateFlux::Copy (

    const OrdinateFlux & source,
    const OpDomain op_domain        // = OpDomain::All
) {

    PRINT_STATUS( "Executing Abstract::OrdinateFlux::%s.\n", __func__ )

    this->DensityFunction::Copy( source, op_domain );

    /*  Update status of halo cells in object.
     *
     *  Only need to update status for interior copy: Abstract::DensityFunction::Copy takes care of the
     *  rest.
     */
    if ( BitmaskHasAll( op_domain, OpDomain::Interior | OpDomain::Halo ) )
        this->upwind_halo_cells_dirty = source.upwind_halo_cells_dirty;
}
