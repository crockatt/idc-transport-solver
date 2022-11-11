//------------------------------------------------------------------------------------------------------------
//! \file   objects/STDG/OrdinateFlux.cpp
//! \brief  Implementation of STDG::OrdinateFlux class.
//!
//! \author Michael M. Crockatt
//! \date   October 2017
//------------------------------------------------------------------------------------------------------------


# include <algorithm>
# include <cinttypes>
# include <cmath>
# include <cstring>
# include <limits>
# include <stdexcept>
# include <string>

# include "objects/RKDG/OrdinateFlux.hpp"
# include "objects/STDG/OrdinateFlux.hpp"
# include "utils/CLog.hpp"
# include "utils/Quadrule/Quadrule.hpp"


using namespace Quadrule;


//============================================================================================================
//=== CONSTRUCTORS, DESTRUCTORS, AND ASSOCIATED HELPER ROUTINES ==============================================
//============================================================================================================


//------------------------------------------------------------------------------------------------------------
//! \brief  Constructs an empty STDG::OrdinateFlux object with zero-initialized parameters.
//------------------------------------------------------------------------------------------------------------
STDG::OrdinateFlux::OrdinateFlux ( void ) :

    DomainDecomposition(),
    Abstract::DensityFunction(),
    Abstract::OrdinateFlux()
{}


//------------------------------------------------------------------------------------------------------------
//! \brief  Constructs an STDG::OrdinateFlux object with given parameters.
//!
//! \param[in]  spatial_params      Contains the parameters of the spatial discretization.
//! \param[in]  DG_degree           Maximum degree of DG basis polynomials.
//! \param[in]  ang_order           Order of the discrete ordinate set of type \pp{ordinate_set} to use.
//! \param[in]  symmetric_reduce    (optional) <br/>
//!                                 Specifies whether or not the number of ordinates should be reduced by
//!                                 symmetry. Only used for 2D problems. Defaults to false.
//! \param[in]  ordinate_type       (optional) <br/>
//!                                 OrdinateType specifying the type of discrete ordinate quadrature to use.
//!                                 Defaults to Gauss-Legendre quadrature in 1D and Chebyshev-Legendre
//!                                 quadrature in multi-D.
//!
//! \see    STDG::OrdinateFlux::Allocate()
//------------------------------------------------------------------------------------------------------------
STDG::OrdinateFlux::OrdinateFlux (

    const DomainDecomposition & spatial_params,
    const int64_t DG_degree,
    const int64_t ang_order,
    const bool symmetric_reduce,        // = false
    const OrdinateType ordinate_type    // Defaulted.

) : // STDG::OrdinateFlux( spatial_params, DG_degree, DG_degree, ang_order, symmetric_reduce, ordinate_type ) {}
    DomainDecomposition( spatial_params ),
    Abstract::DensityFunction( spatial_params ),
    Abstract::OrdinateFlux( ang_order, symmetric_reduce, ordinate_type )
{
    this->DG_degree_x = DG_degree;
    this->DG_degree_t = DG_degree;

    OrdinateFlux::Allocate();
}


//------------------------------------------------------------------------------------------------------------
//! \brief  Constructs an STDG::OrdinateFlux object with given parameters.
//!
//! \param[in]  spatial_params      Contains the parameters of the spatial discretization.
//! \param[in]  DG_degree_x_in      Maximum degree of DG basis polynomials in spatial dimensions.
//! \param[in]  DG_degree_t_in      Maximum degree of DG basis polynomials in temporal dimension.
//! \param[in]  ang_order           Order of the discrete ordinate set of type \pp{ordinate_set} to use.
//! \param[in]  symmetric_reduce    (optional) <br/>
//!                                 Specifies whether or not the number of ordinates should be reduced by
//!                                 symmetry. Only used for 2D problems. Defaults to false.
//! \param[in]  ordinate_type       (optional) <br/>
//!                                 OrdinateType specifying the type of discrete ordinate quadrature to use.
//!                                 Defaults to Gauss-Legendre quadrature in 1D and Chebyshev-Legendre
//!                                 quadrature in multi-D.
//!
//! \see    STDG::OrdinateFlux::Allocate()
//------------------------------------------------------------------------------------------------------------
STDG::OrdinateFlux::OrdinateFlux (

    const DomainDecomposition & spatial_params,
    const int64_t DG_degree_x_in,
    const int64_t DG_degree_t_in,
    const int64_t ang_order,
    const bool symmetric_reduce,        // = false
    const OrdinateType ordinate_type    // Defaulted.
) :
    DomainDecomposition( spatial_params ),
    Abstract::DensityFunction( spatial_params ),
    Abstract::OrdinateFlux( ang_order, symmetric_reduce, ordinate_type )
{
    this->DG_degree_x = DG_degree_x_in;
    this->DG_degree_t = DG_degree_t_in;

    OrdinateFlux::Allocate();
}


//------------------------------------------------------------------------------------------------------------
//! \brief  Constructs an STDG::OrdinateFlux object with the given set of parameters.
//!
//! Delegates to STDG::OrdinateFlux::OrdinateFlux( const DomainDecomposition &, const int64_t, const int64_t, const OrdinateType ).
//!
//! \param[in]      spatial_params  Contains the parameters of the spatial discretization.
//! \param[in]      DG_degree       Maximum degree of DG basis polynomials.
//! \param[in]      ordinate_set    OrdinateSet containing the discrete ordinate quadrature to use.
//!
//! \see    STDG::OrdinateFlux::Allocate()
//------------------------------------------------------------------------------------------------------------
STDG::OrdinateFlux::OrdinateFlux (

    const DomainDecomposition & spatial_params,
    const OrdinateSet & ordinate_set,
    const int64_t DG_degree

) : // STDG::OrdinateFlux( spatial_params, ordinate_set, DG_degree, DG_degree ) {}
    DomainDecomposition( spatial_params ),
    Abstract::DensityFunction( spatial_params ),
    Abstract::OrdinateFlux( ordinate_set.GetAngOrder(), ordinate_set.GetOrdinateSymmetry(),
                            ordinate_set.GetOrdinateType() )
{
    this->DG_degree_x = DG_degree;
    this->DG_degree_t = DG_degree;

    OrdinateFlux::Allocate();
}


//------------------------------------------------------------------------------------------------------------
//! \brief  Constructs an STDG::OrdinateFlux object with the given set of parameters.
//!
//! Delegates to STDG::OrdinateFlux::OrdinateFlux( const DomainDecomposition &, const int64_t, const int64_t, const OrdinateType ).
//!
//! \param[in]      spatial_params  Contains the parameters of the spatial discretization.
//! \param[in]      DG_degree_x_in  Maximum degree of DG basis polynomials in spatial dimensions.
//! \param[in]      DG_degree_t_in  Maximum degree of DG basis polynomials in temporal dimension.
//! \param[in]      ordinate_set    OrdinateSet containing the discrete ordinate quadrature to use.
//!
//! \see    STDG::OrdinateFlux::Allocate()
//------------------------------------------------------------------------------------------------------------
STDG::OrdinateFlux::OrdinateFlux (

    const DomainDecomposition & spatial_params,
    const OrdinateSet & ordinate_set,
    const int64_t DG_degree_x_in,
    const int64_t DG_degree_t_in

) : /* STDG::OrdinateFlux( spatial_params, DG_degree_x_in, DG_degree_t_in, ordinate_set.GetAngOrder(),
                        ordinate_set.GetOrdinateSymmetry(), ordinate_set.GetOrdinateType() )
{} */
    DomainDecomposition( spatial_params ),
    Abstract::DensityFunction( spatial_params ),
    Abstract::OrdinateFlux( ordinate_set.GetAngOrder(), ordinate_set.GetOrdinateSymmetry(),
                            ordinate_set.GetOrdinateType() )
{
    this->DG_degree_x = DG_degree_x_in;
    this->DG_degree_t = DG_degree_t_in;

    OrdinateFlux::Allocate();
}


//------------------------------------------------------------------------------------------------------------
//! \brief  Sets appropriate values for object parameters
//------------------------------------------------------------------------------------------------------------
void STDG::OrdinateFlux::SetDensityDimensions ( void ) {

    // Set value for DOF per cell.
    this->DOF_per_cell = this->nq() * IntPow( this->DG_degree_x + 1, SPACE_DIMS ) * (this->DG_degree_t + 1);

    // Call up class hierarchy.
    this->Abstract::DensityFunction::SetDensityDimensions();
}


//------------------------------------------------------------------------------------------------------------
//! \brief  Reconfigures an STDG::OrdinateFlux object with the given set of parameters.
//!
//! \param[in]  ang_order           Order of the discrete ordinate set of type \pp{ordinate_set} to use.
//! \param[in]  symmetric_reduce    (optional) <br/>
//!                                 Specifies whether or not the number of ordinates should be reduced by
//!                                 symmetry. Only used for 2D problems. Defaults to false.
//! \param[in]  ordinate_type       (optional) <br/>
//!                                 OrdinateType specifying the type of discrete ordinate quadrature to use.
//!                                 Defaults to Gauss-Legendre quadrature in 1D and Chebyshev-Legendre
//!                                 quadrature in multi-D.
//!
//! \see    STDG::OrdinateFlux::Allocate()
//! \see    STDG::OrdinateFlux::Deallocate()
//------------------------------------------------------------------------------------------------------------
STDG::OrdinateFlux & STDG::OrdinateFlux::Reconfigure (

    const int64_t ang_order,
    const bool symmetric_reduce,        // = false
    const OrdinateType ordinate_type    // Defaulted.
) {

    Deallocate();

    this->OrdinateSet::Reconfigure( ang_order, symmetric_reduce, ordinate_type );

    Allocate();

    return *this;
}


//------------------------------------------------------------------------------------------------------------
//! \brief  Reconfigures an STDG::OrdinateFlux object with the parameters of a given OrdinateSet object.
//!
//! \param[in]      ordinate_set    OrdinateSet object containing parameters to use.
//!
//! \see    STDG::OrdinateFlux::Allocate()
//! \see    STDG::OrdinateFlux::Deallocate()
//------------------------------------------------------------------------------------------------------------
STDG::OrdinateFlux & STDG::OrdinateFlux::Reconfigure (

    const OrdinateSet & ordinate_set
) {
    return Reconfigure( ordinate_set.GetAngOrder(), ordinate_set.GetOrdinateSymmetry(),
                        ordinate_set.GetOrdinateType() );
}


//------------------------------------------------------------------------------------------------------------
//! \brief  Reconfigures an STDG::OrdinateFlux object with the given set of parameters.
//!
//! \param[in]  spatial_params      Contains the parameters of the spatial discretization.
//! \param[in]  DG_degree           Maximum degree of DG basis polynomials.
//! \param[in]  ang_order           Order of the discrete ordinate set of type \pp{ordinate_set} to use.
//! \param[in]  symmetric_reduce    (optional) <br/>
//!                                 Specifies whether or not the number of ordinates should be reduced by
//!                                 symmetry. Only used for 2D problems. Defaults to false.
//! \param[in]  ordinate_type       (optional) <br/>
//!                                 OrdinateType specifying the type of discrete ordinate quadrature to use.
//!                                 Defaults to Gauss-Legendre quadrature in 1D and Chebyshev-Legendre
//!                                 quadrature in multi-D.
//!
//! \see    STDG::OrdinateFlux::Allocate()
//! \see    STDG::OrdinateFlux::Deallocate()
//------------------------------------------------------------------------------------------------------------
STDG::OrdinateFlux & STDG::OrdinateFlux::Reconfigure (

    const DomainDecomposition & spatial_params,
    const int64_t DG_degree,
    const int64_t ang_order,
    const bool symmetric_reduce,        // = false
    const OrdinateType ordinate_type    // Defaulted.
) {
    return Reconfigure( spatial_params, DG_degree, DG_degree, ang_order, symmetric_reduce, ordinate_type );
}


//------------------------------------------------------------------------------------------------------------
//! \brief  Reconfigures an STDG::OrdinateFlux object with the given set of parameters.
//!
//! \param[in]  spatial_params      Contains the parameters of the spatial discretization.
//! \param[in]  DG_degree_x_in      Maximum degree of DG basis polynomials in spatial dimensions.
//! \param[in]  DG_degree_t_in      Maximum degree of DG basis polynomials in temporal dimension.
//! \param[in]  ang_order           Order of the discrete ordinate set of type \pp{ordinate_set} to use.
//! \param[in]  symmetric_reduce    (optional) <br/>
//!                                 Specifies whether or not the number of ordinates should be reduced by
//!                                 symmetry. Only used for 2D problems. Defaults to false.
//! \param[in]  ordinate_type       (optional) <br/>
//!                                 OrdinateType specifying the type of discrete ordinate quadrature to use.
//!                                 Defaults to Gauss-Legendre quadrature in 1D and Chebyshev-Legendre
//!                                 quadrature in multi-D.
//!
//! \see    STDG::OrdinateFlux::Allocate()
//! \see    STDG::OrdinateFlux::Deallocate()
//------------------------------------------------------------------------------------------------------------
STDG::OrdinateFlux & STDG::OrdinateFlux::Reconfigure (

    const DomainDecomposition & spatial_params,
    const int64_t DG_degree_x_in,
    const int64_t DG_degree_t_in,
    const int64_t ang_order,
    const bool symmetric_reduce,        // = false
    const OrdinateType ordinate_type    // Defaulted.
) {

    Deallocate();

    DomainDecomposition::operator=( spatial_params );
    this->DG_degree_x = DG_degree_x_in;
    this->DG_degree_t = DG_degree_t_in;
    this->OrdinateSet::Reconfigure( ang_order, symmetric_reduce, ordinate_type );

    Allocate();

    return *this;
}


//------------------------------------------------------------------------------------------------------------
//! \brief  Reconfigures an STDG::OrdinateFlux object with the given set of parameters.
//!
//! \param[in]      spatial_params  Contains the parameters of the spatial discretization.
//! \param[in]      DG_degree       Maximum degree of DG basis polynomials.
//! \param[in]      ordinate_set    OrdinateSet object containing parameters to use.
//!
//! \see    STDG::OrdinateFlux::Allocate()
//! \see    STDG::OrdinateFlux::Deallocate()
//------------------------------------------------------------------------------------------------------------
STDG::OrdinateFlux & STDG::OrdinateFlux::Reconfigure (

    const DomainDecomposition & spatial_params,
    const OrdinateSet & ordinate_set,
    const int64_t DG_degree
) {
    return Reconfigure( spatial_params, ordinate_set, DG_degree, DG_degree );
}


//------------------------------------------------------------------------------------------------------------
//! \brief  Reconfigures an STDG::OrdinateFlux object with the given set of parameters.
//!
//! \param[in]      spatial_params  Contains the parameters of the spatial discretization.
//! \param[in]      DG_degree_x_in  Maximum degree of DG basis polynomials in spatial dimensions.
//! \param[in]      DG_degree_t_in  Maximum degree of DG basis polynomials in temporal dimension.
//! \param[in]      ordinate_set    OrdinateSet object containing parameters to use.
//!
//! \see    STDG::OrdinateFlux::Allocate()
//! \see    STDG::OrdinateFlux::Deallocate()
//------------------------------------------------------------------------------------------------------------
STDG::OrdinateFlux & STDG::OrdinateFlux::Reconfigure (

    const DomainDecomposition & spatial_params,
    const OrdinateSet & ordinate_set,
    const int64_t DG_degree_x_in,
    const int64_t DG_degree_t_in
) {
    return Reconfigure( spatial_params, DG_degree_x_in, DG_degree_t_in, ordinate_set.GetAngOrder(),
                        ordinate_set.GetOrdinateSymmetry(), ordinate_set.GetOrdinateType() );
}


//============================================================================================================
//=== STATIC MEMBER FUNCTIONS ================================================================================
//============================================================================================================


//------------------------------------------------------------------------------------------------------------
//! \brief  Creates and returns an OrdinateFlux object.
//!
//! \param[in]  domain_decomposition    Contains parameters of spatial discretization.
//! \param[in]  ordinate_set            Contains parameters of angular discretization.
//! \param[in]  input_list              List of additional parameters.
//!
//! \return     Returns the newly created object.
//------------------------------------------------------------------------------------------------------------
STDG::OrdinateFlux * STDG::OrdinateFlux::Create (

    const DomainDecomposition & domain_decomposition,
    const OrdinateSet & ordinate_set,
    const ParameterList & input_list
) {

    int64_t DG_degree_x = GetDGDegreeX( input_list );
    int64_t DG_degree_t = GetDGDegreeT( input_list );

    return new STDG::OrdinateFlux( domain_decomposition, ordinate_set, DG_degree_x, DG_degree_t );
}


//------------------------------------------------------------------------------------------------------------
//! \brief  Two STDG::OrdinateFlux objects "are matching" if all of their parameters are equal.
//!
//! \param[in]      first       The first of the two STDG::OrdinateFlux objects to compare.
//! \param[in]      second      The second of the two STDG::OrdinateFlux objects to compare.
//!
//! \return     Returns true if the STDG::OrdinateFlux objects "are matching" and false otherwise.
//------------------------------------------------------------------------------------------------------------
bool STDG::OrdinateFlux::AreMatching (

    const STDG::OrdinateFlux & first,
    const STDG::OrdinateFlux & second
) {

    return      DensityFunction::AreMatching( first, second )
            &&  OrdinateSet::AreMatching( first, second )
            &&  first.DG_degree_x == second.DG_degree_x
            &&  first.DG_degree_t == second.DG_degree_t;
}


//------------------------------------------------------------------------------------------------------------
//! \brief  Swaps the contents of two STDG::OrdinateFlux objects.
//!
//! Throws an error if the two objects are not matching.
//!
//! \todo   Implement STDG::OrdinateFlux::swap for case when objects !AreMatching. (Print warning when this
//!         happens.)
//------------------------------------------------------------------------------------------------------------
void STDG::OrdinateFlux::swap (

    STDG::OrdinateFlux & thing1,
    STDG::OrdinateFlux & thing2
) {

    if ( AreMatching( thing1, thing2 ) ) {

        std::swap( thing1.density,                  thing2.density                  );
        std::swap( thing1.halo_cells_dirty,         thing2.halo_cells_dirty         );
        std::swap( thing1.upwind_halo_cells_dirty,  thing2.upwind_halo_cells_dirty  );

    } else {

        std::string error_message = "Parameter mismatch between STDG::OrdinateFlux objects in '"
                                    + std::string(__func__) + "'.\n";

        PRINT_ERROR( error_message.c_str() )
        throw std::invalid_argument( error_message );
    }
}


//------------------------------------------------------------------------------------------------------------
//! \brief  \f$ \vec{y} \gets \beta \vec{y} + \alpha \vec{x} \f$.
//!
//! Evaluates the STDG::OrdinateFlux object \pp{x} at the endpoint \f$ t_{n+1} \f$ of the timestep interval
//! and performs an operation on the resulting RKDG::OrdinateFlux density similar to the BLAS routines
//! \c ?axpy with the destination \pp{y}.
//!
//! \param[in]      alpha       Value to scale \pp{x} by.
//! \param[in]      beta        Value to scale \pp{y} by.
//! \param[in]      x           STDG::OrdinateFlux to be evaluated, scaled, and added to \pp{y}.
//! \param[in,out]  y           RKDG::OrdinateFlux to be scaled and into which the result stored.
//------------------------------------------------------------------------------------------------------------
void STDG::OrdinateFlux::AXPY (

    const double alpha,
    const double beta,
    const STDG::OrdinateFlux & x,
    RKDG::OrdinateFlux & y
) {

    PRINT_STATUS( "Applying AXPY operation from STDG::OrdinateFlux object to RKDG::OrdinateFlux object.\n" )

# if defined (STRICT_CHECK)

    if (    !DomainDecomposition::AreMatching( x, y )
         || !OrdinateSet::AreMatching( x, y )
         || x.DG_degree_x != y.DG_degree
    ) {
        std::string error_message
            = "Parameter mismatch between STDG::OrdinateFlux  and RKDG::OrdinateFlux objects in '"
              + std::string(__func__) + "'.\n";

        PRINT_ERROR( error_message.c_str() )
        throw std::invalid_argument( error_message );
    }

# endif // if defined (STRICT_CHECK)

    Global::TMR_AF_axpy.Start();

    const int64_t (& nx) [SPACE_DIMS] = y.nx();

# if SPACE_DIMS == 1

    # pragma omp parallel for schedule(static)
    for ( int64_t i = 0; i <= nx[0] + 1;     ++i ) {
    for ( int64_t q = 0; q <  y.nq();        ++q ) {
    for ( int64_t d = 0; d <= x.DG_degree_x; ++d ) {

        y(q,i,d) *= beta;

        for ( int64_t s = 0; s <= x.DG_degree_t; ++s )
            y(q,i,d) += alpha * x(q,i,d,s);
    }}}

# elif SPACE_DIMS == 2

    # pragma omp parallel for collapse(2) schedule(static)
    for ( int64_t i = 0; i <= nx[0] + 1;     ++i ) {
    for ( int64_t j = 0; j <= nx[1] + 1;     ++j ) {
    for ( int64_t q = 0; q <  y.nq();        ++q ) {
    for ( int64_t d = 0; d <= x.DG_degree_x; ++d ) {
    for ( int64_t e = 0; e <= x.DG_degree_x; ++e ) {

        y(q,i,j,d,e) *= beta;

        for ( int64_t s = 0; s <= x.DG_degree_t; ++s )
            y(q,i,j,d,e) += alpha * x(q,i,j,d,e,s);
    }}}}}

# elif SPACE_DIMS == 3

    # pragma omp parallel for collapse(3) schedule(static)
    for ( int64_t i = 0; i <= nx[0] + 1;     ++i ) {
    for ( int64_t j = 0; j <= nx[1] + 1;     ++j ) {
    for ( int64_t k = 0; k <= nx[2] + 1;     ++k ) {
    for ( int64_t q = 0; q <  y.nq();        ++q ) {
    for ( int64_t d = 0; d <= x.DG_degree_x; ++d ) {
    for ( int64_t e = 0; e <= x.DG_degree_x; ++e ) {
    for ( int64_t f = 0; f <= x.DG_degree_x; ++f ) {

        y(q,i,j,k,d,e,f) *= beta;

        for ( int64_t s = 0; s <= x.DG_degree_t; ++s )
            y(q,i,j,k,d,e,f) += alpha * x(q,i,j,k,d,e,f,s);
    }}}}}}}

# endif // if SPACE_DIMS == ?

    y.halo_cells_dirty |= x.halo_cells_dirty;
    y.upwind_halo_cells_dirty |= x.upwind_halo_cells_dirty;

    Global::TMR_AF_axpy.Stop();
}


//------------------------------------------------------------------------------------------------------------
//! \brief  \f$ \vec{y} \gets \beta \vec{y} + alpha \vec{x} \f$.
//!
//! Performs an operation on the RKDG::OrdinateFlux \pp{x} and the space-time moments of the
//! STDG::OrdinateFlux \pp{y} which have degree zero with respect to time similar to the BLAS routines
//! \c ?axpy. The scalar \pp{beta} is applied to \e all space-time moments of \pp{y}.
//!
//! \param[in]      alpha       Value to scale \pp{x} by.
//! \param[in]      beta        Value to scale \pp{y} by.
//! \param[in]      x           RKDG::OrdinateFlux to be scaled and added to \pp{y}.
//! \param[in,out]  y           STDG::OrdinateFlux to be scaled and into which the result stored.
//------------------------------------------------------------------------------------------------------------
void STDG::OrdinateFlux::AXPY (

    const double alpha,
    const double beta,
    const RKDG::OrdinateFlux & x,
    STDG::OrdinateFlux & y
) {

    PRINT_STATUS( "Applying AXPY operation from RKDG::OrdinateFlux object to STDG::OrdinateFlux object.\n" )

# if defined (STRICT_CHECK)

    if (    !DomainDecomposition::AreMatching( x, y )
         || !OrdinateSet::AreMatching( x, y )
         || x.DG_degree != y.DG_degree_x
    ) {
        std::string error_message
            = "Parameter mismatch between STDG::OrdinateFlux  and RKDG::OrdinateFlux objects in '"
              + std::string(__func__) + "'.\n";

        PRINT_ERROR( error_message.c_str() )
        throw std::invalid_argument( error_message );
    }

# endif // if defined (STRICT_CHECK)

    Global::TMR_AF_axpy.Start();

    const int64_t (& nx) [SPACE_DIMS] = y.nx();

# if SPACE_DIMS == 1

    # pragma omp parallel for schedule(static)
    for ( int64_t i = 0; i <= nx[0] + 1;     ++i ) {
    for ( int64_t q = 0; q <  y.nq();        ++q ) {
    for ( int64_t d = 0; d <= y.DG_degree_x; ++d ) {

        for ( int64_t s = 0; s <= y.DG_degree_t; ++s )
            y(q,i,d,s) *= beta;

        y(q,i,d,0) += alpha * x(q,i,d);
    }}}

# elif SPACE_DIMS == 2

    # pragma omp parallel for collapse(2) schedule(static)
    for ( int64_t i = 0; i <= nx[0] + 1;     ++i ) {
    for ( int64_t j = 0; j <= nx[1] + 1;     ++j ) {
    for ( int64_t q = 0; q <  y.nq();        ++q ) {
    for ( int64_t d = 0; d <= y.DG_degree_x; ++d ) {
    for ( int64_t e = 0; e <= y.DG_degree_x; ++e ) {

        for ( int64_t s = 0; s <= y.DG_degree_t; ++s )
            y(q,i,j,d,e,s) *= beta;

        y(q,i,j,d,e,0) += alpha * x(q,i,j,d,e);
    }}}}}

# elif SPACE_DIMS == 3

    # pragma omp parallel for collapse(3) schedule(static)
    for ( int64_t i = 0; i <= nx[0] + 1;     ++i ) {
    for ( int64_t j = 0; j <= nx[1] + 1;     ++j ) {
    for ( int64_t k = 0; k <= nx[2] + 1;     ++k ) {
    for ( int64_t q = 0; q <  y.nq();        ++q ) {
    for ( int64_t d = 0; d <= y.DG_degree_x; ++d ) {
    for ( int64_t e = 0; e <= y.DG_degree_x; ++e ) {
    for ( int64_t f = 0; f <= y.DG_degree_x; ++f ) {

        for ( int64_t s = 0; s <= y.DG_degree_t; ++s )
            y(q,i,j,k,d,e,f,s) *= beta;

        y(q,i,j,k,d,e,f,0) += alpha * x(q,i,j,k,d,e,f);
    }}}}}}}

# endif // if SPACE_DIMS == ?

    y.halo_cells_dirty |= x.halo_cells_dirty;
    y.upwind_halo_cells_dirty |= x.upwind_halo_cells_dirty;

    Global::TMR_AF_axpy.Stop();
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
void STDG::OrdinateFlux::Print (

    const std::string & prefix // = "  "

) const {

    PRINT_STATUS( "Executing STDG::OrdinateFlux::%s.\n", __func__ )

    PRINT_LOG( "\n" )

    PRINT_LOG( "%s%-*s % " PRId64 "\n", prefix.c_str(),
               Global::col_width, "DG degree x:", this->DG_degree_x )

    PRINT_LOG( "%s%-*s % " PRId64 "\n", prefix.c_str(),
               Global::col_width, "DG degree t:", this->DG_degree_t )

    this->Abstract::OrdinateFlux::Print( prefix );
}


//------------------------------------------------------------------------------------------------------------
//! \brief  Updates ghost cells on reflecting boundaries.
//!
//! \todo   3D Implementation.
//------------------------------------------------------------------------------------------------------------
STDG::OrdinateFlux & STDG::OrdinateFlux::ReflectBoundaries ( void ) {

    PRINT_STATUS( "Updating ghost cells with reflecting boundary conditions.\n" )

    # pragma omp master
    Global::TMR_AF_boundary.Start();

    /*  Get thread count and thread IDs before entering parallel region below.
     *
     *  If this function is called from within a parallel region (as is the case for the mesh-sweeping
     *  routines) then each thread executes the function call in its own thread-local scope. Because of this,
     *  the loops below cannot be divided amongst the threads using standard worksharing constructs such as
     *  "omp for". Hence this implementation manually divides the loop iteration amongst the threads.
     *
     *  The manual division is performed using the thread count and thread IDs. There are two cases:
     *
     *      (i) When this function is called from a sequential portion of the code, the parallel region
     *          below is used to create a team of threads to execute the reflect operation. In this case
     *          the thread count and thread IDs must be determined from within the parallel region.
     *
     *     (ii) When this function is called from within a parallel region (the "outer" parallel region)
     *          the parallel region below (the "inner" parallel region) does not spawn additional threads.
     *          However, function calls to the OpenMP runtime library are always bound to the innermost
     *          parallel region, regardless of whether additional threads are spawned there. In order to
     *          obtain the correct thread counts and thread IDs, OpenMP library calls must be made before
     *          entering the innermost parallel region.
     */
    const bool in_parallel =
        # if defined (_OPENMP)
            omp_in_parallel();
        # else
            false;
        # endif

    const int64_t num_threads_outer =
        # if defined (_OPENMP)
            omp_get_num_threads();
        # else
            1;
        # endif

    const int64_t tid_outer =
        # if defined (_OPENMP)
            omp_get_thread_num();
        # else
            0;
        # endif

    // Spawn parallel threads here only if this function was called from a sequential section of the code.
    # pragma omp parallel if( !in_parallel )
    {
        /*  Get thread count and thread IDs for innermost parallel region. */
        const int64_t num_threads_inner =
            # if defined (_OPENMP)
                omp_get_num_threads();
            # else
                1;
            # endif

        const int64_t tid_inner =
            # if defined (_OPENMP)
                omp_get_thread_num();
            # else
                0;
            # endif

        /*  Determine "real" thread count and thread IDs.
         *
         *  If this function was called from within a parallel region, then the thread count and thread IDs
         *  that we want are num_threads_outer and tid_outer, respectively. If this function was called from
         *  a sequential portion of the code, then the thread count and thread IDs that we want are
         *  num_threads_inner and tid_inner, respectively.
         */
        const int64_t num_threads = ( in_parallel ? num_threads_outer : num_threads_inner );
        const int64_t tid         = ( in_parallel ? tid_outer         : tid_inner         );

    # if SPACE_DIMS == 1
        ((void)(num_threads));  // To avoid compiler warning in 1D code.
    # endif


        // --- X_Min ---------------------------------------------------------------------------------- //

        if (    BitmaskHasAll( Global::reflecting_boundaries, BoundaryEdge::X_Min )
        # if defined (ENABLE_MPI)
             && this->MPI_block_coords(0) == 0
        # endif
        ) {

        # if SPACE_DIMS == 1

            if ( tid == 0 ) {

                for ( int64_t q = 0; q <  this->nq();        ++q ) {
                for ( int64_t d = 0; d <= this->DG_degree_x; ++d ) {
                for ( int64_t s = 0; s <= this->DG_degree_t; ++s ) {

                    (*this)(q,0,d,s) = neg1pow(d) * (*this)( this->nq() - q, 1,d,s);
                }}}
            }

        # elif SPACE_DIMS == 2

            const int64_t Q = this->nx(1) / num_threads;    // Quotient and remainder, respectively.
            const int64_t R = this->nx(1) % num_threads;    //

            const int64_t j_start = 1 + tid * Q + std::min( tid, R );
            const int64_t j_stop  = j_start + Q - ( tid < R ? 0 : 1 );

            for ( int64_t quad : { 0, 1 } ) {
            for ( int64_t q = this->Quadrants(quad); q < this->Quadrants(quad + 1); ++q ) {

                for ( int64_t j = j_start; j <= j_stop; ++j ) {
                for ( int64_t d = 0; d <= this->DG_degree_x; ++d ) {
                for ( int64_t e = 0; e <= this->DG_degree_x; ++e ) {
                for ( int64_t s = 0; s <= this->DG_degree_t; ++s ) {

                    (*this)(q, 0, j,d,e,s) = neg1pow(d) * (*this)( ReflectMap(0,q), 1, j,d,e,s);
                }}}}
            }}

        # elif SPACE_DIMS == 3

            # warning "Implementation of STDG::OrdinateFlux>::ReflectBoundaries incomplete."

        # endif // if SPACE_DIMS == ?
        }

        // --- X_Max ---------------------------------------------------------------------------------- //

        if (    BitmaskHasAll( Global::reflecting_boundaries, BoundaryEdge::X_Max )
        # if defined (ENABLE_MPI)
             && this->MPI_block_coords(0) == this->MPI_num_blocks(0) - 1
        # endif
        ) {

        # if SPACE_DIMS == 1

            if ( tid == 0 ) {

                for ( int64_t q = 0; q <  this->nq();        ++q ) {
                for ( int64_t d = 0; d <= this->DG_degree_x; ++d ) {
                for ( int64_t s = 0; s <= this->DG_degree_t; ++s ) {

                    (*this)(q, this->nx(0) + 1, d,s) = neg1pow(d) * (*this)( this->nq() - q, this->nx(0), d,s);
                }}}
            }

        # elif SPACE_DIMS == 2

            const int64_t Q = this->nx(1) / num_threads;    // Quotient and remainder, respectively.
            const int64_t R = this->nx(1) % num_threads;    //

            const int64_t j_start = 1 + tid * Q + std::min( tid, R );
            const int64_t j_stop  = j_start + Q - ( tid < R ? 0 : 1 );

            for ( int64_t quad : { 2, 3 } ) {
            for ( int64_t q = this->Quadrants(quad); q < this->Quadrants(quad + 1); ++q ) {

                for ( int64_t j = j_start; j <= j_stop; ++j ) {
                for ( int64_t d = 0; d <= this->DG_degree_x; ++d ) {
                for ( int64_t e = 0; e <= this->DG_degree_x; ++e ) {
                for ( int64_t s = 0; s <= this->DG_degree_t; ++s ) {

                    (*this)(q, this->nx(0) + 1, j,d,e,s)
                        = neg1pow(d) * (*this)( ReflectMap(0,q), this->nx(0), j,d,e,s);
                }}}}
            }}

        # elif SPACE_DIMS == 3

            # warning "Implementation of STDG::OrdinateFlux>::ReflectBoundaries incomplete."

        # endif // if SPACE_DIMS == ?
        }

    # if SPACE_DIMS >= 2

        // --- Y_Min ---------------------------------------------------------------------------------- //

        if (    BitmaskHasAll( Global::reflecting_boundaries, BoundaryEdge::Y_Min )
        # if defined (ENABLE_MPI)
             && this->MPI_block_coords(1) == 0
        # endif
        ) {

        # if SPACE_DIMS == 2

            const int64_t Q = this->nx(0) / num_threads;    // Quotient and remainder, respectively.
            const int64_t R = this->nx(0) % num_threads;    //

            const int64_t i_start = 1 + tid * Q + std::min( tid, R );
            const int64_t i_stop  = i_start + Q - ( tid < R ? 0 : 1 );

            for ( int64_t quad : { 0, 2 } ) {
            for ( int64_t q = this->Quadrants(quad); q < this->Quadrants(quad + 1); ++q ) {

                for ( int64_t i = i_start; i <= i_stop; ++i ) {
                for ( int64_t d = 0; d <= this->DG_degree_x; ++d ) {
                for ( int64_t e = 0; e <= this->DG_degree_x; ++e ) {
                for ( int64_t s = 0; s <= this->DG_degree_t; ++s ) {

                    (*this)(q,i, 0, d,e,s) = neg1pow(e) * (*this)( ReflectMap(1,q), i, 1, d,e,s);
                }}}}
            }}

        # elif SPACE_DIMS == 3

            # warning "Implementation of STDG::OrdinateFlux>::ReflectBoundaries incomplete."

        # endif // if SPACE_DIMS == ?
        }

        // --- Y_Max ---------------------------------------------------------------------------------- //

        if (    BitmaskHasAll( Global::reflecting_boundaries, BoundaryEdge::Y_Max )
        # if defined (ENABLE_MPI)
             && this->MPI_block_coords(1) == this->MPI_num_blocks(1) - 1
        # endif
        ) {

        # if SPACE_DIMS == 2

            const int64_t Q = this->nx(0) / num_threads;    // Quotient and remainder, respectively.
            const int64_t R = this->nx(0) % num_threads;    //

            const int64_t i_start = 1 + tid * Q + std::min( tid, R );
            const int64_t i_stop  = i_start + Q - ( tid < R ? 0 : 1 );

            for ( int64_t quad : { 0, 1, 2, 3 } ) {
            for ( int64_t q = this->Quadrants(quad); q < this->Quadrants(quad + 1); ++q ) {

                for ( int64_t i = i_start; i <= i_stop; ++i ) {
                for ( int64_t d = 0; d <= this->DG_degree_x; ++d ) {
                for ( int64_t e = 0; e <= this->DG_degree_x; ++e ) {
                for ( int64_t s = 0; s <= this->DG_degree_t; ++s ) {

                    (*this)(q,i, this->nx(1) + 1, d,e,s)
                        = neg1pow(e) * (*this)( ReflectMap(1,q), i, this->nx(1), d,e,s);
                }}}}
            }}

        # elif SPACE_DIMS == 3

            # warning "Implementation of STDG::OrdinateFlux>::ReflectBoundaries incomplete."

        # endif // if SPACE_DIMS == ?
        }

    # endif // if SPACE_DIMS >= 2
    } // end parallel region.

    # pragma omp master
    Global::TMR_AF_boundary.Stop();

    return *this;
}


//------------------------------------------------------------------------------------------------------------
//! \brief  Returns the stride in the global array between subsequent spatial cells, moving in the specified
//!         dimension.
//!
//! \param[in]  dim     Dimension to compute stride for.
//!
//! \return     Returns the stride in the global array for the specified dimension.
//------------------------------------------------------------------------------------------------------------
int64_t STDG::OrdinateFlux::CellStride (

    const int64_t dim

) const {

    switch ( dim ) {

        case 0:
            return
            # if SPACE_DIMS == 1
                Index(0,1,0,0) - Index(0,0,0,0);
            # elif SPACE_DIMS == 2
                Index(0,1,0,0,0,0) - Index(0,0,0,0,0,0);
            # elif SPACE_DIMS == 3
                Index(0,1,0,0,0,0,0,0) - Index(0,0,0,0,0,0,0,0);
            # endif

    # if SPACE_DIMS >= 2

        case 1:
            return
            # if SPACE_DIMS == 2
                Index(0,0,1,0,0,0) - Index(0,0,0,0,0,0);
            # elif SPACE_DIMS == 3
                Index(0,0,1,0,0,0,0,0) - Index(0,0,0,0,0,0,0,0);
            # endif

    # endif // if SPACE_DIMS >= 2

    # if SPACE_DIMS == 3

        case 2:
            return Index(0,0,0,1,0,0,0,0) - Index(0,0,0,0,0,0,0,0);

    # endif // if SPACE_DIMS == 3

        default:
        {   std::string error_message =   "Spatial dimension "
                                        + std::to_string(dim)
                                        + " not in valid range [0,"
                                        + std::to_string( SPACE_DIMS )
                                        + ") for "
                                        + std::string(__func__)
                                        + ".\n";

            PRINT_ERROR( error_message.c_str() )
            throw std::out_of_range( error_message );
        }
    }
}


//------------------------------------------------------------------------------------------------------------
//! \brief  Returns a pointer to the location within DensityFunction::density where the coefficients for the
//!         given spatial cell begin.
//------------------------------------------------------------------------------------------------------------
double * STDG::OrdinateFlux::PointerAtCell (

    const int64_t i
# if SPACE_DIMS >= 2
,   const int64_t j
# endif
# if SPACE_DIMS == 3
,   const int64_t k
# endif

) const {

    return this->density + Index
                            # if SPACE_DIMS == 1
                                (0,i,0,0)
                            # elif SPACE_DIMS == 2
                                (0,i,j,0,0,0)
                            # elif SPACE_DIMS == 3
                                (0,i,j,k,0,0,0,0)
                            # endif
                            ;
}


//------------------------------------------------------------------------------------------------------------
//! \brief  Returns a pointer to the location within DensityFunction::density where the coefficients for the
//!         given spatial cell begin.
//------------------------------------------------------------------------------------------------------------
double * STDG::OrdinateFlux::PointerAtOrdinateCell (

    const int64_t q
,   const int64_t i
# if SPACE_DIMS >= 2
,   const int64_t j
# endif
# if SPACE_DIMS == 3
,   const int64_t k
# endif

) const {

    return this->density + Index
                            # if SPACE_DIMS == 1
                                (q,i,0,0)
                            # elif SPACE_DIMS == 2
                                (q,i,j,0,0,0)
                            # elif SPACE_DIMS == 3
                                (q,i,j,k,0,0,0,0)
                            # endif
                            ;
}


# if SPACE_DIMS == 2

//------------------------------------------------------------------------------------------------------------
//! \brief  Returns a pointer to the location within DensityFunction::density where the coefficients for the
//!         given spatial cell begin.
//------------------------------------------------------------------------------------------------------------
double * STDG::OrdinateFlux::PointerAtQuadrantCell (

    const int64_t quad,
    const int64_t i,
    const int64_t j

) const {

    return this->density + Index( Quadrants(quad), i,j,0,0,0);
}

# elif SPACE_DIMS == 3

//------------------------------------------------------------------------------------------------------------
//! \brief  Returns a pointer to the location within DensityFunction::density where the coefficients for the
//!         given spatial cell begin.
//------------------------------------------------------------------------------------------------------------
double * STDG::OrdinateFlux::PointerAtOctantCell (

    const int64_t oct,
    const int64_t i,
    const int64_t j

) const {

    return this->density + Index( Octants(oct), i,j,k,0,0,0,0);
}

# endif // if SPACE_DIMS == ?
