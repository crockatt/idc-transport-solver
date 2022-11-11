//------------------------------------------------------------------------------------------------------------
//! \file   objects/RKDG/CrossSection.cpp
//! \file   CrossSection.cpp
//! \brief  Implementation of RKDG::CrossSection class.
//!
//! \author Michael Crockatt
//! \date   February 2018
//------------------------------------------------------------------------------------------------------------


# include <cinttypes>
# include <cstdlib>

# include "utils/CLog.hpp"
# include "utils/global.hpp"
# include "objects/RKDG/CrossSection.hpp"


using namespace RKDG;


//============================================================================================================
//=== CONSTRUCTORS, DESTRUCTORS, AND ASSOCIATED HELPER ROUTINES ==============================================
//============================================================================================================


//------------------------------------------------------------------------------------------------------------
//! \brief  Constructs an empty RKDG::CrossSection object with zero-initialized parameters.
//------------------------------------------------------------------------------------------------------------
CrossSection::CrossSection ( void ) :

    Abstract::DensityFunction(),
    RKDG::DensityFunction(),
    tensor { nullptr },
    sizeof_tensor {},
    dimof_tensor {}
{}


//------------------------------------------------------------------------------------------------------------
//! \brief  Destructor for class RKDG::CrossSection.
//------------------------------------------------------------------------------------------------------------
CrossSection::~CrossSection ( void ) {

    Deallocate();
}


//------------------------------------------------------------------------------------------------------------
//! \brief  Sets parameters and allocates memory at internal pointers.
//------------------------------------------------------------------------------------------------------------
void CrossSection::Allocate ( void ) {

    DensityFunction::Allocate();

    this->dimof_tensor = this->dimof_density
        # if SPACE_DIMS == 1
            * (this->DG_degree + 1);
        # elif SPACE_DIMS == 2
            * (this->DG_degree + 1)*(this->DG_degree + 1);
        # elif SPACE_DIMS == 3
            * (this->DG_degree + 1)*(this->DG_degree + 1)*(this->DG_degree + 1);
        # endif

    this->sizeof_tensor = this->dimof_tensor * sizeof(double);

    // Allocate memory for tensor array.
    this->tensor = (double *)
        # if defined (ENABLE_HWLOC) && defined (USE_HWLOC_ALLOC)

            hwloc_alloc_membind( Global::machine_topology, this->sizeof_tensor, Global::active_core_mask,
                                 HWLOC_MEMBIND_FIRSTTOUCH, HWLOC_MEMBIND_PROCESS );

        # elif defined (USE_ALIGNED_ALLOC)
            aligned_alloc( ALIGNED_ALLOC_ALIGNMENT, this->sizeof_tensor );
        # else
            std::malloc( this->sizeof_tensor );
        # endif

    ComputeTensor();
}


//------------------------------------------------------------------------------------------------------------
//! \brief  Deallocates memory at all internal pointers.
//!
//! \see    DensityFunction::Deallocate()
//------------------------------------------------------------------------------------------------------------
void CrossSection::Deallocate ( void ) {

    if ( this->tensor != nullptr ) {

        # if defined (ENABLE_HWLOC) && defined (USE_HWLOC_ALLOC)
            hwloc_free( Global::machine_topology, this->tensor, this->sizeof_tensor );
        # else
            std::free( this->tensor );
        # endif

        this->tensor = nullptr;
    }

    this->sizeof_tensor = 0;
    this->dimof_tensor = 0;

    DensityFunction::Deallocate();
}


//============================================================================================================
//=== ADDITIONAL MEMBER FUNCTIONS ============================================================================
//============================================================================================================


//------------------------------------------------------------------------------------------------------------
//! \brief  Deallocates memory at all internal pointers and zeros all parameters of an RKDG::CrossSection
//!         object.
//!
//! \see    DensityFunction::Zero()
//------------------------------------------------------------------------------------------------------------
CrossSection & CrossSection::Zero ( void ) {

    this->sizeof_tensor = 0;
    this->dimof_tensor = 0;

    DensityFunction::Zero();

    return *this;
}


//------------------------------------------------------------------------------------------------------------
//! \brief  Computes the coefficients of the tensor contraction between the coefficients at
//!         DensityFunction::density and the 3-tensor of Legendre triple product integrals, storing the
//!         resulting 2-tensor at RKDG::CrossSection::tensor.
//!
//! \attention  Currently this function requires that the Global::TPI array be properly allocated and setup
//!             prior to calling this function.
//!
//! \todo   Make the computation of Global::TPI array adaptive -- i.e., compute it from here only if it is
//!         necessary.
//------------------------------------------------------------------------------------------------------------
void CrossSection::ComputeTensor ( void ) {

    const int64_t (& nx) [SPACE_DIMS] = this->nx();

# if SPACE_DIMS == 1

    #pragma omp parallel for schedule(static)
    for ( int64_t i = 0; i <= nx[0] + 1; ++i ) {

        for ( int64_t a = 0; a <= this->DG_degree; ++a ) {
        for ( int64_t d = 0; d <= this->DG_degree; ++d ) {

            (*this)(i,a,d) = 0.0;
        }}

        for ( int64_t a = 0; a <= this->DG_degree; ++a ) {
        for ( int64_t d = 0; d <= this->DG_degree; ++d ) {
        for ( int64_t u = 0; u <= this->DG_degree; ++u ) {

            (*this)(i,a,d) += (*this)(i,u) * Global::TPI[ ITPI(a,d,u) ];
        }}}
    }

# elif SPACE_DIMS == 2

    #pragma omp parallel for collapse(2) schedule(static)
    for ( int64_t i = 0; i <= nx[0] + 1; ++i ) {
    for ( int64_t j = 0; j <= nx[1] + 1; ++j ) {

        for ( int64_t a = 0; a <= this->DG_degree; ++a ) {
        for ( int64_t b = 0; b <= this->DG_degree; ++b ) {
        for ( int64_t d = 0; d <= this->DG_degree; ++d ) {
        for ( int64_t e = 0; e <= this->DG_degree; ++e ) {

            (*this)(i,j,a,b,d,e) = 0.0;
        }}}}

        for ( int64_t a = 0; a <= this->DG_degree; ++a ) {
        for ( int64_t b = 0; b <= this->DG_degree; ++b ) {
        for ( int64_t d = 0; d <= this->DG_degree; ++d ) {
        for ( int64_t e = 0; e <= this->DG_degree; ++e ) {
        for ( int64_t u = 0; u <= this->DG_degree; ++u ) {
        for ( int64_t v = 0; v <= this->DG_degree; ++v ) {

            (*this)(i,j,a,b,d,e) += (*this)(i,j,u,v) * Global::TPI[ ITPI(a,d,u) ]
                                                     * Global::TPI[ ITPI(b,e,v) ];
        }}}}}}
    }}

# elif SPACE_DIMS == 3

    #pragma omp parallel for collapse(3) schedule(static)
    for ( int64_t i = 0; i <= nx[0] + 1; ++i ) {
    for ( int64_t j = 0; j <= nx[1] + 1; ++j ) {
    for ( int64_t k = 0; k <= nx[2] + 1; ++k ) {

        for ( int64_t a = 0; a <= this->DG_degree; ++a ) {
        for ( int64_t b = 0; b <= this->DG_degree; ++b ) {
        for ( int64_t c = 0; c <= this->DG_degree; ++c ) {
        for ( int64_t d = 0; d <= this->DG_degree; ++d ) {
        for ( int64_t e = 0; e <= this->DG_degree; ++e ) {
        for ( int64_t f = 0; f <= this->DG_degree; ++f ) {

            (*this)(i,j,k,a,b,c,d,e,f) = 0.0;
        }}}}}}

        for ( int64_t a = 0; a <= this->DG_degree; ++a ) {
        for ( int64_t b = 0; b <= this->DG_degree; ++b ) {
        for ( int64_t c = 0; c <= this->DG_degree; ++c ) {
        for ( int64_t d = 0; d <= this->DG_degree; ++d ) {
        for ( int64_t e = 0; e <= this->DG_degree; ++e ) {
        for ( int64_t f = 0; f <= this->DG_degree; ++f ) {
        for ( int64_t u = 0; u <= this->DG_degree; ++u ) {
        for ( int64_t v = 0; v <= this->DG_degree; ++v ) {
        for ( int64_t w = 0; w <= this->DG_degree; ++w ) {

            (*this)(i,j,k,a,b,c,d,e,f) += (*this)(i,j,k,u,v,w) * Global::TPI[ ITPI(a,d,u) ]
                                                               * Global::TPI[ ITPI(b,e,v) ]
                                                               * Global::TPI[ ITPI(c,f,w) ];
        }}}}}}}}}
    }}}

# endif // if SPACE_DIMS == ?
}
