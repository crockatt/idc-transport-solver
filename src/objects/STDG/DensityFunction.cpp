//------------------------------------------------------------------------------------------------------------
//! \file   objects/STDG/DensityFunction.cpp
//! \brief  Implementation of STDG::DensityFunction class.
//!
//! \author Michael Crockatt
//! \date   February 2018
//------------------------------------------------------------------------------------------------------------


# include <cinttypes>
# include <cmath>
# include <cstdint>
# include <cstdlib>
# include <cstring>
# include <limits>

# include "objects/STDG/DensityFunction.hpp"
# include "utils/CLog.hpp"
# include "utils/global.hpp"
# include "utils/Quadrule/Quadrule.hpp"


using namespace STDG;


//============================================================================================================
//=== CONSTRUCTORS, DESTRUCTORS, AND ASSOCIATED HELPER ROUTINES ==============================================
//============================================================================================================


//------------------------------------------------------------------------------------------------------------
//! \brief  Constructs an empty STDG::DensityFunction object with zero-initialized parameters.
//------------------------------------------------------------------------------------------------------------
DensityFunction::DensityFunction ( void ) :

    DomainDecomposition(),
    Abstract::DensityFunction(),
    DG_degree_x {},
    DG_degree_t {}
{}


//------------------------------------------------------------------------------------------------------------
//! \brief  Constructs an STDG::DensityFunction object with given parameters.
//!
//! The elements of internal data vectors are initialized to zero.
//!
//! \param[in]      spatial_params  Contains the parameters of the spatial discretization.
//! \param[in]      DG_degree       Maximum degree of DG basis polynomials (with respect to both space and
//!                                 time).
//!
//! \see    STDG::DensityFunction::Allocate()
//------------------------------------------------------------------------------------------------------------
DensityFunction::DensityFunction (

    const DomainDecomposition & spatial_params,
    const int64_t DG_degree

) : // DensityFunction( spatial_params, DG_degree, DG_degree ) {}
    DomainDecomposition( spatial_params ),
    Abstract::DensityFunction( spatial_params ),
    DG_degree_x( DG_degree ),
    DG_degree_t( DG_degree )
{
    Allocate();
}


//------------------------------------------------------------------------------------------------------------
//! \brief  Constructs an STDG::DensityFunction object with given parameters.
//!
//! The elements of internal data vectors are initialized to zero.
//!
//! \param[in]      spatial_params  Contains the parameters of the spatial discretization.
//! \param[in]      DG_degree_x_in  Maximum degree of DG basis polynomials in spatial dimensions.
//! \param[in]      DG_degree_t_in  Maximum degree of DG basis polynomials in temporal dimension.
//!
//! \see    STDG::DensityFunction::Allocate()
//------------------------------------------------------------------------------------------------------------
DensityFunction::DensityFunction (

    const DomainDecomposition & spatial_params,
    const int64_t DG_degree_x_in,
    const int64_t DG_degree_t_in
) :
    DomainDecomposition( spatial_params ),
    Abstract::DensityFunction( spatial_params ),
    DG_degree_x( DG_degree_x_in ),
    DG_degree_t( DG_degree_t_in )
{
    Allocate();
}


//------------------------------------------------------------------------------------------------------------
//! \brief  Sets appropriate values for DensityFunction::DOF_per_cell, DensityFunction::dimof_density and
//!         DensityFunction::sizeof_density.
//------------------------------------------------------------------------------------------------------------
void DensityFunction::SetDensityDimensions ( void ) {

    // Set value for DOF per cell.
    this->DOF_per_cell = IntPow( this->DG_degree_x + 1, SPACE_DIMS ) * (this->DG_degree_t + 1);

    // Call up class hierarchy.
    this->Abstract::DensityFunction::SetDensityDimensions();
}


//------------------------------------------------------------------------------------------------------------
//! \brief  Reconfigures an STDG::DensityFunction object with a different set of parameters.
//!
//! \param[in]      DG_degree       Maximum degree of DG basis polynomials.
//!
//! \see    STDG::DensityFunction::Allocate()
//! \see    STDG::DensityFunction::Deallocate()
//------------------------------------------------------------------------------------------------------------
DensityFunction & DensityFunction::Reconfigure (

    const int64_t DG_degree
) {
    return Reconfigure( DG_degree, DG_degree );
}


//------------------------------------------------------------------------------------------------------------
//! \brief  Reconfigures an STDG::DensityFunction object with a different set of parameters.
//!
//! \param[in]      DG_degree_x_in  Maximum degree of DG basis polynomials in spatial dimensions.
//! \param[in]      DG_degree_t_in  Maximum degree of DG basis polynomials in temporal dimension.
//!
//! \see    STDG::DensityFunction::Allocate()
//! \see    STDG::DensityFunction::Deallocate()
//------------------------------------------------------------------------------------------------------------
DensityFunction & DensityFunction::Reconfigure (

    const int64_t DG_degree_x_in,
    const int64_t DG_degree_t_in
) {

    Deallocate();

    this->DG_degree_x = DG_degree_x_in;
    this->DG_degree_t = DG_degree_t_in;

    Allocate();

    return *this;
}


//------------------------------------------------------------------------------------------------------------
//! \brief  Reconfigures an STDG::DensityFunction object with a different set of parameters.
//!
//! \param[in]      spatial_params  Contains the parameters of the spatial discretization.
//! \param[in]      DG_degree       Maximum degree of DG basis polynomials.
//!
//! \see    STDG::DensityFunction::Allocate()
//! \see    STDG::DensityFunction::Deallocate()
//------------------------------------------------------------------------------------------------------------
DensityFunction & DensityFunction::Reconfigure (

    const DomainDecomposition & spatial_params,
    const int64_t DG_degree
) {
    return Reconfigure( spatial_params, DG_degree, DG_degree );
}


//------------------------------------------------------------------------------------------------------------
//! \brief  Reconfigures an STDG::DensityFunction object with a different set of parameters.
//!
//! \param[in]      spatial_params  Contains the parameters of the spatial discretization.
//! \param[in]      DG_degree_x_in  Maximum degree of DG basis polynomials in spatial dimensions.
//! \param[in]      DG_degree_t_in  Maximum degree of DG basis polynomials in temporal dimension.
//!
//! \see    STDG::DensityFunction::Allocate()
//! \see    STDG::DensityFunction::Deallocate()
//------------------------------------------------------------------------------------------------------------
DensityFunction & DensityFunction::Reconfigure (

    const DomainDecomposition & spatial_params,
    const int64_t DG_degree_x_in,
    const int64_t DG_degree_t_in
) {

    Deallocate();

    DomainDecomposition::operator=( spatial_params );
    this->DG_degree_x = DG_degree_x_in;
    this->DG_degree_t = DG_degree_t_in;

    Allocate();

    return *this;
}


//============================================================================================================
//=== STATIC MEMBER FUNCTIONS ================================================================================
//============================================================================================================


//------------------------------------------------------------------------------------------------------------
//! \brief  Creates and returns a DensityFunction object.
//!
//! \param[in]  domain_decomposition    Contains parameters of spatial discretization.
//! \param[in]  input_list              List of additional parameters.
//!
//! \return     Returns the newly created object.
//------------------------------------------------------------------------------------------------------------
DensityFunction * DensityFunction::Create (

    const DomainDecomposition & domain_decomposition,
    const ParameterList & input_list
) {

    int64_t DG_degree_x = GetDGDegreeX( input_list );
    int64_t DG_degree_t = GetDGDegreeT( input_list );

    return new DensityFunction( domain_decomposition, DG_degree_x, DG_degree_t );
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

    PRINT_STATUS( "Executing STDG::DensityFunction::%s.\n", __func__ )

    PRINT_LOG( "\n" )

    PRINT_LOG( "%s%-*s % " PRId64 "\n", prefix.c_str(),
               Global::col_width, "DG degree x:", this->DG_degree_x )

    PRINT_LOG( "%s%-*s % " PRId64 "\n", prefix.c_str(),
               Global::col_width, "DG degree t:", this->DG_degree_t )

    this->Abstract::DensityFunction::Print( prefix );
}


//------------------------------------------------------------------------------------------------------------
//! \brief  Deallocates memory at all internal pointers and zeros all parameters of an STDG::DensityFunction
//!         object.
//!
//! \see    Abstract::DensityFunction::Zero()
//------------------------------------------------------------------------------------------------------------
DensityFunction & DensityFunction::Zero ( void ) {

    this->DG_degree_x = 0;
    this->DG_degree_t = 0;
    Abstract::DensityFunction::Zero();

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
int64_t DensityFunction::CellStride (

    const int64_t dim

) const {

    switch ( dim ) {

        case 0:
            return
            # if SPACE_DIMS == 1
                Index(1,0,0) - Index(0,0,0);
            # elif SPACE_DIMS == 2
                Index(1,0,0,0,0) - Index(0,0,0,0,0);
            # elif SPACE_DIMS == 3
                Index(1,0,0,0,0,0,0) - Index(0,0,0,0,0,0,0);
            # endif

    # if SPACE_DIMS >= 2

        case 1:
            return
            # if SPACE_DIMS == 2
                Index(0,1,0,0,0) - Index(0,0,0,0,0);
            # elif SPACE_DIMS == 3
                Index(0,1,0,0,0,0,0) - Index(0,0,0,0,0,0,0);
            # endif

    # endif // if SPACE_DIMS >= 2

    # if SPACE_DIMS == 3

        case 2:
            return Index(0,0,1,0,0,0,0) - Index(0,0,0,0,0,0,0);

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
double * DensityFunction::PointerAtCell (

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
                                (i,0,0)
                            # elif SPACE_DIMS == 2
                                (i,j,0,0,0)
                            # elif SPACE_DIMS == 3
                                (i,j,k,0,0,0,0)
                            # endif
                            ;
}
