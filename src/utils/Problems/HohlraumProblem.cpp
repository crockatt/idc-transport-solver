//------------------------------------------------------------------------------------------------------------
//! \file   utils/Problems/HohlraumProblem.cpp
//! \brief  Contains implementation of the HohlraumProblem class.
//!
//! \author Michael M. Crockatt
//! \date   June 2018
//------------------------------------------------------------------------------------------------------------

# include "utils/CLog.hpp"
# include "utils/global.hpp"
# include "utils/Problems/HohlraumProblem.hpp"


# if SPACE_DIMS == 2 || defined (DOXYCOMPILE)


//============================================================================================================
//=== CONSTRUCTORS AND DESTRUCTOR ============================================================================
//============================================================================================================


//------------------------------------------------------------------------------------------------------------
//! \brief  Construct object from given parameters.
//!
//! \param[in]  params  Contains parameters used to construct object.
//------------------------------------------------------------------------------------------------------------
HohlraumProblem::HohlraumProblem (

    const ParameterList & params
) :
    Problem( params )
{
    PRINT_STATUS( "Executing HohlraumProblem::%s.\n", __func__ )

    /* empty */
}


//------------------------------------------------------------------------------------------------------------
//! \brief  Destructor for HohlraumProblem class.
//------------------------------------------------------------------------------------------------------------
HohlraumProblem::~HohlraumProblem ( void ) {

    PRINT_STATUS( "Executing HohlraumProblem::%s.\n", __func__ )

    /* empty */
}


//============================================================================================================
//=== PUBLIC INTERFACE ROUTINES ==============================================================================
//============================================================================================================


//------------------------------------------------------------------------------------------------------------
//! \brief  Returns the string descriptor of the class type.
//------------------------------------------------------------------------------------------------------------
const std::string HohlraumProblem::Descriptor( void ) {

    PRINT_STATUS( "Executing HohlraumProblem::%s.\n", __func__ )

    return "hohlraum";
}


//------------------------------------------------------------------------------------------------------------
//! \brief  Returns the string descriptor for an object's type.
//------------------------------------------------------------------------------------------------------------
const std::string HohlraumProblem::GetDescriptor( void ) const {

    PRINT_STATUS( "Executing HohlraumProblem::%s.\n", __func__ )

    return Descriptor();
}


//------------------------------------------------------------------------------------------------------------
//! \brief  Prints a summary of the problem configuration to the logging interface.
//------------------------------------------------------------------------------------------------------------
void HohlraumProblem::Print (

    const std::string & prefix // = "  "

) const {

    PRINT_STATUS( "Executing HohlraumProblem::%s.\n", __func__ )

    this->Problem::Print( prefix );
}


//------------------------------------------------------------------------------------------------------------
//! \brief  Sets inflow boundary conditions for the test problem.
//!
//! \attention  This routine assumes that the boundary cells of the given object have been zero-initialized.
//!
//! \param[in]  source  Object in which to assign boundary conditions.
//------------------------------------------------------------------------------------------------------------
void HohlraumProblem::SetInflowBoundaries (

    RKDG::DensityFunction & source

) const {

    PRINT_STATUS( "Executing HohlraumProblem::%s.\n", __func__ )

    if ( source.MPI_block_coords(0) == 0 ) {

        for ( int64_t j = 0; j <= source.nx(1) + 1; ++j )
            source(0,j,0,0) = 1.0;
    }
}


//------------------------------------------------------------------------------------------------------------
//! \brief  Overrides values in given ParameterList with requirements imposed by test problem.
//!
//! \param[in]  params  ParameterList to apply overrides to.
//------------------------------------------------------------------------------------------------------------
void HohlraumProblem::OverrideOptions (

    ParameterList & params

) const {

    PRINT_STATUS( "Executing HohlraumProblem::%s.\n", __func__ )

    // Override limits of spatial domain.
    try {
        if ( params.GetValue<double>( "ax" ) != 0.0 ) {

            PRINT_WARNING( "Applying override to parameter 'ax' in ParameterList passed to "
                           "HohlraumProblem::%s.\n", __func__ )

            params.SetValue( "ax", "0.0" );
        }
    } catch (...) {/* empty */}

    try {
        if ( params.GetValue<double>( "bx" ) != 1.3 ) {

            PRINT_WARNING( "Applying override to parameter 'bx' in ParameterList passed to "
                           "HohlraumProblem::%s.\n", __func__ )

            params.SetValue( "bx", "1.3" );
        }
    } catch (...) {/* empty */}

    try {
        if ( params.GetValue<double>( "ay" ) != 0.0 ) {

            PRINT_WARNING( "Applying override to parameter 'ay' in ParameterList passed to "
                           "HohlraumProblem::%s.\n", __func__ )

            params.SetValue( "ay", "0.0" );
        }
    } catch (...) {/* empty */}

    try {
        if ( params.GetValue<double>( "by" ) != 1.3 ) {

            PRINT_WARNING( "Applying override to parameter 'by' in ParameterList passed to "
                           "HohlraumProblem::%s.\n", __func__ )

            params.SetValue( "by", "1.3" );
        }
    } catch (...) {/* empty */}

    // Override periodic boundary conditions.
    try {
        if ( params.GetValue<bool>( "periodic" ) == true ) {

            PRINT_WARNING( "Applying override to parameter 'periodic' in ParameterList passed to "
                           "HohlraumProblem::%s. Periodic boundary conditions cannot be used for "
                           "2D hohlraum test problems.\n", __func__ )

            params.SetValue( "periodic", "false" );
        }
    } catch (...) {/* empty */}

    // Override reflecting condition where inflow is defined.
    try {
        if ( params.GetValue<bool>( "reflect_x_min" ) == true ) {

            PRINT_WARNING( "Applying override to parameter 'reflect_x_min' in ParameterList passed to "
                           "HohlraumProblem::%s. Inflow condition is defined on X_Min edge.\n", __func__ )

            params.SetValue( "reflect_x_min", "false" );
        }
    } catch (...) {/* empty */}
}


//============================================================================================================
//=== PROTECTED INTERFACE ROUTINES ===========================================================================
//============================================================================================================


//------------------------------------------------------------------------------------------------------------
//! \brief  Returns the value of the total cross section without mollification at the specified position.
//!
//! \param[in]  x   \f$ x \f$ coordinate of position.
//! \param[in]  y   \f$ y \f$ coordinate of position.
//------------------------------------------------------------------------------------------------------------
double HohlraumProblem::TotalCross (

    const double x,
    const double y

) const {

    // Absorbing cells.
    if (    x > 1.25
         || y <= 0.05
         || y > 1.25
    ) {
        return 100.0;
    }

    if (    ( y > .25 && y <= 1.05 )
         && ( x <= .05 || ( x > .45 && x <= .85 ) )
    ) {
        return 100.0;
    }

    // Void cells.
    return 0.1;
}


//------------------------------------------------------------------------------------------------------------
//! \brief  Returns the value of the scattering cross section without mollification at the specified position.
//!
//! \param[in]  x   \f$ x \f$ coordinate of position.
//! \param[in]  y   \f$ y \f$ coordinate of position.
//------------------------------------------------------------------------------------------------------------
double HohlraumProblem::ScatterCross (

    const double x,
    const double y

) const {

    // Absorbing cells.
    if (    x >  1.25
         || y <= 0.05
         || y >  1.25
    ) {
        return 1.0;
    }

    if (    ( y >  .25 && y <= 1.05 )
         && ( x <= .05 || ( x > .45 && x <= .85 ) )
    ) {
        return 1.0;
    }

    // Void cells.
    return 0.1;
}


//============================================================================================================
//=== EXPLICIT INSTANTIATIONS ================================================================================
//============================================================================================================

# include "utils/Factory/DerivedFactory.hpp"

template class DerivedFactory< Problem, HohlraumProblem >;


# endif // if SPACE_DIMS == 2
