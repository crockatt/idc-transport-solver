//------------------------------------------------------------------------------------------------------------
//! \file   utils/Problems/ReedsProblem.cpp
//! \brief  Contains implementation of the ReedsProblem class.
//!
//! \author Michael M. Crockatt
//! \date   June 2018
//------------------------------------------------------------------------------------------------------------

# include "utils/CLog.hpp"
# include "utils/global.hpp"
# include "utils/Problems/ReedsProblem.hpp"


# if SPACE_DIMS == 1 || defined (DOXYCOMPILE)


//============================================================================================================
//=== CONSTRUCTORS AND DESTRUCTOR ============================================================================
//============================================================================================================


//------------------------------------------------------------------------------------------------------------
//! \brief  Construct object from given parameters.
//!
//! \param[in]  params  Contains parameters used to construct object.
//------------------------------------------------------------------------------------------------------------
ReedsProblem::ReedsProblem (

    const ParameterList & params
) :
    Problem( params )
{
    PRINT_STATUS( "Executing ReedsProblem::%s.\n", __func__ )

    params.GetValue( "alpha", this->alpha );
}


//------------------------------------------------------------------------------------------------------------
//! \brief  Destructor for ReedsProblem class.
//------------------------------------------------------------------------------------------------------------
ReedsProblem::~ReedsProblem ( void ) {

    PRINT_STATUS( "Executing ReedsProblem::%s.\n", __func__ )

    /* empty */
}


//============================================================================================================
//=== PUBLIC INTERFACE ROUTINES ==============================================================================
//============================================================================================================


//------------------------------------------------------------------------------------------------------------
//! \brief  Returns the string descriptor of the class type.
//------------------------------------------------------------------------------------------------------------
const std::string ReedsProblem::Descriptor( void ) {

    PRINT_STATUS( "Executing ReedsProblem::%s.\n", __func__ )

    return "reeds-problem";
}


//------------------------------------------------------------------------------------------------------------
//! \brief  Returns the string descriptor for an object's type.
//------------------------------------------------------------------------------------------------------------
const std::string ReedsProblem::GetDescriptor( void ) const {

    PRINT_STATUS( "Executing ReedsProblem::%s.\n", __func__ )

    return Descriptor();
}


//------------------------------------------------------------------------------------------------------------
//! \brief  Prints a summary of the problem configuration to the logging interface.
//------------------------------------------------------------------------------------------------------------
void ReedsProblem::Print (

    const std::string & prefix // = "  "

) const {

    PRINT_STATUS( "Executing ReedsProblem::%s.\n", __func__ )

    this->Problem::Print( prefix );

    PRINT_LOG( "%s%-*s  %.4e\n", prefix.c_str(), Global::col_width, "Alpha:", this->alpha )
}


//------------------------------------------------------------------------------------------------------------
//! \brief  Overrides values in given ParameterList with requirements imposed by test problem.
//!
//! \param[in]  params  ParameterList to apply overrides to.
//------------------------------------------------------------------------------------------------------------
void ReedsProblem::OverrideOptions (

    ParameterList & params

) const {

    PRINT_STATUS( "Executing ReedsProblem::%s.\n", __func__ )

    try {
        if ( params.GetValue<double>( "ax" ) != 0.0 ) {

            PRINT_WARNING( "Applying override to parameter 'ax' in ParameterList passed to "
                           "ReedsProblem::%s.\n", __func__ )

            params.SetValue( "ax", "0.0" );
        }
    } catch (...) {/* empty */}

    try {
        if ( params.GetValue<double>( "bx" ) != 16.0 ) {

            PRINT_WARNING( "Applying override to parameter 'bx' in ParameterList passed to "
                           "ReedsProblem::%s.\n", __func__ )

            params.SetValue( "bx", "16.0" );
        }
    } catch (...) {/* empty */}
}


//============================================================================================================
//=== PROTECTED INTERFACE ROUTINES ===========================================================================
//============================================================================================================


//------------------------------------------------------------------------------------------------------------
//! \brief  Returns the value of the source term without mollification at the specified position.
//!
//!  \param[in]  x   \f$ x \f$ coordinate of position.
//------------------------------------------------------------------------------------------------------------
double ReedsProblem::Source (

    const double x

) const {

    // Region 2.
    if (    ( x >  2.0 && x <  3.0 )
         || ( x < 14.0 && x > 13.0 )
    ) {
        return this->alpha;
    }

    // Region 5.
    if ( x > 6.0 && x < 10.0 )
        return 50.0 * this->alpha;

    return 0.0;
}


//------------------------------------------------------------------------------------------------------------
//! \brief  Returns the value of the total cross section without mollification at the specified position.
//!
//! \param[in]  x   \f$ x \f$ coordinate of position.
//------------------------------------------------------------------------------------------------------------
double ReedsProblem::TotalCross (

    const double x

) const {

    // Region 1.
    if (    ( x <  2.0 )
         || ( x > 14.0 )
    ) {
        return this->alpha;
    }

    // Region 2.
    if (    ( x >=  2.0 && x <  3.0 )
         || ( x <= 14.0 && x > 13.0 )
    ) {
        return 0.3 * this->alpha;
    }

    // Region 4.
    if (    ( x >=  5.0 && x <  6.0 )
         || ( x <= 11.0 && x > 10.0 )
    ) {
        return 5.0 * this->alpha;
    }

    // Region 5.
    if ( x >= 6.0 && x <= 10.0 ) {

        return 50.0 * this->alpha;
    }

    return 0.0;
}


//------------------------------------------------------------------------------------------------------------
//! \brief  Returns the value of the scattering cross section without mollification at the specified position.
//!
//! \param[in]  x   \f$ x \f$ coordinate of position.
//------------------------------------------------------------------------------------------------------------
double ReedsProblem::ScatterCross (

    const double x

) const {

    // Region 1.
    if (    ( x <  2.0 )
         || ( x > 14.0 )
    ) {
        return 0.9 * this->alpha;
    }

    // Region 2.
    if (    ( x >=  2.0 && x <  3.0 )
         || ( x <= 14.0 && x > 13.0 )
    ) {
        return 0.1 * this->alpha;
    }

    return 0.0;
}


//============================================================================================================
//=== EXPLICIT INSTANTIATIONS ================================================================================
//============================================================================================================

# include "utils/Factory/DerivedFactory.hpp"

template class DerivedFactory< Problem, ReedsProblem >;


# endif // if SPACE_DIMS == 1
