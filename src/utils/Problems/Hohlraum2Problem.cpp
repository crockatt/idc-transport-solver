//------------------------------------------------------------------------------------------------------------
//! \file   utils/Problems/Hohlraum2Problem.cpp
//! \brief  Contains implementation of the Hohlraum2Problem class.
//!
//! \author Michael M. Crockatt
//! \date   June 2018
//------------------------------------------------------------------------------------------------------------

# include "utils/CLog.hpp"
# include "utils/global.hpp"
# include "utils/Problems/Hohlraum2Problem.hpp"


# if SPACE_DIMS == 2 || defined (DOXYCOMPILE)


//============================================================================================================
//=== CONSTRUCTORS AND DESTRUCTOR ============================================================================
//============================================================================================================


//------------------------------------------------------------------------------------------------------------
//! \brief  Construct object from given parameters.
//!
//! \param[in]  params  Contains parameters used to construct object.
//------------------------------------------------------------------------------------------------------------
Hohlraum2Problem::Hohlraum2Problem (

    const ParameterList & params
) :
    HohlraumProblem( params )
{
    PRINT_STATUS( "Executing Hohlraum2Problem::%s.\n", __func__ )

    /* empty */
}


//------------------------------------------------------------------------------------------------------------
//! \brief  Destructor for Hohlraum2Problem class.
//------------------------------------------------------------------------------------------------------------
Hohlraum2Problem::~Hohlraum2Problem ( void ) {

    PRINT_STATUS( "Executing Hohlraum2Problem::%s.\n", __func__ )

    /* empty */
}


//============================================================================================================
//=== PUBLIC INTERFACE ROUTINES ==============================================================================
//============================================================================================================


//------------------------------------------------------------------------------------------------------------
//! \brief  Returns the string descriptor of the class type.
//------------------------------------------------------------------------------------------------------------
const std::string Hohlraum2Problem::Descriptor( void ) {

    PRINT_STATUS( "Executing Hohlraum2Problem::%s.\n", __func__ )

    return "hohlraum2";
}


//------------------------------------------------------------------------------------------------------------
//! \brief  Returns the string descriptor for an object's type.
//------------------------------------------------------------------------------------------------------------
const std::string Hohlraum2Problem::GetDescriptor( void ) const {

    PRINT_STATUS( "Executing Hohlraum2Problem::%s.\n", __func__ )

    return Descriptor();
}


//------------------------------------------------------------------------------------------------------------
//! \brief  Prints a summary of the problem configuration to the logging interface.
//------------------------------------------------------------------------------------------------------------
void Hohlraum2Problem::Print (

    const std::string & prefix // = "  "

) const {

    PRINT_STATUS( "Executing Hohlraum2Problem::%s.\n", __func__ )

    this->HohlraumProblem::Print( prefix );
}


//============================================================================================================
//=== PROTECTED INTERFACE ROUTINES ===========================================================================
//============================================================================================================


//------------------------------------------------------------------------------------------------------------
//! \brief  Returns the value of the scattering cross section without mollification at the specified position.
//!
//! \param[in]  x   \f$ x \f$ coordinate of position.
//! \param[in]  y   \f$ y \f$ coordinate of position.
//------------------------------------------------------------------------------------------------------------
double Hohlraum2Problem::ScatterCross (

    const double x,
    const double y

) const {

    // Top, bottom, and right walls.
    if ( x > 1.25 || y > 1.25 || y < 0.05 )
        return 50.0;

    // Center block (absorp).
    if ( x > 0.50 && x <= 0.85 && y > 0.30 && y <= 1.00 )
        return 0.0;

    // Center block (left).
    if ( x > 0.45 && x <= 0.50 && y > 0.25 && y <= 1.05 )
        return 90.0;

    // Center block (top).
    if ( x > 0.50 && x <= 0.85 && y > 1.00 && y <= 1.05 )
        return 90.0;

    // Center block (bottom).
    if ( x > 0.50 && x <= 0.85 && y > 0.25 && y <= 0.30 )
        return 90.0;

    // Left wall.
    if ( x <= 0.05 && y > 0.25 && y <= 1.05 )
        return 95.0;

    // Void cells.
    return 0.1;
}


//============================================================================================================
//=== EXPLICIT INSTANTIATIONS ================================================================================
//============================================================================================================

# include "utils/Factory/DerivedFactory.hpp"

template class DerivedFactory< Problem, Hohlraum2Problem >;


# endif // if SPACE_DIMS == 2
