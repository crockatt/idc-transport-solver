//------------------------------------------------------------------------------------------------------------
//! \file   utils/Problems/LatticeProblem.cpp
//! \brief  Contains implementation of the LatticeProblem class.
//!
//! \author Michael M. Crockatt
//! \date   June 2018
//------------------------------------------------------------------------------------------------------------

# include "utils/CLog.hpp"
# include "utils/global.hpp"
# include "utils/Problems/LatticeProblem.hpp"


# if SPACE_DIMS == 1 || SPACE_DIMS == 2 || defined (DOXYCOMPILE)


//============================================================================================================
//=== CONSTRUCTORS AND DESTRUCTOR ============================================================================
//============================================================================================================


//------------------------------------------------------------------------------------------------------------
//! \brief  Construct object from given parameters.
//!
//! \param[in]  params  Contains parameters used to construct object.
//------------------------------------------------------------------------------------------------------------
LatticeProblem::LatticeProblem (

    const ParameterList & params
) :
    Problem( params )
{
    PRINT_STATUS( "Executing LatticeProblem::%s.\n", __func__ )

# if SPACE_DIMS == 1

    try         {  params.GetValue( "epsilon", this->epsilon ); }
    catch (...) {  this->epsilon = 1.0;                         }

# endif
}


//------------------------------------------------------------------------------------------------------------
//! \brief  Destructor for LatticeProblem class.
//------------------------------------------------------------------------------------------------------------
LatticeProblem::~LatticeProblem ( void ) {

    PRINT_STATUS( "Executing LatticeProblem::%s.\n", __func__ )

    /* empty */
}


//============================================================================================================
//=== PUBLIC INTERFACE ROUTINES ==============================================================================
//============================================================================================================


//------------------------------------------------------------------------------------------------------------
//! \brief  Returns the string descriptor of the class type.
//------------------------------------------------------------------------------------------------------------
const std::string LatticeProblem::Descriptor( void ) {

    PRINT_STATUS( "Executing LatticeProblem::%s.\n", __func__ )

    return "lattice";
}


//------------------------------------------------------------------------------------------------------------
//! \brief  Returns the string descriptor for an object's type.
//------------------------------------------------------------------------------------------------------------
const std::string LatticeProblem::GetDescriptor( void ) const {

    PRINT_STATUS( "Executing LatticeProblem::%s.\n", __func__ )

    return Descriptor();
}


//------------------------------------------------------------------------------------------------------------
//! \brief  Prints a summary of the problem configuration to the logging interface.
//------------------------------------------------------------------------------------------------------------
void LatticeProblem::Print (

    const std::string & prefix // = "  "

) const {

    PRINT_STATUS( "Executing LatticeProblem::%s.\n", __func__ )

    this->Problem::Print( prefix );

# if SPACE_DIMS == 1

    PRINT_LOG( "%s%-*s  %.4e\n", prefix.c_str(), Global::col_width, "Epsilon:", this->epsilon )

# endif
}


//------------------------------------------------------------------------------------------------------------
//! \brief  Overrides values in given ParameterList with requirements imposed by test problem.
//!
//! \param[in]  params  ParameterList to apply overrides to.
//------------------------------------------------------------------------------------------------------------
void LatticeProblem::OverrideOptions (

    ParameterList & params

) const {

    PRINT_STATUS( "Executing LatticeProblem::%s.\n", __func__ )

    try {
        if ( params.GetValue<double>( "ax" ) != 0.0 ) {

            PRINT_WARNING( "Applying override to parameter 'ax' in ParameterList passed to "
                           "LatticeProblem::%s.\n", __func__ )

            params.SetValue( "ax", "0.0" );
        }
    } catch (...) {/* empty */}

    try {
        if ( params.GetValue<double>( "bx" ) != 7.0 ) {

            PRINT_WARNING( "Applying override to parameter 'bx' in ParameterList passed to "
                           "LatticeProblem::%s.\n", __func__ )

            params.SetValue( "bx", "7.0" );
        }
    } catch (...) {/* empty */}

# if SPACE_DIMS == 2

    try {
        if ( params.GetValue<double>( "ay" ) != 0.0 ) {

            PRINT_WARNING( "Applying override to parameter 'ay' in ParameterList passed to "
                           "LatticeProblem::%s.\n", __func__ )

            params.SetValue( "ay", "0.0" );
        }
    } catch (...) {/* empty */}

    try {
        if ( params.GetValue<double>( "by" ) != 7.0 ) {

            PRINT_WARNING( "Applying override to parameter 'by' in ParameterList passed to "
                           "LatticeProblem::%s.\n", __func__ )

            params.SetValue( "by", "7.0" );
        }
    } catch (...) {/* empty */}

# endif // if SPACE_DIMS == 2
}


//============================================================================================================
//=== PROTECTED INTERFACE ROUTINES ===========================================================================
//============================================================================================================


//------------------------------------------------------------------------------------------------------------
//! \brief  Returns the value of the source term without mollification at the specified position.
//!
//! \param[in]  x   \f$ x \f$ coordinate of position.
//! \param[in]  y   \f$ y \f$ coordinate of position.
//------------------------------------------------------------------------------------------------------------
double LatticeProblem::Source (

    const double x
# if SPACE_DIMS == 2 || defined (DOXYCOMPILE)
  , const double y
# endif

) const {

    if (    x > 3.0 && x < 4.0
    # if SPACE_DIMS == 2
         && y > 3.0 && y < 4.0
    # endif
    ) {

    # if SPACE_DIMS == 1

        return this->epsilon;

    # elif SPACE_DIMS == 2

        return 1.0;

    # endif
    }

    return 0.0;
}


//------------------------------------------------------------------------------------------------------------
//! \brief  Returns the value of the total cross section without mollification at the specified position.
//!
//! \param[in]  x   \f$ x \f$ coordinate of position.
//! \param[in]  y   \f$ y \f$ coordinate of position.
//------------------------------------------------------------------------------------------------------------
double LatticeProblem::TotalCross (

    const double x
# if SPACE_DIMS == 2 || defined (DOXYCOMPILE)
  , const double y
# endif

) const {

# if SPACE_DIMS == 1

    if (    ( x > 1.0 && x < 2.0 )
         || ( x > 5.0 && x < 6.0 )
    ) {
        return 1.0 / this->epsilon;
    }

    return 0.0;

# elif SPACE_DIMS == 2

    // Absorber cells.
    if ( ( x > 1.0 && x <= 2.0 ) && (    ( y > 1.0 && y <= 2.0 )
                                      || ( y > 3.0 && y <= 4.0 )
                                      || ( y > 5.0 && y <= 6.0 ) )
    ) {
        return 10.0;

    } else if ( ( x > 2.0 && x <= 3.0 ) && (    ( y > 2.0 && y <= 3.0 )
                                             || ( y > 4.0 && y <= 5.0 ) )
    ) {
        return 10.0;

    } else if ( ( x > 3.0 && x <= 4.0 ) && ( y > 5.0 && y <= 6.0 ) ) {

        return 10.0;

    } else if ( ( x > 4.0 && x <= 5.0 ) && (    ( y > 2.0 && y <= 3.0 )
                                             || ( y > 4.0 && y <= 5.0 ) )
    ) {
        return 10.0;

    } else if ( ( x > 5.0 && x <= 6.0 ) && (    ( y > 1.0 && y <= 2.0 )
                                             || ( y > 3.0 && y <= 4.0 )
                                             || ( y > 5.0 && y <= 6.0 ) )
    ) {
        return 10.0;
    }

    // Scattering background cells.
    return 1.0;

# endif // if SPACE_DIMS == ?
}


//------------------------------------------------------------------------------------------------------------
//! \brief  Returns the value of the scattering cross section without mollification at the specified position.
//!
//! \param[in]  x   \f$ x \f$ coordinate of position.
//! \param[in]  y   \f$ y \f$ coordinate of position.
//------------------------------------------------------------------------------------------------------------
double LatticeProblem::ScatterCross (

    const double x
# if SPACE_DIMS == 2 || defined (DOXYCOMPILE)
  , const double y
# endif

) const {

# if SPACE_DIMS == 1

    if (    ( x >= 1.0 && x <= 2.0 )
         || ( x >= 5.0 && x <= 6.0 )
    ) {
        return 1.0 / this->epsilon;
    }

    return 0.0;

# elif SPACE_DIMS == 2

    // Absorber cells.
    if ( ( x > 1.0 && x <= 2.0 ) && (    ( y > 1.0 && y <= 2.0 )
                                      || ( y > 3.0 && y <= 4.0 )
                                      || ( y > 5.0 && y <= 6.0 ) )
    ) {
        return 0.0;

    } else if ( ( x > 2.0 && x <= 3.0 ) && (    ( y > 2.0 && y <= 3.0 )
                                             || ( y > 4.0 && y <= 5.0 ) )
    ) {
        return 0.0;

    } else if ( ( x > 3.0 && x <= 4.0 ) && ( y > 5.0 && y <= 6.0 ) ) {

        return 0.0;

    } else if ( ( x > 4.0 && x <= 5.0 ) && (    ( y > 2.0 && y <= 3.0 )
                                             || ( y > 4.0 && y <= 5.0 ) )
    ) {
        return 0.0;

    } else if ( ( x > 5.0 && x <= 6.0 ) && (    ( y > 1.0 && y <= 2.0 )
                                             || ( y > 3.0 && y <= 4.0 )
                                             || ( y > 5.0 && y <= 6.0 ) )
    ) {
        return 0.0;
    }

    // Scattering background cells.
    return 1.0;

# endif // if SPACE_DIMS == ?
}


//============================================================================================================
//=== EXPLICIT INSTANTIATIONS ================================================================================
//============================================================================================================

# include "utils/Factory/DerivedFactory.hpp"

template class DerivedFactory< Problem, LatticeProblem >;


# endif // if SPACE_DIMS == 1 || SPACE_DIMS == 2
