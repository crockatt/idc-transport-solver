//------------------------------------------------------------------------------------------------------------
//! \file   utils/Problems/BumpProblem.cpp
//! \brief  Contains implementation of the BumpProblem class.
//!
//! \author Michael M. Crockatt
//! \date   June 2018
//------------------------------------------------------------------------------------------------------------

# include "utils/CLog.hpp"
# include "utils/global.hpp"
# include "utils/Problems/BumpProblem.hpp"


//============================================================================================================
//=== CONSTRUCTORS AND DESTRUCTOR ============================================================================
//============================================================================================================


//------------------------------------------------------------------------------------------------------------
//! \brief  Construct object from given parameters.
//!
//! \param[in]  params  Contains parameters used to construct object.
//------------------------------------------------------------------------------------------------------------
BumpProblem::BumpProblem (

    const ParameterList & params
) :
    Problem( params )
{
    PRINT_STATUS( "Executing BumpProblem::%s.\n", __func__ )

    params.GetValue( "sigma_val",   this->sigma_val   );
    params.GetValue( "bump_radius", this->bump_radius );

    this->bump_moll = std::shared_ptr<Mollifier>( new Mollifier( this->bump_radius ) );
}


//------------------------------------------------------------------------------------------------------------
//! \brief  Destructor for BumpProblem class.
//------------------------------------------------------------------------------------------------------------
BumpProblem::~BumpProblem ( void ) {

    PRINT_STATUS( "Executing BumpProblem::%s.\n", __func__ )

    /* empty */
}


//============================================================================================================
//=== PUBLIC INTERFACE ROUTINES ==============================================================================
//============================================================================================================


//------------------------------------------------------------------------------------------------------------
//! \brief  Returns the string descriptor of the class type.
//------------------------------------------------------------------------------------------------------------
const std::string BumpProblem::Descriptor( void ) {

    PRINT_STATUS( "Executing BumpProblem::%s.\n", __func__ )

    return "bump";
}


//------------------------------------------------------------------------------------------------------------
//! \brief  Returns the string descriptor for an object's type.
//------------------------------------------------------------------------------------------------------------
const std::string BumpProblem::GetDescriptor( void ) const {

    PRINT_STATUS( "Executing BumpProblem::%s.\n", __func__ )

    return Descriptor();
}


//------------------------------------------------------------------------------------------------------------
//! \brief  Prints a summary of the problem configuration to the logging interface.
//------------------------------------------------------------------------------------------------------------
void BumpProblem::Print (

    const std::string & prefix // = "  "

) const {

    PRINT_STATUS( "Executing BumpProblem::%s.\n", __func__ )

    this->Problem::Print( prefix );

    PRINT_LOG( "%s%-*s  %.4e\n", prefix.c_str(), Global::col_width, "Bump Radius:", this->bump_radius )
    PRINT_LOG( "%s%-*s  %.4e\n", prefix.c_str(), Global::col_width, "Sigma Val:",   this->sigma_val   )
}


//============================================================================================================
//=== PROTECTED INTERFACE ROUTINES ===========================================================================
//============================================================================================================


//------------------------------------------------------------------------------------------------------------
//! \brief  Returns the value of the initial condition without mollification at the specified position.
//!
//!  \param[in]  x   \f$ x \f$ coordinate of position.
//!  \param[in]  y   \f$ y \f$ coordinate of position.
//   \param[in]  z   \f$ z \f$ coordinate of position.
//------------------------------------------------------------------------------------------------------------
double BumpProblem::InitialCondition (

        const double x
# if SPACE_DIMS >= 2 || defined (DOXYCOMPILE)
    ,   const double y
# endif
# if SPACE_DIMS == 3
    ,   const double z
# endif

) const {

    return this->bump_moll->Exact
                            # if SPACE_DIMS == 1
                                 ( x )
                            # elif SPACE_DIMS == 2
                                 ( x, y )
                            # elif SPACE_DIMS == 3
                                 ( x, y, z )
                            # endif
                                 ;
}


//------------------------------------------------------------------------------------------------------------
//! \brief  Returns the value of the total cross section without mollification at the specified position.
//!
//  \param[in]  x   \f$ x \f$ coordinate of position.
//  \param[in]  y   \f$ y \f$ coordinate of position.
//  \param[in]  z   \f$ z \f$ coordinate of position.
//------------------------------------------------------------------------------------------------------------
double BumpProblem::TotalCross (

        const double // x
# if SPACE_DIMS >= 2 || defined (DOXYCOMPILE)
    ,   const double // y
# endif
# if SPACE_DIMS == 3
    ,   const double // z
# endif

) const {

    return this->sigma_val;
}


//------------------------------------------------------------------------------------------------------------
//! \brief  Returns the value of the scattering cross section without mollification at the specified position.
//!
//  \param[in]  x   \f$ x \f$ coordinate of position.
//  \param[in]  y   \f$ y \f$ coordinate of position.
//  \param[in]  z   \f$ z \f$ coordinate of position.
//------------------------------------------------------------------------------------------------------------
double BumpProblem::ScatterCross (

        const double // x
# if SPACE_DIMS >= 2 || defined (DOXYCOMPILE)
    ,   const double // y
# endif
# if SPACE_DIMS == 3
    ,   const double // z
# endif

) const {

    return this->sigma_val;
}


//============================================================================================================
//=== EXPLICIT INSTANTIATIONS ================================================================================
//============================================================================================================

# include "utils/Factory/DerivedFactory.hpp"

template class DerivedFactory< Problem, BumpProblem >;
