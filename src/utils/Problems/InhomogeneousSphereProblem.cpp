//------------------------------------------------------------------------------------------------------------
//! \file   utils/Problems/InhomogeneousSphereProblem.cpp
//! \brief  Contains implementation of the InhomogeneousSphereProblem class.
//!
//! \author Michael M. Crockatt
//! \date   June 2018
//------------------------------------------------------------------------------------------------------------

# include "utils/CLog.hpp"
# include "utils/global.hpp"
# include "utils/Problems/InhomogeneousSphereProblem.hpp"


//============================================================================================================
//=== CONSTRUCTORS AND DESTRUCTOR ============================================================================
//============================================================================================================


//------------------------------------------------------------------------------------------------------------
//! \brief  Construct object from given parameters.
//!
//! \param[in]  params  Contains parameters used to construct object.
//------------------------------------------------------------------------------------------------------------
InhomogeneousSphereProblem::InhomogeneousSphereProblem (

    const ParameterList & params
) :
    Problem( params )
{
    PRINT_STATUS( "Executing InhomogeneousSphereProblem::%s.\n", __func__ )

    params.GetValue( "sphere_inner_radius",         this->inner_radius          );
    params.GetValue( "sphere_outer_radius",         this->outer_radius          );
    params.GetValue( "sphere_inner_trans_radius",   this->inner_trans_radius    );
    params.GetValue( "sphere_outer_trans_radius",   this->outer_trans_radius    );
    params.GetValue( "sphere_inner_sigma_val",      this->inner_sigma_val       );
    params.GetValue( "sphere_interm_sigma_val",     this->interm_sigma_val      );
    params.GetValue( "sphere_outer_sigma_val",      this->outer_sigma_val       );
}


//------------------------------------------------------------------------------------------------------------
//! \brief  Destructor for InhomogeneousSphereProblem class.
//------------------------------------------------------------------------------------------------------------
InhomogeneousSphereProblem::~InhomogeneousSphereProblem ( void ) {

    PRINT_STATUS( "Executing InhomogeneousSphereProblem::%s.\n", __func__ )

    /* empty */
}


//============================================================================================================
//=== PUBLIC INTERFACE ROUTINES ==============================================================================
//============================================================================================================


//------------------------------------------------------------------------------------------------------------
//! \brief  Returns the string descriptor of the class type.
//------------------------------------------------------------------------------------------------------------
const std::string InhomogeneousSphereProblem::Descriptor( void ) {

    PRINT_STATUS( "Executing InhomogeneousSphereProblem::%s.\n", __func__ )

    return "inhomogeneous-sphere";
}


//------------------------------------------------------------------------------------------------------------
//! \brief  Returns the string descriptor for an object's type.
//------------------------------------------------------------------------------------------------------------
const std::string InhomogeneousSphereProblem::GetDescriptor( void ) const {

    PRINT_STATUS( "Executing InhomogeneousSphereProblem::%s.\n", __func__ )

    return Descriptor();
}


//------------------------------------------------------------------------------------------------------------
//! \brief  Prints a summary of the problem configuration to the logging interface.
//------------------------------------------------------------------------------------------------------------
void InhomogeneousSphereProblem::Print (

    const std::string & prefix // = "  "

) const {

    PRINT_STATUS( "Executing InhomogeneousSphereProblem::%s.\n", __func__ )

    this->Problem::Print( prefix );

    PRINT_LOG( "%s%-*s  %.4e\n", prefix.c_str(), Global::col_width,
               "Inner Radius:", this->inner_radius )

    PRINT_LOG( "%s%-*s  %.4e\n", prefix.c_str(), Global::col_width,
               "Outer Radius:", this->outer_radius )

    PRINT_LOG( "%s%-*s  %.4e\n", prefix.c_str(), Global::col_width,
               "Inner Trans. Radius:", this->inner_trans_radius )

    PRINT_LOG( "%s%-*s  %.4e\n", prefix.c_str(), Global::col_width,
               "Outer Trans. Radius:", this->outer_trans_radius )

    PRINT_LOG( "%s%-*s  %.4e\n", prefix.c_str(), Global::col_width,
               "Inner Sigma:", this->inner_sigma_val )

    PRINT_LOG( "%s%-*s  %.4e\n", prefix.c_str(), Global::col_width,
               "Interm. Sigma:", this->interm_sigma_val )

    PRINT_LOG( "%s%-*s  %.4e\n", prefix.c_str(), Global::col_width,
               "Outer Sigma:", this->outer_sigma_val )
}


//============================================================================================================
//=== PROTECTED INTERFACE ROUTINES ===========================================================================
//============================================================================================================


//------------------------------------------------------------------------------------------------------------
//! \brief  Returns the value of the source term without mollification at the specified position.
//!
//!  \param[in]  x   \f$ x \f$ coordinate of position.
//!  \param[in]  y   \f$ y \f$ coordinate of position.
//   \param[in]  z   \f$ z \f$ coordinate of position.
//------------------------------------------------------------------------------------------------------------
double InhomogeneousSphereProblem::Source (

    const double x
# if SPACE_DIMS >= 2 || defined (DOXYCOMPILE)
  , const double y
# endif
# if SPACE_DIMS == 3
  , const double z
# endif

) const {

    const double r =
        # if SPACE_DIMS == 1
            std::abs(x)
        # elif SPACE_DIMS == 2
            std::sqrt( x*x + y*y )
        # elif SPACE_DIMS == 3
            std::sqrt( x*x + y*y + z*z )
        # endif
            ;

    // Inner region.
    if ( r <= this->inner_radius - 0.5 * this->inner_trans_radius ) {

        return 1.0;

    // Transition region.
    } else if ( r < this->inner_radius + 0.5 * this->inner_trans_radius ) {

        const double rho = 0.5 * this->inner_trans_radius;
        const double R = this->inner_radius;

        const double exp_val = std::exp( 2*rho / ((R + rho) - r) - 2*rho / (r - (R - rho)) );

        return 1.0 / (1.0 + exp_val);

    // Outer region.
    } else {

        return 0.0;
    }
}


//------------------------------------------------------------------------------------------------------------
//! \brief  Returns the value of the total cross section without mollification at the specified position.
//!
//! \param[in]  x   \f$ x \f$ coordinate of position.
//! \param[in]  y   \f$ y \f$ coordinate of position.
//  \param[in]  z   \f$ z \f$ coordinate of position.
//------------------------------------------------------------------------------------------------------------
double InhomogeneousSphereProblem::TotalCross (

    const double x
# if SPACE_DIMS >= 2 || defined (DOXYCOMPILE)
  , const double y
# endif
# if SPACE_DIMS == 3
  , const double z
# endif

) const {

    const double r =
        # if SPACE_DIMS == 1
            std::abs(x)
        # elif SPACE_DIMS == 2
            std::sqrt( x*x + y*y )
        # elif SPACE_DIMS == 3
            std::sqrt( x*x + y*y + z*z )
        # endif
            ;

    // Inner region.
    if ( r <= this->inner_radius - 0.5 * this->inner_trans_radius ) {

        return this->inner_sigma_val;

    // Inner transition region.
    } else if ( r < this->inner_radius + 0.5 * this->inner_trans_radius ) {

        const double rho = 0.5 * this->inner_trans_radius;
        const double R = this->inner_radius;

        const double exp_val = std::exp( 2*rho / ((R + rho) - r) - 2*rho / (r - (R - rho)) );

        return this->interm_sigma_val + (this->inner_sigma_val - this->interm_sigma_val) / (1.0 + exp_val);

    // Intermediate region.
    } else if ( r <= this->outer_radius - 0.5 * this->outer_trans_radius ) {

        return this->interm_sigma_val;

    // Outer transition region.
    } else if ( r < this->outer_radius + 0.5 * this->outer_trans_radius ) {

        const double rho = 0.5 * this->outer_trans_radius;
        const double R = this->outer_radius;

        const double exp_val = std::exp( 2*rho / ((R + rho) - r) - 2*rho / (r - (R - rho)) );

        return this->outer_sigma_val + (this->interm_sigma_val - this->outer_sigma_val) / (1.0 + exp_val);

    // Outer region.
    } else {

        return this->outer_sigma_val;
    }
}


//------------------------------------------------------------------------------------------------------------
//! \brief  Returns the value of the scattering cross section without mollification at the specified position.
//!
//! \param[in]  x   \f$ x \f$ coordinate of position.
//! \param[in]  y   \f$ y \f$ coordinate of position.
//  \param[in]  z   \f$ z \f$ coordinate of position.
//------------------------------------------------------------------------------------------------------------
double InhomogeneousSphereProblem::ScatterCross (

    const double x
# if SPACE_DIMS >= 2 || defined (DOXYCOMPILE)
  , const double y
# endif
# if SPACE_DIMS == 3
  , const double z
# endif

) const {

    return this->TotalCross
        # if SPACE_DIMS == 1
            ( x )
        # elif SPACE_DIMS == 2
            ( x, y )
        # elif SPACE_DIMS == 3
            ( x, y, z )
        # endif
            ;
}


//============================================================================================================
//=== EXPLICIT INSTANTIATIONS ================================================================================
//============================================================================================================

# include "utils/Factory/DerivedFactory.hpp"

template class DerivedFactory< Problem, InhomogeneousSphereProblem >;
