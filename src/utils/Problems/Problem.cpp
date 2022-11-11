//------------------------------------------------------------------------------------------------------------
//! \file   utils/Problems/Problem.cpp
//! \brief  Contains implementation of the Problem class.
//!
//! \author Michael M. Crockatt
//! \date   June 2018
//------------------------------------------------------------------------------------------------------------

# include "utils/CLog.hpp"
# include "utils/global.hpp"
# include "utils/Problems/Problem.hpp"
# include "utils/Quadrule/Quadrule.hpp"


using namespace Quadrule;


//============================================================================================================
//=== CONSTRUCTORS AND DESTRUCTOR ============================================================================
//============================================================================================================


//------------------------------------------------------------------------------------------------------------
//! \brief  Construct object from given parameters.
//!
//! \param[in]  params  Contains parameters used to construct object.
//------------------------------------------------------------------------------------------------------------
Problem::Problem (

    const ParameterList & params
) :
    GK_order{ 15 },
    GK_nodes{ nullptr },
    G_weights{ nullptr },
    K_weights{ nullptr }
{
    PRINT_STATUS( "Executing Problem::%s.\n", __func__ )

    params.GetValue( "moll_radius", this->moll_radius );

    this->moll = std::shared_ptr<Mollifier>( new Mollifier( this->moll_radius ) );

    this->GK_nodes  = new double[ this->GK_order ];
    this->G_weights = new double[ this->GK_order/2 ];
    this->K_weights = new double[ this->GK_order ];

    ComputeGKquad( this->GK_order, this->GK_nodes, this->G_weights, this->K_weights );
}


//------------------------------------------------------------------------------------------------------------
//! \brief  Destructor for Problem class.
//------------------------------------------------------------------------------------------------------------
Problem::~Problem ( void ) {

    PRINT_STATUS( "Executing Problem::%s.\n", __func__ )

    delete [] this->GK_nodes;
    delete [] this->G_weights;
    delete [] this->K_weights;
}


//============================================================================================================
//=== PUBLIC INTERFACE ROUTINES ==============================================================================
//============================================================================================================


//------------------------------------------------------------------------------------------------------------
//! \brief  Returns the key used to search input lists for determining the derived type to construct.
//------------------------------------------------------------------------------------------------------------
std::string Problem::GetInputKey( void ) {

    PRINT_STATUS( "Executing Problem::%s.\n", __func__ )

    return "problem";
}


//------------------------------------------------------------------------------------------------------------
//! \brief  Prints a summary of the problem configuration to the logging interface.
//------------------------------------------------------------------------------------------------------------
void Problem::Print (

    const std::string & prefix // = "  "

) const {

    PRINT_STATUS( "Executing Problem::%s.\n", __func__ )

    PRINT_LOG( "%s%-*s  %s\n", prefix.c_str(), Global::col_width, "Initial Condition:",
               this->GetDescriptor().c_str() )

    PRINT_LOG( "%s%-*s  %.4e\n", prefix.c_str(), Global::col_width, "Mollifier Radius:", this->moll_radius )
}


//------------------------------------------------------------------------------------------------------------
//! \brief  Sets inflow boundary conditions for the test problem.
//!
//! \attention  This routine assumes that the boundary cells of the given object have been zero-initialized.
//!
//  \param[in]  source  Object in which to assign boundary conditions.
//------------------------------------------------------------------------------------------------------------
void Problem::SetInflowBoundaries (

    RKDG::DensityFunction & // source

) const {

    PRINT_STATUS( "Executing Problem::%s.\n", __func__ )

    /* empty */
}


//------------------------------------------------------------------------------------------------------------
//! \brief  Returns the value of the initial condition at the specified position.
//!
//! \param[in]  x   \f$ x \f$ coordinate of position.
//! \param[in]  y   \f$ y \f$ coordinate of position.
//  \param[in]  z   \f$ z \f$ coordinate of position.
//------------------------------------------------------------------------------------------------------------
double Problem::EvalInitialCondition (

    const double x
# if SPACE_DIMS >= 2 || defined (DOXYCOMPILE)
  , const double y
# endif
# if SPACE_DIMS == 3
  , const double z
# endif

) const {

    PRINT_STATUS( "Executing Problem::%s.\n", __func__ )

    return this->EvalMoll( &Problem::InitialCondition, x,
                                                    # if SPACE_DIMS >= 2
                                                       y,
                                                    # endif
                                                    # if SPACE_DIMS == 3
                                                       z,
                                                    # endif
                                                       this->InitialConditionFtSize() );
}


//------------------------------------------------------------------------------------------------------------
//! \brief  Returns the value of the source term at the specified position.
//!
//! \param[in]  x   \f$ x \f$ coordinate of position.
//! \param[in]  y   \f$ y \f$ coordinate of position.
//  \param[in]  z   \f$ z \f$ coordinate of position.
//------------------------------------------------------------------------------------------------------------
double Problem::EvalSource (

    const double x
# if SPACE_DIMS >= 2 || defined (DOXYCOMPILE)
  , const double y
# endif
# if SPACE_DIMS == 3
  , const double z
# endif

) const {

    PRINT_STATUS( "Executing Problem::%s.\n", __func__ )

    return this->EvalMoll( &Problem::Source, x,
                                        # if SPACE_DIMS >= 2
                                             y,
                                        # endif
                                        # if SPACE_DIMS == 3
                                             z,
                                        # endif
                                             this->SourceFtSize() );
}


//------------------------------------------------------------------------------------------------------------
//! \brief  Returns the value of the total cross section at the specified position.
//!
//! \param[in]  x   \f$ x \f$ coordinate of position.
//! \param[in]  y   \f$ y \f$ coordinate of position.
//  \param[in]  z   \f$ z \f$ coordinate of position.
//------------------------------------------------------------------------------------------------------------
double Problem::EvalTotalCross (

    const double x
# if SPACE_DIMS >= 2 || defined (DOXYCOMPILE)
  , const double y
# endif
# if SPACE_DIMS == 3
  , const double z
# endif

) const {

    PRINT_STATUS( "Executing Problem::%s.\n", __func__ )

    return this->EvalMoll( &Problem::TotalCross, x,
                                            # if SPACE_DIMS >= 2
                                                 y,
                                            # endif
                                            # if SPACE_DIMS == 3
                                                 z,
                                            # endif
                                                 this->TotalCrossFtSize() );
}


//------------------------------------------------------------------------------------------------------------
//! \brief  Returns the value of the scattering cross section at the specified position.
//!
//! \param[in]  x   \f$ x \f$ coordinate of position.
//! \param[in]  y   \f$ y \f$ coordinate of position.
//  \param[in]  z   \f$ z \f$ coordinate of position.
//------------------------------------------------------------------------------------------------------------
double Problem::EvalScatterCross (

    const double x
# if SPACE_DIMS >= 2 || defined (DOXYCOMPILE)
  , const double y
# endif
# if SPACE_DIMS == 3
  , const double z
# endif

) const {

    PRINT_STATUS( "Executing Problem::%s.\n", __func__ )

    return this->EvalMoll( &Problem::ScatterCross, x,
                                                # if SPACE_DIMS >= 2
                                                   y,
                                                # endif
                                                # if SPACE_DIMS == 3
                                                   z,
                                                # endif
                                                   this->ScatterCrossFtSize() );
}


//============================================================================================================
//=== PROTECTED INTERFACE ROUTINES ===========================================================================
//============================================================================================================


//------------------------------------------------------------------------------------------------------------
//! \brief  Returns the value of the initial condition at the specified position.
//!
//  \param[in]  x   \f$ x \f$ coordinate of position.
//  \param[in]  y   \f$ y \f$ coordinate of position.
//  \param[in]  z   \f$ z \f$ coordinate of position.
//------------------------------------------------------------------------------------------------------------
double Problem::InitialCondition (

    const double // x
# if SPACE_DIMS >= 2 || defined (DOXYCOMPILE)
  , const double // y
# endif
# if SPACE_DIMS == 3
  , const double // z
# endif

) const {

    PRINT_STATUS( "Executing Problem::%s.\n", __func__ )

    return 0.0;
}


//------------------------------------------------------------------------------------------------------------
//! \brief  Returns the value of the source term at the specified position.
//!
//  \param[in]  x   \f$ x \f$ coordinate of position.
//  \param[in]  y   \f$ y \f$ coordinate of position.
//  \param[in]  z   \f$ z \f$ coordinate of position.
//------------------------------------------------------------------------------------------------------------
double Problem::Source (

    const double // x
# if SPACE_DIMS >= 2 || defined (DOXYCOMPILE)
  , const double // y
# endif
# if SPACE_DIMS == 3
  , const double // z
# endif

) const {

    PRINT_STATUS( "Executing Problem::%s.\n", __func__ )

    return 0.0;
}


//------------------------------------------------------------------------------------------------------------
//! \brief  Returns the value of the total cross section at the specified position.
//!
//  \param[in]  x   \f$ x \f$ coordinate of position.
//  \param[in]  y   \f$ y \f$ coordinate of position.
//  \param[in]  z   \f$ z \f$ coordinate of position.
//------------------------------------------------------------------------------------------------------------
double Problem::TotalCross (

    const double // x
# if SPACE_DIMS >= 2 || defined (DOXYCOMPILE)
  , const double // y
# endif
# if SPACE_DIMS == 3
  , const double // z
# endif

) const {

    PRINT_STATUS( "Executing Problem::%s.\n", __func__ )

    return 0.0;
}


//------------------------------------------------------------------------------------------------------------
//! \brief  Returns the value of the scattering cross section at the specified position.
//!
//  \param[in]  x   \f$ x \f$ coordinate of position.
//  \param[in]  y   \f$ y \f$ coordinate of position.
//  \param[in]  z   \f$ z \f$ coordinate of position.
//------------------------------------------------------------------------------------------------------------
double Problem::ScatterCross (

    const double // x
# if SPACE_DIMS >= 2 || defined (DOXYCOMPILE)
  , const double // y
# endif
# if SPACE_DIMS == 3
  , const double // z
# endif

) const {

    PRINT_STATUS( "Executing Problem::%s.\n", __func__ )

    return 0.0;
}


//============================================================================================================
//=== PROTECTED INTEGRATION ROUTINES =========================================================================
//============================================================================================================


//--------------------------------------------------------------------------------------------------------
//! \brief  Main mollifier routine. Evaluates the mollification of \pp eval_fcn at the specified point.
//!
//! \param[in]  eval_fcn        Pointer to the function to apply mollifier to.
//! \param[in]  x_i             \f$ x \f$ coordinate of point to evaluate mollified function at.
//! \param[in]  y_j             \f$ y \f$ coordinate of point to evaluate mollified function at.
//! \param[in]  feature_size    Feature size of \pp eval_fcn.
//!
//! \return     Returns the value of the mollified parameter at the specified point.
//--------------------------------------------------------------------------------------------------------
double Problem::EvalMoll (

    const ProblemEvalFcn eval_fcn,
    const double x_i,
# if SPACE_DIMS >= 2 || defined (DOXYCOMPILE)
    const double y_j,
# endif
# if SPACE_DIMS == 3
    const double z_k,
# endif
    const double feature_size

) const {

    PRINT_STATUS( "Executing Problem::%s.\n", __func__ )

    const double TOL = 1.0e-14;

    // If the mollification radius is essentially zero, then just return the value of the function.
    if ( this->moll_radius < TOL ) {

        return (this->*eval_fcn)
                            # if SPACE_DIMS == 1
                                ( x_i )
                            # elif SPACE_DIMS == 2
                                ( x_i, y_j )
                            # elif SPACE_DIMS == 3
                                ( x_i, y_j, z_k )
                            # endif
                                ;
    }

    // If feature size is zero, just use standard adaptive integration routine.
    if ( feature_size == 0.0 ) {

        return IntegrateMollSmooth( eval_fcn,
                                    x_i,
                                # if SPACE_DIMS >= 2
                                    y_j,
                                # endif
                                # if SPACE_DIMS == 3
                                    z_k,
                                # endif
                                    x_i - this->moll_radius, x_i + this->moll_radius,
                                # if SPACE_DIMS >= 2
                                    y_j - this->moll_radius, y_j + this->moll_radius,
                                # endif
                                # if SPACE_DIMS == 3
                                    z_k - this->moll_radius, z_k + this->moll_radius,
                                # endif
                                    TOL );
    }

    // Break up convolution integral based on the feature size if necessary.
    if ( 2.0 * this->moll_radius >= feature_size ) {

        double result = 0.0;

        const int64_t num_intervals = 2.0 * this->moll_radius / feature_size + 1;
        const double dx = 2.0 * this->moll_radius / num_intervals;

    # if SPACE_DIMS == 1

        for ( int64_t i = 0; i < num_intervals; ++i ) {

            const double ax = x_i - this->moll_radius + dx * i;
            const double bx = ax + dx;

            result += this->IntegrateMollDiscontinuous( eval_fcn, x_i, ax, bx, TOL );
        }

    # elif SPACE_DIMS == 2

        const double dy = dx;

        for ( int64_t i = 0; i < num_intervals; ++i ) {
        for ( int64_t j = 0; j < num_intervals; ++j ) {

            const double ax = x_i - this->moll_radius + dx * i;
            const double bx = ax + dx;

            const double ay = y_j - this->moll_radius + dy * j;
            const double by = ay + dy;

            result += this->IntegrateMollDiscontinuous( eval_fcn, x_i, y_j, ax, bx, ay, by, TOL );
        }}

    # elif SPACE_DIMS == 3

        const double dy = dx;
        const double dz = dx;

        for ( int64_t i = 0; i < num_intervals; ++i ) {
        for ( int64_t j = 0; j < num_intervals; ++j ) {
        for ( int64_t k = 0; k < num_intervals; ++k ) {

            const double ax = x_i - this->moll_radius + dx * i;
            const double bx = ax + dx;

            const double ay = y_j - this->moll_radius + dy * j;
            const double by = ay + dy;

            const double az = z_k - this->moll_radius + dz * k;
            const double bz = az + dz;

            result += this->IntegrateMollDiscontinuous( eval_fcn, x_i, y_j, z_k, ax, bx, ay, by, az, bz, TOL );
        }}}

    # endif // if SPACE_DIMS == ?
    }

    return this->IntegrateMollDiscontinuous( eval_fcn,
                                             x_i,
                                        # if SPACE_DIMS >= 2
                                             y_j,
                                        # endif
                                        # if SPACE_DIMS == 3
                                             z_k,
                                        # endif
                                             x_i - this->moll_radius, x_i + this->moll_radius,
                                        # if SPACE_DIMS >= 2
                                             y_j - this->moll_radius, y_j + this->moll_radius,
                                        # endif
                                        # if SPACE_DIMS == 3
                                             z_k - this->moll_radius, z_k + this->moll_radius,
                                        # endif
                                             TOL );
}


//============================================================================================================
//=== ONE-DIMENSIONAL IMPLEMENTATIONS ========================================================================
//============================================================================================================

# if SPACE_DIMS == 1


//------------------------------------------------------------------------------------------------------------
//! \brief  Used in the computation of a numerical approximation of the convolution integral used to define
//!         the mollification of the given function at the specified point.
//!
//! This function assumes that the integrand has _at most one_ discontinuity within the domain of integration.
//! The purpose of this function is to find the (single) discontinuity (if one exists) and then call the
//! function IntegrateMollSmooth() in a way that integrates around the discontinuity.
//!
//! \param[in]      eval_fcn    Pointer to the function to integrate against. This is the function to which
//!                             the mollifier is applied.
//! \param[in]      x_i         Center point of the mollifier.
//! \param[in]      ax          Integrate over \f$ [a_x,b_x] \f$.
//! \param[in]      bx          Integrate over \f$ [a_x,b_x] \f$.
//! \param[in]      tol         Error tolerance.
//!
//! \return     Returns an approximation of the specified integral.
//------------------------------------------------------------------------------------------------------------
double Problem::IntegrateMollDiscontinuous (

    const ProblemEvalFcn eval_fcn,
    const double x_i,
    const double ax,
    const double bx,
    const double tol

) const {

    PRINT_STATUS( "Executing Problem::%s.\n", __func__ )

    // Locate discontinuity via bisection method.
    double d_left_point  = ax;
    double d_right_point = bx;

    double d_left_value = (this->*eval_fcn)( d_left_point );

    while ( std::fabs( d_right_point - d_left_point ) > tol ) {

        const double d_mid_point = (d_left_point + d_right_point) / 2.0;
        const double d_mid_value = (this->*eval_fcn)( d_mid_point );

        if ( d_mid_value == d_left_value ) {

            d_left_point = d_mid_point;
            d_left_value = d_mid_value;

        } else {

            d_right_point = d_mid_point;
        }
    }

    // Now integrate around the discontinuity using adaptive integration rules.
    return IntegrateMollSmooth( eval_fcn, x_i, ax,            d_left_point, tol )
         + IntegrateMollSmooth( eval_fcn, x_i, d_right_point, bx,           tol );
}


//------------------------------------------------------------------------------------------------------------
//! \brief  Used in the computation of numerical approximations of convolution integrals used to define
//!         mollified problem parameters.
//!
//! A recursive implementation of adaptive Gauss-Kronrod quadrature is used to compute the integral. This
//! function assumes that the integrand is _smooth_, and hence makes no attempt to accurately handle
//! discontinuities or singularities.
//!
//! \param[in]      eval_fcn    Pointer to the function to integrate against. This is the function to which
//!                             the mollifier is applied.
//! \param[in]      x_i         Center point of the mollifier.
//! \param[in]      ax          Integrate over \f$ [a_x,b_x] \f$.
//! \param[in]      bx          Integrate over \f$ [a_x,b_x] \f$.
//! \param[in]      tol         Error tolerance.
//!
//! \return     Returns an approximation of the specified integral.
//------------------------------------------------------------------------------------------------------------
double Problem::IntegrateMollSmooth (

    const ProblemEvalFcn eval_fcn,
    const double x_i,
    const double ax,
    const double bx,
    const double tol

) const {

    PRINT_STATUS( "Executing Problem::%s.\n", __func__ )

    const double hx = (bx - ax) / 2.0;
    const double cx = (bx + ax) / 2.0;

    double G = 0.0;
    double K = 0.0;

    for ( int64_t i = 0; i < this->GK_order; ++i ) {

        const double x_pt = cx + hx * this->GK_nodes[i];
        const double f_val = (this->*eval_fcn)( x_pt ) * (*this->moll)( x_i - x_pt );

        K += this->K_weights[i] * f_val;

        if ( i % 2 ) {  G += this->G_weights[i/2] * f_val;  }
    }

    K *= hx;
    G *= hx;

    if ( std::fabs(K - G) < tol ) {  return K;  }

    return this->IntegrateMollSmooth( eval_fcn, x_i, ax, cx, tol )
         + this->IntegrateMollSmooth( eval_fcn, x_i, cx, bx, tol );
}


# endif


//============================================================================================================
//=== TWO-DIMENSIONAL IMPLEMENTATIONS ========================================================================
//============================================================================================================

# if SPACE_DIMS == 2 || defined (DOXYCOMPILE)


//------------------------------------------------------------------------------------------------------------
//! \brief  Used in the computation of a numerical approximation of the convolution integral used to define
//!         the mollification of the given function at the specified point.
//!
//! This function assumes that the integrand has _at most one_ discontinuity within the domain of integration.
//! The purpose of this function is to find the (single) discontinuity (if one exists) and then call the
//! function IntegrateMollSmooth() in a way that integrates around the discontinuity.
//!
//! \param[in]      eval_fcn    Pointer to the function to integrate against. This is the function to which
//!                             the mollifier is applied.
//! \param[in]      x_i         \f$ x \f$ coordinate of center point of the mollifier.
//! \param[in]      y_j         \f$ y \f$ coordinate of center point of the mollifier.
//! \param[in]      ax          Lower limit of integration wrt \f$ x \f$.
//! \param[in]      bx          Upper limit of integration wrt \f$ x \f$.
//! \param[in]      ay          Lower limit of integration wrt \f$ y \f$.
//! \param[in]      by          Upper limit of integration wrt \f$ y \f$.
//! \param[in]      tol         Error tolerance.
//!
//! \return     Returns an approximation of the specified integral.
//------------------------------------------------------------------------------------------------------------
double Problem::IntegrateMollDiscontinuous (

    const ProblemEvalFcn eval_fcn,
    const double x_i,
    const double y_j,
    const double ax,
    const double bx,
    const double ay,
    const double by,
    const double tol

) const {

    PRINT_STATUS( "Executing Problem::%s.\n", __func__ )

    // Locate corner via bisection method.
    double d_x_left_point  = ax;
    double d_x_right_point = bx;
    double d_y_left_point  = ay;
    double d_y_right_point = by;

    // Apply bisection in the x dimension.
    double y_point = std::nan("0");

    if (    (this->*eval_fcn)( d_x_left_point, d_y_left_point )
         == (this->*eval_fcn)( d_x_right_point, d_y_left_point )
    ) {
        y_point = d_y_right_point;
    } else {
        y_point = d_y_left_point;
    }

    while ( std::abs( d_x_right_point - d_x_left_point ) > tol ) {

        const double d_x_mid_point = (d_x_left_point + d_x_right_point) / 2.0;

        if (    (this->*eval_fcn)( d_x_mid_point, y_point )
             == (this->*eval_fcn)( d_x_left_point, y_point )
        ) {
            d_x_left_point = d_x_mid_point;
        } else {
            d_x_right_point = d_x_mid_point;
        }
    }

    // Apply bisection in the y dimension.
    double x_point = std::nan("0");

    if (    (this->*eval_fcn)( d_x_left_point, d_y_left_point )
         == (this->*eval_fcn)( d_x_left_point, d_y_right_point )
    ) {
        x_point = d_x_right_point;
    } else {
        x_point = d_x_left_point;
    }

    while ( std::abs( d_y_right_point - d_y_left_point ) > tol ) {

        const double d_y_mid_point = (d_y_left_point + d_y_right_point) / 2.0;

        if (    (this->*eval_fcn)( x_point, d_y_mid_point )
             == (this->*eval_fcn)( x_point, d_y_left_point )
        ) {
            d_y_left_point = d_y_mid_point;
        } else {
            d_y_right_point = d_y_mid_point;
        }
    }

    // Now integrate around the corner using adaptive integration rule.
    return IntegrateMollSmooth( eval_fcn, x_i, y_j, ax, d_x_left_point,  ay, d_y_left_point,  tol )
         + IntegrateMollSmooth( eval_fcn, x_i, y_j, ax, d_x_left_point,  d_y_right_point, by, tol )
         + IntegrateMollSmooth( eval_fcn, x_i, y_j, d_x_right_point, bx, d_y_right_point, by, tol )
         + IntegrateMollSmooth( eval_fcn, x_i, y_j, d_x_right_point, bx, ay, d_y_left_point,  tol );
}


//------------------------------------------------------------------------------------------------------------
//! \brief  Used in the computation of numerical approximations of convolution integrals used to define
//!         mollified problem parameters.
//!
//! A recursive implementation of adaptive Gauss-Kronrod quadrature is used to compute the integral. This
//! function assumes that the integrand is _smooth_, and hence makes no attempt to accurately handle
//! discontinuities or singularities.
//!
//! \param[in]      eval_fcn    Pointer to the function to integrate against. This is the function to which
//!                             the mollifier is applied.
//! \param[in]      x_i         \f$ x \f$ coordinate of center point of the mollifier.
//! \param[in]      y_j         \f$ y \f$ coordinate of center point of the mollifier.
//! \param[in]      ax          Lower limit of integration wrt \f$ x \f$.
//! \param[in]      bx          Upper limit of integration wrt \f$ x \f$.
//! \param[in]      ay          Lower limit of integration wrt \f$ y \f$.
//! \param[in]      by          Upper limit of integration wrt \f$ y \f$.
//! \param[in]      tol         Error tolerance.
//!
//! \return     Returns an approximation of the specified integral.
//------------------------------------------------------------------------------------------------------------
double Problem::IntegrateMollSmooth (

    const ProblemEvalFcn eval_fcn,
    const double x_i,
    const double y_j,
    const double ax,
    const double bx,
    const double ay,
    const double by,
    const double tol

) const {

    PRINT_STATUS( "Executing Problem::%s.\n", __func__ )

    const double hx = (bx - ax) / 2.0;
    const double hy = (by - ay) / 2.0;

    const double cx = (bx + ax) / 2.0;
    const double cy = (by + ay) / 2.0;

    double G = 0.0;
    double K = 0.0;

    for ( int64_t i = 0; i < this->GK_order; ++i ) {
    for ( int64_t j = 0; j < this->GK_order; ++j ) {

        const double x_pt = cx + hx * this->GK_nodes[i];
        const double y_pt = cy + hy * this->GK_nodes[j];

        const double f_val = (this->*eval_fcn)( x_pt, y_pt ) * (*this->moll)( x_i - x_pt, y_j - y_pt );

        K += this->K_weights[i] * this->K_weights[j] * f_val;

        if ( i % 2 && j % 2 ) {  G += this->G_weights[i/2] * this->G_weights[j/2] * f_val;  }
    }}

    K *= hx;
    G *= hx;

    if ( std::fabs(K - G) < tol ) {  return K;  }

    return IntegrateMollSmooth( eval_fcn, x_i, y_j, ax, cx, ay, cy, tol )
         + IntegrateMollSmooth( eval_fcn, x_i, y_j, cx, bx, ay, cy, tol )
         + IntegrateMollSmooth( eval_fcn, x_i, y_j, ax, cx, cy, by, tol )
         + IntegrateMollSmooth( eval_fcn, x_i, y_j, cx, bx, cy, by, tol );
}


# endif // if SPACE_DIMS == 2
