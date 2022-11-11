//------------------------------------------------------------------------------------------------------------
//! \file   utils/Mollifier.cpp
//! \brief  Implementation of Mollifier class.
//!
//! \author Michael Crockatt
//! \date   May 2017
//------------------------------------------------------------------------------------------------------------


# include <cinttypes>
# include <cmath>
# include <cstdlib>
# include <cstring>
# include <limits>

# include "utils/CLog.hpp"
# include "utils/global.hpp"
# include "utils/Mollifier.hpp"
# include "utils/Quadrule/Quadrule.hpp"


using namespace Quadrule;


//============================================================================================================
//=== CONSTRUCTORS AND DESTRUCTOR ============================================================================
//============================================================================================================


//------------------------------------------------------------------------------------------------------------
//! \brief  Constructs a Mollifier object with the given parameters.
//!
//! Sets values specifying parameters of the object and computes the projection of the exact mollifier
//! function onto a DG finite element approximation with the specified tolerance \pp{tol}.
//!
//! \todo   Project unscaled bump. Show equation here, reference function ExponentialBump() above.
//!
//! \param[in]      radius_in       Desired radius of mollifier.
//! \param[in]      DG_degree_in    Maximum degree of DG basis polynomials for projection.
//! \param[in]      tol_in          Tolerance for DG approximation of exact mollifier function.
//------------------------------------------------------------------------------------------------------------
Mollifier::Mollifier (

    const double radius_in,
    const int64_t DG_degree_in, // Defaulted
    const double tol_in         // Defaulted.
):
    nx( 1 ),
    radius( radius_in ),
    DG_degree( DG_degree_in ),
    density{ nullptr },
    tol( tol_in ),
    normalization_const( 1.0 )
{
    PRINT_STATUS( "Constructing Mollifier object.\n" )

    // --- Setup quadrature. -------------------------------------------------------------------------- //

    //!
    //! Order of Gauss-Kronrod quadrature rule to use for computing integrals. This should be set to
    //! 15, 21, 31, or 41.
    //!
    const int64_t GK_order = 15;

    double * const GK_nodes  = new double[ GK_order ];
    double * const G_weights = new double[ GK_order/2 ];
    double * const K_weights = new double[ GK_order ];

    ComputeGKquad( GK_order, GK_nodes, G_weights, K_weights );

    this->normalization_const
        = 1.0 / ComputeNorm( 0.0, this->radius, GK_order, GK_nodes, G_weights, K_weights );

    // --- Compute projection onto DG basis. ---------------------------------------------------------- //

    double error = 0.0;
    double last_error = std::numeric_limits<double>::max();
    double current_error = std::numeric_limits<double>::max() / 2.0;

    while (     current_error > this->tol
            &&  last_error > current_error
    ) {

        last_error = current_error;
        current_error = 0.0;

        this->nx *= 2;
        this->dx = 1.0 / this->nx;

        delete [] this->density;

        this->dimof_density = this->nx * (this->DG_degree + 1);
        this->sizeof_density = this->dimof_density * sizeof(double);

        this->density = new double[ this->sizeof_density ];
        std::memset( this->density, 0, this->sizeof_density );

        for ( int64_t i = 0; i < this->nx; ++i ) {

            // Compute coefficients of expansion in each cell using Gaussian nodes.
            for ( int64_t l = 1; l < GK_order; l += 2 ) {

                // Quadrature point scaled to current cell.
                const double x_l = i * this->dx + (GK_nodes[l] + 1.0) * this->dx / 2.0;

                const double f_val = ExponentialBump( x_l );

                for ( int64_t d = 0; d <= this->DG_degree; ++d )
                    (*this)(i,d) += G_weights[l/2] * Legendre( d, GK_nodes[l] ) * (2*d + 1) * f_val / 2.0;
            }

            // Compute maximum pointwise error of approximation in current cell using Kronrod nodes.
            for ( int64_t l = 0; l < GK_order; ++l ) {

                // Quadrature point scaled to current cell.
                const double x_l = i * this->dx + (GK_nodes[l] + 1.0) * this->dx / 2.0;

                // Evaluate current approximation at given point.
                error = 0.0;

                for ( int64_t d = 0; d <= this->DG_degree; ++d )
                    error += (*this)(i,d) * Legendre( d, GK_nodes[l] );

                // Compute error in approximation at current point and store it if necessary.
                error = std::abs( error - ExponentialBump( x_l ) );

                current_error = std::max( current_error, error );
            }
        }

        PRINT_STATUS( "In %s last_error = %.4e  current_error = %.4e\n", __func__, last_error, current_error )
    }

    this->cell_factor = 1.0 / ( this->radius * this->radius * this->dx );

    // Cleanup memory.
    delete [] GK_nodes;
    delete [] G_weights;
    delete [] K_weights;

# if LOGLEVEL >= 3
    Print();
# endif
}


//------------------------------------------------------------------------------------------------------------
//! \brief  Destructor for Mollifier class.
//------------------------------------------------------------------------------------------------------------
Mollifier::~Mollifier ( void ) {

    delete [] this->density;
}


//============================================================================================================
//=== PUBLIC MEMBER FUNCTIONS ================================================================================
//============================================================================================================


//------------------------------------------------------------------------------------------------------------
//! \brief  Prints the parameters of the data stored in a Mollifier object to the global logging interface.
//!
//! \todo   Add prefix to Mollifier::Print statement.
//------------------------------------------------------------------------------------------------------------
void Mollifier::Print ( void ) const {

    PRINT_LOG( "\n" )
    PRINT_LOG( "Contents of mollifier object:\n" )
    PRINT_LOG( "\n" )

    PRINT_LOG( "   %-*s   %-*" PRId64 "\n", Global::col_width, "DG order:", Global::col_width, this->DG_degree )

    PRINT_LOG( "\n" )

    PRINT_LOG( "   %-*s   %-*" PRId64 "\n", Global::col_width, "Num cells:", Global::col_width, this->nx )
    PRINT_LOG( "   %-*s   %-*.2e\n",        Global::col_width, "Radius:",    Global::col_width, this->radius )
    PRINT_LOG( "   %-*s   %-*.2e\n",        Global::col_width, "Delta x:",   Global::col_width, this->dx )

    PRINT_LOG( "\n" )

    PRINT_LOG( "   %-*s   %-*.2e\n", Global::col_width, "Tolerance:",
                                     Global::col_width, this->tol )

    PRINT_LOG( "   %-*s   %-*.2e\n", Global::col_width, "Normalization:",
                                     Global::col_width, this->normalization_const )

    PRINT_LOG( "\n" )

    const size_t total_bytes = sizeof(Mollifier) + this->sizeof_density;

    PRINT_LOG( "   %-*s   %-*zu\n", Global::col_width, "Memory (Bytes):",
                                    Global::col_width, total_bytes )

    PRINT_LOG( "\n" )
}


//============================================================================================================
//=== PRIVATE MEMBER FUNCTIONS ===============================================================================
//============================================================================================================


//------------------------------------------------------------------------------------------------------------
//! \brief  One-dimensional \f$ C^{\infty} \f$ exponential function with compact support in
//!         \f$ |x| \leq 1.0 \f$.
//!
//! \todo   Documentation of mollifier bump.
//!
//! \param[in]      x           Point to evaluate exponential function at.
//!
//! \return     Returns the value of the exponential function with radius \f$ 1.0 \f$ at the specified point.
//------------------------------------------------------------------------------------------------------------
double Mollifier::ExponentialBump (

    const double x

) const {

    if ( x < 1.0 )
        return std::exp( 1.0 / ( x - 1.0 ) );
    else
        return 0.0;
}


//------------------------------------------------------------------------------------------------------------
//! \brief  Evaluates the exact value of the exponential mollifier at the given point.
//!
//! Evaluates the exponential mollifier
//! \f[ \eta_r (x) = \frac{1}{C} \begin{cases} \exp \left[ \left( \frac{ |x|^2 }{ r^2 } - 1 \right)^{-1} \right], & |x|^2 < r^2 \\ 0, & \text{otherwise}, \end{cases} \f]
//! at the given point, where \f$ C \f$ is a normalization constant.
//!
//! \param[in]      x           The point at which to evaluate the mollifier.
//!
//! \return     Returns the value of the exponential mollifier at the specified point.
//!
//! \see    ExponentialBump()
//------------------------------------------------------------------------------------------------------------
double Mollifier::Exact (

    const double x

) const {

    return ExponentialBump( (x*x) / (this->radius * this->radius) ) * this->normalization_const;
}


//------------------------------------------------------------------------------------------------------------
//! \brief  Evaluates the exact value of the exponential mollifier at the given point.
//!
//! Evaluates the exponential mollifier
//! \f[ \eta_r (x) = \frac{1}{C} \begin{cases} \exp \left[ \left( \frac{ |x|^2 }{ r^2 } - 1 \right)^{-1} \right], & |x|^2 < r^2 \\ 0, & \text{otherwise}, \end{cases} \f]
//! at the given point, where \f$ C \f$ is a normalization constant.
//!
//! \param[in]      x   \f$ x \f$ coordinate of the point at which to evaluate the mollifier.
//! \param[in]      y   \f$ y \f$ coordinate of the point at which to evaluate the mollifier.
//!
//! \return     Returns the value of the exponential mollifier at the specified point.
//!
//! \see    ExponentialBump()
//------------------------------------------------------------------------------------------------------------
double Mollifier::Exact (

    const double x,
    const double y

) const {

    return ExponentialBump( (x*x + y*y) / (this->radius * this->radius) ) * this->normalization_const;
}


//------------------------------------------------------------------------------------------------------------
//! \brief  Computes a numerical approximation to the \f$ L^1 \f$ of the exact mollifier.
//!
//! Computes a numerical approximation to the \f$ L^1 \f$ norm of the exact representation of the given
//! mollifier to the error tolerance given by Mollifier::tol using a recursive implementation of an adaptive
//! Gauss-Kronrod quadrature.
//!
//! The mollifier is represented radially, therefore the normalization constant is computed using an integral
//! in polar coordinates.
//!
//! \todo   Explanation of polar integrals for computing norm of mollifier.
//!
//! \param[in]      ax          Lower limit of integration.
//! \param[in]      bx          Upper limit of integration.
//! \param[in]      GK_order    Order of Gauss-Kronrod quadrature given to use.
//! \param[in]      GK_nodes    Nodes of Gauss-Kronrod quadrature given to use.
//! \param[in]      G_weights   Gaussian weights of quadrature given to use.
//! \param[in]      K_weights   Kronrod weights of quadrature given to use.
//!
//! \return     Returns an approximation of
//!             \f$ \displaystyle \int_{a_x}^{b_x} \texttt{ExponentialBump}(x,r) \, dx \f$
//!             to within an error estimated to be no larger than Mollifier::tol.
//------------------------------------------------------------------------------------------------------------
double Mollifier::ComputeNorm (

    const double ax,
    const double bx,
    const int64_t GK_order,
    const double * const GK_nodes,
    const double * const G_weights,
    const double * const K_weights

) const {

    const double hx = (bx - ax) / 2.0;
    const double cx = (bx + ax) / 2.0;

    double G = 0.0;
    double K = 0.0;

    for ( int64_t i = 0; i < GK_order; ++i ) {

        const double x_pt = cx + hx * GK_nodes[i];

    # if SPACE_DIMS == 1

        // Value of exponential bump times two (+/- symmetry).
        const double moll_val = 2.0 * Exact( x_pt );

    # elif SPACE_DIMS == 2

        // Value of exponential bump times circumference of circle.
        const double moll_val = 2.0 * M_PI * Exact( x_pt ) * x_pt;

    # elif SPACE_DIMS == 3

        // Value of exponential bump times surface of sphere.
        const double moll_val = 4.0 * M_PI * Exact( x_pt ) * x_pt * x_pt;

    # endif

        K += K_weights[i] * moll_val;

        if ( i % 2 ) {  G += G_weights[i/2] * moll_val;  }
    }

    K *= hx;
    G *= hx;

    if ( std::abs(K - G) < this->tol ) {  return K;  }

    return  ComputeNorm( ax, cx, GK_order, GK_nodes, G_weights, K_weights )
         +  ComputeNorm( cx, bx, GK_order, GK_nodes, G_weights, K_weights );
}


//--------------------------------------------------------------------------------------------------------
//! \brief  Evaluates the radially-defined mollifier at the specified radius from the origin.
//!
//! \param[in]  r   Radius to evaluate mollifier at.
//!
//! \return     Returns the value of the mollifier at the radius \pp{r}.
//--------------------------------------------------------------------------------------------------------
double Mollifier::EvalAtRadius (

    const double r

) const {

    const int64_t i = (int64_t) r;

    if ( i > this->nx ) {  return 0.0;  }

    // Compute \bar{x} to evaluate Legendre polynomials with.
    const double x_bar = 2.0 * r - (2*i + 1);

    // --- Evaluate density function at given point. -------------------------------------------------- //
    double pn  = 1.0;
    double pn1 = 0.0;
    double pn2;

    double ret_val = (*this)(i,0) * pn;

    for ( int64_t d = 0; d < this->DG_degree; ++d ) {

        pn2 = pn1;
        pn1 = pn;
        pn  = ( (2*d + 1) * x_bar * pn1 - d * pn2 ) / (d + 1);

        ret_val += (*this)(i,d+1) * pn;
    }

    return ret_val * this->normalization_const;
}


//============================================================================================================
//=== ONE SPATIAL DIMENSION ==================================================================================
//============================================================================================================

# if SPACE_DIMS == 1


//------------------------------------------------------------------------------------------------------------
//! \brief  Routine to obtain the value of mollifier at specified point.
//!
//! \param[in]      x           \f$ x \f$ coordinate of point to evaluate mollifier at.
//!
//! \return     Returns the value of the mollifier at the given spatial point.
//------------------------------------------------------------------------------------------------------------
double Mollifier::operator() (

    const double x

) const {

    return EvalAtRadius( ( x*x ) * this->cell_factor );
}


# endif // if SPACE_DIMS == 1


//============================================================================================================
//=== TWO SPATIAL DIMENSIONS =================================================================================
//============================================================================================================

# if SPACE_DIMS == 2


//------------------------------------------------------------------------------------------------------------
//! \brief  Routine to obtain the value of mollifier at specified point.
//!
//! \param[in]      x           \f$ x \f$ coordinate of point to evaluate mollifier at.
//! \param[in]      y           \f$ y \f$ coordinate of point to evaluate mollifier at.
//!
//! \return     Returns the value of the mollifier at the given spatial point.
//------------------------------------------------------------------------------------------------------------
double Mollifier::operator() (

    const double x,
    const double y

) const {

    return EvalAtRadius( ( x*x + y*y ) * this->cell_factor );
}


# endif // if SPACE_DIMS == 2


//============================================================================================================
//=== THREE SPATIAL DIMENSIONS ===============================================================================
//============================================================================================================

# if SPACE_DIMS == 3


//------------------------------------------------------------------------------------------------------------
//! \brief  Routine to obtain the value of mollifier at specified point.
//!
//! \param[in]      x           \f$ x \f$ coordinate of point to evaluate mollifier at.
//! \param[in]      y           \f$ y \f$ coordinate of point to evaluate mollifier at.
//! \param[in]      z           \f$ z \f$ coordinate of point to evaluate mollifier at.
//!
//! \return     Returns the value of the mollifier at the given spatial point.
//------------------------------------------------------------------------------------------------------------
double Mollifier::operator() (

    const double x,
    const double y,
    const double z

) const {

    return EvalAtRadius( ( x*x + y*y + z*z ) * this->cell_factor );
}


# endif // if SPACE_DIMS == 3
