//------------------------------------------------------------------------------------------------------------
//! \file   utils/Quadrule/OrdinateSet.cpp
//! \brief  Implementation of Quadrule::OrdinateSet class.
//!
//! \author Michael M. Crockatt
//! \date   December 2017
//------------------------------------------------------------------------------------------------------------


# include <cinttypes>
# include <cmath>
# include <cstring>
# include <random>

# include "utils/CLog.hpp"
# include "utils/global.hpp"
# include "utils/Quadrule/OrdinateSet.hpp"
# include "utils/Quadrule/Quadrule.hpp"
# include "utils/Quadrule/SphericalTriangle.hpp"

# if defined (ENABLE_LEBEDEV)
# include "utils/Quadrule/Lebedev.hpp"
# endif


using namespace Quadrule;


//============================================================================================================
//=== STATIC DEFINITIONS =====================================================================================
//============================================================================================================


//------------------------------------------------------------------------------------------------------------
//! \brief  Used to convert OrdinateSet::OrdinateType values to their corresponding string descriptors.
//------------------------------------------------------------------------------------------------------------
std::map< OrdinateSet::OrdinateType, std::string > OrdinateSet::OrdinateType_to_String = {

    // Common to all dimensions.
    { OrdinateSet::OrdinateType::None,                  "none"                      },

    // 1D only.
# if SPACE_DIMS == 1 || defined (DOXYCOMPILE)

    { OrdinateSet::OrdinateType::GaussLegendre,         "gauss-legendre"            },
    { OrdinateSet::OrdinateType::GaussLobatto,          "gauss-lobatto"             },
    { OrdinateSet::OrdinateType::DoubleGauss,           "double-gauss"              },
    { OrdinateSet::OrdinateType::DoubleRadau,           "double-radau"              }

# endif

    // 2D & 3D only.
# if SPACE_DIMS >= 2 || defined (DOXYCOMPILE)

    { OrdinateSet::OrdinateType::ChebyshevLegendre,     "chebyshev-legendre"        },
    { OrdinateSet::OrdinateType::SphericalTriangle,     "spherical-triangle"        },

# if defined (ENABLE_LEBEDEV)
    { OrdinateSet::OrdinateType::Lebedev,               "lebedev"                   }
# endif

# endif
};


//------------------------------------------------------------------------------------------------------------
//! \brief  Used to convert string descriptors to OrdinateSet::OrdinateType values.
//------------------------------------------------------------------------------------------------------------
std::map< std::string, OrdinateSet::OrdinateType > OrdinateSet::String_to_OrdinateType = {

    // Common to all dimensions.
    { "none",                       OrdinateSet::OrdinateType::None                 },

    // 1D only.
# if SPACE_DIMS == 1 || defined (DOXYCOMPILE)

    { "gauss-legendre",             OrdinateSet::OrdinateType::GaussLegendre        },
    { "gauss-lobatto",              OrdinateSet::OrdinateType::GaussLobatto         },
    { "double-gauss",               OrdinateSet::OrdinateType::DoubleGauss          },
    { "double-radau",               OrdinateSet::OrdinateType::DoubleRadau          }

# endif

    // 2D & 3D only.
# if SPACE_DIMS >= 2 || defined (DOXYCOMPILE)

    { "chebyshev-legendre",         OrdinateSet::OrdinateType::ChebyshevLegendre    },
    { "spherical-triangle",         OrdinateSet::OrdinateType::SphericalTriangle    },

# if defined (ENABLE_LEBEDEV)
    { "lebedev",                    OrdinateSet::OrdinateType::Lebedev              }
# endif

# endif
};


//------------------------------------------------------------------------------------------------------------
//! \brief  Caches previously computed OrdinateData structs so that they may be reused with new objects.
//------------------------------------------------------------------------------------------------------------
std::map< std::tuple< OrdinateSet::OrdinateType, int64_t, bool >, OrdinateSet::OrdinateData >
    OrdinateSet::ordinate_cache {};


//============================================================================================================
//=== CONSTRUCTORS, DESTRUCTOR, ALLOCATORS, AND ASSOCIATED HELPER ROUTINES ===================================
//============================================================================================================


//------------------------------------------------------------------------------------------------------------
//! \brief  Constructs an empty OrdinateSet object with zero-initialized parameters.
//------------------------------------------------------------------------------------------------------------
OrdinateSet::OrdinateSet ( void ) :

    ordinate_data {}
{}


//------------------------------------------------------------------------------------------------------------
//! \brief  Constructs an OrdinateSet object given an OrdinateType and order.
//!
//! \param[in]  ang_order           Order of discrete ordinate quadrature of specified type.
//! \param[in]  symmetric_reduce    (optional) <br/>
//!                                 Specifies whether or not the number of ordinates should be reduced by
//!                                 symmetry. Only used for 2D problems. Defaults to false.
//! \param[in]  ordinate_type       (optional) <br/>
//!                                 Type of discrete ordinate quadrature to use.
//!                                 Defaults to Gauss-Legendre quadrature in 1D and Chebyshev-Legendre
//!                                 quadrature in multi-D.
//------------------------------------------------------------------------------------------------------------
OrdinateSet::OrdinateSet (

    const int64_t ang_order,
    const bool symmetric_reduce,        // = false
    const OrdinateType ordinate_type    // Defaulted.
) :
    ordinate_data {}
{
    Reconfigure( ang_order, symmetric_reduce, ordinate_type );
}



//------------------------------------------------------------------------------------------------------------
//! \brief  Constructs an OrdinateSet object using the parameters of a given OrdinateSet object.
//!
//! Delegates to OrdinateSet::OrdinateSet( const int64_t, const bool, const OrdinateType ).
//!
//! \param[in]      that        OrdinateSet object containing parameters to use during construction.
//------------------------------------------------------------------------------------------------------------
OrdinateSet::OrdinateSet (

    const OrdinateSet & that

) : OrdinateSet( that.ordinate_data.order, that.ordinate_data.symmetric_reduce, that.ordinate_data.type ) {}


//------------------------------------------------------------------------------------------------------------
//! \brief  Reconfigures an OrdinateSet object using the parameters of a given OrdinateSet object.
//!
//! \param[in]      that        OrdinateSet object containing parameters with which to reconfigure the object.
//------------------------------------------------------------------------------------------------------------
OrdinateSet & OrdinateSet::Reconfigure (

    const OrdinateSet & that
) {
    return Reconfigure( that.GetAngOrder(), that.GetOrdinateSymmetry(), that.GetOrdinateType() );
}


//------------------------------------------------------------------------------------------------------------
//! \brief  Reconfigures an OrdinateSet object using the given parameters.
//!
//! \param[in]  ang_order           Order of discrete ordinates quadrature of specified type.
//! \param[in]  symmetric_reduce    (optional) <br/>
//!                                 Specifies whether or not the number of ordinates should be reduced by
//!                                 symmetry. Only used for 2D problems. Defaults to false.
//! \param[in]  ordinate_type       (optional) <br/>
//!                                 Type of discrete ordinate quadrature to use.
//!                                 Defaults to Gauss-Legendre quadrature in 1D and Chebyshev-Legendre
//!                                 quadrature in multi-D.
//------------------------------------------------------------------------------------------------------------
OrdinateSet & OrdinateSet::Reconfigure (

    const int64_t ang_order,
    const bool symmetric_reduce,        // = false
    const OrdinateType ordinate_type    // Defaulted.
) {

    // Avoid compiler warning with 1D code.
# if SPACE_DIMS == 1
    (void) symmetric_reduce;
# endif

    // Search cache for ordinate data.
    const std::tuple< OrdinateType, int64_t, bool > ordinate_tuple{ ordinate_type, ang_order,
                                                                # if SPACE_DIMS == 2
                                                                    symmetric_reduce
                                                                # else
                                                                    false
                                                                # endif
                                                                  };

    auto search = ordinate_cache.find( ordinate_tuple );

    // If data not found in cache, create it in the cache, and then try again.
    if ( search == ordinate_cache.end() ) {

        switch ( ordinate_type ) {

        # if SPACE_DIMS == 1

            case OrdinateType::GaussLegendre:

                SetupGaussLegendre( ang_order );
                break;

            case OrdinateType::GaussLobatto:

                SetupGaussLobatto( ang_order );
                break;

            case OrdinateType::DoubleGauss:

                SetupDoubleGauss( ang_order );
                break;

            case OrdinateType::DoubleRadau:

                SetupDoubleRadau( ang_order );
                break;

        # endif // if SPACE_DIMS == 1

        # if SPACE_DIMS >= 2

            case OrdinateType::ChebyshevLegendre:

                SetupChebyshevLegendre( ang_order, symmetric_reduce );
                break;

            case OrdinateType::SphericalTriangle:

                SetupSphericalTriangle( ang_order, symmetric_reduce );
                break;

        # if defined (ENABLE_LEBEDEV)

            case OrdinateType::Lebedev:

                SetupLebedev( ang_order, symmetric_reduce );
                break;

        # endif // if defined (ENABLE_LEBEDEV)

        # endif // if SPACE_DIMS >= 2

            default:
            {   std::string error_message = "Invalid OrdinateType '"
                                            + OrdinateType_to_String.at( ordinate_type )
                                            + "' for constructing OrdinateSet.\n";

                PRINT_ERROR( error_message.c_str() );
                throw std::invalid_argument( error_message );
            }
        }

        search = ordinate_cache.find( ordinate_tuple );
    }

    if ( search == ordinate_cache.end() ) {

        std::string error_message = "In '" + std::string(__func__)
                                    + "': OrdinateData for OrdinateSet '"
                                    + OrdinateType_to_String.at( ordinate_type )
                                    + ":" + std::to_string( ang_order )
                                    + "' not found in cache after add. "
                                    + "Something must have gone horribly wrong.\n";

        PRINT_ERROR( error_message.c_str() )
        throw std::runtime_error( error_message );
    }

    this->ordinate_data = std::get<1>( *search );

    return *this;
}


//============================================================================================================
//=== ONE-DIMENSIONAL QUADRATURE INITIALIZATION ROUTINES =====================================================
//============================================================================================================

# if SPACE_DIMS == 1 || defined (DOXYCOMPILE)


//------------------------------------------------------------------------------------------------------------
//! \brief  Adds the Gauss-Legendre quadrature of order \pp{order} to OrdinateSet::cache.
//!
//! \param[in]      order       Order of the Gauss-Legendre quadrature to add to the cache.
//------------------------------------------------------------------------------------------------------------
void OrdinateSet::SetupGaussLegendre (

    const int64_t order
) {

    OrdinateData new_data {};

    new_data.type = OrdinateType::GaussLegendre;
    new_data.order = order;
    new_data.num_nodes = order;

    new_data.xi_ptr  = new double[ new_data.num_nodes ];
    new_data.weights = new double[ new_data.num_nodes ];

    ComputeQuadrature( new_data.num_nodes, new_data.xi_ptr, new_data.weights, NodesType::GaussLegendre );

    ordinate_cache.insert( std::make_pair( std::make_tuple( new_data.type, new_data.order, false ), new_data ) );
}


//------------------------------------------------------------------------------------------------------------
//! \brief  Adds the Gauss-Lobatto quadrature of order \pp{order} to OrdinateSet::cache.
//!
//! \param[in]      order       Order of the Gauss-Lobatto quadrature to add to the cache.
//------------------------------------------------------------------------------------------------------------
void OrdinateSet::SetupGaussLobatto (

    const int64_t order
) {

    OrdinateData new_data {};

    new_data.type = OrdinateType::GaussLobatto;
    new_data.order = order;
    new_data.num_nodes = order;

    new_data.xi_ptr  = new double[ new_data.num_nodes ];
    new_data.weights = new double[ new_data.num_nodes ];

    ComputeQuadrature( new_data.num_nodes, new_data.xi_ptr, new_data.weights, NodesType::GaussLobatto );

    ordinate_cache.insert( std::make_pair( std::make_tuple( new_data.type, new_data.order, false ), new_data ) );
}


//------------------------------------------------------------------------------------------------------------
//! \brief  Adds the double-Gauss quadrature of order \pp{order} to OrdinateSet::cache.
//!
//! \param[in]      order       Order of the double-Gauss quadrature to add to the cache.
//------------------------------------------------------------------------------------------------------------
void OrdinateSet::SetupDoubleGauss (

    const int64_t order
) {

    if ( order % 2 != 0 ) {

        std::string error_message = "Order of double-Gauss quadrature must be even ("
                                    + std::to_string( order ) + " given).\n";

        PRINT_ERROR( error_message.c_str() )
        throw std::invalid_argument( error_message );
    }

    OrdinateData new_data {};

    new_data.type = OrdinateType::DoubleGauss;
    new_data.order = order;
    new_data.num_nodes = order;

    new_data.xi_ptr  = new double[ new_data.num_nodes ];
    new_data.weights = new double[ new_data.num_nodes ];

    // Compute Gauss-Legendre quadrature.
    ComputeQuadrature( new_data.num_nodes /2, new_data.xi_ptr + new_data.order /2,
                       new_data.weights + new_data.order /2, NodesType::GaussLegendre, 0.0, 1.0 );

    // Reflect quadrature across origin.
    for ( int64_t i = 0; i < new_data.order /2; ++i ) {

        new_data.xi_ptr[i] = - new_data.xi_ptr[ new_data.order -i-1 ];
        new_data.weights[i] = new_data.weights[ new_data.order -i-1 ];
    }

    ordinate_cache.insert( std::make_pair( std::make_tuple( new_data.type, new_data.order, false ), new_data ) );
}


//------------------------------------------------------------------------------------------------------------
//! \brief  Adds the double-Radau quadrature of order \pp{order} to OrdinateSet::cache.
//!
//! \param[in]      order       Order of the double-Radu quadrature to add to the cache.
//------------------------------------------------------------------------------------------------------------
void OrdinateSet::SetupDoubleRadau (

    const int64_t order
) {

    if ( order % 2 != 0 ) {

        std::string error_message = "Order of double-Radau quadrature must be even ("
                                    + std::to_string( order ) + " given).\n";

        PRINT_ERROR( error_message.c_str() )
        throw std::invalid_argument( error_message );
    }

    OrdinateData new_data {};

    new_data.type = OrdinateType::DoubleRadau;
    new_data.order = order;
    new_data.num_nodes = order;

    new_data.xi_ptr  = new double[ new_data.num_nodes ];
    new_data.weights = new double[ new_data.num_nodes ];

    // Compute Gauss-Legendre quadrature.
    ComputeQuadrature( new_data.num_nodes /2, new_data.xi_ptr + new_data.order /2,
                       new_data.weights + new_data.order /2, NodesType::GaussRadau, 0.0, 1.0 );

    // Reflect quadrature across origin.
    for ( int64_t i = 0; i < new_data.order /2; ++i ) {

        new_data.xi_ptr[i] = - new_data.xi_ptr[ new_data.order -i-1 ];
        new_data.weights[i] = new_data.weights[ new_data.order -i-1 ];
    }

    ordinate_cache.insert( std::make_pair( std::make_tuple( new_data.type, new_data.order, false ), new_data ) );
}


# endif // if SPACE_DIMS == 1


//============================================================================================================
//=== MULTI-DIMENSIONAL QUADRATURE INITIALIZATION ROUTINES ===================================================
//============================================================================================================

# if SPACE_DIMS >= 2 || defined (DOXYCOMPILE)


//------------------------------------------------------------------------------------------------------------
//! \brief  Adds the Chebyshev-Legendre quadrature of order \pp{order} to OrdinateSet::cache.
//!
//! \param[in]  order               Order of the Legendre-Chebyshev quadrature to add to the cache.
//! \param[in]  symmetric_reduce    Specifies whether or not the number of ordinates should be reduced by
//!                                 symmetry. Only used for 2D problems.
//------------------------------------------------------------------------------------------------------------
void OrdinateSet::SetupChebyshevLegendre (

    const int64_t order,
    const bool symmetric_reduce
) {

    OrdinateData new_data {};

    new_data.type = OrdinateType::ChebyshevLegendre;
    new_data.order = order;
    new_data.num_nodes = 2 * order * order;
    new_data.symmetric_reduce = symmetric_reduce;

    new_data.xi_ptr  = new double[ new_data.num_nodes ];
    new_data.eta_ptr = new double[ new_data.num_nodes ];
    new_data.mu_ptr  = new double[ new_data.num_nodes ];
    new_data.phi_ptr = new double[ new_data.num_nodes ];
    new_data.weights = new double[ new_data.num_nodes ];

    double * const gauss_nodes   = new double[ order ];
    double * const gauss_weights = new double[ order ];

    ComputeQuadrature( order, gauss_nodes, gauss_weights, NodesType::GaussLegendre );

    int64_t q = 0;

    for ( int64_t q1 = 0; q1 < order;   ++q1 ) {
    for ( int64_t q2 = 0; q2 < 2*order; ++q2 ) {

        const double phi = (q2 + 0.5) * M_PI / order;
        const double mu  = gauss_nodes[q1];

        new_data.weights[q] = M_PI / order * gauss_weights[q1];

        new_data.xi_ptr[q]  = std::sqrt( 1.0 - mu*mu ) * std::cos( phi );
        new_data.eta_ptr[q] = std::sqrt( 1.0 - mu*mu ) * std::sin( phi );
        new_data.mu_ptr[q]  = mu;
        new_data.phi_ptr[q] = phi - M_PI;  // Compute values for ϕ ∈ [-π,π].

        ++q;
    }}

    new_data.Sort();
    new_data.LocateOctants();

    delete [] gauss_nodes;
    delete [] gauss_weights;

# if SPACE_DIMS == 2

    if ( new_data.symmetric_reduce )
        new_data.SymmetricReduce();

    if ( BitmaskHasAny( Global::reflecting_boundaries, BoundaryEdge::All ) )
        new_data.ComputeReflectMap();

# endif // if SPACE_DIMS == 2

    const std::tuple< OrdinateType, int64_t, bool > ordinate_tuple{ new_data.type, new_data.order,
                                                                # if SPACE_DIMS == 2
                                                                    new_data.symmetric_reduce
                                                                # else
                                                                    false
                                                                # endif
                                                                  };

    ordinate_cache.insert( std::make_pair( ordinate_tuple, new_data ) );
}


//------------------------------------------------------------------------------------------------------------
//! \brief  Adds the spherical triangle (\f$ T_N \f$) quadrature of order \pp{order} to OrdinateSet::cache.
//!
//! \param[in]  order               Order of quadrature to add to the cache.
//! \param[in]  symmetric_reduce    Specifies whether or not the number of ordinates should be reduced by
//!                                 symmetry. Only used for 2D problems.
//------------------------------------------------------------------------------------------------------------
void OrdinateSet::SetupSphericalTriangle (

    const int64_t order,
    const bool symmetric_reduce
) {

    OrdinateData new_data {};

    new_data.type = OrdinateType::SphericalTriangle;
    new_data.order = order;
    new_data.num_nodes = 8 * new_data.order * new_data.order;
    new_data.symmetric_reduce = symmetric_reduce;

    new_data.xi_ptr  = new double[ new_data.num_nodes ];
    new_data.eta_ptr = new double[ new_data.num_nodes ];
    new_data.mu_ptr  = new double[ new_data.num_nodes ];
    new_data.phi_ptr = new double[ new_data.num_nodes ];
    new_data.weights = new double[ new_data.num_nodes ];

    ComputeSphericalTriangleQuadrature( order, new_data.xi_ptr, new_data.eta_ptr, new_data.mu_ptr, new_data.weights );

    new_data.FillByOctahedralReflection();
    new_data.ComputePhi();
    new_data.LocateOctants();

# if SPACE_DIMS == 2

    if ( new_data.symmetric_reduce )
        new_data.SymmetricReduce();

    if ( BitmaskHasAny( Global::reflecting_boundaries, BoundaryEdge::All ) )
        new_data.ComputeReflectMap();

# endif // if SPACE_DIMS == 2

    const std::tuple< OrdinateType, int64_t, bool > ordinate_tuple{ new_data.type, new_data.order,
                                                                # if SPACE_DIMS == 2
                                                                    new_data.symmetric_reduce
                                                                # else
                                                                    false
                                                                # endif
                                                                  };

    ordinate_cache.insert( std::make_pair( ordinate_tuple, new_data ) );
}


# if defined (ENABLE_LEBEDEV)

//------------------------------------------------------------------------------------------------------------
//! \brief  Adds the Lebedev quadrature of order \pp{order} to OrdinateSet::cache.
//!
//! \param[in]  order               Order of the Lebedev quadrature to add to the cache.
//! \param[in]  symmetric_reduce    Specifies whether or not the number of ordinates should be reduced by
//!                                 symmetry. Only used for 2D problems.
//------------------------------------------------------------------------------------------------------------
void OrdinateSet::SetupLebedev (

    const int64_t order,
    const bool symmetric_reduce
) {

    OrdinateData new_data {};

    new_data.type = OrdinateType::Lebedev;
    new_data.order = order;
    new_data.symmetric_reduce = symmetric_reduce;

    try         {  new_data.num_nodes = LebedevOrderToCount.at( order );  }
    catch (...) {

        PRINT_ERROR( "Invalid order %" PRId64 " for Lebedev quadrature in %s.\n", order, __func__ )
        throw;
    }

    new_data.xi_ptr  = new double[ new_data.num_nodes ];
    new_data.eta_ptr = new double[ new_data.num_nodes ];
    new_data.mu_ptr  = new double[ new_data.num_nodes ];
    new_data.phi_ptr = new double[ new_data.num_nodes ];
    new_data.weights = new double[ new_data.num_nodes ];

    ComputeLebedevQuadrature( order, new_data.xi_ptr, new_data.eta_ptr, new_data.mu_ptr, new_data.weights );

    // Scale quadrature weights from unit normalization to 4π normalization.
    for ( int64_t q = 0; q < new_data.num_nodes; ++q )
        new_data.weights[q] *= 4.0 * M_PI;

    new_data.ComputePhi();
    new_data.Sort();
    new_data.LocateOctants();

# if SPACE_DIMS == 2

    if ( new_data.symmetric_reduce )
        new_data.SymmetricReduce();

    if ( BitmaskHasAny( Global::reflecting_boundaries, BoundaryEdge::All ) )
        new_data.ComputeReflectMap();

# endif // if SPACE_DIMS == 2

    const std::tuple< OrdinateType, int64_t, bool > ordinate_tuple{ new_data.type, new_data.order,
                                                                # if SPACE_DIMS == 2
                                                                    new_data.symmetric_reduce
                                                                # else
                                                                    false
                                                                # endif
                                                                  };

    ordinate_cache.insert( std::make_pair( ordinate_tuple, new_data ) );
}

# endif // if defined (ENABLE_LEBEDEV)


# endif // if SPACE_DIMS >= 2


//============================================================================================================
//=== STATIC MEMBER FUNCTIONS ================================================================================
//============================================================================================================

# if SPACE_DIMS >= 2

//------------------------------------------------------------------------------------------------------------
//! \brief  Sets an array of boolean values specifying the signs on the ordinates in a given quadrant (in 2D)
//!         or octant (in 3D).
//!
//! \param[in]  quad    Index of quadrant (in 2D) or octant (in 3D) to set signs for.
//! \param[out] signs   Array in which to set boolean values.
//------------------------------------------------------------------------------------------------------------
void OrdinateSet::SetSigns (

    const int64_t quad,
    bool (& signs) [SPACE_DIMS]
) {
# if SPACE_DIMS == 2

    switch ( quad ) {

        // First quadrant: ξ > 0, η > 0 -- index: 0b11 - 0b11 = 0.
        case 0:
            signs[0] = false;
            signs[1] = false;
            break;

        // Second quadrant: ξ < 0, η > 0 -- index: 0b11 - 0b01 = 2.
        case 2:
            signs[0] = true;
            signs[1] = false;
            break;

        // Third quadrant: ξ < 0, η < 0 -- index: 0b11 - 0b00 = 3.
        case 3:
            signs[0] = true;
            signs[1] = true;
            break;

        // Fourth quadrant: ξ > 0, η < 0 -- index: 0b11 - 0b10 = 1.
        case 1:
            signs[0] = false;
            signs[1] = true;
            break;

        default:
        {   std::string error_message = "Invalid quadrant index " + std::to_string(quad)
                                        + "' in '" + std::string(__func__) + "'.\n";

            PRINT_ERROR( error_message.c_str() )
            throw std::runtime_error( error_message );
        }
    }

# elif SPACE_DIMS == 3

    # warning "Implementation of OrdinateSet::SetSigns incomplete."

# endif // if SPACE_DIMS == ?
}

# endif // if SPACE_DIMS >= 2


//============================================================================================================
//=== ADDITIONAL MEMBER FUNCTIONS ============================================================================
//============================================================================================================


//------------------------------------------------------------------------------------------------------------
//! \brief  Prints the quadrature rule to the logging interface.
//------------------------------------------------------------------------------------------------------------
void OrdinateSet::Print (

    const std::string & prefix // = "  "

) const {

    char line[100] = "\0";
    strncat( line, "-------------------------------------------------", Global::col_width );

    PRINT_LOG( "%s%-*s  %s\n", prefix.c_str(),
               Global::col_width, "Ordinate type:",
               OrdinateType_to_String.at( this->ordinate_data.type ).c_str() )

    PRINT_LOG( "%s%-*s % " PRId64 "\n", prefix.c_str(),
               Global::col_width, "Angular order:",
               this->GetAngOrder() )

# if SPACE_DIMS >= 2

    PRINT_LOG( "%s%-*s  %s\n", prefix.c_str(),
               Global::col_width, "Reduce ordinates:",
               ( this->GetOrdinateSymmetry() ? "yes" : "no" ) )

# endif // if SPACE_DIMS >= 2

# if LOGLEVEL >= 3

    PRINT_LOG( "%s%-*s % " PRId64 "\n", prefix.c_str(),
               Global::col_width, "Num ordinates:",
               this->nq() )

    PRINT_LOG( "\n" )

    PRINT_LOG( "   %-*s   %-*s   %-*s   %-*s   %-*s\n",
               Global::col_width, "  xi",   Global::col_width, "  eta",     Global::col_width, "  mu",
               Global::col_width, "  phi",  Global::col_width, "  weight" )

    PRINT_LOG( "   %-*s   %-*s   %-*s   %-*s   %-*s\n",
               Global::col_width, line,     Global::col_width, line,    Global::col_width, line,
               Global::col_width, line,     Global::col_width, line )

# if SPACE_DIMS == 1

    for ( int64_t q = 0; q < this->nq(); ++q ) {

        PRINT_LOG( "     % -*.2e     % -*.2e\n",
                   Global::col_width -2, this->xi(q),
                   Global::col_width -2, this->w(q) )
    }

# elif SPACE_DIMS >= 2

    for ( int64_t oct = 0; oct < 8; ++oct ) {

        PRINT_LOG( "Octant %d:  [ %d, %d )\n", (int) oct,
                   (int) this->Octants(oct), (int) this->Octants(oct + 1) )

        PRINT_LOG( "\n" )

        for ( int64_t q = this->Octants(oct); q < this->Octants(oct+1); ++q ) {

            PRINT_LOG( "     % -*.2e     % -*.2e     % -*.2e     % -*.2e     % -*.2e\n",
                       Global::col_width -2, this->xi(q),
                       Global::col_width -2, this->eta(q),
                       Global::col_width -2, this->mu(q),
                       Global::col_width -2, this->OrdinatePhi(q),
                       Global::col_width -2, this->w(q) )
        }

        PRINT_LOG( "\n" )
    }

# endif // if SPACE_DIMS == ?

    PRINT_LOG( "\n" )

# endif // if LOGLEVEL >= 3
}


//============================================================================================================
//=== MEMBER FUNCTIONS OF OrdinateData SUB-CLASS =============================================================
//============================================================================================================

# if SPACE_DIMS >= 2


//------------------------------------------------------------------------------------------------------------
//! \brief  Generates nodes and weights for octants II through VIII from those of octant I by reflections.
//!
//! \attention  Assumes that all ordinates lie on the interior of the octant (i.e., none of the Cartesian
//!             coordinates of any ordinate is zero).
//!
//! Only generates values for OrdinateSet::xi_ptr, OrdinateSet::eta_ptr, OrdinateSet::mu_ptr, and
//! OrdinateSet::weights (\e not OrdinateSet::phi_ptr).
//!
//! \see    OrdinateSet::OrdinateData::ComputePhi()
//------------------------------------------------------------------------------------------------------------
void OrdinateSet::OrdinateData::FillByOctahedralReflection ( void ) const {

    // Reflect octant I to octant V by switching the sign on μ.
    for ( int64_t q = 0; q < num_nodes / 8; ++q ) {

        xi_ptr [ num_nodes / 8 + q ] =   xi_ptr [ q ];
        eta_ptr[ num_nodes / 8 + q ] =   eta_ptr[ q ];
        mu_ptr [ num_nodes / 8 + q ] = - mu_ptr [ q ];
        weights[ num_nodes / 8 + q ] =   weights[ q ];
    }

    // Reflect octants I,V to octants IV,VIII by switching the sign on η.
    for ( int64_t q = 0; q < num_nodes / 4; ++q ) {

        xi_ptr [ num_nodes / 4 + q ] =   xi_ptr [ q ];
        eta_ptr[ num_nodes / 4 + q ] = - eta_ptr[ q ];
        mu_ptr [ num_nodes / 4 + q ] =   mu_ptr [ q ];
        weights[ num_nodes / 4 + q ] =   weights[ q ];
    }

    // Reflect octants I,V,IV,VIII to octants II,VI,III,VII by switching the sign on ξ.
    for ( int64_t q = 0; q < num_nodes / 2; ++q ) {

        xi_ptr [ num_nodes / 2 + q ] = - xi_ptr [ q ];
        eta_ptr[ num_nodes / 2 + q ] =   eta_ptr[ q ];
        mu_ptr [ num_nodes / 2 + q ] =   mu_ptr [ q ];
        weights[ num_nodes / 2 + q ] =   weights[ q ];
    }
}


//------------------------------------------------------------------------------------------------------------
//! \brief  Computes azimuthal angles for each ordinate from the Cartesian coordinates.
//------------------------------------------------------------------------------------------------------------
void OrdinateSet::OrdinateData::ComputePhi ( void ) const {

    for ( int64_t q = 0; q < num_nodes; ++q )
        phi_ptr[q] = std::atan2( eta_ptr[q], xi_ptr[q] );
}


//------------------------------------------------------------------------------------------------------------
//! \brief  Returns an unsigned integer encoding a binary description of the octant in which the ordinate
//!         lies.
//!
//! Uses signs of Cartesian coordinates of nodes to determine the octant in which the ordinate with given
//! index lies. Boundaries between octants are treated as special cases to balance the number of ordinates
//! per octant for quadratures using points on the boundaries between octants.
//!
//! \todo   Detailed documentation of OrdinateSet::OrdinateData::DetermineOctant.
//------------------------------------------------------------------------------------------------------------
uint64_t OrdinateSet::OrdinateData::DetermineOctant (

    const int64_t q

) const {

    uint64_t oct;

    // Poles of the sphere.
    if (    xi_ptr[q] == 0.0
         && eta_ptr[q] == 0.0
    ) {
        if ( mu_ptr[q] > 0.0 ) {  oct = 0b111u;  }
        else                   {  oct = 0b010u;  }

        return oct;
    }

    // Equator of the sphere.
    if ( mu_ptr[q] == 0.0 ) {

        if      ( phi_ptr[q] < -0.75*M_PI ) {  oct = 0b000u;  }
        else if ( phi_ptr[q] < -0.50*M_PI ) {  oct = 0b001u;  }
        else if ( phi_ptr[q] < -0.25*M_PI ) {  oct = 0b101u;  }
        else if ( phi_ptr[q] <  0.0       ) {  oct = 0b100u;  }
        else if ( phi_ptr[q] <  0.25*M_PI ) {  oct = 0b110u;  }
        else if ( phi_ptr[q] <  0.50*M_PI ) {  oct = 0b111u;  }
        else if ( phi_ptr[q] <  0.75*M_PI ) {  oct = 0b011u;  }
        else if ( phi_ptr[q] <       M_PI ) {  oct = 0b010u;  }
        else                                {  oct = 0b000u;  }

        return oct;
    }

    // Other boundaries between octants.
    if ( eta_ptr[q] == 0.0 ) {

        if ( xi_ptr[q] > 0.0 ) {  oct = 0b110u;  }
        else                   {  oct = 0b000u;  }

        return oct | ((uint64_t)(mu_ptr[q] > 0.0));
    }

    if ( xi_ptr[q] == 0.0 ) {

        if ( eta_ptr[q] > 0.0 ) {  oct = 0b010;  }
        else                    {  oct = 0b100;  }

        return oct | ((uint64_t)(mu_ptr[q] > 0.0));
    }

    // Formula valid for points on the interior of an octant.
    oct =   ( ((uint64_t)(xi_ptr [q] > 0.0)) << 2 )
          | ( ((uint64_t)(eta_ptr[q] > 0.0)) << 1 )
          | ( ((uint64_t)(mu_ptr [q] > 0.0)) << 0 );

    return oct;
}



//------------------------------------------------------------------------------------------------------------
//! \brief  Locates starting indices of each octant and stores the values in OrdinateData::octants.
//!
//! \attention  This routine should always be called immediately after an OrdinateData object is filled with
//!             quadrature values otherwise the values stored in OrdinateData::octants are undefined.
//!
//! Assumes that the ordinates are sorted by octant in the order: I, V, IV, VIII, II, VI, III, VII
//! (this is the same ordering produced by OrdinateSet::OrdinateData::Sort).
//------------------------------------------------------------------------------------------------------------
void OrdinateSet::OrdinateData::LocateOctants ( void ) {

    int64_t q = 0;
    octants[0] = q;

    while ( q < num_nodes  &&  DetermineOctant(q) == 0b111u ) {  ++q;  }    octants[1] = q;
    while ( q < num_nodes  &&  DetermineOctant(q) == 0b110u ) {  ++q;  }    octants[2] = q;
    while ( q < num_nodes  &&  DetermineOctant(q) == 0b101u ) {  ++q;  }    octants[3] = q;
    while ( q < num_nodes  &&  DetermineOctant(q) == 0b100u ) {  ++q;  }    octants[4] = q;
    while ( q < num_nodes  &&  DetermineOctant(q) == 0b011u ) {  ++q;  }    octants[5] = q;
    while ( q < num_nodes  &&  DetermineOctant(q) == 0b010u ) {  ++q;  }    octants[6] = q;
    while ( q < num_nodes  &&  DetermineOctant(q) == 0b001u ) {  ++q;  }    octants[7] = q;
    while ( q < num_nodes  &&  DetermineOctant(q) == 0b000u ) {  ++q;  }    octants[8] = q;
}


//------------------------------------------------------------------------------------------------------------
//! \brief  Compares the coordinates of the ordinates with given indices to define an ordering on the
//!         ordinates suitable to sort them by octant.
//!
//! Returns +1 if \pp{i} is "greater" than \pp{j}, -1 if \pp{i} is "less" than \pp{j}, and 0 if they are
//! equal. If the ordinates are sorted from "largest" to "smallest", then they are grouped into octants in the
//! order: I, V, IV, VIII, II, VI, III, VII.
//!
//! \attention  Requires that the following coordinates for each point be defined:
//!             \f$ \xi, \eta, \mu, \phi \f$.
//------------------------------------------------------------------------------------------------------------
int64_t OrdinateSet::OrdinateData::Compare (

    const int64_t i,
    const int64_t j

) const {

    const uint64_t i_val = DetermineOctant(i);
    const uint64_t j_val = DetermineOctant(j);

    // First compare octants.
    if ( i_val > j_val ) {  return +1;  }
    if ( i_val < j_val ) {  return -1;  }

    // If ordinates lie in the same octant, then compare the values of the Cartesian coordinates.
    if ( xi_ptr[i] > xi_ptr[j] ) {  return +1;  }
    if ( xi_ptr[i] < xi_ptr[j] ) {  return -1;  }

    if ( eta_ptr[i] > eta_ptr[j] ) {  return +1;  }
    if ( eta_ptr[i] < eta_ptr[j] ) {  return -1;  }

    if ( mu_ptr[i] > mu_ptr[j] ) {  return +1;  }
    if ( mu_ptr[i] < mu_ptr[j] ) {  return -1;  }

    return 0;
}


//------------------------------------------------------------------------------------------------------------
//! \brief  Swaps the position of two ordinates with given indices.
//------------------------------------------------------------------------------------------------------------
void OrdinateSet::OrdinateData::Swap (

    const int64_t i,
    const int64_t j

) const {

    if ( i == j ) {  return;  }

    double temp;

    temp = xi_ptr [i];      xi_ptr [i] = xi_ptr [j];    xi_ptr [j] = temp;
    temp = eta_ptr[i];      eta_ptr[i] = eta_ptr[j];    eta_ptr[j] = temp;
    temp = mu_ptr [i];      mu_ptr [i] = mu_ptr [j];    mu_ptr [j] = temp;
    temp = phi_ptr[i];      phi_ptr[i] = phi_ptr[j];    phi_ptr[j] = temp;
    temp = weights[i];      weights[i] = weights[j];    weights[j] = temp;
}


//------------------------------------------------------------------------------------------------------------
// \brief  Returns the index of the parent of the node \pp{i}.
//------------------------------------------------------------------------------------------------------------
// static int64_t Parent( const int64_t i ) {  return (i - 1) / 2;  }


//------------------------------------------------------------------------------------------------------------
//! \brief  Returns the index of the left child of the node \pp{i}.
//------------------------------------------------------------------------------------------------------------
static int64_t Left( const int64_t i ) {  return 2*i + 1;  }


//------------------------------------------------------------------------------------------------------------
//! \brief  Returns the index of the right child of the node \pp{i}.
//------------------------------------------------------------------------------------------------------------
static int64_t Right( const int64_t i ) {  return 2*(i + 1);  }


//------------------------------------------------------------------------------------------------------------
//! \brief  Enforces the min-heap property on the subtree with root index \pp{i}, considering only the
//!         sub-array of length \pp{n}.
//------------------------------------------------------------------------------------------------------------
void OrdinateSet::OrdinateData::Heapify (

    const int64_t i,
    const int64_t n

) const {

    const int64_t l = Left(i);
    const int64_t r = Right(i);

    int64_t smallest = i;

    if ( ( l < n ) && ( Compare( l, smallest ) < 0 ) ) {  smallest = l;  }
    if ( ( r < n ) && ( Compare( r, smallest ) < 0 ) ) {  smallest = r;  }

    if ( smallest != i ) {

        Swap( i, smallest );
        Heapify( smallest, n );
    }
}


//------------------------------------------------------------------------------------------------------------
//! \brief  Sorts ordinates by octant.
//!
//! Upon return, ordinates are sorted so that the octants are grouped in the order:
//! I, V, IV, VIII, II, VI, III, VII.
//------------------------------------------------------------------------------------------------------------
void OrdinateSet::OrdinateData::Sort ( void ) const {

    // Build min heap.
    for ( int64_t i = num_nodes / 2 - 1; i >= 0; --i ) {  Heapify( i, num_nodes );  }

    // Sort from "largest" to "smallest" using heap property.
    int64_t n = num_nodes;

    for ( int64_t i = num_nodes - 1; i > 0; --i ) {

        Swap(0,i);
        Heapify( 0, --n );
    }
}


//------------------------------------------------------------------------------------------------------------
//! \brief  Allocates memory for and fills the array specifying mappings between ordinate indices at
//!         reflecting boundary conditions.
//!
//! \attention  When initializing a new quadrature set, this routine should be called _after_
//!             OrdinateSet::OrdinateData::SymmetricReduce().
//!
//! \attention  Implementation supports only two or more dimensions. Reflecting conditions for one-dimensional
//!             slab geometry are straightforward.
//!
//! \see    OrdinateSet::OrdinateData::Reflect()
//------------------------------------------------------------------------------------------------------------
void OrdinateSet::OrdinateData::ComputeReflectMap ( void ) {

    const double TOL = 1.0e-14;

    this->reflect_map = new int64_t[ this->num_nodes * SPACE_DIMS ];

    // Initialize with invalid value.
    for ( int64_t dim = 0; dim < SPACE_DIMS; ++dim ) {
    for ( int64_t q = 0; q < this->num_nodes; ++q ) {

        Reflect(dim,q) = -1;
    }}

    // Construct reflection map for each spatial dimension.
    for ( int64_t dim = 0; dim < SPACE_DIMS; ++dim ) {
    for ( int64_t q = 0; q < this->num_nodes; ++q ) {

        // Search for matching ordinate for reflection.
        for ( int64_t k = 0; k < this->num_nodes; ++k ) {

            // Check whether ordinate k satisfies the reflection conditions for ordinate q.
            bool matches = true;

            for ( int64_t d = 0; d < 3; ++d )
                matches &= std::abs( Coordinate(d,q) - ( d == dim ? -1 : +1 ) * Coordinate(d,k) ) < TOL;

            if ( matches ) {

                Reflect(dim,q) = k;

                goto next_ordinate;
            }
        }

        // If none found, then angular quadrature cannot be used with reflecting boundary conditions.
        {
            std::string error_message =   "Unable to find ordinate satisfying reflection condition in dimension "
                                        + std::to_string(dim)
                                        + " for ordinate with index "
                                        + std::to_string(q)
                                        + ".\n";

            throw std::runtime_error( error_message );
        }

        next_ordinate:;
    }}

# if defined STRICT_CHECK

    // Perform sanity checks.
    for ( int64_t dim = 0; dim < SPACE_DIMS; ++dim ) {
    for ( int64_t q = 0; q < this->num_nodes; ++q ) {

        if ( Reflect(dim,Reflect(dim,q)) != q ) {

            std::string error_message =   "Ordinate with index "
                                        + std::to_string(q)
                                        + " fails sanity check across dimension "
                                        + std::to_string(dim)
                                        + " for reflection map in "
                                        + std::string(__func__)
                                        + ".\n";

            throw std::runtime_error( error_message );
        }
    }}

# endif // if defined STRICT_CHECK
}


# if SPACE_DIMS == 2


//------------------------------------------------------------------------------------------------------------
//! \brief  Reduces the number of ordinates for 2D problems using symmetry in the \f$ z \f$ dimension.
//------------------------------------------------------------------------------------------------------------
void OrdinateSet::OrdinateData::SymmetricReduce ( void ) {

    for ( int64_t i = 0;   i < this->num_nodes; ++i ) {
    for ( int64_t j = i+1; j < this->num_nodes; ++j ) {

        // If matching ordinates found...
        if (    this->xi_ptr[i]  == this->xi_ptr[j]
             && this->eta_ptr[i] == this->eta_ptr[j]
             && this->mu_ptr[i]  == - this->mu_ptr[j]
        ) {
            // ...combine their quadrature weights...
            this->weights[i] += this->weights[j];

            // ...and remove the second.
            Remove(j);

            break;
        }
    }}
}


//------------------------------------------------------------------------------------------------------------
//! \brief  Removes the ordinate with specified index.
//!
//! The ordinate is removed by shifting all ordinates with indices greater than \pp{i} forward by one spot and
//! decrementing OrdinateData::num_nodes. Assumes that the ordinates were sorted before calling this function
//! and hence remain sorted once the specified ordinate has been removed. OrdinateData::LocateOctants is
//! called to update OrdinateData::octants after the specified ordinate has been removed. The arrays used to
//! store the ordinate values are not reallocated.
//!
//! \param[in]  i   Index of ordinate to remove.
//------------------------------------------------------------------------------------------------------------
void OrdinateSet::OrdinateData::Remove (

    const int64_t i
) {

    if ( i < 0 || i >= this->num_nodes ) {

        std::string error_message = "Index " + std::to_string(i)
                                    + " not in range of valid ordinate indices [0,"
                                    + std::to_string( this->num_nodes ) + ").\n";

        PRINT_ERROR( error_message.c_str() )
        throw std::out_of_range( error_message );
    }

    // Shift ordinates with indices greater than i forward by one.
    for ( int64_t j = i+1; j < this->num_nodes; ++j ) {

        this->xi_ptr [j-1] = this->xi_ptr [j];
        this->eta_ptr[j-1] = this->eta_ptr[j];
        this->mu_ptr [j-1] = this->mu_ptr [j];
        this->phi_ptr[j-1] = this->phi_ptr[j];
        this->weights[j-1] = this->weights[j];
    }

    // Decrement number of ordinates and set unused array elements to NAN.
    this->num_nodes--;

      this->xi_ptr [ this->num_nodes ]
    = this->eta_ptr[ this->num_nodes ]
    = this->mu_ptr [ this->num_nodes ]
    = this->phi_ptr[ this->num_nodes ]
    = this->weights[ this->num_nodes ]
    = std::nan("0");

    // Update locations of octants.
    LocateOctants();
}


# endif // if SPACE_DIMS == 2
# endif // if SPACE_DIMS >= 2
