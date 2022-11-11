//------------------------------------------------------------------------------------------------------------
//! \file   utils/Quadrule/Quadrule.hpp
//!
//! \author Michael Crockatt
//! \date   August 2017
//------------------------------------------------------------------------------------------------------------

# ifndef __QUADRULE_HPP__
# define __QUADRULE_HPP__


# include <cstdint>
# include <map>
# include <string>


# ifndef M_PI
    # define M_PI 3.141592653589793238462643383279502884197169399375105820974944592307816406286208998628034825
# endif


namespace Quadrule {


//------------------------------------------------------------------------------------------------------------
//! \brief  Specifies type of quadrature nodes for finite 1D intervals.
//!
//! \see    Quadrule::NodesType_to_String
//------------------------------------------------------------------------------------------------------------
enum class NodesType {

    None,               //!< Null value.
    GaussLegendre,      //!< Gauss-Legendre quadrature.
    GaussRadau,         //!< Gauss-Radau quadrature.
    GaussLobatto,       //!< Gauss-Lobatto quadrature.
    Chebyshev,          //!< Chebyshev quadrature.
    ChebyshevRadau,     //!< Chebyshev-Radau quadrature.
    EquispacedBoth,     //!< Equispaced quadrature including both left and right endpoints.
    EquispacedRight     //!< Equispaced quadrature including the right endpoint but not the left.
};


//------------------------------------------------------------------------------------------------------------
//! \brief  Used to convert Quadrule::NodesType values to their corresponding string descriptors.
//------------------------------------------------------------------------------------------------------------
const std::map< NodesType, std::string > NodesType_to_String = {

    { NodesType::None,              "none"              },
    { NodesType::GaussLegendre,     "gauss-legendre"    },
    { NodesType::GaussRadau,        "gauss-radau"       },
    { NodesType::GaussLobatto,      "gauss-lobatto"     },
    { NodesType::Chebyshev,         "chebyshev"         },
    { NodesType::ChebyshevRadau,    "chebyshev-radau"   },
    { NodesType::EquispacedBoth,    "equispaced-both"   },
    { NodesType::EquispacedRight,   "equispaced-right"  }
};


//------------------------------------------------------------------------------------------------------------
//! \brief  Used to convert string descriptors to Quadrule::NodesType values.
//------------------------------------------------------------------------------------------------------------
const std::map< std::string, NodesType > String_to_NodesType = {

    { "none",               NodesType::None             },
    { "gauss-legendre",     NodesType::GaussLegendre    },
    { "gauss-radau",        NodesType::GaussRadau       },
    { "gauss-lobatto",      NodesType::GaussLobatto     },
    { "chebyshev",          NodesType::Chebyshev        },
    { "chebyshev-radau",    NodesType::ChebyshevRadau   },
    { "equispaced-both",    NodesType::EquispacedBoth   },
    { "equispaced-right",   NodesType::EquispacedRight  }
};


int64_t gcd( int64_t a, int64_t b );
int64_t lcm( int64_t a, int64_t b );

bool IsPowerOfTwo( int64_t x );
bool IsPowerOfTwo( uint64_t x );

double Jacobi( const int n, const double x, const double alpha, const double beta );
double JacobiPrime( const int n, const double x, const double alpha, const double beta );

double Legendre( const int n, const double x );
double LegendrePrime( const int n, const double x );
double Radau( const int n, const double x );
double RadauPrime( const int n, const double x );

void testQuad( const double * const nodes, const double * const weights, const int num );

void ComputeQuadrature( const int64_t n, double * const nodes, double * const weights, const NodesType type,
                        const double a = -1.0, const double b = 1.0 );

void ComputeQuadratureNodes( const int64_t n, double * const nodes, const NodesType type,
                             const double a = -1.0, const double b = 1.0 );

void ComputeQuadratureWeights( const int64_t n, const double * const nodes, double * const weights,
                               const double a = -1.0, const double b = 1.0 );

void ComputeGKquad( const int n, double * const nodes, double * const G_weights, double * const K_weights );

double AssociatedLegendre( const int l, const int m, const double x );
double SphericalHarmonic( const int n, const int l, const double z, const double phi );


} // namespace Quadrule


# endif // ifndef __QUADRULE_HPP__
