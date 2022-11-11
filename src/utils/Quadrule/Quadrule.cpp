//------------------------------------------------------------------------------------------------------------
//! \file   utils/Quadrule/Quadrule.cpp
//! \brief  Implements routines to construct quadrature rules from orthogonal polynomials.
//!
//! \author Michael Crockatt
//! \date   August 2017
//------------------------------------------------------------------------------------------------------------


# include <cmath>
# include <cstdint>
# include <cstdlib>
# include <cstring>
# include <stdexcept>
# include <limits>
# include <set>
# include <string>
# include <vector>

# include "utils/CLog.hpp"
# include "utils/Quadrule/Quadrule.hpp"
# include "utils/Quadrule/Tabulated.hpp"


using namespace Quadrule;


# if defined (__cplusplus)
    extern "C"
# else
    extern
# endif
int dgesv_ ( int *, int *, double *, int *, int *, double *, int *, int * );


//============================================================================================================
//=== NUMBER THEORETIC INTEGER ALGORITHMS ====================================================================
//============================================================================================================


//------------------------------------------------------------------------------------------------------------
//! \brief  Routine to compute the greatest common divisor of two positive integers.
//!
//! The GCD is computed using
//! <a href="https://en.wikipedia.org/wiki/Binary_GCD_algorithm">Stein's algorithm</a>.
//!
//! \param[in]      a           First integer.
//! \param[in]      b           Second integer.
//!
//! \return     Returns the greatest common divisor of \f$ a \f$ and \f$ b \f$.
//------------------------------------------------------------------------------------------------------------
int64_t Quadrule::gcd (

    int64_t a,
    int64_t b
) {

    // Check to make sure inputs are nonnegative.
    if ( a < 0 || b < 0 ) {

        std::string error_message = "Arguments to function '" + std::string(__func__)
                                    + "' must be nonnegative ("
                                    + std::to_string(a) + " and " + std::to_string(b)
                                    + " provided).\n";

        PRINT_ERROR( error_message.c_str() )
        throw std::invalid_argument( error_message );
    }

    int64_t s,t;

    if ( a == 0 ) {  return b;  }
    if ( b == 0 ) {  return a;  }

    // Extract all factors of 2 common to a and b.
    for ( s = 0; !((a | b) & 1); ++s ) {

        a >>= 1;
        b >>= 1;
    }

    // Extract all additional factors of 2 from a. After this a is odd.
    while ( !(a & 1) ) {  a >>= 1;  }

    do {

        while ( !(b & 1) ) {  b >>= 1;  }

        if ( a > b ) {  t = b; b = a; a = t;  } // Swap a and b.

        b -= a;

    } while ( b );

    return a << s;
}

//------------------------------------------------------------------------------------------------------------
//! \brief  Routine to compute the least common multiple of two positive integers.
//!
//! \param[in]      a           First integer.
//! \param[in]      b           Second integer.
//!
//! \return     Returns the least common multiple of \f$ a \f$ and \f$ b \f$.
//------------------------------------------------------------------------------------------------------------
int64_t Quadrule::lcm (

    int64_t a,
    int64_t b
) {
    return a * b / gcd( a, b );
}


//------------------------------------------------------------------------------------------------------------
//! \brief  Determines whether or not the given integer is a power of two.
//!
//! An exception is thrown if a negative integer is provided.
//!
//! \param[in]  x   Integer value to examine.
//!
//! \return     Returns true if \pp{x} is a power of two, and false otherwise.
//------------------------------------------------------------------------------------------------------------
bool Quadrule::IsPowerOfTwo (

    int64_t x
) {

    if ( x < 0 ) {

        std::string error_message = "Argument to function '" + std::string(__func__)
                                    + "' must be nonnegative (" + std::to_string(x) + " provided).\n";

        PRINT_ERROR( error_message.c_str() )
        throw std::invalid_argument( error_message );
    }

    return IsPowerOfTwo( (uint64_t) x );
}


//------------------------------------------------------------------------------------------------------------
//! \brief  Determines whether or not the given integer is a power of two.
//!
//! \param[in]  x   Integer value to examine.
//!
//! \return     Returns true if \pp{x} is a power of two, and false otherwise.
//------------------------------------------------------------------------------------------------------------
bool Quadrule::IsPowerOfTwo (

    uint64_t x
) {
    return ((x != 0) && !(x & (x - 1)));
}


//============================================================================================================
//=== ROUTINES FOR EVALUATING ORTHOGONAL POLYNOMIALS =========================================================
//============================================================================================================


//------------------------------------------------------------------------------------------------------------
//! \brief  Routine to evaluate
//!         <a href="https://en.wikipedia.org/wiki/Jacobi_polynomials">
//!             Jacobi polynomials
//!         </a> using a recurrence relation.
//!
//! Evaluates the \f$ (\alpha,\beta) \f$-Jacobi polynomial of degree \f$ n \f$, where the
//! \f$ (\alpha,\beta) \f$-Jacobi polynomials are the polynomials orthogonal over \f$ [-1,1] \f$ with respect
//! to the weight function \f[ (1-x)^{\alpha} (1+x)^{\beta} \f]
//! The recurrence relation used here can be found on page 74 of \cite Shen2011.
//!
//! \param[in]      n           Evaluate the Jacobi polynomial of degree \c n.
//! \param[in]      x           Point to evaluate the Jacobi polynomial at.
//! \param[in]      alpha       Alpha value specifying class of Jacobi polynomials.
//! \param[in]      beta        Beta value specifying class of Jacobi polynomials.
//!
//! \return     Returns the value of the \f$ (\alpha,\beta) \f$-Jacobi polynomial of degree \f$ n \f$ at the
//!             point \f$ x \f$.
//!
//! \see    JacobiPrime()
//------------------------------------------------------------------------------------------------------------
double Quadrule::Jacobi (

    const int n,
    const double x,
    const double alpha,
    const double beta
) {

    // Base case.
    if ( n == 0 ) {

        return 1.0;

    // Apply recurrence relation.
    } else {

        double a,b,c;
        double pn  = 0.5 * ( ( alpha + beta + 2.0 ) * x + alpha - beta );
        double pn1 = 1.0;
        double pn2;

        for ( int j = 1; j < n; ++j ) {

            pn2 = pn1;
            pn1 = pn;

            a = (2*j + alpha + beta + 1) * (2*j + alpha + beta + 2)
                    / ( (2*j + 2) * (j + alpha + beta + 1) );

            b = (beta*beta - alpha*alpha) * (2*j + alpha + beta + 1)
                    / ( (2*j + 2) * (j + alpha + beta + 1) * (2*j + alpha + beta) );

            c = (j + alpha) * (j + beta) * (2*j + alpha + beta + 2)
                    / ( (j + 1) * (j + alpha + beta + 1) * (2*j + alpha + beta) );

            pn = (a*x - b) * pn1 - c * pn2;
        }

        return pn;
    }
}


//------------------------------------------------------------------------------------------------------------
//! \brief  Routine to evaluate derivatives of Jacobi polynomials.
//!
//! Uses the formula
//! \f[ \frac{d}{dx} P_n^{(\alpha,\beta)}(x) = \frac{ \Gamma(\alpha+\beta+n+2) }{ 2\Gamma(\alpha+\beta+n+2) } P_{n-1}^{(\alpha+1,\beta+1)} (x). \f]
//!
//! \param[in]      n           Evaluate the derivative of the \f$n\f$th Jacobi polynomial.
//! \param[in]      x           Point to evaluate derivative at.
//! \param[in]      alpha       Alpha value specifying class of Jacobi polynomials.
//! \param[in]      beta        Beta value specifying class of Jacobi polynomials.
//!
//! \return     Returns the value of the first derivative of the \f$ (\alpha,\beta) \f$-Jacobi polynomial of
//!             degree \f$ n \f$ at the point \f$ x \f$.
//!
//! \see    Jacobi()
//------------------------------------------------------------------------------------------------------------
double Quadrule::JacobiPrime (

    const int n,
    const double x,
    const double alpha,
    const double beta
) {

    return std::tgamma( alpha + beta + n + 2.0 ) * Jacobi( n-1, x, alpha+1, beta+1 )
            / ( 2.0 * tgamma( alpha + beta + n + 1 ) );
}


//------------------------------------------------------------------------------------------------------------
//! \brief  Routine to evaluate
//!         <a href="https://en.wikipedia.org/wiki/Legendre_polynomials">
//!             Legendre polynomials
//!         </a>.
//!
//! NOTE: The Legendre polynomials are equivalent to the \f$ (0,0) \f$-Jacobi polynomials.
//!
//! Performs the evaluation using Bonnet's recursion formula:
//! \f[ (n+1) P_{n+1}(x) = (2n+1) x P_n(x) - n P_{n-1} (x). \f]
//!
//! \param[in]      n           Degree of Legendre polynomial to evaluate.
//! \param[in]      x           Point to evaluate the Legendre polynomial at.
//!
//! \return     Returns the value of the \f$ n \f$th Legendre polynomial at the point \f$ x \f$.
//!
//! \see    LegendrePrime()
//! \see    Jacobi()
//------------------------------------------------------------------------------------------------------------
double Quadrule::Legendre (

    const int n,
    const double x
) {

    // --- Compute as Jacobi polynomial. -------------------------------------------------------------- //

//     return Jacobi(n,x,0,0);

    // --- Compute by recurrence relation. ------------------------------------------------------------ //

    // Base case.
    if ( n == 0 ) {

        return 1.0;

    // Apply recurrence relation.
    } else {

        double pn  = x;
        double pn1 = 1.0;
        double pn2;

        for ( int j = 1; j < n; ++j ) {

            pn2 = pn1;
            pn1 = pn;
            pn  = ( (2*j + 1) * x * pn1 - j * pn2 ) / (j + 1);
        }

        return pn;
    }
}


//------------------------------------------------------------------------------------------------------------
//! \brief  Routine to evaluate derivatives of Legendre polynomials.
//!
//! \attention  Some methods used by this function are not valid for \f$ x = ±1 \f$. If an invalid situation
//!             is encountered with these methods, an error message is printed to \c stderr and a \c NAN value
//!             is returned.
//!
//! The default method used by this function is to apply the first derivative of Bonnet's recursion formula:
//! \f[ (n+1) P'_{n+1} (x) = (2n+1) \left[ P_n (x) + x P'_n (x) \right] - n P'_{n-1} (x), \f]
//! which _is_ valid at the points \f$ x = \pm 1 \f$.
//!
//! \param[in]      n           Evaluate the first derivative of the Legendre polynomial of degree \f$ n \f$.
//! \param[in]      x           Point to evaluate derivative at.
//!
//! \return     Returns the value of the derivative of the \f$ n \f$th Legendre polynomial at the point
//!             \f$ x \f$.
//!
//! \see    Legendre()
//------------------------------------------------------------------------------------------------------------
double Quadrule::LegendrePrime (

    const int n,
    const double x
) {

    // --- Compute as Jacobi polynomial. -------------------------------------------------------------- //

//     return JacobiPrime(n,x,0,0);

    // --- Compute by formula. ------------------------------------------------------------------------ //

//     if ( fabs(x) == 1.0 )
//         fprintf( stderr,
//                  __RED "ERROR: " __RESET "[%s:%d]  "
//                  "Division by zero in %s.\n",
//                  __FILE__, __LINE__, __func__ );
//
//     return  n * ( Legendre( n-1, x ) - x * Legendre( n, x )  ) / (1.0 - x*x);

    // --- Compute by derivative of recurrence relation. ---------------------------------------------- //

    // Base case.
    if ( n == 0 ) {

        return 0.0;

    // Apply derivative of recurrence relation.
    } else {

        double pn  = 1.0;
        double pn1 = 0.0;
        double pn2;

        for ( int j = 1; j < n; ++j ) {

            pn2 = pn1;
            pn1 = pn;
            pn  = ( (2*j + 1) * ( Legendre(j,x) + x * pn1 ) - j * pn2 ) / (j + 1);
        }

        return pn;
    }
}


//------------------------------------------------------------------------------------------------------------
//! \brief  Routine to evaluate second derivatives of Legendre polynomials.
//!
//! The method used here is tailored to the computation of Gauss-Lobatto quadrature as described on page 99
//! of \cite Shen2011.
//!
//! \param[in]      n           Evaluate the first derivative of the Legendre polynomial of degree \f$ n \f$.
//! \param[in]      x           Point to evaluate derivative at.
//!
//! \return     Returns the value of the second derivative of the \f$ n \f$th Legendre polynomial at the point
//!             \f$ x \f$.
//------------------------------------------------------------------------------------------------------------
double LegendreDoublePrime (

    const int n,
    const double x
) {

    return ( 2.0 * x * LegendrePrime(n,x) - n*(n + 1) * Legendre(n,x) ) / ( 1.0 - x*x );
}


//------------------------------------------------------------------------------------------------------------
//! \brief  Routine to evaluate Radau-right (IIA) polynomials.
//!
//! NOTE: The Radau-right polynomials are equivalent to the \f$ (-1,0) \f$-Jacobi polynomials.
//!
//! \param[in]      n           Evaluate the Radau polynomial of degree \f$ n \f$.
//! \param[in]      x           Point to evaluate polynomial at.
//!
//! \return     Returns the value of the \f$ n \f$th Radau polynomial at the point \f$ x \f$.
//!
//! \see    RadauPrime()
//! \see    Jacobi()
//!
//! \todo   Description of recurrence relation here.
//------------------------------------------------------------------------------------------------------------
double Quadrule::Radau (

    const int n,
    const double x
) {

    // --- Compute as Jacobi polynomial. -------------------------------------------------------------- //

//     return Jacobi(n,x,-1,0);

    // --- Compute by relationship to Legendre polynomials. ------------------------------------------- //

//     return Legendre( n, x ) - Legendre( n-1, x );

    // --- Compute by recurrence relation. ------------------------------------------------------------ //

    // Base case.
    if ( n == 0 ) {

        return 1.0;

    } else {

        double pn  = 0.5 * (x - 1.0);
        double pn1 = 1.0;
        double pn2;

        for ( int j = 1 ; j < n ; ++j ) {

            pn2 = pn1;
            pn1 = pn;
            pn  = ( ( (2*j + 1) * (2*j - 1) * x + 1.0 ) * pn1
                    - (j - 1) * (2*j + 1) * pn2 ) / ( (j + 1) * (2*j - 1) );
        }

        return pn;
    }
}


//------------------------------------------------------------------------------------------------------------
//! \brief  Routine to evaluate derivatives of Radau-right (IIA) polynomials.
//!
//! \attention  Depending on the method used in this function and also the LegendrePrime() function, the
//!             associated formulas may not be valid for \f$ x = \pm 1 \f$. If an invalid situation is
//!             encountered with these methods, an error message is printed to \c stderr and a \c NAN value is
//!             returned.
//!
//! The default method used by this function is to apply the derivative of the recurrence formula for Radau
//! polynomials:
//! \todo   Insert recurrence relation for Radau derivatives here.
//!
//! which _is_ valid at the points \f$ x = \pm 1 \f$.
//!
//! \param[in]      n           Evaluate the derivative of the Radau polynomial of degree \f$ n \f$.
//! \param[in]      x           Point to evaluate the derivative at.
//!
//! \return     Returns the value of the derivative of the \f$ n \f$th Radau polynomial at the point
//!             \f$ x \f$.
//!
//! \see    Radau()
//------------------------------------------------------------------------------------------------------------
double Quadrule::RadauPrime (

    const int n,
    const double x
) {

    // --- Compute as Jacobi polynomial. -------------------------------------------------------------- //

//     return JacobiPrime(n,x,-1,0);

    // --- Compute by relationship to Legendre polynomials. ------------------------------------------- //

//     return LegendrePrime( n, x ) - LegendrePrime( n-1, x );

    // --- Compute by derivative of recurrence relation. ---------------------------------------------- //

    // Base case.
    if ( n == 0 ) {

        return 0.0;

    } else {

        double pn  = 0.5;
        double pn1 = 0.0;
        double pn2;

        for ( int j = 1 ; j < n ; ++j ) {

            pn2 = pn1;
            pn1 = pn;
            pn  = ( (2*j + 1) * (2*j - 1) * Radau(j,x)
                    + ( (2*j + 1) * (2*j - 1) * x + 1.0 ) * pn1
                    - (j - 1) * (2*j + 1) * pn2 ) / ( (j + 1) * (2*j - 1) );
        }

        return pn;
    }
}


//------------------------------------------------------------------------------------------------------------
//! \brief  Wrapper function for computing interior nodes of Gauss-Lobatto quadrature.
//!
//! \param[in]      n           Degree of polynomial to evaluate.
//! \param[in]      x           Point to evaluate polynomial at.
//!
//! \return     Returns the value of the first derivative of the \f$ (n+1) \f$st Legendre polynomial at the
//!             point \f$ x \f$.
//!
//! \see    LobattoPrime()
//------------------------------------------------------------------------------------------------------------
double Lobatto (

    const int n,
    const double x
) {

    return LegendrePrime( n+1, x );
}


//------------------------------------------------------------------------------------------------------------
//! \brief  Wrapper function for computing interior nodes of Gauss-Lobatto quadrature.
//!
//! \param[in]      n           Degree of polynomial to evaluate derivative of.
//! \param[in]      x           Point to evaluate the derivative at.
//!
//! \return     Returns the value of the second derivative of the \f$ (n+1) \f$st Legendre polynomial at the
//!             point \f$ x \f$.
//!
//! \see    Lobatto()
//------------------------------------------------------------------------------------------------------------
double LobattoPrime (

    const int n,
    const double x
) {

    return LegendreDoublePrime( n+1, x );
}


//============================================================================================================
//=== ROUTINES FOR EVALUATING SPHERICAL HARMONICS ON 𝕊²  ======================================================
//============================================================================================================


//------------------------------------------------------------------------------------------------------------
//! \brief  Routine to evaluate
//!         <a href="https://en.wikipedia.org/wiki/Associated_Legendre_polynomials">
//!             associated Legendre functions
//!         </a>.
//!
//! Uses the recurrence relation defined on page 247 of \cite Press1992
//! \f[ (\ell - m) P^m_{\ell} = x (2\ell - 1) P^m_{\ell - 1} + (\ell + m - 1) P^m)_{\ell-2} \f]
//! with starting values
//! \f[ P^m_m = (-1)^m (2m - 1)!! (1-x^2)^{m/2} \f]
//! \f[ P^m_{m+1} = x (2m + 1) P^m_m \f]
//! where \f$ n!! \f$ denotes the product of all odd integers less than or equal to \f$ n \f$.
//!
//! \param[in]      l           Degree of associated Legendre function.
//! \param[in]      m           Order of associated Legendre function.
//! \param[in]      x           Point to evaluate Legendre function at.
//!
//! \return     Returns the value of the associated Legendre function \f$ P_{\ell}^m \f$ at the point
//!             \f$ x \f$.
//------------------------------------------------------------------------------------------------------------
double Quadrule::AssociatedLegendre (

    const int l,
    const int m,
    const double x
) {

    // Test to make sure -l <= m <= l.
    if ( std::abs(m) > l || l < 0 ) {

        std::string error_message = "Values l = " + std::to_string(l) + " and m = " + std::to_string(m)
                                    + " are invalid input for function '" + std::string(__func__) + "'.\n";

        PRINT_ERROR( error_message.c_str() )
        throw std::invalid_argument( error_message );
    }

    // --- If m is zero the associated Legendre function is the Legendre polynomial of degree l. ------ //

    if ( m == 0 ) {  return Legendre( l, x );  }

    // --- If m is negative, use proportionality relation and re-call function with positive m. ------- //

    if ( m < 0 ) {

        double factor = ( m % 2 ? -1.0 : 1.0 );

        for ( int k = l+m+1; k <= l-m; ++k ) {  factor /= k;  }

        return factor * AssociatedLegendre( l, std::abs(m), x );
    }

    // --- Otherwise use the recurrence relation on l. ------------------------------------------------ //

    // Set base case P_m^m.
    double pl = ( m % 2 ? -1.0 : 1.0 );

    for ( int k = 2*m - 1; k > 1; k -= 2 ) {  pl *= k;  }

    pl *= pow( 1.0 - x*x, m / 2.0 );

    // Return base case if appropriate.
    if ( l == m ) {  return pl;  }

    // Otherwise, use recurrence to compute P_l^m of point.
    double pl1 = 0;
    double pl2;

    for ( int j = m+1; j <= l; ++j ) {

        pl2 = pl1;
        pl1 = pl;
        pl  = ( x * (2.0*j - 1.0) * pl1 - (j + m - 1.0) * pl2 ) / (j - m);
    }

    return pl;
}


//------------------------------------------------------------------------------------------------------------
//! \brief  Routine to evaluate normalized spherical harmonics.
//!
//! The spherical harmonics are indexed as given by equation (4.7) on page 133 of \cite Atkinson2012
//! \f{align*}{ Y_{n,1} (\phi,\theta) & = c_n P_n (\cos\theta) && n \geq 0 \\ Y_{n,2\ell} (\phi,\theta) & = c_{n,\ell} P_n^{\ell} (\cos\theta) \cos(\ell\phi) && n \geq 1; \; \ell = 1, \ldots, n \\ Y_{n,2\ell+1} (\phi,\theta) & = c_{n,\ell} P_n^{\ell} (\cos\theta) \sin (\ell\phi) && n \geq 1; \; \ell = 1, \ldots, n \f}
//! with normalization constants
//! \f[ c_n = \sqrt{\frac{2n+1}{4\pi}}, \qquad c_{n,\ell} = \sqrt{\frac{2n+1}{2\pi}\frac{(n-\ell)!}{(n+\ell)!}}. \f]
//!
//! The harmonics are indexed by \f$ (n,\ell) \f$ where \f$ n \f$ is the degree of the harmonic and
//! \f$ m = \ell /2 \f$ is the order of the harmonic satisfying \f$ 1 <= \ell <= 2n+1 \f$. This function
//! evaluates the spherical harmonics in terms of the variables \f$ (z,\phi) \f$ where
//! \f$ \phi \in [0,2\pi] \f$ is the azimuthal angle (around the latitudes of the sphere) and
//! \f$ z = \cos(θ) \f$ where \f$ θ \in [0,\pi] \f$ is the polar angle (from the top of the sphere).
//!
//! \param[in]      n           Degree of spherical harmonic to evaluate.
//! \param[in]      l           Index for order of spherical harmonic to evaluate.
//! \param[in]      z           Cosine of the polar angle \f$ \theta \f$ in spherical coordinates to evaluate
//!                             spherical harmonic at.
//! \param[in]      phi         Azimuthal angle in spherical coordinates to evaluate spherical harmonic at.
//!
//! \return     Returns the value obtained from evaluating the specified normalized spherical harmonic at the
//!             point \f$ ( z, \phi ) \f$.
//!
//! \see    Legendre()
//! \see    AssociatedLegendre()
//------------------------------------------------------------------------------------------------------------
double Quadrule::SphericalHarmonic (

    const int n,
    const int l,
    const double z,
    const double phi
) {

    double factor;

    const int m   = l / 2;
    const int par = l % 2;

    // --- Error checking on l and n. ----------------------------------------------------------------- //

    if (    ( n < 0 )
         || ( l < 1 )
         || ( m > n )
    ) {
        std::string error_message = "Values n = " + std::to_string(n) + " and l = " + std::to_string(l)
                                    + " are invalid input for function '" + std::string(__func__) + "'.\n";

        PRINT_ERROR( error_message.c_str() )
        throw std::invalid_argument( error_message );
     }

    // --- Zeroth order harmonic. --------------------------------------------------------------------- //

    if ( m == 0 ) {

        factor = std::sqrt( (2.0*n + 1.0) / (4.0 * M_PI) );

        return Legendre( n, z ) * factor;
    }

    // --- Parity of l determines whether azimuthal component uses sine or cosine. -------------------- //

    factor = (2.0*n + 1.0) / (2.0 * M_PI);

    for ( int k = n-m+1; k <= n+m; ++k ) { factor /= k; }

    factor = std::sqrt( factor );

    // Odd parity, use sine function.
    if ( par ) {  return factor * AssociatedLegendre( n, m, z ) * std::sin( m * phi );  }

    // Even parity, use cosine function.
    else {  return factor * AssociatedLegendre( n, m, z ) * std::cos( m * phi );  }
}


//============================================================================================================
//=== MISCELLANEOUS HELPER ROUTINES ==========================================================================
//============================================================================================================


//------------------------------------------------------------------------------------------------------------
//! \brief  Performs accuracy tests on the given quadrature rule.
//!
//! Integrates products of Legendre polynomials, printing results to \c stdout.
//!
//! \param[in]      nodes       Pointer to array of quadrature nodes.
//! \param[in]      weights     Pointer to array of quadrature weights.
//! \param[in]      num         Number of nodes in quadrature set.
//------------------------------------------------------------------------------------------------------------
void Quadrule::testQuad (

    const double * const nodes,
    const double * const weights,
    const int num
) {

    double quadValue;
    char buffer[100];

    int colWidth = 14;
    double tol   = 1e-15;

    char line[100] = {};
    for ( int k = 0; k < colWidth; ++k ) {  line[k] = '-';  }
    line[colWidth] = '\0';

    // Print table header.
    printf( "\n %-*s ", colWidth, " " );

    for ( int n = 0 ; n < num ; ++n ) {

        sprintf( buffer, "  L%d", n );
        printf( " %-*s ", colWidth, buffer );
    }
    printf( "\n %-*s ", colWidth, " " );

    for ( int n = 0 ; n < num ; ++n )
        printf( " %-*s ", colWidth, line );

    // Compute tests.
    for ( int m = 0 ; m < num ; ++m ) {

        // Print row label.
        sprintf( buffer, "  L%d", m );
        printf( "\n %-*s:", colWidth, buffer );

        for ( int n = 0 ; n < num ; ++n ) {

            // Compute test.
            quadValue = 0.0;

            for ( int k = 0 ; k < num ; ++k )
                quadValue += weights[k] * Legendre( m, nodes[k] ) * Legendre( n, nodes[k] );

            // Print result of test.
            if ( fabs( quadValue ) > tol )
                printf( __RED " % *.2e   " __RESET, colWidth-2, quadValue );
            else
                printf( " % *.2e   ", colWidth-2, quadValue );
        }
    }

    printf( "\n" );
}


//============================================================================================================
//=== ROUTINES FOR COMPUTING QUADRATURES =====================================================================
//============================================================================================================


//!
//! \brief Tolerance for floating point equality.
//!
//! Two floating point values a and b are determined to be equal if std::abs( b-a ) <= tol;
//!
static const double s_double_tol = 5e-16;


//------------------------------------------------------------------------------------------------------------
//! \brief  Functor used to compare double-precision floating-point values for near-equality.
//!
//! Checks if the first argument is more than #s_double_tol less than the second.
//------------------------------------------------------------------------------------------------------------
struct CompareDoubles {

    bool operator() (

        const double & lhs,
        const double & rhs
    ) {
        return lhs < (rhs - s_double_tol);
    }
};


//------------------------------------------------------------------------------------------------------------
//! \brief  Computes \c n distinct roots of a function via Newton iterations.
//!
//! Uses Newton iteration to find \c n distinct roots of the function \c p with derivative \c pPrime. This
//! function is used to compute quadrature nodes as roots of families of orthogonal polynomials. The resulting
//! set of quadrature nodes can be passed to computeIntervalQuadWeights() or computeStandardQuadWeights() to
//! obtain optimal quadrature weights for the given node set.
//!
//! The function \c p and its derivative \c pPrime are passed as function pointers, and take two arguments.
//! The first argument is the degree of polynomial (which here is \c n) and the second the point to evaluate
//! the polynomial at.
//!
//! \attention  This function assumes that the memory address pointed to by \c roots has sufficient memory
//!             allocated to it to hold all \c n roots which are found.
//!
//! \param[in]      n           Number of roots desired.
//! \param[in]      p           Pointer to function which evaluates the function to find roots of.
//! \param[in]      pPrime      Pointer to function to evaluate the derivative of the function to find roots
//!                             of.
//! \param[out]     roots       Pointer to array in which to store the resulting roots. The values are sorted
//!                             into increasing order before return.
//!
//! \see    computeIntervalQuadWeights()
//! \see    computeStandardQuadWeights()
//! \see    Legendre()
//! \see    LegendrePrime()
//! \see    Radau()
//! \see    RadauPrime()
//------------------------------------------------------------------------------------------------------------
static void ComputeRoots (

    const int n,
    double (*p)( const int, const double ),
    double (*pPrime)( const int, const double ),
    double * const roots
) {

    if ( n <= 0 )
        throw std::range_error( "Number of roots must be nonnegative in" + std::string(__func__) + "." );

    int64_t N = n;  // Number of test points to iterate on.

    // Contains set of (unique) roots that have been found.
    std::set< double, CompareDoubles > roots_found;

    while ( roots_found.size() < (uint64_t) n ) {

        N *= 2;

        std::vector<double> points(N);
        std::vector<double> errors(N);

        // Uniformly spaced initial points on [-1,1] for Newton iteration.
        double dp = 2.0 / N;

        for ( int64_t i = 0; i < N; ++i ) {

            points[i] = -1.0 + dp * (i + 0.5);
            errors[i] = std::numeric_limits<double>::max();
        }

        // Iterate on each initial point until convergence.
        # pragma omp parallel for schedule(dynamic,1)
        for ( int64_t i = 0; i < N; ++i ) {

            double temp;

            while ( errors[i] > s_double_tol ) {

                temp = points[i] - p( n, points[i] ) / pPrime( n, points[i] );

                errors[i] = std::abs( points[i] - temp );
                points[i] = temp;
            }
        }

        // Add new points to set of found roots.
        for ( auto root : points )
            roots_found.insert( root );
    }

    // Return found roots at pointer.
    int64_t i = 0;

    for ( auto root : roots_found ) {

        roots[i] = root;
        ++i;
    }
}


//------------------------------------------------------------------------------------------------------------
//! \brief  Wrapper function to compute the Gauss-Radau IIA quadrature nodes on the interval \f$ [-1,1] \f$.
//!
//! \attention  This function assumes that the memory address pointed to by \c nodes has sufficient memory
//!             allocated to it to hold all \f$ n \f$ quadrature weights.
//!
//! \param[in]      n           Size of quadrature rule desired.
//! \param[in]      nodes       Upon return contains the \f$ n \f$th set of Gauss-Radau IIA quadrature nodes.
//!
//! \see    ComputeRoots()
//! \see    Radau()
//! \see    RadauPrime()
//------------------------------------------------------------------------------------------------------------
void ComputeRadauNodes (

    const int n,
    double * const nodes
) {
    ComputeRoots( n, Radau, RadauPrime, nodes );
}


//------------------------------------------------------------------------------------------------------------
//! \brief  Wrapper function to compute the Gauss-Lobatto quadrature nodes on the interval \f$ [-1,1] \f$.
//!
//! \attention  This function assumes that the memory address pointed to by \c nodes has sufficient memory
//!             allocated to it to hold all \f$ n \f$ quadrature weights.
//!
//! \param[in]      n           Size of quadrature rule desired.
//! \param[in]      nodes       Upon return contains the \f$ n \f$th set of Gauss-Lobatto quadrature nodes.
//!
//! \see    ComputeRoots()
//! \see    Lobatto()
//! \see    LobattoPrime()
//------------------------------------------------------------------------------------------------------------
void ComputeLobattoNodes (

    const int n,
    double * const nodes
) {

    nodes[0]   = -1.0;
    nodes[n-1] =  1.0;

    if ( n > 2 ) {  ComputeRoots( n-2, Lobatto, LobattoPrime, nodes+1 );  }
}


//------------------------------------------------------------------------------------------------------------
//! \brief  Wrapper function to compute the Gauss-Legendre quadrature nodes on the interval \f$ [-1,1] \f$.
//!
//! \attention  This function assumes that the memory address pointed to by \c nodes has sufficient memory
//!             allocated to it to hold all \f$ n \f$ quadrature weights.
//!
//! \todo   Precomputed high-resolution nodes.
//!
//! \param[in]      n           Size of quadrature rule desired.
//! \param[in]      nodes       Upon return contains the \f$ n \f$th set of Gauss-Legendre quadrature nodes.
//!
//! \see    ComputeRoots()
//! \see    Quadrule::Legendre()
//! \see    Quadrule::LegendrePrime()
//------------------------------------------------------------------------------------------------------------
void ComputeGLnodes (

    const int n,
    double * const nodes
) {
    try         {  GLTable( n, nodes );                                 }
    catch (...) {  ComputeRoots( n, Legendre, LegendrePrime, nodes );   }
}


//------------------------------------------------------------------------------------------------------------
//! \brief  Computes quadrature nodes and weights of the specified type.
//!
//! \param[in]      n           Number of quadrature nodes of the specified type.
//! \param[out]     nodes       Pointer at which to store the computed nodes. <br/>
//!                             It is assumed that the necessary memory has been properly allocated at this
//!                             pointer prior to any calls to this function.
//! \param[in]      weights     Upon return contains the optimal quadrature weights for integrating across
//!                             \f$ [a,b] \f$.
//! \param[in]      type        Quadrule::NodesType specifying the type of quadrature nodes to compute.
//! \param[in]      a           (optional) <br/>
//!                             Left endpoint of interval to scale quadrature nodes to.
//!                             Default value is -1.0.
//! \param[in]      b           (optional) <br/>
//!                             Right endpoint of interval to scale quadrature nodes to.
//!                             Default value is 1.0.
//!
//! \see    Quadrule::ComputeQuadratureNodes()
//! \see    Quadrule::ComputeQuadratureWeights()
//------------------------------------------------------------------------------------------------------------
void Quadrule::ComputeQuadrature (

    const int64_t n,
    double * const nodes,
    double * const weights,
    const NodesType type,
    const double a,         // = -1.0
    const double b          // = 1.0
) {

    switch ( type ) {

        // Override for tabulated Gauss-Legendre quadratures.
        case NodesType::GaussLegendre:
        {
            // Try to obtain tabulated values.
            try {

                // Compute quadrature on reference interval [-1,1].
                GLTable( n, nodes, weights );

                // Map to desired interval [a,b].
                if ( a != -1.0 || b != 1.0 ) {

                    // Shift and scale nodes.
                    for ( int i = 0; i < n; ++i )
                        nodes[i] = a + (b - a)/2.0 * ( nodes[i] + 1.0 );

                    // Scale weights.
                    const double scale = (b - a) / 2.0;

                    for ( int i = 0; i < n; ++i )
                        weights[i] *= scale;
                }

                return;

            // Otherwise compute nodes and weights manually.
            } catch (...) {/* empty */}

        } /* fall through */

        default:
        {
            // Compute quadrature nodes.
            ComputeQuadratureNodes( n, nodes, type, a, b );

            // Compute quadrature weights.
            ComputeQuadratureWeights( n, nodes, weights, a, b );

        } return;

    }
}


//------------------------------------------------------------------------------------------------------------
//! \brief  Computes quadrature nodes of the specified type.
//!
//! \param[in]      n       Number of quadrature nodes of the specified type.
//! \param[out]     nodes   Pointer at which to store the computed nodes. <br/>
//!                         It is assumed that the necessary memory has been properly allocated at this
//!                         pointer prior to any calls to this function.
//! \param[in]      type    Quadrule::NodesType specifying the type of quadrature nodes to compute.
//! \param[in]      a       (optional) <br/>
//!                         Left endpoint of interval to scale quadrature nodes to.
//!                         Default value is -1.0.
//! \param[in]      b       (optional) <br/>
//!                         Right endpoint of interval to scale quadrature nodes to.
//!                         Default value is 1.0.
//!
//! \see    Quadrule::ComputeQuadratureWeights()
//------------------------------------------------------------------------------------------------------------
void Quadrule::ComputeQuadratureNodes (

    const int64_t n,
    double * const nodes,
    const NodesType type,
    const double a,         // = -1.0
    const double b          // = 1.0
) {

    // Compute standard quadrature nodes on reference interval [-1,1].
    switch ( type ) {

        case NodesType::GaussLegendre:

            ComputeGLnodes( n, nodes );
            break;

        case NodesType::GaussRadau:

            ComputeRadauNodes( n, nodes );
            break;

        case NodesType::GaussLobatto:

            ComputeLobattoNodes( n, nodes );
            break;

        case NodesType::Chebyshev:

            for ( int64_t i = 0; i < n; ++i )
                nodes[n-i-1] = std::cos( (double) (2*i + 1) * M_PI / (2*n) );

            break;

        case NodesType::ChebyshevRadau:

            for ( int64_t i = 0; i < n; ++i )
                nodes[n-i-1] = std::cos( 2*i * M_PI / (2*n - 1) );

            break;

        case NodesType::EquispacedBoth:

            for ( int64_t i = 0; i < n; ++i )
                nodes[i] = 2.0*i/(n - 1.0) - 1.0;

            break;

        case NodesType::EquispacedRight:

            for ( int64_t i = 0; i < n; ++i )
                nodes[i] = 2.0*(i + 1.0)/n - 1.0;

            break;

        default:
        {   std::string error_message = "Invalid NodesType '" + NodesType_to_String.at( type )
                                        + "' in '" + std::string(__func__) + "'.\n";

            PRINT_ERROR( error_message.c_str() )
            throw std::invalid_argument( error_message );
        }
    }

    // Scale quadrature nodes as needed.
    if ( a != -1.0 || b != 1.0 ) {

        for ( int i = 0; i < n; ++i )
            nodes[i] = a + (b - a)/2.0 * ( nodes[i] + 1.0 );
    }
}


//------------------------------------------------------------------------------------------------------------
//! \brief  Computes the quadrature weights for the given nodes which are optimal for integration across the
//!         interval \f$ [a,b] \f$.
//!
//! Computes the quadrature weights for the set of \f$ n \f$ nodes (stored in \pp{nodes}) which are optimal
//! for integrating polynomials across the interval \f$ [a,b] \f$ in the sense that the maximum number of
//! lowest degree polynomials are integrated exactly with the computed quadrature rule.
//!
//! The quadrature weights are computed using a Vandermonde-type system as described on page 143 of
//! \cite Press1992. The Vandermonde-type system is computed using a basis of Legendre polynomials, due to the
//! relatively good conditioning of resulting systems \cite Gautschi1990.
//!
//! \attention  It is assumed that the necessary memory has been properly allocated at the pointer
//!             \pp{weights} prior to any calls to this function.
//!
//! \param[in]      n           Number of quadrature nodes.
//! \param[in]      nodes       Pointer to array containing quadrature nodes to compute weights for.
//! \param[in]      weights     Upon return contains the optimal quadrature weights for integrating across
//!                             \f$ [a,b] \f$.
//! \param[in]      a           (optional) <br/>
//!                             Lower limit of integration.
//!                             Default value is -1.0.
//! \param[in]      b           (optional)
//!                             Upper limit of integration.
//!                             Default value is 1.0.
//!
//! \see    #Quadrule::ComputeQuadratureNodes()
//------------------------------------------------------------------------------------------------------------
void Quadrule::ComputeQuadratureWeights (

    const int64_t n,
    const double * const nodes,
    double * const weights,
    const double a,             // = -1.0
    const double b              // = 1.0
) {

    PRINT_NOTE( "Computing quadrature weights for %d nodes across the interval [ %.2e, %.2e ] "
                "by inverting Vandermonde matrix.\n", n, a, b )

    // Variables for LAPACK solve.
    int N      = n;
    int NRHS   = 1;
    double * A = new double[ N*N ];
    int LDA    = N;
    int * IPIV = new int[N];
    double * B = weights;
    int LDB    = N;
    int INFO   = 0;

    // --- Compute coefficients of Vandermonde-type matrix A. ----------------------------------------- //

    for ( int k = 0; k < n; ++k ) {
    for ( int j = 0; j < n; ++j ) {

        A[ k + j*N ] = Legendre( k, nodes[j] );
    }}

    // --- Compute coefficients of right-hand side of system. ----------------------------------------- //

    std::memset( B, 0, n * sizeof(double) );

    if ( a == -1.0  &&  b == 1.0 ) {

        B[0] = 2.0;

    } else {

        // Construct Gauss-Legendre quadrature.
        double * gauss_nodes   = new double[n];
        double * gauss_weights = new double[n];

        ComputeQuadrature( n, gauss_nodes, gauss_weights, NodesType::GaussLegendre );

        // Use Gauss-Legendre quadrature scaled to [a,b] to compute integrals.
        for ( int k = 0; k < n; ++k ) {
        for ( int l = 0; l < n; ++l ) {

            B[k] += gauss_weights[l] * Legendre( k, 0.5*( b*(1 + gauss_nodes[l]) + a*(1 - gauss_nodes[l]) ) );
        }}

        for ( int k = 0 ; k < n ; ++k )
            B[k] *= 0.5*(b - a);

        delete [] gauss_nodes;
        delete [] gauss_weights;
    }

    // --- Solve Vandermonde-type system for quadrature weights. -------------------------------------- //

    dgesv_( &N, &NRHS, A, &LDA, IPIV, B, &LDB, &INFO );

    // Check for error in solve.
    if ( INFO ) {

        std::string error_message = "dgesv_ returned INFO = " + std::to_string(INFO) + " inside '"
                                    + std::string(__func__) + "'.\n";

        PRINT_ERROR( error_message.c_str() )
        throw std::runtime_error( error_message );
    }

    // Cleanup auxiliary variables.
    delete [] A;
    delete [] IPIV;
}


//------------------------------------------------------------------------------------------------------------
//! \brief  Routine used to obtain precomputed Gauss-Kronrod quadrature rules for integration across the
//!         interval \f$ [-1,1] \f$.
//!
//! This routine only provides precomputed rules with \f$ n = 15, 21, 31, 41 \f$ computed with high precision
//! libraries.
//!
//! \attention  The size of the quadrature rule \f$ n \f$ corresponds to the Kronrod quadrature rule, which
//!             contains \f$ 2k + 1 \f$ nodes, where \f$ k \f$ is the number of nodes in the embedded Gauss
//!             quadrature rule.
//!
//! \param[in]      n           Size of quadrature rule to compute.
//! \param[out]     nodes       Upon return contains the set of Kronrod quadrature nodes. The gauss nodes are
//!                             stored at the odd indices.
//! \param[out]     G_weights   Upon return contains the set of quadrature weights for the Gauss nodes.
//!                             The weight for the Gauss node stored at index \f$ i \f$ in \c nodes is stored
//!                             at the index \f$ i/2 \f$.
//! \param[out]     K_weights   Upon return contains the set of quadrature weights for the Kronrod nodes.
//------------------------------------------------------------------------------------------------------------
void Quadrule::ComputeGKquad (

    const int n,
    double * const nodes,
    double * const G_weights,
    double * const K_weights
) {
    try         {  GKTable( n, nodes, G_weights, K_weights );  }
    catch (...) {  throw;                                      }
}
