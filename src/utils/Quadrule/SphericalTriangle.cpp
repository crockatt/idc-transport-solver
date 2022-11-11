//------------------------------------------------------------------------------------------------------------
//! \file   utils/Quadrule/SphericalTriangle.cpp
//! \brief  Implementation of routines to generate \f$ T_N \f$ quadratures based on a spherical triangle
//!         decomposition of \f$ \mathbb{S}^2 \f$.
//!
//! Uses a tessellation of \f$ \mathbb{S}^2 \f$ by spherical triangles to construct the nodes and weights of
//! the \f$ T_N \f$ quadrature \cite Thurgood1992 \cite Thurgood1995 .
//!
//! \author Michael M. Crockatt
//! \date   August 2017
//------------------------------------------------------------------------------------------------------------


# include "utils/global.hpp"
# include "utils/Quadrule/Quadrule.hpp"
# include "utils/Quadrule/SphericalTriangle.hpp"

# include <list>


using namespace Quadrule;


//============================================================================================================
//=== STATIC HELPER ROUTINES =================================================================================
//============================================================================================================


//------------------------------------------------------------------------------------------------------------
//! \brief  Converts from \f$ (u,v) \f$ coordinates in the basal triangle to \f$ (x,y,z) \f$ coordinates in
//!         the basal plane \f$ x + y + z - 1 = 0 \f$.
//!
//! Cf. equation (3.10) in \cite Thurgood1992.
//!
//! \param[in]  uv      \f$ (u,v) \f$ coordinates (in that order).
//! \param[out] xyz     Upon return, contains the coordinates \f$ (x,y,z) \f$ (in that order) corresponding to
//!                     the given coordinates \f$ (u,v) \f$.
//------------------------------------------------------------------------------------------------------------
static void uv2xyz (

    const double uv[2],
    double xyz[3]
) {

    const double & u = uv[0];
    const double & v = uv[1];

    double & x = xyz[0];
    double & y = xyz[1];
    double & z = xyz[2];

    x = 0.5 - u - 0.5*v;
    y = 0.5 + u - 0.5*v;
    z = v;
}


//------------------------------------------------------------------------------------------------------------
//! \brief  Projects a point with coordinates \f$ (x,y,z) \f$ in the basal plane onto the first octant of
//!         \f$ \mathbb{S}^2 \f$.
//!
//! Cf. equation (3.8) in \cite Thurgood1992.
//!
//! \param[in]  xyz     \f$ (x,y,z) \f$ coordinates (in that order).
//! \param[out] sph     Upon return, contains the coordinates of the projection of the points \f$ (x,y,z) \f$
//!                     onto the sphere.
//------------------------------------------------------------------------------------------------------------
static void xyz2sph (

    const double xyz[3],
    double sph[3]
) {

    const double & x = xyz[0];
    const double & y = xyz[1];
    const double & z = xyz[2];

    double & a = sph[0];
    double & b = sph[1];
    double & c = sph[2];

    const double norm = std::sqrt( x*x + y*y + z*z );

    a = x / norm;
    b = y / norm;
    c = z / norm;
}


//============================================================================================================
//=== TRIANGLE CLASS DEFINITIONS FOR RECURSIVE IMPLEMENTATION ================================================
//============================================================================================================


//------------------------------------------------------------------------------------------------------------
//! \brief  Triangle struct used for recursive computation of \f$ T_N \f$ triangular quadrature.
//------------------------------------------------------------------------------------------------------------
struct Triangle {

    // (u,v) coordinates of the vertices of the triangle wrt. the basal triangle.
    double vertex1_uv[2];
    double vertex2_uv[2];
    double vertex3_uv[2];

    // (x,y,z) coordinates of the vertices of the triangle wrt. the basal plane.
    double vertex1_xyz[3];
    double vertex2_xyz[3];
    double vertex3_xyz[3];

    // (ξ,η,μ) coordinates of the vertices of the triangle on 𝕊².
    double vertex1_sph[3];
    double vertex2_sph[3];
    double vertex3_sph[3];

    // Coordinates of the centroid of the triangle.
    double centroid_uv[2];      // (u,v) coordinates wrt. the basal triangle.
    double centroid_xyz[3];     // (x,y,z) coordinates wrt. the basal plane.
    double centroid_sph[3];     // (ξ,η,μ) coordinates on 𝕊².


    //! \brief  Delete default constructor.
    Triangle( void ) = delete;

    //! \brief  Use default copy constructor.
    Triangle( const Triangle & ) = default;

    //! \brief  Construct a triangle given the \f$ (u,v) \f$ coordinates of the triangles vertices.
    Triangle(
        const double vertex1_uv_in[2],
        const double vertex2_uv_in[2],
        const double vertex3_uv_in[2]
    );


    //! \brief  Computes the quadrature weight associated with the given triangle.
    double ComputeWeight( void ) const;

    //! \brief  Subdivides a triangle object by powers of two as determined by \pp{order}.
    std::list<Triangle> Subdivide( const uint64_t order ) const;
};


//------------------------------------------------------------------------------------------------------------
//! \brief  Construct a triangle given the \f$ (u,v) \f$ coordinates of the triangles vertices.
//------------------------------------------------------------------------------------------------------------
Triangle::Triangle (

    const double vertex1_uv_in[2],
    const double vertex2_uv_in[2],
    const double vertex3_uv_in[2]
) {

    // Copy input values.
    this->vertex1_uv[0] = vertex1_uv_in[0];
    this->vertex1_uv[1] = vertex1_uv_in[1];

    this->vertex2_uv[0] = vertex2_uv_in[0];
    this->vertex2_uv[1] = vertex2_uv_in[1];

    this->vertex3_uv[0] = vertex3_uv_in[0];
    this->vertex3_uv[1] = vertex3_uv_in[1];

    // Compute (u,v) coordinates of centroid.
    centroid_uv[0] = ( vertex1_uv[0] + vertex2_uv[0] + vertex3_uv[0] ) / 3.0;
    centroid_uv[1] = ( vertex1_uv[1] + vertex2_uv[1] + vertex3_uv[1] ) / 3.0;

    // Compute (x,y,z) coordinates of all points wrt. the basal plane.
    uv2xyz( vertex1_uv, vertex1_xyz );
    uv2xyz( vertex2_uv, vertex2_xyz );
    uv2xyz( vertex3_uv, vertex3_xyz );
    uv2xyz( centroid_uv, centroid_xyz );

    // Compute (ξ,η,μ) coordinates of projection of all points onto 𝕊².
    xyz2sph( vertex1_xyz, vertex1_sph );
    xyz2sph( vertex2_xyz, vertex2_sph );
    xyz2sph( vertex3_xyz, vertex3_sph );
    xyz2sph( centroid_xyz, centroid_sph );
}


//------------------------------------------------------------------------------------------------------------
//! \brief  Computes the quadrature weight associated with the given triangle.
//!
//! \return     Returns the quadrature weight of the associated projected triangle.
//------------------------------------------------------------------------------------------------------------
double Triangle::ComputeWeight ( void ) const {

    // Compute arc lengths of sides of spherical triangle.
    const double a = std::acos( vertex2_sph[0]*vertex3_sph[0] + vertex2_sph[1]*vertex3_sph[1]
                                                              + vertex2_sph[2]*vertex3_sph[2] );

    const double b = std::acos( vertex1_sph[0]*vertex3_sph[0] + vertex1_sph[1]*vertex3_sph[1]
                                                              + vertex1_sph[2]*vertex3_sph[2] );

    const double c = std::acos( vertex1_sph[0]*vertex2_sph[0] + vertex1_sph[1]*vertex2_sph[1]
                                                              + vertex1_sph[2]*vertex2_sph[2] );

    // Compute semi-perimeter of spherical triangle.
    const double s = ( a + b + c ) / 2.0;

    // Compute interior angles of spherical triangle.
    const double tanr = std::sqrt( std::sin(s-a) * std::sin(s-b) * std::sin(s-c) / std::sin(s) );

    const double A = 2.0 * std::atan( tanr / std::sin(s-a) );
    const double B = 2.0 * std::atan( tanr / std::sin(s-b) );
    const double C = 2.0 * std::atan( tanr / std::sin(s-c) );

    // Associated quadrature weight is the area of the spherical triangle in 𝕊².
    return A + B + C - M_PI;
}


//------------------------------------------------------------------------------------------------------------
//! \brief
//------------------------------------------------------------------------------------------------------------
std::list<Triangle> Triangle::Subdivide (

    const uint64_t order

) const {

    // Base case for recursion: return the triangle itself.
    if ( order == 0 ) {  return std::list<Triangle>{ *this };  }

    // Compute midpoints on edges of triangle.
    const double midpoint1_uv[2] = { ( vertex1_uv[0] + vertex2_uv[0] ) / 2.0,
                                     ( vertex1_uv[1] + vertex2_uv[1] ) / 2.0 };

    const double midpoint2_uv[2] = { ( vertex1_uv[0] + vertex3_uv[0] ) / 2.0,
                                     ( vertex1_uv[1] + vertex3_uv[1] ) / 2.0 };

    const double midpoint3_uv[2] = { ( vertex2_uv[0] + vertex3_uv[0] ) / 2.0,
                                     ( vertex2_uv[1] + vertex3_uv[1] ) / 2.0 };

    // Compute ordered list of sub-triangles.
    std::list<Triangle> subtriangles {};

    subtriangles.splice( subtriangles.end(),
                         Triangle( vertex1_uv, midpoint1_uv, midpoint2_uv ).Subdivide( order >> 1 ) );

    subtriangles.splice( subtriangles.end(),
                         Triangle( midpoint1_uv, vertex2_uv, midpoint3_uv ).Subdivide( order >> 1 ) );

    subtriangles.splice( subtriangles.end(),
                         Triangle( midpoint3_uv, midpoint1_uv, midpoint2_uv ).Subdivide( order >> 1 ) );

    subtriangles.splice( subtriangles.end(),
                         Triangle( midpoint2_uv, midpoint3_uv, vertex3_uv ).Subdivide( order >> 1 ) );

    return subtriangles;
}


//============================================================================================================
//=== INTERFACE ROUTINES TO OBTAIN QUADRATURE ================================================================
//============================================================================================================


//------------------------------------------------------------------------------------------------------------
//! \brief  Computes the \f$ T_N \f$ quadrature rule of the specified order.
//!
//! \param[in]  order   Order of quadrature rule to compute.
//! \param[out] x       Pointer to memory location at which to store \f$ x \f$ coordinates of computed nodes.
//! \param[out] y       Pointer to memory location at which to store \f$ y \f$ coordinates of computed nodes.
//! \param[out] z       Pointer to memory location at which to store \f$ z \f$ coordinates of computed nodes.
//! \param[out] w       Pointer to memory location at which to store quadrature weights.
//------------------------------------------------------------------------------------------------------------
void Quadrule::ComputeSphericalTriangleQuadrature (

    const int64_t order,
    double * const x,
    double * const y,
    double * const z,
    double * const w
) {

    std::list<Triangle> subtriangles {};

    //
    // If order is a power of two, generate subtriangles using a recursive decomposition. This ordering of the
    // ordinates permits a hierarchical relabeling procedure from the resulting quadrature set to T_N
    // quadratures of higher order.
    //
    if ( IsPowerOfTwo( order ) ) {

        // Vertices of basal triangle.
        double vertex1_uv[2] = {  0.0, 1.0 };
        double vertex2_uv[2] = { -0.5, 0.0 };
        double vertex3_uv[2] = {  0.5, 0.0 };

        subtriangles = Triangle( vertex1_uv, vertex2_uv, vertex3_uv ).Subdivide( order >> 1 );

    //
    // If order is not a power of two, then compute subtriangles using a generic row-based decomposition. This
    // ordering of the ordinates does not permit relabeling from the associated quadrature set.
    //
    } else {

        for ( int64_t n = 0; n < order;   ++n ) {   // n: Indexes row of triangles in tessellation.
        for ( int64_t k = 0; k < 2*n + 1; ++k ) {   // k: Indexes triangle across row in tessellation.

            // (u,v) coordinates of the vertices of the triangle wrt. the basal triangle.
            double vertex1_uv[2] = {0,0};
            double vertex2_uv[2] = {0,0};
            double vertex3_uv[2] = {0,0};

            // Compute vertices of triangle in (u,v) coordinates.
            if ( k % 2 ) { // k odd: Triangle points down.

                vertex1_uv[0] = ((double)(k-n+1)) / (2*order);
                vertex1_uv[1] = 1.0 - ((double)(n)) / order;

                vertex2_uv[0] = ((double)(k-n-1)) / (2*order);
                vertex2_uv[1] = 1.0 - ((double)(n)) / order;

                vertex3_uv[0] = ((double)(k-n)) / (2*order);
                vertex3_uv[1] = 1.0 - ((double)(n+1)) / order;

            } else { // k even: Triangle points up.

                vertex1_uv[0] = ((double)(k-n)) / (2*order);
                vertex1_uv[1] = 1.0 - ((double)(n)) / order;

                vertex2_uv[0] = ((double)(k-n-1)) / (2*order);
                vertex2_uv[1] = 1.0 - ((double)(n+1)) / order;

                vertex3_uv[0] = ((double)(k-n+1)) / (2*order);
                vertex3_uv[1] = 1.0 - ((double)(n+1)) / order;
            }

            // Add triangle object to list.
            subtriangles.emplace_back( Triangle( vertex1_uv, vertex2_uv, vertex3_uv ) );
        }}
    }

    // Compute quadrature from list of triangles.
    int64_t q = 0;

    for ( auto & triangle : subtriangles ) {

        // Compute coordinates of quadrature node.
        x[q] = triangle.centroid_sph[0];
        y[q] = triangle.centroid_sph[1];
        z[q] = triangle.centroid_sph[2];

        // Compute quadrature weight.
        w[q] = triangle.ComputeWeight();

        ++q;
    }
}
