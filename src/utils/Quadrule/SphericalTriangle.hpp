//------------------------------------------------------------------------------------------------------------
//! \file   utils/Quadrule/SphericalTriangle.hpp
//! \brief  Header file for implementation of routines to construct \f$ T_N \f$ or spherical triangle
//!         quadrature over \f$ \mathbb{S}^2 \f$R.
//!
//! \author Michael M. Crockatt
//! \date   July 2017
//------------------------------------------------------------------------------------------------------------

# ifndef __QUADRULE__SPHERICAL_TRIANGLE_HPP__
# define __QUADRULE__SPHERICAL_TRIANGLE_HPP__


namespace Quadrule {


void ComputeSphericalTriangleQuadrature( const int64_t order,
                                         double * const x, double * const y, double * const z,
                                         double * const w );


} // namespace Quadrule


# endif // ifndef __QUADRULE__SPHERICAL_TRIANGLE_HPP__
