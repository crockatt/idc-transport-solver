//------------------------------------------------------------------------------------------------------------
//! \file   utils/Quadrule/Tabulated.hpp
//!
//! \author Michael Crockatt
//! \date   June 2017
//------------------------------------------------------------------------------------------------------------

# ifndef __QUADRULE__TABULATED_HPP__
# define __QUADRULE__TABULATED_HPP__


# include "utils/Quadrule/Quadrule.hpp"


namespace Quadrule {


void GLTable( const int64_t n, double * const nodes, double * const weights = nullptr );

void GKTable( const int64_t n, double * const nodes, double * const G_weights, double * const K_weights );


} // namespace Quadrule


# endif // ifndef __QUADRULE__TABULATED_HPP__
