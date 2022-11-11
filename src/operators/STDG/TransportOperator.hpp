//------------------------------------------------------------------------------------------------------------
//! \file   operators/STDG/TransportOperator.hpp
//! \brief  Header for STDG TransportOperator routines.
//!
//! \author Michael M. Crockatt
//! \date   October 2017
//------------------------------------------------------------------------------------------------------------

# ifndef __STDG__TRANSPORT_OPERATOR_HPP__
# define __STDG__TRANSPORT_OPERATOR_HPP__


# include "objects/RKDG/CrossSection.hpp"
# include "objects/STDG/OrdinateFlux.hpp"
# include "objects/STDG/ScalarFlux.hpp"


namespace TransportOperator {


//------------------------------------------------------------------------------------------------------------
//! \brief  Computes \f$ \vec{y} \gets \beta \vec{y} + \alpha \mathcal{P} \vec{x} \f$ where
//!         \f$ \mathcal{P} \f$ is the operator that integrates over the angular dimension.
//------------------------------------------------------------------------------------------------------------
void Pmv( const double alpha, const double beta, const STDG::OrdinateFlux & x, STDG::ScalarFlux & y,
          const OpDomain op_domain = OpDomain::All );


//------------------------------------------------------------------------------------------------------------
//! \brief  Computes \f$ \vec{y} \gets \beta \vec{y} + \alpha \mathcal{S} \vec{x} \f$ where
//!         \f$ \mathcal{S} \f$ is the angular redistribution operator.
//------------------------------------------------------------------------------------------------------------
void Smv( const double alpha, const double beta, const RKDG::CrossSection & sigma, const STDG::ScalarFlux & x,
          STDG::OrdinateFlux & y, const OpDomain op_domain = OpDomain::All );


} // namespace TransportOperator


# endif // ifndef __STDG__TRANSPORT_OPERATOR_HPP__
