//------------------------------------------------------------------------------------------------------------
//! \file   operators/RKDG/TransportOperator.hpp
//! \brief  Header for RKDG TransportOperator routines.
//!
//! \author Michael M. Crockatt
//! \date   October 2017
//------------------------------------------------------------------------------------------------------------

# ifndef __RKDG__TRANSPORT_OPERATOR_HPP__
# define __RKDG__TRANSPORT_OPERATOR_HPP__


# include "objects/RKDG/CrossSection.hpp"
# include "objects/RKDG/OrdinateFlux.hpp"
# include "objects/RKDG/ScalarFlux.hpp"


namespace TransportOperator {


//------------------------------------------------------------------------------------------------------------
//! \brief  Computes \f$ \vec{y} \gets \beta \vec{y} + \alpha \mathcal{P} \vec{x} \f$ where
//!         \f$ \mathcal{P} \f$ is the operator that integrates over the angular dimension.
//------------------------------------------------------------------------------------------------------------
void Pmv( const double alpha, const double beta, const RKDG::OrdinateFlux & x, RKDG::ScalarFlux & y,
          const OpDomain op_domain = OpDomain::All );


//------------------------------------------------------------------------------------------------------------
//! \brief  Computes \f$ \vec{y} \gets \beta \vec{y} + \alpha \mathcal{S} \vec{x} \f$ where
//!         \f$ \mathcal{S} \f$ is the angular redistribution operator.
//------------------------------------------------------------------------------------------------------------
void Smv( const double alpha, const double beta, const RKDG::CrossSection & sigma, const RKDG::ScalarFlux & x,
          RKDG::OrdinateFlux & y, const OpDomain op_domain = OpDomain::All );


//------------------------------------------------------------------------------------------------------------
//! \brief  Computes \f$ \vec{y} \gets \beta \vec{y} + \alpha \mathcal{L} \vec{x} \f$ where
//!         \f$ \mathcal{L} \f$ is the advection-absorption operator with total cross-section given by
//!         \pp{sigma}.
//------------------------------------------------------------------------------------------------------------
void Lmv( const double alpha, const double beta, const RKDG::CrossSection & sigma,
          RKDG::OrdinateFlux & x, RKDG::OrdinateFlux & y );


} // namespace TransportOperator


# endif // ifndef __RKDG__TRANSPORT_OPERATOR_HPP__
