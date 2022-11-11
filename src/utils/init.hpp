//------------------------------------------------------------------------------------------------------------
//! \file   utils/init.hpp
//! \brief  Header for initialization routine.
//!
//! \author Michael M. Crockatt
//! \date   June 2018
//------------------------------------------------------------------------------------------------------------

# ifndef __INIT_HPP__
# define __INIT_HPP__

# include "objects/RKDG/OrdinateFlux.hpp"
# include "objects/RKDG/CrossSection.hpp"
# include "utils/global.hpp"


//------------------------------------------------------------------------------------------------------------
//! \brief  Performs initialization of cross sections, initial condition, and sources.
//------------------------------------------------------------------------------------------------------------
void Init(
    RKDG::OrdinateFlux & psi,
    RKDG::OrdinateFlux & source,
    RKDG::CrossSection & sigma_t,
    RKDG::CrossSection & sigma_s,
    DomainDecomposition & spatial_params
);


# endif // ifndef __INIT_HPP__
