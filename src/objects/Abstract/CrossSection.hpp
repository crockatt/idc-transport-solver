//------------------------------------------------------------------------------------------------------------
//! \file   objects/Abstract/CrossSection.hpp
//! \brief  Header for Abstract::CrossSection typedef declaration.
//!
//! \author Michael Crockatt
//! \date   April 2017
//------------------------------------------------------------------------------------------------------------

# ifndef __ABSTRACT__CROSS_SECTION_HPP__
# define __ABSTRACT__CROSS_SECTION_HPP__


# include "objects/Abstract/DensityFunction.hpp"


namespace Abstract {

//------------------------------------------------------------------------------------------------------------
//! \brief  Class for storing discrete scalar flux densities.
//!
//!
//------------------------------------------------------------------------------------------------------------
typedef DensityFunction CrossSection;


} // namespace Abstract


# endif // ifndef __ABSTRACT__CROSS_SECTION_HPP__
