//------------------------------------------------------------------------------------------------------------
//! \file   linear_solvers/Abstract/SweepPattern/SweepPatternUniCartWavefront.hpp
//! \brief  Header file containing declaration of SweepPatternUniCartWavefront class.
//!
//! \author Michael M. Crockatt
//! \date   March 2018
//------------------------------------------------------------------------------------------------------------

# ifndef __ABSTRACT__SWEEP_PATTERN_UNI_CART_WAVEFRONT_HPP__
# define __ABSTRACT__SWEEP_PATTERN_UNI_CART_WAVEFRONT_HPP__

# if SPACE_DIMS == 2


# include "linear_solvers/Abstract/SweepOperator.hpp"
# include "linear_solvers/Abstract/SweepPattern/SweepPatternUniCart.hpp"


//------------------------------------------------------------------------------------------------------------
//! \brief  Declaration of SweepPatternUniCartWavefront class.
//!
//! Spatial mesh is traversed using a wavefront pattern. All MPI ranks update simultaneously using existing
//! halo information.
//!
//! This class implements:
//!     - SweepPattern::Print
//!     - SweepPattern::SetParameters
//!     - SweepPattern::Sweep
//!
//------------------------------------------------------------------------------------------------------------
template< class OrdinateFlux >
class SweepOperator<OrdinateFlux>::SweepPatternUniCartWavefront :
    public SweepOperator<OrdinateFlux>::SweepPatternUniCart
{

    //!
    //! \brief  Friend instantiation of derived factory class.
    //!
    friend typename SweepPattern::template SweepPatternDerivedFactory<SweepPatternUniCartWavefront>;


    //========================================================================================================
    //=== CONSTRUCTORS AND DESTRUCTOR ========================================================================
    //========================================================================================================

protected:

    //!
    //! \brief  Use protected constructors of SweepPatternUniCart.
    //!
    using SweepPatternUniCart::SweepPatternUniCart;

public:

    //--------------------------------------------------------------------------------------------------------
    //! \brief  Virtual destructor.
    //--------------------------------------------------------------------------------------------------------
    virtual ~SweepPatternUniCartWavefront( void );


    //========================================================================================================
    //=== OPERATOR OVERLOADS =================================================================================
    //========================================================================================================

    //--------------------------------------------------------------------------------------------------------
    //! \brief  Delete copy assignment operator.
    //--------------------------------------------------------------------------------------------------------
    SweepPatternUniCartWavefront & operator=( const SweepPatternUniCartWavefront & ) = delete;


    //========================================================================================================
    //=== INTERFACE OVERRIDES ================================================================================
    //========================================================================================================

    //--------------------------------------------------------------------------------------------------------
    //! \brief  Returns the string descriptor of the class type.
    //--------------------------------------------------------------------------------------------------------
    static const std::string Descriptor( void );

    //--------------------------------------------------------------------------------------------------------
    //! \brief  Returns the string descriptor for an objects type.
    //--------------------------------------------------------------------------------------------------------
    const std::string GetDescriptor( void ) const override;

    //--------------------------------------------------------------------------------------------------------
    //! \brief  Prints a summary of the object's configuration to the logging interface.
    //--------------------------------------------------------------------------------------------------------
    void Print( const std::string & = "  " ) const override;

    //--------------------------------------------------------------------------------------------------------
    //! \brief  Sets parameters given in a ParameterList.
    //--------------------------------------------------------------------------------------------------------
    void SetParameters( const ParameterList & ) override;

    //--------------------------------------------------------------------------------------------------------
    //! \brief  Performs a transport sweep.
    //--------------------------------------------------------------------------------------------------------
    void Sweep (
        const OrdinateFlux & source,
        OrdinateFlux & result,
        const RKDG::CrossSection & sigma,
        const RKDG::OrdinateFlux * const initial = nullptr
    ) const override;

};


# endif // if SPACE_DIMS == 2
# endif // ifndef __ABSTRACT__SWEEP_PATTERN_UNI_CART_WAVEFRONT_HPP__
