//------------------------------------------------------------------------------------------------------------
//! \file   linear_solvers/Abstract/SweepPattern/SweepPatternUniCartSequential.hpp
//! \brief  Header file containing declaration of SweepPatternUniCartSequential class.
//!
//! \author Michael M. Crockatt
//! \date   March 2018
//------------------------------------------------------------------------------------------------------------

# ifndef __ABSTRACT__SWEEP_PATTERN_UNI_CART_SEQUENTIAL_HPP__
# define __ABSTRACT__SWEEP_PATTERN_UNI_CART_SEQUENTIAL_HPP__


# include "linear_solvers/Abstract/SweepOperator.hpp"
# include "linear_solvers/Abstract/SweepPattern/SweepPatternUniCart.hpp"


//------------------------------------------------------------------------------------------------------------
//! \brief  Declaration of SweepPatternUniCartSequential class.
//!
//! Spatial mesh is traversed cell-by-cell. All MPI ranks update simultaneously using existing halo
//! information.
//!
//! This class implements:
//!     - SweepPattern::Print
//!     - SweepPattern::SetParameters
//!     - SweepPattern::Sweep
//!
//------------------------------------------------------------------------------------------------------------
template< class OrdinateFlux >
class SweepOperator<OrdinateFlux>::SweepPatternUniCartSequential :
    public SweepOperator<OrdinateFlux>::SweepPatternUniCart
{

    //!
    //! \brief  Friend instantiation of derived factory class.
    //!
    friend typename SweepPattern::template SweepPatternDerivedFactory<SweepPatternUniCartSequential>;


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
    virtual ~SweepPatternUniCartSequential( void );


    //========================================================================================================
    //=== OPERATOR OVERLOADS =================================================================================
    //========================================================================================================

    //--------------------------------------------------------------------------------------------------------
    //! \brief  Delete copy assignment operator.
    //--------------------------------------------------------------------------------------------------------
    SweepPatternUniCartSequential & operator=( const SweepPatternUniCartSequential & ) = delete;


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


# endif // ifndef __ABSTRACT__SWEEP_PATTERN_UNI_CART_SEQUENTIAL_HPP__
