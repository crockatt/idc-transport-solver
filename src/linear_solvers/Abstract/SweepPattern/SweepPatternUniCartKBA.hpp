//------------------------------------------------------------------------------------------------------------
//! \file   linear_solvers/Abstract/SweepPattern/SweepPatternUniCartKBA.hpp
//! \brief  Header file containing declaration of SweepPatternUniCartKBA class.
//!
//! \author Michael M. Crockatt
//! \date   March 2018
//------------------------------------------------------------------------------------------------------------

# ifndef __ABSTRACT__SWEEP_PATTERN_UNI_CART_KBA_HPP__
# define __ABSTRACT__SWEEP_PATTERN_UNI_CART_KBA_HPP__

# if defined (ENABLE_MPI)
# if SPACE_DIMS == 2


# include "linear_solvers/Abstract/SweepOperator.hpp"
# include "linear_solvers/Abstract/SweepPattern/SweepPatternUniCart.hpp"


//------------------------------------------------------------------------------------------------------------
//! \brief  Declaration of SweepPatternUniCartKBA class.
//!
//! \attention  This implementation assumes that the code has been compiled with OpenMP enabled and that the
//!             OpenMP environment has been configured to allow at least two threads.
//!
//! This implementation is motivated by the observation that most MPI implementations are unable to
//! progress background communications well enough to enable effective overlapping of communication and
//! computation \cite Denis2016. This implementation seeks to alleviate this shortcoming through the use of
//! a dedicated communication thread.
//!
//! This class implements:
//!     - SweepPattern::Print
//!     - SweepPattern::SetParameters
//!     - SweepPattern::Sweep
//!
//------------------------------------------------------------------------------------------------------------
template< class OrdinateFlux >
class SweepOperator<OrdinateFlux>::SweepPatternUniCartKBA :
    public SweepOperator<OrdinateFlux>::SweepPatternUniCart
{

    //!
    //! \brief  Friend instantiation of derived factory class.
    //!
    friend typename SweepPattern::template SweepPatternDerivedFactory<SweepPatternUniCartKBA>;


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
    virtual ~SweepPatternUniCartKBA( void );


    //========================================================================================================
    //=== OPERATOR OVERLOADS =================================================================================
    //========================================================================================================

    //--------------------------------------------------------------------------------------------------------
    //! \brief  Delete copy assignment operator.
    //--------------------------------------------------------------------------------------------------------
    SweepPatternUniCartKBA & operator=( const SweepPatternUniCartKBA & ) = delete;


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
# endif // if defined (ENABLE_MPI)
# endif // ifndef __ABSTRACT__SWEEP_PATTERN_UNI_CART_KBA_HPP__
