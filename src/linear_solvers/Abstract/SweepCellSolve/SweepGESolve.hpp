//------------------------------------------------------------------------------------------------------------
//! \file   linear_solvers/Abstract/SweepCellSolve/SweepGESolve.hpp
//! \brief  Header file containing declaration of SweepGESolve class.
//!
//! \author Michael M. Crockatt
//! \date   March 2018
//------------------------------------------------------------------------------------------------------------

# ifndef __ABSTRACT__SWEEP_GE_SOLVE_HPP__
# define __ABSTRACT__SWEEP_GE_SOLVE_HPP__


# include "linear_solvers/Abstract/SweepOperator.hpp"
# include "linear_solvers/Abstract/SweepCellSolve/SweepLinearSolveSpec.hpp"


//------------------------------------------------------------------------------------------------------------
//! \brief  Declaration of SweepGESolve class.
//!
//! Solves local linear systems in each mesh element using Gaussian elimination.
//!
//! This class template implements:
//!     - SweepCellSolveBase::Print
//!     - SweepCellSolveBase::SolveSystem
//!
//------------------------------------------------------------------------------------------------------------
template< class OrdinateFlux, int64_t SIMD_length >
class SweepGESolve :
    public SweepLinearSolveSpec< OrdinateFlux, SIMD_length >
{

public:

    //========================================================================================================
    //=== CONSTRUCTORS AND DESTRUCTOR ========================================================================
    //========================================================================================================

    //!
    //! \brief  Use constructors of SweepLinearSolveSpec.
    //!
    using SweepLinearSolveSpec< OrdinateFlux, SIMD_length >::SweepLinearSolveSpec;

    //--------------------------------------------------------------------------------------------------------
    //! \brief  Virtual destructor.
    //--------------------------------------------------------------------------------------------------------
    virtual ~SweepGESolve( void );


    //========================================================================================================
    //=== OPERATOR OVERLOADS =================================================================================
    //========================================================================================================

    //--------------------------------------------------------------------------------------------------------
    //! \brief  Delete copy assignment operator.
    //--------------------------------------------------------------------------------------------------------
    SweepGESolve< OrdinateFlux, SIMD_length > &
        operator=( const SweepGESolve< OrdinateFlux, SIMD_length > & ) = delete;


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

protected:

    //--------------------------------------------------------------------------------------------------------
    //! \brief  Solves the system for the updated coefficients.
    //--------------------------------------------------------------------------------------------------------
    void SolveSystem (
        double * const A,
        double * const B,
        void * const // work_ptr
    ) const override;

};


# endif // ifndef __ABSTRACT__SWEEP_GE_SOLVE_HPP__
