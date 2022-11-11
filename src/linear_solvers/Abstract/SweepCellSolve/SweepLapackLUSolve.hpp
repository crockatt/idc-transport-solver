//------------------------------------------------------------------------------------------------------------
//! \file   linear_solvers/Abstract/SweepCellSolve/SweepLapackLUSolve.hpp
//! \brief  Header file containing declaration of SweepLapackLUSolve class.
//!
//! \author Michael M. Crockatt
//! \date   March 2018
//------------------------------------------------------------------------------------------------------------

# ifndef __ABSTRACT__SWEEP_LAPACK_LU_SOLVE_HPP__
# define __ABSTRACT__SWEEP_LAPACK_LU_SOLVE_HPP__


# include "linear_solvers/Abstract/SweepOperator.hpp"
# include "linear_solvers/Abstract/SweepCellSolve/SweepLinearSolveSpec.hpp"


//------------------------------------------------------------------------------------------------------------
//! \brief  Declaration of SweepLapackLUSolve class.
//!
//! Solves local linear systems in each mesh element using Gaussian elimination.
//!
//! This class template implements:
//!     - SweepCellSolveBase::Print
//!     - SweepCellSolveBase::SolveSystem
//!
//------------------------------------------------------------------------------------------------------------
template< class OrdinateFlux, int64_t SIMD_length >
class SweepLapackLUSolve :
    public SweepLinearSolveSpec< OrdinateFlux, SIMD_length >
{

public:

    //========================================================================================================
    //=== CONSTRUCTORS AND DESTRUCTOR ========================================================================
    //========================================================================================================

    //!
    //! \brief  Use protected constructors of SweepLinearSolveSpec.
    //!
    using SweepLinearSolveSpec< OrdinateFlux, SIMD_length >::SweepLinearSolveSpec;

    //--------------------------------------------------------------------------------------------------------
    //! \brief  Virtual destructor.
    //--------------------------------------------------------------------------------------------------------
    virtual ~SweepLapackLUSolve( void );


    //========================================================================================================
    //=== OPERATOR OVERLOADS =================================================================================
    //========================================================================================================

    //--------------------------------------------------------------------------------------------------------
    //! \brief  Delete copy assignment operator.
    //--------------------------------------------------------------------------------------------------------
    SweepLapackLUSolve< OrdinateFlux, SIMD_length > &
        operator=( const SweepLapackLUSolve< OrdinateFlux, SIMD_length > & ) = delete;


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
        void * const work_ptr
    ) const override;

};


# endif // ifndef __ABSTRACT__SWEEP_LAPACK_LU_SOLVE_HPP__
