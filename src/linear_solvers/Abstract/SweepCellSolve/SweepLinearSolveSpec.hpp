//------------------------------------------------------------------------------------------------------------
//! \file   linear_solvers/Abstract/SweepCellSolve/SweepLinearSolveSpec.hpp
//! \brief  Header file containing declaration of SweepLinearSolveSpec abstract class.
//!
//! \author Michael M. Crockatt
//! \date   March 2018
//------------------------------------------------------------------------------------------------------------

# ifndef __ABSTRACT__SWEEP_LINEAR_SOLVE_SPEC_HPP__
# define __ABSTRACT__SWEEP_LINEAR_SOLVE_SPEC_HPP__


# include "linear_solvers/Abstract/SweepCellSolve/SweepLinearSolveBase.hpp"


//------------------------------------------------------------------------------------------------------------
//! \brief  Declaration of SweepLinearSolveSpec abstract class.
//!
//! This class provides no implementations.
//! Specializations implement:
//!     - SweepCellSolveBase::CalcA
//!     - SweepCellSolveBase::CalcB
//!
//------------------------------------------------------------------------------------------------------------
template< class OrdinateFlux, int64_t SIMD_length >
class SweepLinearSolveSpec :
    public SweepLinearSolveBase< OrdinateFlux, SIMD_length >
{

public:

    //========================================================================================================
    //=== CONSTRUCTORS AND DESTRUCTOR ========================================================================
    //========================================================================================================

    //--------------------------------------------------------------------------------------------------------
    //! \brief  Delete default constructor.
    //--------------------------------------------------------------------------------------------------------
    SweepLinearSolveSpec( void ) = delete;

    //--------------------------------------------------------------------------------------------------------
    //! \brief  Delete copy constructor.
    //--------------------------------------------------------------------------------------------------------
    SweepLinearSolveSpec( const SweepLinearSolveSpec< OrdinateFlux, SIMD_length > & ) = delete;

    //--------------------------------------------------------------------------------------------------------
    //! \brief  Construct solver object from given parameters.
    //--------------------------------------------------------------------------------------------------------
    SweepLinearSolveSpec (
        const SweepOperator<OrdinateFlux> & enclosing,
        const ParameterList & input_list,
        const double dt_in = std::numeric_limits<double>::infinity()
    );

    //--------------------------------------------------------------------------------------------------------
    //! \brief  Virtual destructor.
    //--------------------------------------------------------------------------------------------------------
    virtual ~SweepLinearSolveSpec( void );


    //========================================================================================================
    //=== OPERATOR OVERLOADS =================================================================================
    //========================================================================================================

    //--------------------------------------------------------------------------------------------------------
    //! \brief  Delete copy assignment operator.
    //--------------------------------------------------------------------------------------------------------
    SweepLinearSolveSpec< OrdinateFlux, SIMD_length > &
        operator=( const SweepLinearSolveSpec< OrdinateFlux, SIMD_length > & ) = delete;


protected:

    //========================================================================================================
    //=== PROTECTED HELPER FUNCTIONS =========================================================================
    //========================================================================================================

    //--------------------------------------------------------------------------------------------------------
    //! \brief  Computes the elements of the matrix/matrices for the linear system(s) to be solved at the
    //!         given step of a sweep.
    //--------------------------------------------------------------------------------------------------------
    void CalcA (
        double * const A,
        const SIMD_BlkIdx<0> & idx,
        const RKDG::CrossSection & sigma,
        void * const work_ptr
    ) const override;

    //--------------------------------------------------------------------------------------------------------
    //! \brief  Computes the elements of the right-hand vectors \f$ b \f$ to be used for solving the linear
    //!         systems at the given step of a sweep.
    //--------------------------------------------------------------------------------------------------------
    void CalcB (
        double * const B,
        const SIMD_BlkIdx<0> & idx,
        const OrdinateFlux & source,
        const OrdinateFlux & result,
        void * const work_ptr,
        const RKDG::OrdinateFlux * const initial
    ) const override;

    //--------------------------------------------------------------------------------------------------------
    //! \brief  Computes the angular component of the matrices that define the linear systems.
    //--------------------------------------------------------------------------------------------------------
    void ComputeMatrices( void ) const override;

};


//
// Include headers containing specializations.
//
# include "linear_solvers/RKDG/SweepCellSolve/SweepLinearSolveSpec.hpp"
# include "linear_solvers/STDG/SweepCellSolve/SweepLinearSolveSpec.hpp"


# endif // ifndef __ABSTRACT__SWEEP_LINEAR_SOLVE_SPEC_HPP__
