//------------------------------------------------------------------------------------------------------------
//! \file   linear_solvers/Abstract/SweepCellSolve/SweepCellSolveBase.cpp
//! \brief  Contains implementations and instantiations of methods from SweepCellSolveBase class template.
//!
//! \author Michael M. Crockatt
//! \date   March 2018
//------------------------------------------------------------------------------------------------------------

# include <limits>
# include <map>
# include <string>
# include <typeinfo>

# include "linear_solvers/Abstract/SweepCellSolve/SweepCellSolveBase.hpp"


//============================================================================================================
//=== MEMBER DEFINITIONS FOR SweepCellSolveBase CLASS ========================================================
//============================================================================================================


//------------------------------------------------------------------------------------------------------------
//! \brief  Initializes a SweepCellSolveBase object.
//!
//! \param[in]  enclosing   Reference to instance of enclosing SweepOperator class.
//! \param[in]  dt_in       Initial timestep size.
//------------------------------------------------------------------------------------------------------------
template< class OrdinateFlux >
SweepCellSolveBase<OrdinateFlux>::SweepCellSolveBase (

    const SweepOperator<OrdinateFlux> & enclosing,
    const double dt_in // = std::numeric_limits<double>::infinity()
) :
    dt{ dt_in },
    system_dim{ 0 },
    sw_op{ enclosing }
{
    PRINT_STATUS( "Executing SweepCellSolveBase<%s>::%s.\n", typeid(OrdinateFlux).name(), __func__ )

    /* empty */
}


//------------------------------------------------------------------------------------------------------------
//! \brief  Virtual destructor.
//------------------------------------------------------------------------------------------------------------
template< class OrdinateFlux >
SweepCellSolveBase<OrdinateFlux>::~SweepCellSolveBase ( void ) {

    PRINT_STATUS( "Executing SweepCellSolveBase<%s>::%s.\n", typeid(OrdinateFlux).name(), __func__ )

    /* empty */
}


//------------------------------------------------------------------------------------------------------------
//! \brief  Sets the timestep size for the sweep, recomputing values if needed.
//------------------------------------------------------------------------------------------------------------
template< class OrdinateFlux >
void SweepCellSolveBase<OrdinateFlux>::SetDt (

    const double dt_in
) {
    PRINT_STATUS( "Executing SweepCellSolveBase<%s>::%s.\n", typeid(OrdinateFlux).name(), __func__ )

    this->dt = dt_in;
}


//------------------------------------------------------------------------------------------------------------
//! \brief  Constructs and solves the local system for a given mesh element.
//!
//! Used to enforce ordering on CalcA, CalcB, etc.
//!
//! \param[in]      A           Pointer to memory location containing matrix for linear system to be solved.
//! \param[in,out]  B           Initially contains the RHS vector for the linear system to be solved.
//!                             Upon return, contains the solution of the system.
//! \param[in]      idx         Contains indices for system(s) to construct and solve.
//! \param[in]      sigma       Total cross section \f$ \sigma_{\mathrm{t}} \f$.
//! \param[in]      source      Contains the source term \f$ Q \f$ of the transport equation.
//! \param[in]      result      Contains the solution of the transport equation.
//! \param[in]      work_ptr    Pointer to temporary workspace used by solver.
//! \param[in]      initial     Pointer to object containing initial condition.
//------------------------------------------------------------------------------------------------------------
template< class OrdinateFlux >
void SweepCellSolveBase<OrdinateFlux>::ConstructAndSolve (

    double * const A,
    double * const B,
    const SIMD_BlkIdx<0> & idx,
    const RKDG::CrossSection & sigma,
    const OrdinateFlux & source,
    const OrdinateFlux & result,
    void * const work_ptr,
    const RKDG::OrdinateFlux * const initial

) const {

    this->CalcB( B, idx, source, result, work_ptr, initial );
    this->CalcA( A, idx, sigma, work_ptr );
    this->SolveSystem( A, B, work_ptr );
    this->StoreResult( B, idx, result );
}


//============================================================================================================
//=== EXPLICIT INSTANTIATIONS ================================================================================
//============================================================================================================

# include "objects/RKDG/OrdinateFlux.hpp"
# include "objects/STDG/OrdinateFlux.hpp"


template class SweepCellSolveBase<RKDG::OrdinateFlux>;
template class SweepCellSolveBase<STDG::OrdinateFlux>;
