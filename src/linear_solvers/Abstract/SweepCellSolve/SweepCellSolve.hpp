//------------------------------------------------------------------------------------------------------------
//! \file   linear_solvers/Abstract/SweepCellSolve/SweepCellSolve.hpp
//! \brief  Header file containing declaration of SweepCellSolve abstract class.
//!
//! \author Michael M. Crockatt
//! \date   February 2018
//------------------------------------------------------------------------------------------------------------

# ifndef __ABSTRACT__SWEEP_CELL_SOLVE_HPP__
# define __ABSTRACT__SWEEP_CELL_SOLVE_HPP__


# include <cstring>
# include <limits>

# include "linear_solvers/Abstract/SweepOperator.hpp"
# include "linear_solvers/Abstract/SweepCellSolve/SweepCellSolveBase.hpp"
# include "linear_solvers/Abstract/SweepCellSolve/SweepCellSolveFactory.hpp"


//------------------------------------------------------------------------------------------------------------
//! \brief  Declaration of SweepCellSolve abstract class.
//!
//! This class is not intended to be instantiated implicitly.
//!
//! This class provides no implementations.
//! Specializations of this class implement:
//!     - Indexing functions for local systems.
//!     - SweepCellSolveBase::StoreResult.
//!
//------------------------------------------------------------------------------------------------------------
template< class OrdinateFlux, int64_t SIMD_length >
class SweepCellSolve :
    public SweepCellSolveBase<OrdinateFlux>
{

public:

    //========================================================================================================
    //=== CONSTRUCTORS AND DESTRUCTOR ========================================================================
    //========================================================================================================

    //--------------------------------------------------------------------------------------------------------
    //! \brief  Delete default constructor.
    //--------------------------------------------------------------------------------------------------------
    SweepCellSolve( void ) = delete;

    //--------------------------------------------------------------------------------------------------------
    //! \brief  Delete copy constructor.
    //--------------------------------------------------------------------------------------------------------
    SweepCellSolve( const SweepCellSolve< OrdinateFlux, SIMD_length > & ) = delete;

    //--------------------------------------------------------------------------------------------------------
    //! \brief  Construct solver object from given parameters.
    //--------------------------------------------------------------------------------------------------------
    SweepCellSolve (
        const SweepOperator<OrdinateFlux> & enclosing,
        const ParameterList & input_list,
        const double dt = std::numeric_limits<double>::infinity()
    );

    //--------------------------------------------------------------------------------------------------------
    //! \brief  Virtual destructor.
    //--------------------------------------------------------------------------------------------------------
    virtual ~SweepCellSolve( void );


    //========================================================================================================
    //=== OPERATOR OVERLOADS =================================================================================
    //========================================================================================================

    //--------------------------------------------------------------------------------------------------------
    //! \brief  Delete copy assignment operator.
    //--------------------------------------------------------------------------------------------------------
    SweepCellSolve< OrdinateFlux, SIMD_length > &
        operator=( const SweepCellSolve< OrdinateFlux, SIMD_length > & ) = delete;


    //========================================================================================================
    //=== OVERRIDES OF INTERFACE FUNCTIONS ===================================================================
    //========================================================================================================

protected:

    //--------------------------------------------------------------------------------------------------------
    //! \brief  Stores the coefficients obtained from the solve into the output object.
    //--------------------------------------------------------------------------------------------------------
    void StoreResult (
        const double * const B,
        const SIMD_BlkIdx<0> & idx,
        const OrdinateFlux & result
    ) const override;

};


//
// Include headers containing specializations.
//
# include "linear_solvers/RKDG/SweepCellSolve/SweepCellSolve.hpp"
# include "linear_solvers/STDG/SweepCellSolve/SweepCellSolve.hpp"


# endif // ifndef __ABSTRACT__SWEEP_CELL_SOLVE_HPP__
