//------------------------------------------------------------------------------------------------------------
//! \file   linear_solvers/Abstract/SweepCellSolve/SweepCellSolveBase.hpp
//! \brief  Header file containing declaration of SweepCellSolveBase class template.
//!
//! \author Michael M. Crockatt
//! \date   March 2018
//------------------------------------------------------------------------------------------------------------

# ifndef __ABSTRACT__SWEEP_CELL_SOLVE_BASE_HPP__
# define __ABSTRACT__SWEEP_CELL_SOLVE_BASE_HPP__


# include <cstddef>
# include <limits>
# include <string>

# include "linear_solvers/Abstract/SweepOperator.hpp"
# include "objects/RKDG/CrossSection.hpp"
# include "utils/SIMD.hpp"


//------------------------------------------------------------------------------------------------------------
//! \brief  Abstract class outlining local solve methods for discrete ordinates sweep solvers.
//!
//! This class defines the interface used by the SweepCellSolveManager class for solving the local systems
//! within each mesh element during transport sweeps.
//!
//! Constructors are declared protected and nested class SweepCellSolveFactory is declared public. This is
//! done to force construction of SweepCellSolve objects through the factory interface.
//!
//! \attention  The interface defined by this class is intended to be executed in the following order:
//!                 - SweepCellSolveBase::CalcB
//!                 - SweepCellSolveBase::CalcA
//!                 - SweepCellSolveBase::SolveSystem
//!                 - SweepCellSolveBase::StoreResult
//!
//! \see    SweepCellSolve
//! \see    SweepLinearSolve
//! \see    SweepLapackLUSolve
//! \see    SweepGESolve
//! \see    SweepHessSolve
//------------------------------------------------------------------------------------------------------------
template< class OrdinateFlux >
class SweepCellSolveBase {

public:

    //========================================================================================================
    //=== CONSTRUCTORS AND DESTRUCTOR ========================================================================
    //========================================================================================================

    //--------------------------------------------------------------------------------------------------------
    //! \brief  Delete default constructor.
    //--------------------------------------------------------------------------------------------------------
    SweepCellSolveBase( void ) = delete;

    //--------------------------------------------------------------------------------------------------------
    //! \brief  Delete copy constructor.
    //--------------------------------------------------------------------------------------------------------
    SweepCellSolveBase( const SweepCellSolveBase & ) = delete;

    //--------------------------------------------------------------------------------------------------------
    //! \brief  Initializes a SweepCellSolveBase object.
    //--------------------------------------------------------------------------------------------------------
    SweepCellSolveBase (
        const SweepOperator<OrdinateFlux> & enclosing,
        const double dt_in = std::numeric_limits<double>::infinity()
    );

    //--------------------------------------------------------------------------------------------------------
    //! \brief  Destructor for SweepCellSolveBase classes.
    //--------------------------------------------------------------------------------------------------------
    virtual ~SweepCellSolveBase( void );


    //========================================================================================================
    //=== OPERATOR OVERLOADS =================================================================================
    //========================================================================================================

    //--------------------------------------------------------------------------------------------------------
    //! \brief  Delete copy assignment operator.
    //--------------------------------------------------------------------------------------------------------
    SweepCellSolveBase & operator=( const SweepCellSolveBase & ) = delete;


    //========================================================================================================
    //=== PUBLIC INTERFACE ROUTINES ==========================================================================
    //========================================================================================================

    //--------------------------------------------------------------------------------------------------------
    //! \brief  Returns the string descriptor for an objects type.
    //--------------------------------------------------------------------------------------------------------
    virtual const std::string GetDescriptor( void ) const = 0;

    //--------------------------------------------------------------------------------------------------------
    //! \brief  Prints a summary of the object's configuration to the logging interface.
    //--------------------------------------------------------------------------------------------------------
    virtual void Print( const std::string & = "  " ) const = 0;

    //--------------------------------------------------------------------------------------------------------
    //! \brief  Returns the dimension of the array to be allocated for A.
    //--------------------------------------------------------------------------------------------------------
    virtual size_t GetDimofA( void ) const = 0;

    //--------------------------------------------------------------------------------------------------------
    //! \brief  Returns the dimension of the array to be allocated for B.
    //--------------------------------------------------------------------------------------------------------
    virtual size_t GetDimofB( void ) const = 0;

    //--------------------------------------------------------------------------------------------------------
    //! \brief  Returns the dimension of the array to be allocated for extra work space.
    //--------------------------------------------------------------------------------------------------------
    virtual size_t GetDimofWork( void ) const = 0;

    //--------------------------------------------------------------------------------------------------------
    //! \brief  Sets the timestep size for the sweep, recomputing values if needed.
    //--------------------------------------------------------------------------------------------------------
    virtual void SetDt( const double dt_in );

    //--------------------------------------------------------------------------------------------------------
    //! \brief  Constructs and solves the local system for a given mesh element.
    //--------------------------------------------------------------------------------------------------------
    virtual void ConstructAndSolve (
        double * const A,
        double * const B,
        const SIMD_BlkIdx<0> & idx,
        const RKDG::CrossSection & sigma,
        const OrdinateFlux & source,
        const OrdinateFlux & result,
        void * const work_ptr,
        const RKDG::OrdinateFlux * const initial
    ) const;


protected:

    //========================================================================================================
    //=== PROTECTED HELPER FUNCTIONS =========================================================================
    //========================================================================================================

    //--------------------------------------------------------------------------------------------------------
    //! \brief  Computes the elements of the systems to be solved at the given step of a sweep.
    //--------------------------------------------------------------------------------------------------------
    virtual void CalcA (
        double * const A,
        const SIMD_BlkIdx<0> & idx,
        const RKDG::CrossSection & sigma,
        void * const work_ptr
    ) const = 0;

    //--------------------------------------------------------------------------------------------------------
    //! \brief  Computes the elements of the right-hand vectors of the systems to be solved at the given step
    //!         of a sweep.
    //--------------------------------------------------------------------------------------------------------
    virtual void CalcB (
        double * const B,
        const SIMD_BlkIdx<0> & idx,
        const OrdinateFlux & source,
        const OrdinateFlux & result,
        void * const work_ptr,
        const RKDG::OrdinateFlux * const initial = nullptr
    ) const = 0;

    //--------------------------------------------------------------------------------------------------------
    //! \brief  Solves the system for the updated coefficients.
    //--------------------------------------------------------------------------------------------------------
    virtual void SolveSystem (
        double * const A,
        double * const B,
        void * const work_ptr
    ) const = 0;

    //--------------------------------------------------------------------------------------------------------
    //! \brief  Stores the coefficients obtained from the solve into the output object.
    //--------------------------------------------------------------------------------------------------------
    virtual void StoreResult (
        const double * const B,
        const SIMD_BlkIdx<0> & idx,
        const OrdinateFlux & result
    ) const = 0;


    //========================================================================================================
    //=== INTERNAL VARIABLES =================================================================================
    //========================================================================================================

    double dt;                  //!< Timestep size for current implicit solve.
    int64_t system_dim;         //!< Dimension of the local linear systems solved during the sweep algorithm.

    //!
    //! \brief  Reference to instance of enclosing SweepOperator class.
    //!
    const SweepOperator<OrdinateFlux> & sw_op;

};


# endif // ifndef __ABSTRACT__SWEEP_CELL_SOLVE_BASE_HPP__
