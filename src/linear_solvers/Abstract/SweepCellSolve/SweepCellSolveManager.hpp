//------------------------------------------------------------------------------------------------------------
//! \file   linear_solvers/Abstract/SweepCellSolve/SweepCellSolveManager.hpp
//! \brief  Header file containing declaration of SweepCellSolveManager class template.
//!
//! \author Michael M. Crockatt
//! \date   March 2018
//------------------------------------------------------------------------------------------------------------

# ifndef __ABSTRACT__SWEEP_CELL_SOLVE_MANAGER_HPP__
# define __ABSTRACT__SWEEP_CELL_SOLVE_MANAGER_HPP__


# include <array>
# include <memory>
# include <utility>

# include "linear_solvers/Abstract/SweepOperator.hpp"
# include "linear_solvers/Abstract/SweepCellSolve/SweepCellSolveBase.hpp"
# include "objects/RKDG/CrossSection.hpp"
# include "utils/ParameterList.hpp"


//------------------------------------------------------------------------------------------------------------
//! \brief  Abstract class used for managing local solve methods for discrete ordinates sweep solvers.
//------------------------------------------------------------------------------------------------------------
template< class OrdinateFlux >
class SweepOperator<OrdinateFlux>::SweepCellSolveManager {

public:

    //========================================================================================================
    //=== CONSTRUCTORS AND DESTRUCTOR ========================================================================
    //========================================================================================================

    //--------------------------------------------------------------------------------------------------------
    //! \brief  Delete default constructor.
    //--------------------------------------------------------------------------------------------------------
    SweepCellSolveManager( void ) = delete;

    //--------------------------------------------------------------------------------------------------------
    //! \brief  Delete copy constructor.
    //--------------------------------------------------------------------------------------------------------
    SweepCellSolveManager( const SweepCellSolveManager & ) = delete;

    //--------------------------------------------------------------------------------------------------------
    //! \brief  Initializes a SweepCellSolveManager object.
    //--------------------------------------------------------------------------------------------------------
    SweepCellSolveManager (
        const SweepOperator<OrdinateFlux> & enclosing,
        const ParameterList & input_list,
        const double dt = std::numeric_limits<double>::infinity()
    );

    //--------------------------------------------------------------------------------------------------------
    //! \brief  Destructor for SweepCellSolveManager classes.
    //--------------------------------------------------------------------------------------------------------
    ~SweepCellSolveManager( void );


    //========================================================================================================
    //=== OPERATOR OVERLOADS =================================================================================
    //========================================================================================================

    //--------------------------------------------------------------------------------------------------------
    //! \brief  Delete copy assignment operator.
    //--------------------------------------------------------------------------------------------------------
    SweepCellSolveManager & operator=( const SweepCellSolveManager & ) = delete;


    //========================================================================================================
    //=== PUBLIC INTERFACE ROUTINES ==========================================================================
    //========================================================================================================

    //--------------------------------------------------------------------------------------------------------
    //! \brief  Prints a summary of the object's configuration to the logging interface.
    //--------------------------------------------------------------------------------------------------------
    void Print( const std::string & = "  " ) const;

    //--------------------------------------------------------------------------------------------------------
    //! \brief  Returns the dimension of the array to be allocated for A.
    //--------------------------------------------------------------------------------------------------------
    size_t GetDimofA( void ) const;

    //--------------------------------------------------------------------------------------------------------
    //! \brief  Returns the dimension of the array to be allocated for B.
    //--------------------------------------------------------------------------------------------------------
    size_t GetDimofB( void ) const;

    //--------------------------------------------------------------------------------------------------------
    //! \brief  Returns the dimension of the array to be allocated for extra work space.
    //--------------------------------------------------------------------------------------------------------
    size_t GetDimofWork( void ) const;

    //--------------------------------------------------------------------------------------------------------
    //! \brief  Sets the timestep size for the sweep, recomputing values if needed.
    //--------------------------------------------------------------------------------------------------------
    void SetDt( const double dt );

    //--------------------------------------------------------------------------------------------------------
    //! \brief  Constructs and solves the local system for a given mesh element.
    //--------------------------------------------------------------------------------------------------------
    void ConstructAndSolve (
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
    //! \brief  Template used to initialize SweepCellSolveManager::cell_solvers.
    //--------------------------------------------------------------------------------------------------------
    template< int64_t... max_length >
    static constexpr auto SetupCellSolvers (
        const SweepOperator<OrdinateFlux> & enclosing,
        const ParameterList & input_list,
        std::integer_sequence< int64_t, max_length... >,
        const double dt = std::numeric_limits<double>::infinity()
    );


    //========================================================================================================
    //=== INTERNAL VARIABLES =================================================================================
    //========================================================================================================

    //!
    //! \brief  Reference to instance of enclosing SweepOperator class.
    //!
    const SweepOperator<OrdinateFlux> & sw_op;

    //!
    //! \brief  Array containing pointers to Create static member functions of concrete implementations
    //!         derived from SweepCellSolveBase.
    //!
    std::array< std::shared_ptr< SweepCellSolveBase<OrdinateFlux> >, SWEEP_SIMD_LEN +1 > cell_solvers;

};


# endif // ifndef __ABSTRACT__SWEEP_CELL_SOLVE_MANAGER_HPP__
