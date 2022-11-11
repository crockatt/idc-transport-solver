//------------------------------------------------------------------------------------------------------------
//! \file   linear_solvers/Abstract/SweepCellSolve/SweepCellSolveFactory.hpp
//! \brief  Header file containing declaration of SweepCellSolveFactory class template.
//!
//! \author Michael M. Crockatt
//! \date   March 2018
//------------------------------------------------------------------------------------------------------------

# ifndef __ABSTRACT__SWEEP_CELL_SOLVE_FACTORY_HPP__
# define __ABSTRACT__SWEEP_CELL_SOLVE_FACTORY_HPP__


# include "linear_solvers/Abstract/SweepCellSolve/SweepCellSolveBase.hpp"


//------------------------------------------------------------------------------------------------------------
//! \brief  Abstract class defining interface for factories used to construct objects derived from
//!         SweepCellSolve.
//!
//! Concrete classes derived from SweepCellSolve should use the class template
//! SweepCellSolveDerivedFactory (which derives from this class) to enable their construction.
//!
//------------------------------------------------------------------------------------------------------------
template< class OrdinateFlux, int64_t SIMD_length >
class SweepCellSolveFactory {

    //========================================================================================================
    //=== CONSTRUCTORS AND DESTRUCTOR ========================================================================
    //========================================================================================================

protected:

    SweepCellSolveFactory( void ) = default;

    SweepCellSolveFactory( const SweepCellSolveFactory< OrdinateFlux, SIMD_length > & ) = delete;

    SweepCellSolveFactory< OrdinateFlux, SIMD_length > &
        operator=( const SweepCellSolveFactory< OrdinateFlux, SIMD_length > & ) = delete;

public:

    ~SweepCellSolveFactory( void ) = default;


    //========================================================================================================
    //=== PUBLIC INTERFACE ROUTINES ==========================================================================
    //========================================================================================================

    //--------------------------------------------------------------------------------------------------------
    //! \brief  Creates and returns an object of a type derived from SweepCellSolveBase.
    //--------------------------------------------------------------------------------------------------------
    static SweepCellSolveBase<OrdinateFlux> * CreateCellSolver (
        const SweepOperator<OrdinateFlux> & enclosing,
        const ParameterList & input_list,
        const double dt = std::numeric_limits<double>::infinity()
    );


protected:

    //========================================================================================================
    //=== PROTECTED INTERFACE ROUTINES =======================================================================
    //========================================================================================================

    //--------------------------------------------------------------------------------------------------------
    //! \brief  Returns a reference to the map stored at ImplicitSolverFactory::factory_map_ptr.
    //--------------------------------------------------------------------------------------------------------
    static std::map< std::string, SweepCellSolveFactory< OrdinateFlux, SIMD_length > * > &
        GetFactoryMap( void );

    //--------------------------------------------------------------------------------------------------------
    //! \brief  Adds a reference to a factory class to the map SweepCellSolveFactory::factory_map_ptr.
    //--------------------------------------------------------------------------------------------------------
    void AddFactory (
        const std::string & descriptor,
        SweepCellSolveFactory< OrdinateFlux, SIMD_length > & factory
    );

    //--------------------------------------------------------------------------------------------------------
    //! \brief  Virtual interface for derived factory classes to create objects of their respective types.
    //--------------------------------------------------------------------------------------------------------
    virtual SweepCellSolveBase<OrdinateFlux> * Create (
        const SweepOperator<OrdinateFlux> & enclosing,
        const ParameterList & input_list,
        const double dt = std::numeric_limits<double>::infinity()
    ) const = 0;


    //========================================================================================================
    //=== INTERNAL VARIABLES =================================================================================
    //========================================================================================================

    //!
    //! \brief  Pointer to map for mapping between string descriptors and factories for classes deriving from
    //!         SweepCellSolve.
    //!
    //! Factories are implemented as singletons and are responsible for adding themselves to this map
    //! through the SweepCellSolveFactory::AddFactory method.
    //!
    static std::map< std::string, SweepCellSolveFactory< OrdinateFlux, SIMD_length > * > * factory_map_ptr;

};


//
// Include header for SweepCellSolveDerivedFactory ***AFTER*** declaration of SweepCellSolveFactory class.
//
# include "linear_solvers/Abstract/SweepCellSolve/SweepCellSolveDerivedFactory.hpp"


# endif // ifndef __ABSTRACT__SWEEP_CELL_SOLVE_FACTORY_HPP__
