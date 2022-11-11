//------------------------------------------------------------------------------------------------------------
//! \file   linear_solvers/Abstract/SweepCellSolve/SweepCellSolveDerivedFactory.hpp
//! \brief  Header file containing declaration of SweepCellSolveDerivedFactory class template.
//!
//! \author Michael M. Crockatt
//! \date   February 2018
//------------------------------------------------------------------------------------------------------------

# ifndef __ABSTRACT__SWEEP_CELL_SOLVE_DERIVED_FACTORY_HPP__
# define __ABSTRACT__SWEEP_CELL_SOLVE_DERIVED_FACTORY_HPP__


//============================================================================================================
//=== CLASS DECLARATION ======================================================================================
//============================================================================================================


//------------------------------------------------------------------------------------------------------------
//! \brief  Class template for factories of classes derived from SweepCellSolveBase.
//!
//! Implements the interface defined by SweepCellSolveFactory.
//!
//! Let \c DerivedCellSolve denote a class that has been derived from SweepCellSolveBase. Use of this template
//! to register \c DerivedCellSolve in SweepCellSolveFactory::factory_map so that objects of type
//! \c DerivedCellSolve may be created requires two steps:
//!     1. Declare a specialization of this template for \c DerivedCellSolve.
//!     2. Perform an explicit instantiation of this specialization in the implementation file for
//!        \c DerivedCellSolve.
//!
//------------------------------------------------------------------------------------------------------------
template<
    class OrdinateFlux,
    int64_t SIMD_length,
    template< typename, int64_t > class DerivedCellSolve
>
class SweepCellSolveDerivedFactory :
    public SweepCellSolveFactory< OrdinateFlux, SIMD_length >
{

    //========================================================================================================
    //=== CONSTRUCTORS AND DESTRUCTOR ========================================================================
    //========================================================================================================

protected:

    SweepCellSolveDerivedFactory( void );

    SweepCellSolveDerivedFactory(
        const SweepCellSolveDerivedFactory< OrdinateFlux, SIMD_length, DerivedCellSolve > &
    ) = delete;

    SweepCellSolveDerivedFactory< OrdinateFlux, SIMD_length, DerivedCellSolve > &
        operator=(
            const SweepCellSolveDerivedFactory< OrdinateFlux, SIMD_length, DerivedCellSolve > &
        ) = delete;

public:

    ~SweepCellSolveDerivedFactory( void ) = default;


protected:

    //========================================================================================================
    //=== PROTECTED INTERFACE ROUTINES =======================================================================
    //========================================================================================================

    //!
    //! \brief  Inherit SweepCellSolveFactory::AddFactory.
    //!
    using SweepCellSolveFactory< OrdinateFlux, SIMD_length >::AddFactory;

    //--------------------------------------------------------------------------------------------------------
    //! \brief  Creates and returns an object of a type derived from SweepCellSolveBase.
    //--------------------------------------------------------------------------------------------------------
    SweepCellSolveBase<OrdinateFlux> * Create (
        const SweepOperator<OrdinateFlux> & enclosing,
        const ParameterList & input_list,
        const double dt = std::numeric_limits<double>::infinity()
    ) const override;


    //========================================================================================================
    //=== INTERNAL VARIABLES =================================================================================
    //========================================================================================================

    //!
    //! \brief  Singleton instance.
    //!
    static SweepCellSolveDerivedFactory< OrdinateFlux, SIMD_length, DerivedCellSolve > instance;
};


//============================================================================================================
//=== DEFINITION OF MEMBER FUNCTIONS =========================================================================
//============================================================================================================


//------------------------------------------------------------------------------------------------------------
//! \brief  Default constructor.
//!
//! Adds singleton object to SweepCellSolveFactory::factory_map using
//! SweepCellSolveDerivedFactory::Descriptor().
//!
//------------------------------------------------------------------------------------------------------------
template<
    class OrdinateFlux,
    int64_t SIMD_length,
    template< typename, int64_t > class DerivedCellSolve
>
SweepCellSolveDerivedFactory< OrdinateFlux, SIMD_length, DerivedCellSolve >
::SweepCellSolveDerivedFactory ( void ) {

    PRINT_STATUS( "Executing SweepCellSolveFactory<%s,%" PRId64 ",%s>::%s.\n",
                  typeid(OrdinateFlux).name(), SIMD_length,
                  typeid( DerivedCellSolve< OrdinateFlux, SIMD_length > ).name(), __func__ )

    AddFactory( DerivedCellSolve< OrdinateFlux, SIMD_length >::Descriptor(), this->instance );
}


//------------------------------------------------------------------------------------------------------------
//! \brief  Creates and returns an object of type derived from SweepCellSolveBase.
//!
//! \param[in]  enclosing   Reference to instance of enclosing SweepOperator class.
//! \param[in]  input_list  Contains parameters for initialization of object.
//! \param[in]  dt          Initial timestep size.
//!
//! \return     Returns a pointer to the created object.
//------------------------------------------------------------------------------------------------------------
template<
    class OrdinateFlux,
    int64_t SIMD_length,
    template< typename, int64_t > class DerivedCellSolve
>
SweepCellSolveBase<OrdinateFlux> *
SweepCellSolveDerivedFactory< OrdinateFlux, SIMD_length, DerivedCellSolve >
::Create (

    const SweepOperator<OrdinateFlux> & enclosing,
    const ParameterList & input_list,
    const double dt                 // = std::numeric_limits<double>::infinity()

) const {

    PRINT_STATUS( "Executing SweepCellSolveFactory<%s,%" PRId64 ",%s>::%s.\n",
                  typeid(OrdinateFlux).name(), SIMD_length,
                  typeid( DerivedCellSolve< OrdinateFlux, SIMD_length > ).name(), __func__ )

    return new DerivedCellSolve< OrdinateFlux, SIMD_length >( enclosing, input_list, dt );
}


//============================================================================================================
//=== INSTANTIATION OF STATIC MEMBERS ========================================================================
//============================================================================================================


//!
//! \brief  Define static singleton instance SweepCellSolveDerivedFactory::instance.
//!
template<
    class OrdinateFlux,
    int64_t SIMD_length,
    template< typename, int64_t > class DerivedCellSolve
>
SweepCellSolveDerivedFactory< OrdinateFlux, SIMD_length, DerivedCellSolve >
SweepCellSolveDerivedFactory< OrdinateFlux, SIMD_length, DerivedCellSolve >::instance;


# endif // ifndef __ABSTRACT__SWEEP_CELL_SOLVE_DERIVED_FACTORY_HPP__
