//------------------------------------------------------------------------------------------------------------
//! \file   linear_solvers/Abstract/ImplicitSolver/ImplicitSolverDerivedFactory.hpp
//! \brief  Header file containing declaration of ImplicitSolverDerivedFactory class template.
//!
//! \author Michael M. Crockatt
//! \date   February 2018
//------------------------------------------------------------------------------------------------------------

# ifndef __ABSTRACT__IMPLICIT_SOLVER_DERIVED_FACTORY_HPP__
# define __ABSTRACT__IMPLICIT_SOLVER_DERIVED_FACTORY_HPP__


//============================================================================================================
//=== CLASS DECLARATION ======================================================================================
//============================================================================================================


//------------------------------------------------------------------------------------------------------------
//! \brief  Class template for factories of classes derived from ImplicitSolver.
//!
//! Implements the interface defined by ImplicitSolverFactory.
//!
//! Let \c DerivedSolver denote a class that has been derived from ImplicitSolver. Use of this template
//! to register \c DerivedSolver in ImplicitSolverFactory::factory_map so that objects of type
//! \c DerivedSolver may be created requires two steps:
//!     1. Declare a specialization of this template for \c DerivedSolver.
//!     2. Perform an explicit instantiation of this specialization in the implementation file for
//!        \c DerivedSolver.
//!
//------------------------------------------------------------------------------------------------------------
template< class AngularFlux >
template< class DerivedSolver >
class ImplicitSolver<AngularFlux>::ImplicitSolverDerivedFactory :
    public ImplicitSolver<AngularFlux>::ImplicitSolverFactory
{

    //========================================================================================================
    //=== CONSTRUCTORS AND DESTRUCTOR ========================================================================
    //========================================================================================================

protected:

    ImplicitSolverDerivedFactory( void );

    ImplicitSolverDerivedFactory( const ImplicitSolverDerivedFactory<DerivedSolver> & ) = delete;

    ImplicitSolverDerivedFactory<DerivedSolver> &
        operator=( const ImplicitSolverDerivedFactory<DerivedSolver> & ) = delete;

public:

    ~ImplicitSolverDerivedFactory( void ) = default;


protected:

    //========================================================================================================
    //=== PROTECTED INTERFACE ROUTINES =======================================================================
    //========================================================================================================

    //!
    //! \brief  Inherit ImplicitSolverFactory::AddFactory.
    //!
    using ImplicitSolver<AngularFlux>::ImplicitSolverFactory::AddFactory;

    //--------------------------------------------------------------------------------------------------------
    //! \brief  Creates and returns an object of type SweepGESolve.
    //--------------------------------------------------------------------------------------------------------
    ImplicitSolver * Create (
        const DomainDecomposition & domain_decomposition,
        const Quadrule::OrdinateSet & ordinate_set,
        const ParameterList & input_list,
        const double dt = std::numeric_limits<double>::infinity()
    ) const override;


    //========================================================================================================
    //=== INTERNAL VARIABLES =================================================================================
    //========================================================================================================

    //!
    //! \brief  Singleton instance.
    //!
    static ImplicitSolverDerivedFactory<DerivedSolver> instance;
};


//============================================================================================================
//=== DEFINITION OF MEMBER FUNCTIONS =========================================================================
//============================================================================================================


//------------------------------------------------------------------------------------------------------------
//! \brief  Default constructor.
//!
//! Adds singleton object to ImplicitSolverFactory::factory_map using
//! ImplicitSolverDerivedFactory::Descriptor().
//!
//------------------------------------------------------------------------------------------------------------
template< class AngularFlux >
template< class DerivedSolver >
ImplicitSolver<AngularFlux>::ImplicitSolverDerivedFactory<DerivedSolver>::
ImplicitSolverDerivedFactory ( void ) {

    PRINT_STATUS( "Executing ImplicitSolver<%s>::ImplicitSolverDerivedFactory<%s>::%s.\n",
                  typeid(AngularFlux).name(), typeid(DerivedSolver).name(), __func__ )

    AddFactory( DerivedSolver::Descriptor(), this->instance );
}


//------------------------------------------------------------------------------------------------------------
//! \brief  Creates and returns an object of type derived from ImplicitSolver.
//!
//! \param[in]  domain_decomposition    Contains parameters of spatial discretization.
//! \param[in]  ordinate_set            Contains parameters of angular discretization.
//! \param[in]  input_list              Contains additional parameters.
//! \param[in]  dt                      Initial timestep size.
//!
//! \return     Returns a pointer to the created object.
//------------------------------------------------------------------------------------------------------------
template< class AngularFlux >
template< class DerivedSolver >
ImplicitSolver<AngularFlux> *
ImplicitSolver<AngularFlux>::ImplicitSolverDerivedFactory<DerivedSolver>::Create (

    const DomainDecomposition & domain_decomposition,
    const Quadrule::OrdinateSet & ordinate_set,
    const ParameterList & input_list,
    const double dt                 // = std::numeric_limits<double>::infinity()

) const {

    PRINT_STATUS( "Executing ImplicitSolver<%s>::ImplicitSolverDerivedFactory<%s>::%s.\n",
                  typeid(AngularFlux).name(), typeid(DerivedSolver).name(), __func__ )

    return new DerivedSolver( domain_decomposition, ordinate_set, input_list, dt );
}


//============================================================================================================
//=== INSTANTIATION OF STATIC MEMBERS ========================================================================
//============================================================================================================


//!
//! \brief  Define static singleton instance ImplicitSolverDerivedFactory::instance.
//!
template< class AngularFlux >
template< class DerivedSolver >
typename ImplicitSolver<AngularFlux>::template ImplicitSolverDerivedFactory<DerivedSolver>
ImplicitSolver<AngularFlux>::ImplicitSolverDerivedFactory<DerivedSolver>::instance;


# endif // ifndef __ABSTRACT__IMPLICIT_SOLVER_DERIVED_FACTORY_HPP__
