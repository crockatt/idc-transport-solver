//------------------------------------------------------------------------------------------------------------
//! \file   linear_solvers/Abstract/ImplicitSolver.cpp
//! \brief  Contains implementations and instantiations of methods from ImplicitSolver class template.
//!
//! \author Michael M. Crockatt
//! \date   January 2018
//------------------------------------------------------------------------------------------------------------

# include <typeinfo>

# include "linear_solvers/Abstract/ImplicitSolver.hpp"
# include "linear_solvers/Abstract/ImplicitSolver/SweepSolver.hpp"


//============================================================================================================
//=== MEMBER DEFINITIONS FOR ImplicitSolver CLASS ============================================================
//============================================================================================================


//------------------------------------------------------------------------------------------------------------
//! Initializes an ImplicitSolver object.
//------------------------------------------------------------------------------------------------------------
template< class AngularFlux >
ImplicitSolver<AngularFlux>::ImplicitSolver (

    const DomainDecomposition & domain_decomposition,
    const Quadrule::OrdinateSet & ordinate_set,
    const ParameterList &,
    const double                    // = std::numeric_limits<double>::infinity()
) :
    DomainDecomposition( domain_decomposition ),
    Quadrule::OrdinateSet( ordinate_set )
{
    PRINT_STATUS( "Executing ImplicitSolver<%s>::%s.\n", typeid(AngularFlux).name(), __func__ )

    /* empty */
}


//------------------------------------------------------------------------------------------------------------
//! \brief  Virtual destructor.
//------------------------------------------------------------------------------------------------------------
template< class AngularFlux >
ImplicitSolver<AngularFlux>::~ImplicitSolver ( void ) {

    PRINT_STATUS( "Executing ImplicitSolver<%s>::%s.\n", typeid(AngularFlux).name(), __func__ )

    /* empty */
}


//------------------------------------------------------------------------------------------------------------
//! \brief  Prints a summary of the object's configuration to the logging interface.
//!
//! \param[in]  prefix  String prepended to each line of output, intended to be used for controlling
//!                     indentation of output.
//------------------------------------------------------------------------------------------------------------
template< class AngularFlux >
void ImplicitSolver<AngularFlux>::Print (

    const std::string & prefix // = "  "

) const {

    PRINT_STATUS( "Executing ImplicitSolver<%s>::%s.\n", typeid(AngularFlux).name(), __func__ )

    this->Quadrule::OrdinateSet::Print( prefix );
}


//------------------------------------------------------------------------------------------------------------
//! \brief  Sets parameters given a ParameterList.
//------------------------------------------------------------------------------------------------------------
template< class AngularFlux >
void ImplicitSolver<AngularFlux>::SetParameters (

    const ParameterList &
) {
    PRINT_STATUS( "Executing ImplicitSolver<%s>::%s.\n", typeid(AngularFlux).name(), __func__ )

    /* empty */
}


//------------------------------------------------------------------------------------------------------------
//! \brief  Sets the initial guess for the next solve.
//!
//! Base implementation takes no action (suitable for non-iterative solvers).
//------------------------------------------------------------------------------------------------------------
template< class AngularFlux >
void ImplicitSolver<AngularFlux>::SetInitialGuess (

    const AngularFlux * const  // = nullptr
) {
    PRINT_STATUS( "Executing ImplicitSolver<%s>::%s.\n", typeid(AngularFlux).name(), __func__ )

    /* empty */
}


//============================================================================================================
//=== MEMBER DEFINITIONS FOR ImplicitSolverFactory CLASS =====================================================
//============================================================================================================


//!
//! \brief  Definition of ImplicitSolverFactory::factory_map_ptr.
//!
template< class AngularFlux >
std::map< std::string, typename ImplicitSolver<AngularFlux>::ImplicitSolverFactory * > *
ImplicitSolver<AngularFlux>::ImplicitSolverFactory::factory_map_ptr{ nullptr };


//------------------------------------------------------------------------------------------------------------
//! \brief  Returns a reference to the map stored at ImplicitSolverFactory::factory_map_ptr.
//!
//! Responsible for allocating the map if ImplicitSolverFactory::factory_map_ptr is null.
//------------------------------------------------------------------------------------------------------------
template< class AngularFlux >
std::map< std::string, typename ImplicitSolver<AngularFlux>::ImplicitSolverFactory * > &
ImplicitSolver<AngularFlux>::ImplicitSolverFactory::GetFactoryMap ( void ) {

    if ( factory_map_ptr == nullptr )
        factory_map_ptr = new std::map< std::string, ImplicitSolverFactory * >();

    return *factory_map_ptr;
}


//------------------------------------------------------------------------------------------------------------
//! \brief  Adds a reference to a factory class to the map ImplicitSolverFactory::factory_map.
//!
//! Factory classes deriving from ImplicitSolverFactory should call this routine to add a mapping to their
//! singleton instance from the string descriptor of the ImplicitSolver subclass which they construct.
//!
//! \note   Because this function uses the \c insert routine of std::map each factory may only be added once.
//!
//! \param[in]  descriptor  String descriptor of the ImplicitSolver associated with the factory \pp{factory}.
//! \param[in]  factory     Factory object for creating objects of types deriving from ImplicitSolver.
//------------------------------------------------------------------------------------------------------------
template< class AngularFlux >
void ImplicitSolver<AngularFlux>::ImplicitSolverFactory::AddFactory (

    const std::string & descriptor,
    ImplicitSolverFactory & factory
) {

    PRINT_STATUS( "Executing ImplicitSolver<%s>::ImplicitSolverFactory::%s.\n",
                  typeid(AngularFlux).name(), __func__ )

    auto & factory_map = this->GetFactoryMap();

    if ( factory_map.count( descriptor ) ) {

        std::string error_message =   "Factory for ImplicitSolver '"
                                    + descriptor
                                    + "' already present in ImplicitSolver<"
                                    + std::string( typeid(AngularFlux).name() )
                                    + ">::ImplicitSolverFactory::factory_map.\n";

        PRINT_ERROR( error_message.c_str() )
        throw std::runtime_error( error_message );
    }

    factory_map.insert( std::pair< std::string, ImplicitSolverFactory * >( descriptor, &factory ) );
}


//------------------------------------------------------------------------------------------------------------
//! \brief  Calls creator for desired ImplicitSolver implementation.
//!
//! Uses the map ImplicitSolverFactory::factory_map to determine which factory to call to create object.
//!
//! \param[in]  domain_decomposition    Contains parameters of spatial discretization.
//! \param[in]  ordinate_set            Contains parameters of angular discretization.
//! \param[in]  input_list              Contains additional parameters.
//! \param[in]  dt                      Initial timestep size.
//!
//! \return     Returns a pointer to the created object.
//------------------------------------------------------------------------------------------------------------
template< class AngularFlux >
ImplicitSolver<AngularFlux> * ImplicitSolver<AngularFlux>::ImplicitSolverFactory::CreateSolver (

    const DomainDecomposition & domain_decomposition,
    const Quadrule::OrdinateSet & ordinate_set,
    const ParameterList & input_list,
    const double dt                 // = std::numeric_limits<double>::infinity()
) {

    PRINT_STATUS( "Executing ImplicitSolver<%s>::ImplicitSolverFactory::%s.\n",
                  typeid(AngularFlux).name(), __func__ )

    std::string solver_str;

    // Get string descriptor from input list.
    try {

        solver_str = input_list.GetValue<std::string>( "solve_type" );

    } catch (...) {

        PRINT_WARNING( "Using default solve type.\n" )
        solver_str = "default";
    }

    ImplicitSolverFactory * solver_factory = nullptr;

    // Retrieve factory pointer with given string descriptor.
    try {

        auto & factory_map = GetFactoryMap();
        solver_factory = factory_map.at( solver_str );

    } catch (...) {

        std::string error_message =   "Failed to access factory for ImplicitSolver '"
                                    + solver_str
                                    + "' in ImplicitSolver<"
                                    + std::string( typeid(AngularFlux).name() )
                                    + ">::ImplicitSolverFactory::"
                                    + std::string(__func__)
                                    + ".\n";

        throw std::out_of_range( error_message );
    }

    // Create desired object using factory.
    return solver_factory->Create( domain_decomposition, ordinate_set, input_list, dt );
}


//============================================================================================================
//=== EXPLICIT INSTANTIATIONS ================================================================================
//============================================================================================================

# include "objects/RKDG/OrdinateFlux.hpp"
# include "objects/STDG/OrdinateFlux.hpp"


//!
//! \brief  Instantiation of ImplicitSolver<RKDG::OrdinateFlux> class.
//!
template class ImplicitSolver<RKDG::OrdinateFlux>;


//!
//! \brief  Instantiation of ImplicitSolver<STDG::OrdinateFlux> class.
//!
template class ImplicitSolver<STDG::OrdinateFlux>;
