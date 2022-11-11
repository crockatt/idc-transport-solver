//------------------------------------------------------------------------------------------------------------
//! \file   linear_solvers/Abstract/SweepCellSolve/SweepCellSolveFactory.cpp
//! \brief  Contains implementations and instantiations of methods from SweepCellSolveFactory class template.
//!
//! \author Michael M. Crockatt
//! \date   March 2018
//------------------------------------------------------------------------------------------------------------

# include "linear_solvers/Abstract/SweepCellSolve/SweepCellSolveFactory.hpp"


//!
//! \brief  Definition of SweepCellSolveFactory::factory_map_ptr.
//!
template< class OrdinateFlux, int64_t SIMD_length >
std::map<
    std::string,
    SweepCellSolveFactory< OrdinateFlux, SIMD_length > *
> *
SweepCellSolveFactory< OrdinateFlux, SIMD_length >::factory_map_ptr{ nullptr };


//------------------------------------------------------------------------------------------------------------
//! \brief  Returns a reference to the map stored at SweepCellSolveFactory::factory_map_ptr.
//!
//! Responsible for allocating the map if SweepCellSolveFactory::factory_map_ptr is null.
//------------------------------------------------------------------------------------------------------------
template< class OrdinateFlux, int64_t SIMD_length >
std::map<
    std::string,
    SweepCellSolveFactory< OrdinateFlux, SIMD_length > *
> &
SweepCellSolveFactory< OrdinateFlux, SIMD_length >::GetFactoryMap ( void ) {

    if ( factory_map_ptr == nullptr )
        factory_map_ptr = new std::map< std::string, SweepCellSolveFactory< OrdinateFlux, SIMD_length > * >();

    return *factory_map_ptr;
}


//------------------------------------------------------------------------------------------------------------
//! \brief  Adds a reference to a factory class to the map SweepCellSolveFactory::factory_map_ptr.
//!
//! Factory classes deriving from SweepCellSolveFactory should call this routine to add a mapping to their
//! singleton instance from the string descriptor of the SweepCellSolve subclass which they construct.
//!
//! \note   Because this function uses the \c insert routine of std::map each factory may only be added once.
//!
//! \param[in]  descriptor  String descriptor of the SweepCellSolve associated with the factory \pp{factory}.
//! \param[in]  factory     Factory object for creating objects of types deriving from SweepCellSolve.
//------------------------------------------------------------------------------------------------------------
template< class OrdinateFlux, int64_t SIMD_length >
void SweepCellSolveFactory< OrdinateFlux, SIMD_length >::AddFactory (

    const std::string & descriptor,
    SweepCellSolveFactory< OrdinateFlux, SIMD_length > & factory
) {

    PRINT_STATUS( "Executing SweepCellSolveFactory<%s,% " PRId64 ">::%s.\n",
                  typeid(OrdinateFlux).name(), SIMD_length, __func__ )

    auto & factory_map = this->GetFactoryMap();

    if ( factory_map.count( descriptor ) ) {

        std::string error_message =   "Factory for SweepCellSolve '"
                                    + descriptor
                                    + "' already present in SweepCellSolveFactory<"
                                    + std::string( typeid(OrdinateFlux).name() )
                                    + ","
                                    + std::to_string( SIMD_length )
                                    + ">::factory_map.\n";

        PRINT_ERROR( error_message.c_str() )
        throw std::runtime_error( error_message );
    }

    factory_map.insert( std::pair< std::string,
                                   SweepCellSolveFactory< OrdinateFlux, SIMD_length > *
                                 >( descriptor, &factory ) );
}


//------------------------------------------------------------------------------------------------------------
//! \brief  Calls creator for desired SweepCellSolve implementation.
//!
//! Uses the map SweepCellSolveFactory::factory_map to determine which factory to call to create object.
//!
//! \param[in]  enclosing   Reference to instance of enclosing SweepOperator class.
//! \param[in]  input_list  Contains parameters for initialization of object.
//! \param[in]  dt          Timestep size.
//!
//! \return     Returns a pointer to the created object.
//------------------------------------------------------------------------------------------------------------
template< class OrdinateFlux, int64_t SIMD_length >
SweepCellSolveBase<OrdinateFlux> *
SweepCellSolveFactory< OrdinateFlux, SIMD_length >::CreateCellSolver (

    const SweepOperator<OrdinateFlux> & enclosing,
    const ParameterList & input_list,
    const double dt                 // = std::numeric_limits<double>::infinity()
) {

    PRINT_STATUS( "Executing SweepCellSolveFactory<%s,% " PRId64 ">::%s.\n",
                  typeid(OrdinateFlux).name(), SIMD_length, __func__ )

    std::string sweep_solve_str;

    // Get string descriptor from input list.
    try {

        sweep_solve_str = input_list.GetValue<std::string>( "sweep_solve" );

    } catch (...) {

        PRINT_WARNING( "Using default sweep solve.\n" )
        sweep_solve_str = "default";
    }

    SweepCellSolveFactory< OrdinateFlux, SIMD_length > * sweep_solve_factory = nullptr;

    // Retrieve factory pointer for given string descriptor.
    try {

        auto & factory_map = GetFactoryMap();
        sweep_solve_factory = factory_map.at( sweep_solve_str );

    } catch (...) {

        std::string error_message =   "Failed to access factory for SweepCellSolve '"
                                    + sweep_solve_str
                                    + "' in SweepCellSolveFactory<"
                                    + std::string( typeid(OrdinateFlux).name() )
                                    + ","
                                    + std::to_string( SIMD_length )
                                    + ">::"
                                    + std::string(__func__)
                                    + ".\n";

        throw std::out_of_range( error_message );
    }

    // Create desired object using factory.
    return sweep_solve_factory->Create( enclosing, input_list, dt );
}


//============================================================================================================
//=== EXPLICIT INSTANTIATIONS ================================================================================
//============================================================================================================

# include "objects/RKDG/OrdinateFlux.hpp"
# include "objects/STDG/OrdinateFlux.hpp"


# if SWEEP_SIMD_LEN >= (1 << 0)
    template class SweepCellSolveFactory< RKDG::OrdinateFlux, (1 << 0) >;
    template class SweepCellSolveFactory< STDG::OrdinateFlux, (1 << 0) >;
# endif
# if SWEEP_SIMD_LEN >= (1 << 1)
    template class SweepCellSolveFactory< RKDG::OrdinateFlux, (1 << 1) >;
    template class SweepCellSolveFactory< STDG::OrdinateFlux, (1 << 1) >;
# endif
# if SWEEP_SIMD_LEN >= (1 << 2)
    template class SweepCellSolveFactory< RKDG::OrdinateFlux, (1 << 2) >;
    template class SweepCellSolveFactory< STDG::OrdinateFlux, (1 << 2) >;
# endif
# if SWEEP_SIMD_LEN >= (1 << 3)
    template class SweepCellSolveFactory< RKDG::OrdinateFlux, (1 << 3) >;
    template class SweepCellSolveFactory< STDG::OrdinateFlux, (1 << 3) >;
# endif
# if SWEEP_SIMD_LEN >= (1 << 4)
    template class SweepCellSolveFactory< RKDG::OrdinateFlux, (1 << 4) >;
    template class SweepCellSolveFactory< STDG::OrdinateFlux, (1 << 4) >;
# endif
# if SWEEP_SIMD_LEN >= (1 << 5)
    template class SweepCellSolveFactory< RKDG::OrdinateFlux, (1 << 5) >;
    template class SweepCellSolveFactory< STDG::OrdinateFlux, (1 << 5) >;
# endif
# if SWEEP_SIMD_LEN > (1 << 5)
    # error "Maximum SIMD vector length supported is (1 << 5)."
# endif

