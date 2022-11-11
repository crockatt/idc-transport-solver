//------------------------------------------------------------------------------------------------------------
//! \file   linear_solvers/Abstract/SweepCellSolve/SweepCellSolveManager.cpp
//! \brief  Contains implementations and instantiations of methods from SweepCellSolveManager class template.
//!
//! \author Michael M. Crockatt
//! \date   March 2018
//------------------------------------------------------------------------------------------------------------

# include <string>
# include <typeinfo>

# include "linear_solvers/Abstract/SweepCellSolve/SweepCellSolveFactory.hpp"
# include "linear_solvers/Abstract/SweepCellSolve/SweepCellSolveManager.hpp"



//============================================================================================================
//=== CONSTRUCTORS AND DESTRUCTOR ============================================================================
//============================================================================================================


//------------------------------------------------------------------------------------------------------------
//! \brief  Template used to initialize SweepCellSolveManager::cell_solvers.
//------------------------------------------------------------------------------------------------------------
template< class OrdinateFlux >
template< int64_t... max_length >
constexpr auto SweepOperator<OrdinateFlux>::SweepCellSolveManager::SetupCellSolvers (

    const SweepOperator<OrdinateFlux> & enclosing,
    const ParameterList & input_list,
    std::integer_sequence< int64_t, max_length... >,
    const double dt
) {

    return std::array< std::shared_ptr< SweepCellSolveBase<OrdinateFlux> >, sizeof...(max_length) > {{
            (      max_length > 0
                && 1 << IntLog2(max_length) == max_length
              ? std::shared_ptr< SweepCellSolveBase<OrdinateFlux> >(
                    SweepCellSolveFactory< OrdinateFlux, max_length >
                        ::CreateCellSolver( enclosing, input_list, dt )
                )
              : nullptr
            ) ...
        }};
}


//------------------------------------------------------------------------------------------------------------
//! \brief  Initializes a SweepCellSolveManager object.
//!
//! \param[in]  enclosing   Reference to instance of enclosing SweepOperator class.
//! \param[in]  input_list  Contains parameters for initialization of object.
//! \param[in]  dt          Timestep size.
//------------------------------------------------------------------------------------------------------------
template< class OrdinateFlux >
SweepOperator<OrdinateFlux>::SweepCellSolveManager::SweepCellSolveManager (

    const SweepOperator<OrdinateFlux> & enclosing,
    const ParameterList & input_list,
    const double dt                                 // = std::numeric_limits<double>::infinity()
) :
    sw_op{ enclosing },
    cell_solvers{ SetupCellSolvers( enclosing, input_list,
                                    std::make_integer_sequence< int64_t, SWEEP_SIMD_LEN + 1 >{}, dt
                ) }
{
    PRINT_STATUS( "Executing SweepOperator<%s>::SweepCellSolveManager::%s.\n",
                  typeid(OrdinateFlux).name(), __func__ )
}


//------------------------------------------------------------------------------------------------------------
//! \brief  Destructor for SweepCellSolveManager classes.
//------------------------------------------------------------------------------------------------------------
template< class OrdinateFlux >
SweepOperator<OrdinateFlux>::SweepCellSolveManager::~SweepCellSolveManager ( void ) {

    PRINT_STATUS( "Executing SweepOperator<%s>::SweepCellSolveManager::%s.\n",
                  typeid(OrdinateFlux).name(), __func__ )

    /* empty */
}


//========================================================================================================
//=== PUBLIC INTERFACE ROUTINES ==========================================================================
//========================================================================================================


//------------------------------------------------------------------------------------------------------------
//! \brief  Prints a summary of the object's configuration to the logging interface.
//------------------------------------------------------------------------------------------------------------
template< class OrdinateFlux >
void SweepOperator<OrdinateFlux>::SweepCellSolveManager::Print (

    const std::string & prefix // = "  "

) const {

    PRINT_STATUS( "Executing SweepOperator<%s>::SweepCellSolveManager::%s.\n",
                  typeid(OrdinateFlux).name(), __func__ )

    // Solver with vector length one always exists.
    this->cell_solvers[1]->Print( prefix );
}


//------------------------------------------------------------------------------------------------------------
//! \brief  Returns the dimension of the array to be allocated for A.
//------------------------------------------------------------------------------------------------------------
template< class OrdinateFlux >
size_t SweepOperator<OrdinateFlux>::SweepCellSolveManager::GetDimofA ( void ) const {

    PRINT_STATUS( "Executing SweepOperator<%s>::SweepCellSolveManager::%s.\n",
                  typeid(OrdinateFlux).name(), __func__ )

    // Solver with vector length one always exists.
    return SWEEP_SIMD_LEN * this->cell_solvers[1]->GetDimofA();
}


//------------------------------------------------------------------------------------------------------------
//! \brief  Returns the dimension of the array to be allocated for B.
//------------------------------------------------------------------------------------------------------------
template< class OrdinateFlux >
size_t SweepOperator<OrdinateFlux>::SweepCellSolveManager::GetDimofB ( void ) const {

    PRINT_STATUS( "Executing SweepOperator<%s>::SweepCellSolveManager::%s.\n",
                  typeid(OrdinateFlux).name(), __func__ )

    // Solver with vector length one always exists.
    return SWEEP_SIMD_LEN * this->cell_solvers[1]->GetDimofB();
}


//------------------------------------------------------------------------------------------------------------
//! \brief  Returns the dimension of the array to be allocated for extra work space.
//------------------------------------------------------------------------------------------------------------
template< class OrdinateFlux >
size_t SweepOperator<OrdinateFlux>::SweepCellSolveManager::GetDimofWork ( void ) const {

    PRINT_STATUS( "Executing SweepOperator<%s>::SweepCellSolveManager::%s.\n",
                  typeid(OrdinateFlux).name(), __func__ )

    // Solver with vector length one always exists.
    return SWEEP_SIMD_LEN * this->cell_solvers[1]->GetDimofWork();
}


//------------------------------------------------------------------------------------------------------------
//! \brief  Sets the timestep size for the sweep, recomputing values if needed.
//------------------------------------------------------------------------------------------------------------
template< class OrdinateFlux >
void SweepOperator<OrdinateFlux>::SweepCellSolveManager::SetDt (

    const double dt
){

    PRINT_STATUS( "Executing SweepOperator<%s>::SweepCellSolveManager::%s.\n",
                  typeid(OrdinateFlux).name(), __func__ )

    // Update all solvers with lengths equal to a power of two.
    for ( int64_t i = 1; i <= SWEEP_SIMD_LEN; i <<= 1 )
        this->cell_solvers[i]->SetDt( dt );
}


//------------------------------------------------------------------------------------------------------------
//! \brief  Constructs and solves the local system for a given mesh element.
//------------------------------------------------------------------------------------------------------------
template< class OrdinateFlux >
void SweepOperator<OrdinateFlux>::SweepCellSolveManager::ConstructAndSolve (

    double * const A,
    double * const B,
    const SIMD_BlkIdx<0> & idx,
    const RKDG::CrossSection & sigma,
    const OrdinateFlux & source,
    const OrdinateFlux & result,
    void * const work_ptr,
    const RKDG::OrdinateFlux * const initial

) const {

    PRINT_STATUS( "Executing SweepOperator<%s>::SweepCellSolveManager::%s.\n",
                  typeid(OrdinateFlux).name(), __func__ )

    const SweepCellSolveBase<OrdinateFlux> & solver = *this->cell_solvers[ idx.len ];

    solver.ConstructAndSolve( A, B, idx, sigma, source, result, work_ptr, initial );
}


//============================================================================================================
//=== EXPLICIT INSTANTIATIONS ================================================================================
//============================================================================================================

# include "objects/RKDG/OrdinateFlux.hpp"
# include "objects/STDG/OrdinateFlux.hpp"


template class SweepOperator<RKDG::OrdinateFlux>::SweepCellSolveManager;
template class SweepOperator<STDG::OrdinateFlux>::SweepCellSolveManager;
