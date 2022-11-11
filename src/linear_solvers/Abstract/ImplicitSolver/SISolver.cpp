//------------------------------------------------------------------------------------------------------------
//! \file   linear_solvers/Abstract/ImplicitSolver/SISolver.cpp
//! \brief  Contains implementations and instantiations of methods from SISolver class template.
//!
//! \author Michael M. Crockatt
//! \date   March 2018
//------------------------------------------------------------------------------------------------------------

# include <limits>
# include <typeinfo>

# include "linear_solvers/Abstract/ImplicitSolver/SISolver.hpp"
# include "operators/TransportOperator.hpp"


//============================================================================================================
//=== MEMBER DEFINITIONS FOR SISolver CLASS ==================================================================
//============================================================================================================


//--------------------------------------------------------------------------------------------------------
//! \brief  Returns the string descriptor of the class type.
//--------------------------------------------------------------------------------------------------------
template< class AngularFlux >
const std::string SISolver<AngularFlux>::Descriptor( void ) {

    return "source-iteration";
}


//--------------------------------------------------------------------------------------------------------
//! \brief  Returns the string descriptor for an objects type.
//--------------------------------------------------------------------------------------------------------
template< class AngularFlux >
const std::string SISolver<AngularFlux>::GetDescriptor( void ) const {

    return Descriptor();
}


//------------------------------------------------------------------------------------------------------------
//! \brief  Initializes an SISolver object.
//!
//! \param[in]  domain_decomposition    Contains parameters of spatial discretization.
//! \param[in]  ordinate_set            Contains parameters of angular discretization.
//! \param[in]  input_list              Contains additional parameters.
//! \param[in]  dt                      Initial timestep size.
//------------------------------------------------------------------------------------------------------------
template< class AngularFlux >
SISolver<AngularFlux>::SISolver (

    const DomainDecomposition & domain_decomposition,
    const Quadrule::OrdinateSet & ordinate_set,
    const ParameterList & input_list,
    const double dt                 // = std::numeric_limits<double>::infinity()
) :
    ImplicitSolver<AngularFlux>( domain_decomposition, ordinate_set, input_list, dt ),

    abs_tol{ 1.0e-8 },
    max_its{ 500 }
{
    PRINT_STATUS( "Executing SISolver<%s>::%s.\n", typeid(AngularFlux).name(), __func__ )

    this->sweep_operator = std::shared_ptr<SweepOperator<AngularFlux>>(
            new SweepOperator<AngularFlux>( domain_decomposition, ordinate_set, input_list, dt )
        );

    SISolver<AngularFlux>::SetParameters( input_list );

    this->psi_ptr = std::shared_ptr<AngularFlux>(
            AngularFlux::Create( domain_decomposition, ordinate_set, input_list )
        );

    this->phi.push_back( std::shared_ptr<ScalarFlux>(
            ScalarFlux::Create( domain_decomposition, input_list )
        ) );

    this->phi.push_back( std::shared_ptr<ScalarFlux>(
            ScalarFlux::Create( domain_decomposition, input_list )
        ) );
}


//------------------------------------------------------------------------------------------------------------
//! \brief  Virtual destructor.
//------------------------------------------------------------------------------------------------------------
template< class AngularFlux >
SISolver<AngularFlux>::~SISolver ( void ) {

    PRINT_STATUS( "Executing SISolver<%s>::%s.\n", typeid(AngularFlux).name(), __func__ )

    /* empty */
}


//------------------------------------------------------------------------------------------------------------
//! \brief  Prints a summary of the object's configuration to the logging interface.
//!
//! \param[in]  prefix  String prepended to each line of output, intended to be used for controlling
//!                     indentation of output.
//------------------------------------------------------------------------------------------------------------
template< class AngularFlux >
void SISolver<AngularFlux>::Print (

    const std::string & prefix // = "  "

) const {

    PRINT_STATUS( "Executing SISolver<%s>::%s.\n", typeid(AngularFlux).name(), __func__ )

    this->ImplicitSolver<AngularFlux>::Print( prefix );

    PRINT_LOG( "%s%-*s  %s\n", prefix.c_str(), Global::col_width, "Solve Type:", this->Descriptor().c_str() )

    PRINT_LOG( "%s%-*s % .4e\n",        prefix.c_str(), Global::col_width, "SI abstol:",  this->abs_tol )
    PRINT_LOG( "%s%-*s % " PRId64 "\n", prefix.c_str(), Global::col_width, "SI max its:", this->max_its )

    this->sweep_operator->Print( prefix );
}


//------------------------------------------------------------------------------------------------------------
//! \brief  Sets parameters given a ParameterList.
//!
//! \param[in]  input_list  Contains parameters to set in the object.
//------------------------------------------------------------------------------------------------------------
template< class AngularFlux >
void SISolver<AngularFlux>::SetParameters (

    const ParameterList & input_list
) {
    PRINT_STATUS( "Executing SISolver<%s>::%s.\n", typeid(AngularFlux).name(), __func__ )

    try { this->abs_tol = input_list.GetValue<double> ( "abs_tol" );  } catch (...) {/* empty */}
    try { this->max_its = input_list.GetValue<int64_t>( "max_its" );  } catch (...) {/* empty */}

    ImplicitSolver<AngularFlux>::SetParameters( input_list );

    this->sweep_operator->SetParameters( input_list );
}


//------------------------------------------------------------------------------------------------------------
//! \brief  Sets the initial guess for the next solve.
//!
//! \param[in]  guess   (optional) <br>
//!                     Pointer to object from which to set the initial guess. If this is null (default) then
//!                     the initial guess is set to the zero vector and this pointer is not dereferenced.
//------------------------------------------------------------------------------------------------------------
template< class AngularFlux >
void SISolver<AngularFlux>::SetInitialGuess (

    const AngularFlux * const guess    // = nullptr
) {
    PRINT_STATUS( "Executing SISolver<%s>::%s.\n", typeid(AngularFlux).name(), __func__ )

    if ( guess == nullptr ) {

        PRINT_STATUS( "Setting initial guess to zero.\n" )

        this->psi_ptr->ZeroDensity();

    // Otherwise use guess.
    } else {

        PRINT_STATUS( "Setting initial guess from given object.\n" )

        this->psi_ptr->Copy( *guess );
    }
}


//------------------------------------------------------------------------------------------------------------
//! \brief  Solves a system of equations using the specified parameters.
//!
//! \param[in,out]  source      Initially contains the coefficients of the source term \f$ Q \f$ for
//!                             the transport system. Upon return, contains the coefficients \f$ \Psi \f$ that
//!                             satisfy the transport equation.
//! \param[in]      sigma_t     Total cross section \f$ \sigma_{\mathrm{t}} \f$.
//! \param[in]      sigma_s     Scattering cross section \f$ \sigma_{\mathrm{s}} \f$.
//! \param[in]      dt          (optional) <br>
//!                             Timestep size. <br>
//!                             For solving steady-state problems on should use
//!                             <tt>\pp{dt} = inf</tt> (default value).
//! \param[in]      initial     (optional) <br>
//!                             Pointer to object containing the initial condition for the timestep.
//!                             This pointer is not dereferenced if it is null (default value).
//------------------------------------------------------------------------------------------------------------
template< class AngularFlux >
void SISolver<AngularFlux>::Solve (

    AngularFlux & source,
    const CrossSection & sigma_t,
    const CrossSection & sigma_s,
    const double dt,                        // = std::numeric_limits<double>::infinity(),
    const AngularFluxInit * const initial   // = nullptr
) {

    PRINT_STATUS( "Executing SISolver<%s>::%s.\n", typeid(AngularFlux).name(), __func__ )

    AngularFlux & psi = *this->psi_ptr;

    ScalarFlux * phi_current = this->phi.at(0).get();
    ScalarFlux * phi_last    = this->phi.at(1).get();

    int64_t iter = 0;
    double error = std::numeric_limits<double>::max();

    phi_current->ZeroDensity();
    phi_last->ZeroDensity();

    // Set initial scalar flux.
    TransportOperator::Pmv( 1.0, 0.0, psi, *phi_current, OpDomain::Interior );

    while (     iter < this->max_its
            &&  error > this->abs_tol
    ) {
        ++iter;
        std::swap( phi_current, phi_last );

        // --- Compute next iterate. ------------------------------------------------------------------ //

        if ( psi.UpwindHalosAreDirty() ) {  psi.SynchronizeHalos();  }

        // Q.
        AngularFlux::Copy( psi, source, OpDomain::Interior | OpDomain::Boundary );

        // SPΨⁿ + Q
        TransportOperator::Smv( 1.0, 1.0, sigma_s, *phi_last, psi, OpDomain::Interior );

        // L⁻¹ ( SPΨⁿ + Q ).
        this->sweep_operator->Apply( psi, sigma_t, dt, initial );

        // --- Measure error between successive iterates. --------------------------------------------- //

        TransportOperator::Pmv( 1.0, 0.0, psi, *phi_current, OpDomain::Interior );

        Global::TMR_PETSc.Start();
        error = ScalarFlux::l2Dist( *phi_current, *phi_last );
        Global::TMR_PETSc.Stop();

    # if LOGLEVEL >= 2
        PRINT_LOG( "In source iteration solve, iteration = %" PRId64 "  error = % .4e\n", iter, error )
    # endif
    }

    // --- Source iterations finished. ---------------------------------------------------------------- //

    AngularFlux::Copy( source, psi );

# if LOGLEVEL >= 1
    PRINT_LOG( "Iterations: %5" PRId64"       Final Error: %9.4e\n", iter, error )
# endif
}


//============================================================================================================
//=== EXPLICIT INSTANTIATIONS ================================================================================
//============================================================================================================

# include "linear_solvers/Abstract/ImplicitSolver/ImplicitSolverDerivedFactory.hpp"
# include "objects/RKDG/OrdinateFlux.hpp"
# include "objects/STDG/OrdinateFlux.hpp"


//!
//! \brief  Instantiation of SISolver<RKDG::OrdinateFlux> class.
//!
template class SISolver<RKDG::OrdinateFlux>;

//!
//! \brief  Instantiation of factory for SISolver<RKDG::OrdinateFlux>.
//!
template class
ImplicitSolver<RKDG::OrdinateFlux>::ImplicitSolverDerivedFactory<
        SISolver<RKDG::OrdinateFlux>
    >;


//!
//! \brief  Instantiation of SISolver<STDG::OrdinateFlux> class.
//!
template class SISolver<STDG::OrdinateFlux>;

//!
//! \brief  Instantiation of factory for SISolver<STDG::OrdinateFlux>.
//!
template class
ImplicitSolver<STDG::OrdinateFlux>::ImplicitSolverDerivedFactory<
        SISolver<STDG::OrdinateFlux>
    >;
