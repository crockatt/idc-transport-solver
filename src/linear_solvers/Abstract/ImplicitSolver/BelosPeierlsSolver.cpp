//------------------------------------------------------------------------------------------------------------
//! \file   linear_solvers/Abstract/ImplicitSolver/BelosPeierlsSolver.cpp
//! \brief  Contains implementations and instantiations of methods from BelosPeierlsSolver class template.
//!
//! \author Michael M. Crockatt
//! \date   March 2018
//------------------------------------------------------------------------------------------------------------

# if defined (ENABLE_BELOS_TPETRA) || defined (DOXYCOMPILE)


# include <limits>
# include <typeinfo>

# include <BelosTpetraAdapter.hpp>
# include <BelosPseudoBlockGmresSolMgr.hpp>

# include "linear_solvers/Abstract/ImplicitSolver/BelosPeierlsSolver.hpp"
# include "operators/TransportOperator.hpp"
# include "utils/BelosTpetraInterface.hpp"


//============================================================================================================
//=== MEMBER DEFINITIONS FOR BelosPeierlsSolver CLASS =====================================================
//============================================================================================================


//--------------------------------------------------------------------------------------------------------
//! \brief  Returns the string descriptor of the class type.
//--------------------------------------------------------------------------------------------------------
template< class AngularFlux >
const std::string BelosPeierlsSolver<AngularFlux>::Descriptor( void ) {

    return "belos-peierls-gmres";
}


//--------------------------------------------------------------------------------------------------------
//! \brief  Returns the string descriptor for an objects type.
//--------------------------------------------------------------------------------------------------------
template< class AngularFlux >
const std::string BelosPeierlsSolver<AngularFlux>::GetDescriptor( void ) const {

    return Descriptor();
}


//------------------------------------------------------------------------------------------------------------
//! \brief  Initializes a BelosPeierlsSolver object.
//!
//! \param[in]  domain_decomposition    Contains parameters of spatial discretization.
//! \param[in]  ordinate_set            Contains parameters of angular discretization.
//! \param[in]  input_list              Contains additional parameters.
//! \param[in]  dt                      Initial timestep size.
//------------------------------------------------------------------------------------------------------------
template< class AngularFlux >
BelosPeierlsSolver<AngularFlux>::BelosPeierlsSolver (

    const DomainDecomposition & domain_decomposition,
    const Quadrule::OrdinateSet & ordinate_set,
    const ParameterList & input_list,
    const double dt                 // = std::numeric_limits<double>::infinity()
) :
    ImplicitSolver<AngularFlux>( domain_decomposition, ordinate_set, input_list, dt ),

    dt{ std::nan("0") },
    sigma_t{ nullptr },
    sigma_s{ nullptr }
{
    PRINT_STATUS( "Executing BelosPeierlsSolver<%s>::%s.\n", typeid(AngularFlux).name(), __func__ )

    // Setup sweep operator.
    this->sweep_operator = std::shared_ptr<SweepOperator<AngularFlux>>(
            new SweepOperator<AngularFlux>( domain_decomposition, ordinate_set, input_list, dt )
        );

    // Setup temporary objects.
    this->psi = std::shared_ptr<AngularFlux>(
            AngularFlux::Create( domain_decomposition, ordinate_set, input_list )
        );

    this->phi.push_back( std::shared_ptr<ScalarFlux>(
            ScalarFlux::Create( domain_decomposition, input_list )
        ) );

    // Setup Trilinos operator and vectors.
    this->trilinos_mat = Teuchos::rcp( new IPLS( *this ) );

    this->trilinos_vec_x = Teuchos::rcp( new TpetraVec( this->trilinos_mat->getDomainMap(), 1 ) );
    this->trilinos_vec_b = Teuchos::rcp( new TpetraVec( this->trilinos_mat->getRangeMap(),  1 ) );

    this->trilinos_vec_x->putScalar( 0.0 );
    this->trilinos_vec_b->putScalar( 0.0 );

    // Setup Belos::LinearProblem.
    this->belos_problem = Teuchos::rcp(
            new Belos::LinearProblem< TpetraScalar, TpetraVec, TpetraOperator >()
        );

    this->belos_problem->setOperator( this->trilinos_mat );

    // Set default solver options.
    this->solver_params = Teuchos::parameterList();

//     this->solver_params->set( "Verbosity", Belos::Errors | Belos::Warnings | Belos::IterationDetails );
//     this->solver_params->set( "Output Frequency", 1 );

    this->solver_params->set( "Implicit Residual Scaling", "None" );
    this->solver_params->set( "Explicit Residual Scaling", "None" );

    this->solver_params->set( "Convergence Tolerance", 1.0e-8 );
    this->solver_params->set( "Maximum Iterations",    2000   );
    this->solver_params->set( "Num Blocks",            30     );

    // Create Belos solver object.
    this->belos_solver = Teuchos::rcp(
            new Belos::PseudoBlockGmresSolMgr
                < TpetraScalar, TpetraVec, TpetraOperator >
                ( this->belos_problem, this->solver_params )
        );

    BelosPeierlsSolver<AngularFlux>::SetParameters( input_list );
}


//------------------------------------------------------------------------------------------------------------
//! \brief  Virtual destructor.
//------------------------------------------------------------------------------------------------------------
template< class AngularFlux >
BelosPeierlsSolver<AngularFlux>::~BelosPeierlsSolver ( void ) {

    PRINT_STATUS( "Executing BelosPeierlsSolver<%s>::%s.\n", typeid(AngularFlux).name(), __func__ )
}


//------------------------------------------------------------------------------------------------------------
//! \brief  Prints a summary of the object's configuration to the logging interface.
//!
//! \param[in]  prefix  String prepended to each line of output, intended to be used for controlling
//!                     indentation of output.
//------------------------------------------------------------------------------------------------------------
template< class AngularFlux >
void BelosPeierlsSolver<AngularFlux>::Print (

    const std::string & prefix // = "  "

) const {

    PRINT_STATUS( "Executing BelosPeierlsSolver<%s>::%s.\n", typeid(AngularFlux).name(), __func__ )

    this->ImplicitSolver<AngularFlux>::Print( prefix );

    PRINT_LOG( "%s%-*s  %s\n", prefix.c_str(), Global::col_width, "Solve Type:", this->Descriptor().c_str() )

    const double abstol  = this->solver_params->template get<double>( "Convergence Tolerance" );
    const int    maxits  = this->solver_params->template get<int>   ( "Maximum Iterations"    );
    const int    restart = this->solver_params->template get<int>   ( "Num Blocks"            );

    PRINT_LOG( "%s%-*s % .4e\n", prefix.c_str(), Global::col_width, "GMRES tol:",     abstol  )
    PRINT_LOG( "%s%-*s % d\n",   prefix.c_str(), Global::col_width, "GMRES maxits:",  maxits  )
    PRINT_LOG( "%s%-*s % d\n",   prefix.c_str(), Global::col_width, "GMRES restart:", restart )

    PRINT_LOG( "%s%-*s  %s\n", prefix.c_str(), Global::col_width, "Preconditioner:", "none" )

    this->sweep_operator->Print( prefix );
}


//------------------------------------------------------------------------------------------------------------
//! \brief  Sets parameters given a ParameterList.
//!
//! \param[in]  input_list  Contains parameters to set in the object.
//------------------------------------------------------------------------------------------------------------
template< class AngularFlux >
void BelosPeierlsSolver<AngularFlux>::SetParameters (

    const ParameterList & input_list
) {
    PRINT_STATUS( "Executing BelosPeierlsSolver<%s>::%s.\n", typeid(AngularFlux).name(), __func__ )

    try {  this->solver_params->set( "Convergence Tolerance",
                                     input_list.GetValue<double>( "abs_tol" ) );
    } catch (...) {/* empty */}

    try {  this->solver_params->set( "Maximum Iterations",
                                     input_list.GetValue<int>( "max_its" ) );
    } catch (...) {/* empty */}

    try {  this->solver_params->set( "Num Blocks",
                                     input_list.GetValue<int>( "gmres_restart" ) );
    } catch (...) {/* empty */}

    this->belos_solver->setParameters( this->solver_params );

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
void BelosPeierlsSolver<AngularFlux>::SetInitialGuess (

    const AngularFlux * const guess    // = nullptr
) {
    PRINT_STATUS( "Executing BelosPeierlsSolver<%s>::%s.\n", typeid(AngularFlux).name(), __func__ )

    if ( guess == nullptr ) {

        PRINT_STATUS( "Setting initial guess to zero.\n" )

        this->trilinos_vec_x->putScalar( 0.0 );

    // Otherwise use guess.
    } else {

        PRINT_STATUS( "Setting initial guess from given object.\n" )

        ScalarFlux & phi = *this->phi.at(0);

        TransportOperator::Pmv( 1.0, 0.0, *guess, phi );
        phi.PackExternalVector( *this->trilinos_vec_x, 1.0, 0.0 );
    }
}


//------------------------------------------------------------------------------------------------------------
//! \brief  Solves a system of equations using the specified parameters.
//!
//! \param[in,out]  source              Initially contains the coefficients of the source term \f$ Q \f$ for
//!                                     the transport system. Upon return, contains the coefficients
//!                                     \f$ \Psi \f$ that satisfy the transport equation.
//! \param[in]      sigma_t             Total cross section \f$ \sigma_{\mathrm{t}} \f$.
//! \param[in]      sigma_s             Scattering cross section \f$ \sigma_{\mathrm{s}} \f$.
//! \param[in]      dt                  (optional) <br>
//!                                     Timestep size. <br>
//!                                     For solving steady-state problems on should use
//!                                     <tt>\pp{dt} = inf</tt> (default value).
//! \param[in]      initial_condition   (optional) <br>
//!                                     Pointer to object containing the initial condition for the timestep.
//!                                     This pointer is not dereferenced if it is null (default value).
//------------------------------------------------------------------------------------------------------------
template< class AngularFlux >
void BelosPeierlsSolver<AngularFlux>::Solve (

    AngularFlux & source,
    const CrossSection & sigma_t,
    const CrossSection & sigma_s,
    const double dt,                                    // = std::numeric_limits<double>::infinity(),
    const AngularFluxInit * const initial_condition     // = nullptr
) {

    PRINT_STATUS( "Executing BelosPeierlsSolver<%s>::%s.\n", typeid(AngularFlux).name(), __func__ )

    // --- Initialization. ---------------------------------------------------------------------------- //

    // Alias solve parameters for simplicity.
    AngularFlux & psi = *this->psi;
    ScalarFlux & phi = *this->phi.at(0);

    // Save parameters for solve.
    this->dt = dt;
    this->sigma_t = &sigma_t;
    this->sigma_s = &sigma_s;

    // --- Compute right hand side of system. --------------------------------------------------------- //

    // L⁻¹Q
    psi.Copy( source );
    this->sweep_operator->Apply( psi, sigma_t, dt, initial_condition );

    // PL⁻¹Q
    TransportOperator::Pmv( 1.0, 0.0, psi, phi );
    phi.PackExternalVector( *this->trilinos_vec_b, 1.0, 0.0 );

    // --- Solve for the scalar flux. ----------------------------------------------------------------- //

    Global::TMR_PETSc.Start();

    // Solve system.
    this->belos_problem->setProblem( this->trilinos_vec_x, this->trilinos_vec_b );
    const Belos::ReturnType result = this->belos_solver->solve();

    // Get information about the solve.
# if LOGLEVEL >= 1

    const int num_its = this->belos_solver->getNumIters();
    const double final_norm = this->belos_solver->achievedTol();

    PRINT_LOG( "Iterations = %5d \tFinal Residual = % 10.4e\n", num_its, final_norm )

# endif // if LOGLEVEL >= 1

    Global::TMR_PETSc.Stop();

    if ( result != Belos::Converged ) {

        std::string error_message = "Trilinos GMRES transport solve failed to converge.\n";

        PRINT_ERROR( error_message.c_str() )
        throw std::runtime_error( error_message );
    }

    // --- Compute angular flux from scalar flux. ----------------------------------------------------- //

    phi.UnpackExternalVector( *this->trilinos_vec_x );

    // L⁻¹(Q + SPΨ).
    TransportOperator::Smv( 1.0, 1.0, sigma_s, phi, source );
    this->sweep_operator->Apply( source, sigma_t, dt, initial_condition );

    // --- Cleanup. ----------------------------------------------------------------------------------- //

    // Reset solve parameters.
    this->dt = std::nan("0");
    this->sigma_t = nullptr;
    this->sigma_s = nullptr;
}


//============================================================================================================
//=== MEMBER DEFINITIONS FOR BelosPeierlsSolver::IPLS CLASS ===============================================
//============================================================================================================


//------------------------------------------------------------------------------------------------------------
//! \brief  Constructs an IPLS object from the given parameters.
//!
//! \note   This routine assumes that the enclosing object \pp{solver} has already initialized
//!         BelosPeierlsSolver::phi with the appropriate objects.
//!
//! \param[in]  solver          Reference to enclosing solver object.
//------------------------------------------------------------------------------------------------------------
template< class AngularFlux >
BelosPeierlsSolver<AngularFlux>::IPLS::IPLS (

    const BelosPeierlsSolver<AngularFlux> & solver
) :
    enc{ solver }
{

    PRINT_STATUS( "Executing BelosPeierlsSolver<%s>::IPLS::%s.\n", typeid(AngularFlux).name(), __func__ )

    Teuchos::RCP<const Teuchos::Comm<int>> comm = Teuchos::rcp(
        # if defined (ENABLE_MPI)
            new Teuchos::MpiComm<int>( Global::MPI_cart_comm )
        # else
            new Teuchos::SerialComm<int>()
        # endif
        );

    domain_decomp_map = Teuchos::rcp(
            new TpetraMap( solver.phi.at(0)->GetLocalArrayDim(),
                           solver.phi.at(0)->GetGlobalArrayDim(),
                           0, comm
        ) );
}


//------------------------------------------------------------------------------------------------------------
//! \brief  Computes
//!         \f$ \vec{y} \gets \beta \vec{y} + \alpha (\mathcal{I} - \mathcal{PL}^{-1}\mathcal{S}) \vec{x} \f$.
//!
//! \param[in]      X       Input vector \f$ x \f$.
//! \param[in,out]  Y       Output vector \f$ y \f$.
//! \param[in]      mode    Transpose mode applied to operator. Only Teuchos::NO_TRANS is supported.
//! \param[in]      alpha   The scalar \f$ \alpha \f$.
//! \param[in]      beta    The scalar \f$ \beta \f$.
//!
//------------------------------------------------------------------------------------------------------------
template< class AngularFlux >
void BelosPeierlsSolver<AngularFlux>::IPLS::apply (

    const TpetraVec & X,
    TpetraVec & Y,
    Teuchos::ETransp mode,  // = Teuchos::NO_TRANS
    TpetraScalar alpha,     // = Teuchos::ScalarTraits<TpetraScalar>::one()
    TpetraScalar beta       // = Teuchos::ScalarTraits<TpetraScalar>::zero()

) const {

    PRINT_STATUS( "Executing BelosPeierlsSolver<%s>::IPLS::%s.\n", typeid(AngularFlux).name(), __func__ )

# if defined (STRICT_CHECK)

    if ( mode != Teuchos::NO_TRANS ) {

        std::string error_message =   "Unsupported transpose mode '"
                                    + std::to_string( mode )
                                    + "' in '"
                                    + std::string(__func__)
                                    + "'.\n";

        PRINT_ERROR( error_message.c_str() )
        throw std::invalid_argument( error_message );
    }

# else // if defined (STRICT_CHECK)
    (void) mode;
# endif // if defined (STRICT_CHECK)

    Global::TMR_PETSc.Stop();

    const double dt = this->enc.dt;
    const CrossSection & sigma_t = *this->enc.sigma_t;
    const CrossSection & sigma_s = *this->enc.sigma_s;

    AngularFlux & psi = *this->enc.psi;
    ScalarFlux & phi = *this->enc.phi.at(0);

    phi.UnpackExternalVector( X );

    // (I - PL⁻¹S)(PΨ)
    TransportOperator::Smv( 1.0, 0.0, sigma_s, phi, psi );
    this->enc.sweep_operator->Apply( psi, sigma_t, dt );
    TransportOperator::Pmv( -1.0, 1.0, psi, phi );

    phi.PackExternalVector( Y, alpha, beta );

    Global::TMR_PETSc.Start();
}


//============================================================================================================
//=== EXPLICIT INSTANTIATIONS ================================================================================
//============================================================================================================

# include "linear_solvers/Abstract/ImplicitSolver/ImplicitSolverDerivedFactory.hpp"
# include "objects/RKDG/OrdinateFlux.hpp"
# include "objects/STDG/OrdinateFlux.hpp"


//!
//! \brief  Instantiation of BelosPeierlsSolver<RKDG::OrdinateFlux> class.
//!
template class BelosPeierlsSolver<RKDG::OrdinateFlux>;

//!
//! \brief  Instantiation of factory for BelosPeierlsSolver<RKDG::OrdinateFlux>.
//!
template class
ImplicitSolver<RKDG::OrdinateFlux>::ImplicitSolverDerivedFactory<
        BelosPeierlsSolver<RKDG::OrdinateFlux>
    >;


//!
//! \brief  Instantiation of BelosPeierlsSolver<STDG::OrdinateFlux> class.
//!
template class BelosPeierlsSolver<STDG::OrdinateFlux>;

//!
//! \brief  Instantiation of factory for BelosPeierlsSolver<STDG::OrdinateFlux>.
//!
template class
ImplicitSolver<STDG::OrdinateFlux>::ImplicitSolverDerivedFactory<
        BelosPeierlsSolver<STDG::OrdinateFlux>
    >;


# endif // if defined (ENABLE_PETSC)
