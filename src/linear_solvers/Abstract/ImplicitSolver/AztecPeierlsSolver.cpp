//------------------------------------------------------------------------------------------------------------
//! \file   linear_solvers/Abstract/ImplicitSolver/AztecPeierlsSolver.cpp
//! \brief  Contains implementations and instantiations of methods from AztecPeierlsSolver class template.
//!
//! \author Michael M. Crockatt
//! \date   March 2018
//------------------------------------------------------------------------------------------------------------

# if defined (ENABLE_AZTEC_EPETRA) || defined (DOXYCOMPILE)


# include <limits>
# include <typeinfo>

# include <Epetra_Comm.h>
# include <Epetra_Map.h>
# include <Teuchos_ParameterList.hpp>

# if defined (ENABLE_MPI)
    # include <Epetra_MpiComm.h>
# else
    # include <Epetra_SerialComm.h>
# endif

# include "linear_solvers/Abstract/ImplicitSolver/AztecPeierlsSolver.hpp"
# include "operators/TransportOperator.hpp"
# include "utils/BelosTpetraInterface.hpp"


//============================================================================================================
//=== MEMBER DEFINITIONS FOR AztecPeierlsSolver CLASS =======================================================
//============================================================================================================


//--------------------------------------------------------------------------------------------------------
//! \brief  Returns the string descriptor of the class type.
//--------------------------------------------------------------------------------------------------------
template< class AngularFlux >
const std::string AztecPeierlsSolver<AngularFlux>::Descriptor( void ) {

    return "aztec-peierls-gmres";
}


//--------------------------------------------------------------------------------------------------------
//! \brief  Returns the string descriptor for an objects type.
//--------------------------------------------------------------------------------------------------------
template< class AngularFlux >
const std::string AztecPeierlsSolver<AngularFlux>::GetDescriptor( void ) const {

    return Descriptor();
}



//------------------------------------------------------------------------------------------------------------
//! \brief  Initializes a AztecPeierlsSolver object.
//!
//! \param[in]  domain_decomposition    Contains parameters of spatial discretization.
//! \param[in]  ordinate_set            Contains parameters of angular discretization.
//! \param[in]  input_list              Contains additional parameters.
//! \param[in]  dt                      Initial timestep size.
//------------------------------------------------------------------------------------------------------------
template< class AngularFlux >
AztecPeierlsSolver<AngularFlux>::AztecPeierlsSolver (

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
    PRINT_STATUS( "Executing AztecPeierlsSolver<%s>::%s.\n", typeid(AngularFlux).name(), __func__ )

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
    this->ipls = std::shared_ptr<Epetra_Operator>( new IPLS( *this ) );

    this->epetra_vec_x = std::shared_ptr<Epetra_MultiVector>(
            new Epetra_MultiVector( this->ipls->OperatorDomainMap(), 1 )
        );

    this->epetra_vec_b = std::shared_ptr<Epetra_MultiVector>(
            new Epetra_MultiVector( this->ipls->OperatorRangeMap(), 1 )
        );

    // Setup Epetra_LinearProblem and AztecOO solver.
    this->transport_prob = std::shared_ptr<Epetra_LinearProblem>(
            new Epetra_LinearProblem( this->ipls.get(), this->epetra_vec_x.get(), this->epetra_vec_b.get() )
        );

    this->transport_solver = std::shared_ptr<AztecOO>(
            new AztecOO( *this->transport_prob )
        );

    // Set default solver options.
    this->transport_params = std::shared_ptr<Teuchos::ParameterList>( new Teuchos::ParameterList() );

    this->transport_params->set( "AZ_solver",  "AZ_gmres" );  // Default, but set anyway.
    this->transport_params->set( "AZ_precond", "AZ_none"  );  // Default, but set anyway.

    // Determine convergence based on 2-norm of the residual.
    this->transport_params->set( "AZ_conv", "AZ_noscaled" );

    this->transport_params->set( "AZ_output", "AZ_none" );   // Don't print things, I'll do it myself.

    // Defaults for convergence.
    this->transport_params->set( "AZ_max_iter", 2000   ); // Max iterations.
    this->transport_params->set( "AZ_kspace",   30     ); // Restart.
    this->transport_params->set( "AZ_tol",      1.0e-8 ); // Tolerance.

    AztecPeierlsSolver<AngularFlux>::SetParameters( input_list );
}


//------------------------------------------------------------------------------------------------------------
//! \brief  Virtual destructor.
//------------------------------------------------------------------------------------------------------------
template< class AngularFlux >
AztecPeierlsSolver<AngularFlux>::~AztecPeierlsSolver ( void ) {

    PRINT_STATUS( "Executing AztecPeierlsSolver<%s>::%s.\n", typeid(AngularFlux).name(), __func__ )
}


//------------------------------------------------------------------------------------------------------------
//! \brief  Prints a summary of the object's configuration to the logging interface.
//!
//! \param[in]  prefix  String prepended to each line of output, intended to be used for controlling
//!                     indentation of output.
//------------------------------------------------------------------------------------------------------------
template< class AngularFlux >
void AztecPeierlsSolver<AngularFlux>::Print (

    const std::string & prefix // = "  "

) const {

    PRINT_STATUS( "Executing AztecPeierlsSolver<%s>::%s.\n", typeid(AngularFlux).name(), __func__ )

    this->ImplicitSolver<AngularFlux>::Print( prefix );

    PRINT_LOG( "%s%-*s  %s\n", prefix.c_str(), Global::col_width, "Solve Type:", this->Descriptor().c_str() )

    const double abstol  = this->transport_params->template get<double>( "AZ_tol"      );
    const int    maxits  = this->transport_params->template get<int>   ( "AZ_max_iter" );
    const int    restart = this->transport_params->template get<int>   ( "AZ_kspace"   );

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
void AztecPeierlsSolver<AngularFlux>::SetParameters (

    const ParameterList & input_list
) {
    PRINT_STATUS( "Executing AztecPeierlsSolver<%s>::%s.\n", typeid(AngularFlux).name(), __func__ )

    try {  this->transport_params->set( "AZ_tol",
                                        input_list.GetValue<double>( "abs_tol" ) );
    } catch (...) {/* empty */}

    try {  this->transport_params->set( "AZ_max_iter",
                                        input_list.GetValue<int>( "max_its" ) );
    } catch (...) {/* empty */}

    try {  this->transport_params->set( "AZ_kspace",
                                        input_list.GetValue<int>( "gmres_restart" ) );
    } catch (...) {/* empty */}

    this->transport_solver->SetParameters( *this->transport_params, true );

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
void AztecPeierlsSolver<AngularFlux>::SetInitialGuess (

    const AngularFlux * const guess    // = nullptr
) {
    PRINT_STATUS( "Executing AztecPeierlsSolver<%s>::%s.\n", typeid(AngularFlux).name(), __func__ )

    if ( guess == nullptr ) {

        PRINT_STATUS( "Setting initial guess to zero.\n" )

        this->epetra_vec_x->PutScalar( 0.0 );

    // Otherwise use guess.
    } else {

        PRINT_STATUS( "Setting initial guess from given object.\n" )

        ScalarFlux & phi = *this->phi.at(0);

        TransportOperator::Pmv( 1.0, 0.0, *guess, phi );
        phi.PackExternalVector( *this->epetra_vec_x, 1.0, 0.0 );
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
void AztecPeierlsSolver<AngularFlux>::Solve (

    AngularFlux & source,
    const CrossSection & sigma_t,
    const CrossSection & sigma_s,
    const double dt,                                    // = std::numeric_limits<double>::infinity(),
    const AngularFluxInit * const initial_condition     // = nullptr
) {

    PRINT_STATUS( "Executing AztecPeierlsSolver<%s>::%s.\n", typeid(AngularFlux).name(), __func__ )

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
    phi.PackExternalVector( *this->epetra_vec_b, 1.0, 0.0 );

    // --- Solve for the scalar flux. ----------------------------------------------------------------- //

    Global::TMR_PETSc.Start();

    // Solve system.
    int result = this->transport_solver->Iterate(
            this->transport_params->template get<int>( "AZ_max_iter" ),
            this->transport_params->template get<double>( "AZ_tol" )
        );

    // Get information about the solve.
# if LOGLEVEL >= 1

    const int num_its = this->transport_solver->NumIters();
    const double final_norm = this->transport_solver->TrueResidual();

    PRINT_LOG( "Iterations = %5d \tFinal Residual = % 10.4e\n", num_its, final_norm )

# endif // if LOGLEVEL >= 1

    Global::TMR_PETSc.Stop();

    if ( result != 0 ) {

        std::string error_message = "AztecOO GMRES transport solve failed to converge.\n";

        PRINT_ERROR( error_message.c_str() )
        throw std::runtime_error( error_message );
    }

    // --- Compute angular flux from scalar flux. ----------------------------------------------------- //

    phi.UnpackExternalVector( *this->epetra_vec_x );

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
//=== MEMBER DEFINITIONS FOR AztecPeierlsSolver::IPLS CLASS ===============================================
//============================================================================================================


//------------------------------------------------------------------------------------------------------------
//! \brief  Constructs an IPLS object from the given parameters.
//!
//! \note   This routine assumes that the enclosing object \pp{solver} has already initialized
//!         AztecPeierlsSolver::phi with the appropriate objects.
//!
//! \param[in]  solver          Reference to enclosing solver object.
//------------------------------------------------------------------------------------------------------------
template< class AngularFlux >
AztecPeierlsSolver<AngularFlux>::IPLS::IPLS (

    const AztecPeierlsSolver<AngularFlux> & solver
) :
    enc{ solver }
{

    PRINT_STATUS( "Executing AztecPeierlsSolver<%s>::IPLS::%s.\n", typeid(AngularFlux).name(), __func__ )

    this->comm = std::shared_ptr<Epetra_Comm>(
        # if defined (ENABLE_MPI)
            new Epetra_MpiComm( Global::MPI_cart_comm )
        # else
            new Epetra_SerialComm()
        # endif
        );

    this->domain_decomp_map = std::shared_ptr<Epetra_Map>(
            new Epetra_Map( (long long int) enc.phi.at(0)->GetGlobalArrayDim(),
                            (long long int) enc.phi.at(0)->GetLocalArrayDim(),
                            (int) 0, *this->comm )
        );
}


//------------------------------------------------------------------------------------------------------------
//! \brief  Returns -1 because this implementation of the Epetra_Operator interface does not support
//!         transpose use.
//------------------------------------------------------------------------------------------------------------
template< class AngularFlux >
int AztecPeierlsSolver<AngularFlux>::IPLS::SetUseTranspose ( bool ) {

    return -1;
}


//------------------------------------------------------------------------------------------------------------
//! \brief  Computes
//!         \f$ \vec{y} \gets \beta \vec{y} + \alpha (\mathcal{I} - \mathcal{PL}^{-1}\mathcal{S}) \vec{x} \f$.
//!
//! \param[in]      X       Input vector \f$ x \f$.
//! \param[in,out]  Y       Output vector \f$ y \f$.
//!
//------------------------------------------------------------------------------------------------------------
template< class AngularFlux >
int AztecPeierlsSolver<AngularFlux>::IPLS::Apply (

    const Epetra_MultiVector & X,
    Epetra_MultiVector & Y

) const {

    PRINT_STATUS( "Executing AztecPeierlsSolver<%s>::IPLS::%s.\n", typeid(AngularFlux).name(), __func__ )

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

    phi.PackExternalVector( Y, 1.0, 0.0 );

    Global::TMR_PETSc.Start();

    return 0;
}


//------------------------------------------------------------------------------------------------------------
//! \brief  Returns -1 because this implementation does not support inverse use.
//------------------------------------------------------------------------------------------------------------
template< class AngularFlux >
int AztecPeierlsSolver<AngularFlux>::IPLS::ApplyInverse (

    const Epetra_MultiVector &,
    Epetra_MultiVector &

) const {

    return -1;
}


//--------------------------------------------------------------------------------------------------------
//! \brief  Returns -1 because this implementation does not support use of the infinity norm.
//--------------------------------------------------------------------------------------------------------
template< class AngularFlux >
double AztecPeierlsSolver<AngularFlux>::IPLS::NormInf ( void ) const {

    return -1.0;
}


//--------------------------------------------------------------------------------------------------------
//! \brief  Returns a character string describing the operator.
//--------------------------------------------------------------------------------------------------------
template< class AngularFlux >
const char * AztecPeierlsSolver<AngularFlux>::IPLS::Label ( void ) const {

    return this->label;
}


//--------------------------------------------------------------------------------------------------------
//! \brief  Returns false because this implementation does not support transpose use.
//--------------------------------------------------------------------------------------------------------
template< class AngularFlux >
bool AztecPeierlsSolver<AngularFlux>::IPLS::UseTranspose ( void ) const {

    return false;
}


//--------------------------------------------------------------------------------------------------------
//! \brief  Returns false because this object cannot provide an approximate Inf-norm.
//--------------------------------------------------------------------------------------------------------
template< class AngularFlux >
bool AztecPeierlsSolver<AngularFlux>::IPLS::HasNormInf ( void ) const {

    return false;
}


//--------------------------------------------------------------------------------------------------------
//! \brief  Returns a reference to the Epetra_Comm communicator associated with this operator.
//--------------------------------------------------------------------------------------------------------
template< class AngularFlux >
const Epetra_Comm & AztecPeierlsSolver<AngularFlux>::IPLS::Comm ( void ) const {

    return *this->comm;
}


//--------------------------------------------------------------------------------------------------------
//! \brief  Returns a reference to the Epetra_Map object associated with the domain of this operator.
//--------------------------------------------------------------------------------------------------------
template< class AngularFlux >
const Epetra_Map & AztecPeierlsSolver<AngularFlux>::IPLS::OperatorDomainMap ( void ) const {

    return *this->domain_decomp_map;
}


//--------------------------------------------------------------------------------------------------------
//! \brief  Returns a reference to the Epetra_Map object associated with the range of this operator.
//--------------------------------------------------------------------------------------------------------
template< class AngularFlux >
const Epetra_Map & AztecPeierlsSolver<AngularFlux>::IPLS::OperatorRangeMap ( void ) const {

    return *this->domain_decomp_map;
}


//============================================================================================================
//=== EXPLICIT INSTANTIATIONS ================================================================================
//============================================================================================================

# include "linear_solvers/Abstract/ImplicitSolver/ImplicitSolverDerivedFactory.hpp"
# include "objects/RKDG/OrdinateFlux.hpp"
# include "objects/STDG/OrdinateFlux.hpp"


//!
//! \brief  Instantiation of AztecPeierlsSolver<RKDG::OrdinateFlux> class.
//!
template class AztecPeierlsSolver<RKDG::OrdinateFlux>;

//!
//! \brief  Instantiation of factory for AztecPeierlsSolver<RKDG::OrdinateFlux>.
//!
template class
ImplicitSolver<RKDG::OrdinateFlux>::ImplicitSolverDerivedFactory<
        AztecPeierlsSolver<RKDG::OrdinateFlux>
    >;


//!
//! \brief  Instantiation of AztecPeierlsSolver<STDG::OrdinateFlux> class.
//!
template class AztecPeierlsSolver<STDG::OrdinateFlux>;

//!
//! \brief  Instantiation of factory for AztecPeierlsSolver<STDG::OrdinateFlux>.
//!
template class
ImplicitSolver<STDG::OrdinateFlux>::ImplicitSolverDerivedFactory<
        AztecPeierlsSolver<STDG::OrdinateFlux>
    >;


# endif // if defined (ENABLE_PETSC)
