//------------------------------------------------------------------------------------------------------------
//! \file   linear_solvers/Abstract/ImplicitSolver/PETScGMRESSolver.cpp
//! \brief  Contains implementations and instantiations of methods from PETScGMRESSolver class template.
//!
//! \author Michael M. Crockatt
//! \date   March 2018
//------------------------------------------------------------------------------------------------------------

# if defined (ENABLE_PETSC) || defined (DOXYCOMPILE)


# include <limits>
# include <typeinfo>

# include "linear_solvers/Abstract/ImplicitSolver/PETScGMRESSolver.hpp"
# include "operators/TransportOperator.hpp"


//============================================================================================================
//=== MEMBER DEFINITIONS FOR PETScGMRESSolver CLASS ========================================================
//============================================================================================================


//!
//! \brief  Initialization of instance pointer.
//!
template< class AngularFlux >
PETScGMRESSolver<AngularFlux> * PETScGMRESSolver<AngularFlux>::instance{ nullptr };


//--------------------------------------------------------------------------------------------------------
//! \brief  Returns the string descriptor of the class type.
//--------------------------------------------------------------------------------------------------------
template< class AngularFlux >
const std::string PETScGMRESSolver<AngularFlux>::Descriptor( void ) {

    return "petsc-gmres";
}


//--------------------------------------------------------------------------------------------------------
//! \brief  Returns the string descriptor for an objects type.
//--------------------------------------------------------------------------------------------------------
template< class AngularFlux >
const std::string PETScGMRESSolver<AngularFlux>::GetDescriptor( void ) const {

    return Descriptor();
}


//------------------------------------------------------------------------------------------------------------
//! \brief  Initializes an PETScGMRESSolver object.
//!
//! \param[in]  domain_decomposition    Contains parameters of spatial discretization.
//! \param[in]  ordinate_set            Contains parameters of angular discretization.
//! \param[in]  input_list              Contains additional parameters.
//! \param[in]  dt                      Initial timestep size.
//------------------------------------------------------------------------------------------------------------
template< class AngularFlux >
PETScGMRESSolver<AngularFlux>::PETScGMRESSolver (

    const DomainDecomposition & domain_decomposition,
    const Quadrule::OrdinateSet & ordinate_set,
    const ParameterList & input_list,
    const double dt                 // = std::numeric_limits<double>::infinity()
) :
    ImplicitSolver<AngularFlux>( domain_decomposition, ordinate_set, input_list, dt ),

    dt_save{ std::nan("0") },
    sigma_t_save{ nullptr },
    sigma_s_save{ nullptr }
{
    PRINT_STATUS( "Executing PETScGMRESSolver<%s>::%s.\n", typeid(AngularFlux).name(), __func__ )

    // Setup sweep.
    this->sweep_operator = std::shared_ptr<SweepOperator<AngularFlux>>(
            new SweepOperator<AngularFlux>( domain_decomposition, ordinate_set, input_list, dt )
        );

    // Setup temporary objects.
    this->psi = std::shared_ptr<AngularFlux>(
            AngularFlux::Create( domain_decomposition, ordinate_set, input_list )
        );

    this->phi = std::shared_ptr<ScalarFlux>(
            ScalarFlux::Create( domain_decomposition, input_list )
        );

    // Setup PETSc variables.
    const PetscInt N = this->psi->GetLocalArrayDim();

    // Setup vectors.
# if defined (ENABLE_MPI)
    VecCreateMPI( Global::MPI_cart_comm, N, PETSC_DETERMINE, &this->transp_vec_sol );
# else
    VecCreate( PETSC_COMM_WORLD, &this->transp_vec_sol );
    VecSetSizes( this->transp_vec_sol, PETSC_DECIDE, N );
# endif

    VecSetFromOptions( this->transp_vec_sol );
    VecDuplicate( this->transp_vec_sol, &this->transp_vec_rhs );

    // Setup shell matrix.
    MatCreateShell(
        # if defined (ENABLE_MPI)
            Global::MPI_cart_comm,
        # else
            PETSC_COMM_WORLD,
        # endif
            N, N, PETSC_DETERMINE, PETSC_DETERMINE, NULL, &this->transp_mat );

    MatShellSetOperation( this->transp_mat, MATOP_MULT,
                          (void(*)(void)) PETScGMRESSolver<AngularFlux>::ILLSP );

    // Setup KSP object for transport solver.
    KSPCreate(
        # if defined (ENABLE_MPI)
            Global::MPI_cart_comm,
        # else
            PETSC_COMM_WORLD,
        # endif
            &this->transp_solver );

    KSPSetOperators( this->transp_solver, this->transp_mat, this->transp_mat );
    KSPSetType( this->transp_solver, KSPLGMRES );
    KSPSetInitialGuessNonzero( this->transp_solver, PETSC_TRUE );
    KSPSetFromOptions( this->transp_solver );

    /* Setup parameter list with default values.
     *
     * The default PETSc behavior is to have a reasonable relative tolerance (1.0e-6) and a very small
     * absolute tolerance (1.0e-50) to force convergence with respect to the relative tolerance except in
     * cases where the absolute tolerance is exceptionally low. We want the opposite behavior, so default
     * values different from the PETSc defaults are provided here.
     */
    this->transp_params = std::shared_ptr<ParameterList>( new ParameterList() );

    this->transp_params->SetValue( "rtol",           "1.0e-60"   );
    this->transp_params->SetValue( "abs_tol",        "1.0e-8"    );
    this->transp_params->SetValue( "dtol",           "1.0e+8"    );
    this->transp_params->SetValue( "max_its",        "3000"      );
    this->transp_params->SetValue( "gmres_restart",  "30"        );

    PETScGMRESSolver<AngularFlux>::SetParameters( input_list );
}


//------------------------------------------------------------------------------------------------------------
//! \brief  Virtual destructor.
//------------------------------------------------------------------------------------------------------------
template< class AngularFlux >
PETScGMRESSolver<AngularFlux>::~PETScGMRESSolver ( void ) {

    PRINT_STATUS( "Executing PETScGMRESSolver<%s>::%s.\n", typeid(AngularFlux).name(), __func__ )

    KSPDestroy( &this->transp_solver );
    MatDestroy( &this->transp_mat );
    VecDestroy( &this->transp_vec_sol );
    VecDestroy( &this->transp_vec_rhs );
}


//------------------------------------------------------------------------------------------------------------
//! \brief  Prints a summary of the object's configuration to the logging interface.
//------------------------------------------------------------------------------------------------------------
template< class AngularFlux >
void PETScGMRESSolver<AngularFlux>::Print (

    const std::string & prefix // = "  "

) const {

    PRINT_STATUS( "Executing PETScGMRESSolver<%s>::%s.\n", typeid(AngularFlux).name(), __func__ )

    this->ImplicitSolver<AngularFlux>::Print( prefix );

    // Print tolerances, etc. for GMRES solver.
    PRINT_LOG( "%s%-*s  %s\n", prefix.c_str(), Global::col_width, "Solve Type:", this->Descriptor().c_str() )

    PRINT_LOG( "%s%-*s % .4e\n", prefix.c_str(), Global::col_width, "GMRES rtol:",
               this->transp_params->template GetValue<double>( "rtol" ) )

    PRINT_LOG( "%s%-*s % .4e\n", prefix.c_str(), Global::col_width, "GMRES abstol:",
               this->transp_params->template GetValue<double>( "abs_tol" ) )

    PRINT_LOG( "%s%-*s % .4e\n", prefix.c_str(), Global::col_width, "GMRES dtol:",
               this->transp_params->template GetValue<double>( "dtol" ) )

    PRINT_LOG( "%s%-*s % d\n", prefix.c_str(), Global::col_width, "GMRES maxits:",
               this->transp_params->template GetValue<int>( "max_its" ) )

    PRINT_LOG( "%s%-*s % d\n", prefix.c_str(), Global::col_width, "GMRES restart:",
               this->transp_params->template GetValue<int>( "gmres_restart" ) )

    this->sweep_operator->Print( prefix );
}


//------------------------------------------------------------------------------------------------------------
//! \brief  Sets parameters given a ParameterList.
//!
//! \param[in]  input_list  Contains parameters to set in the object.
//------------------------------------------------------------------------------------------------------------
template< class AngularFlux >
void PETScGMRESSolver<AngularFlux>::SetParameters (

    const ParameterList & input_list
) {
    PRINT_STATUS( "Executing PETScGMRESSolver<%s>::%s.\n", typeid(AngularFlux).name(), __func__ )

    // Update values in this->transp_params.
    for ( const std::string & key : { "rtol", "abs_tol", "dtol", "max_its", "gmres_restart" } ) {

        try {
            this->transp_params->SetValue( key,
                                           input_list.template GetValue<std::string>( key ) );
        } catch (...) {/* empty */}
    }

    // Set values in solver(s) from this->transp_params.

    PetscReal rtol      = this->transp_params->template GetValue<PetscReal>( "rtol" ),
              abs_tol   = this->transp_params->template GetValue<PetscReal>( "abs_tol" ),
              dtol      = this->transp_params->template GetValue<PetscReal>( "dtol" );

    PetscInt max_its    = this->transp_params->template GetValue<PetscInt>( "max_its" ),
             restart    = this->transp_params->template GetValue<PetscInt>( "gmres_restart" );

    KSPSetTolerances( this->transp_solver, rtol, abs_tol, dtol, max_its );
    KSPGMRESSetRestart( this->transp_solver, restart );

    // Propagate call up and down solver hierarchy.
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
void PETScGMRESSolver<AngularFlux>::SetInitialGuess (

    const AngularFlux * const guess    // = nullptr
) {
    PRINT_STATUS( "Executing PETScGMRESSolver<%s>::%s.\n", typeid(AngularFlux).name(), __func__ )

    if ( guess == nullptr ) {

        PRINT_STATUS( "Setting initial guess to zero.\n" )

        VecSet( this->transp_vec_sol, 0.0 );

    // Otherwise use guess.
    } else {

        PRINT_STATUS( "Setting initial guess from given object.\n" )

        guess->PackExternalVector( this->transp_vec_sol, 1.0, 0.0 );
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
void PETScGMRESSolver<AngularFlux>::Solve (

    AngularFlux & source,
    const CrossSection & sigma_t,
    const CrossSection & sigma_s,
    const double dt,                                    // = std::numeric_limits<double>::infinity(),
    const AngularFluxInit * const initial_condition     // = nullptr
) {

    PRINT_STATUS( "Executing PETScGMRESSolver<%s>::%s.\n", typeid(AngularFlux).name(), __func__ )

    PetscInt its;
    PetscReal rnorm;
    KSPConvergedReason reason;

    // --- Initialization. ---------------------------------------------------------------------------- //

    // Alias solve parameters for simplicity.
    AngularFlux & psi_ref = *this->psi;

    // Save parameters for solve.
    this->dt_save = dt;
    this->sigma_t_save = &sigma_t;
    this->sigma_s_save = &sigma_s;
    instance = this;

    // --- Compute right hand side of system. --------------------------------------------------------- //

    // L⁻¹Q
    psi_ref.Copy( source );
    psi_ref.ZeroDensity( OpDomain::Halo );
    this->sweep_operator->Apply( psi_ref, sigma_t, dt, initial_condition );
    psi_ref.PackExternalVector( this->transp_vec_rhs, 1.0, 0.0 );

    // --- Perform solve. ----------------------------------------------------------------------------- //

    Global::TMR_PETSc.Start();

    // Solve system.
    KSPSolve( this->transp_solver, this->transp_vec_rhs, this->transp_vec_sol );

    // Get information about the solve.
# if LOGLEVEL >= 1

    KSPGetIterationNumber( this->transp_solver, &its );
    KSPGetResidualNorm( this->transp_solver, &rnorm );
    KSPGetConvergedReason( this->transp_solver, &reason );

    PRINT_LOG( "Iterations = %5d \tFinal Residual = % 10.4e \tStop Reason = %s\n",
               its, rnorm, KSPConvergedReasons[reason] )

# endif // if LOGLEVEL >= 1

    Global::TMR_PETSc.Stop();

    if ( reason < 0 ) {

        std::string error_message = "GMRES transport solve failed to converge with reason '"
                                    + std::string( KSPConvergedReasons[reason] )
                                    + "'.\n";

        PRINT_ERROR( error_message.c_str() )
        throw std::runtime_error( error_message );
    }

    // Get solution.
    source.UnpackExternalVector( this->transp_vec_sol );

    // --- Cleanup. ----------------------------------------------------------------------------------- //

    // Reset solve parameters.
    this->dt_save = std::nan("0");
    this->sigma_t_save = nullptr;
    this->sigma_s_save = nullptr;
    instance = nullptr;
}


//------------------------------------------------------------------------------------------------------------
//! \brief  Applies the operator \f$ (\mathcal{I} - \mathcal{PL}^{-1}\mathcal{S}) \f$ to a given vector.
//!
//! \f$ \vec{y} \gets (\mathcal{I} - \mathcal{PL}^{-1}\mathcal{S}) \vec{x} \f$.
//!
// \param[in]  A   PETSc shell matrix object that is to be applied.
//! \param[in]  x   PETSc vector to apply operator to.
//! \param[in]  y   PETSc vector in which to store result.
//!
//! \return     Returns a \c PetscErrorCode of 0 upon successful execution.
//------------------------------------------------------------------------------------------------------------
template< class AngularFlux >
PetscErrorCode PETScGMRESSolver<AngularFlux>::ILLSP (

    PETScMat,   // A,
    PETScVec x,
    PETScVec y
) {

    PRINT_STATUS( "Executing PETScGMRESSolver<%s>::%s.\n", typeid(AngularFlux).name(), __func__ )

    Global::TMR_PETSc.Stop();

    // Alias solve parameters for simplicity.
    AngularFlux & psi = *instance->psi;
    ScalarFlux & phi = *instance->phi;

    const double dt = instance->dt_save;
    const CrossSection & sigma_t = *instance->sigma_t_save;
    const CrossSection & sigma_s = *instance->sigma_s_save;

    psi.UnpackExternalVector( x );
    VecCopy( x, y );

    // (L₁ + SP)Ψ
    TransportOperator::Pmv( 1.0, 0.0, psi, phi, OpDomain::Interior );

    if ( psi.UpwindHalosAreDirty() ) {  psi.SynchronizeHalos();  }

    TransportOperator::Smv( 1.0, 0.0, sigma_s, phi, psi, OpDomain::Interior );

    // L₀⁻¹(L₁ + SP)Ψ
    instance->sweep_operator->Apply( psi, sigma_t, dt );

    // (I - L₀⁻¹(L₁ + SP))Ψ
    psi.PackExternalVector( y, -1.0, 1.0 );

    // Print iteration status.
# if LOGLEVEL >= 2
    PetscInt its;
    PetscReal rnorm;

    Global::TMR_PETSc.Start();
    KSPGetResidualNorm( instance->transp_solver, &rnorm );
    KSPGetIterationNumber( instance->transp_solver, &its );
    Global::TMR_PETSc.Stop();

    PRINT_LOG( "In GMRES transport solve: iteration = %4d \tresidual = % 10.4e\n", its, rnorm )
# endif

    Global::TMR_PETSc.Start();
    return 0;
}


//============================================================================================================
//=== EXPLICIT INSTANTIATIONS ================================================================================
//============================================================================================================

# include "linear_solvers/Abstract/ImplicitSolver/ImplicitSolverDerivedFactory.hpp"
# include "objects/RKDG/OrdinateFlux.hpp"
# include "objects/STDG/OrdinateFlux.hpp"


//!
//! \brief  Instantiation of PETScGMRESSolver<RKDG::OrdinateFlux> class.
//!
template class PETScGMRESSolver<RKDG::OrdinateFlux>;

//!
//! \brief  Instantiation of factory for PETScGMRESSolver<RKDG::OrdinateFlux>.
//!
template class
ImplicitSolver<RKDG::OrdinateFlux>::ImplicitSolverDerivedFactory<
        PETScGMRESSolver<RKDG::OrdinateFlux>
    >;


//!
//! \brief  Instantiation of PETScGMRESSolver<STDG::OrdinateFlux> class.
//!
template class PETScGMRESSolver<STDG::OrdinateFlux>;

//!
//! \brief  Instantiation of factory for PETScGMRESSolver<STDG::OrdinateFlux>.
//!
template class
ImplicitSolver<STDG::OrdinateFlux>::ImplicitSolverDerivedFactory<
        PETScGMRESSolver<STDG::OrdinateFlux>
    >;


# endif // if defined (ENABLE_PETSC)
