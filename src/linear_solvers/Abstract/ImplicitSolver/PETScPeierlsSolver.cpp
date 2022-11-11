//------------------------------------------------------------------------------------------------------------
//! \file   linear_solvers/Abstract/ImplicitSolver/PETScPeierlsSolver.cpp
//! \brief  Contains implementations and instantiations of methods from PETScPeierlsSolver class template.
//!
//! \author Michael M. Crockatt
//! \date   March 2018
//------------------------------------------------------------------------------------------------------------

# if defined (ENABLE_PETSC) || defined (DOXYCOMPILE)


# include <limits>
# include <typeinfo>

# include "linear_solvers/Abstract/ImplicitSolver/PETScPeierlsSolver.hpp"
# include "operators/TransportOperator.hpp"


//============================================================================================================
//=== MEMBER DEFINITIONS FOR PETScPeierlsSolver CLASS ========================================================
//============================================================================================================


//!
//! \brief  Initialization of instance pointer.
//!
template< class AngularFlux >
PETScPeierlsSolver<AngularFlux> * PETScPeierlsSolver<AngularFlux>::instance{ nullptr };


//--------------------------------------------------------------------------------------------------------
//! \brief  Returns the string descriptor of the class type.
//--------------------------------------------------------------------------------------------------------
template< class AngularFlux >
const std::string PETScPeierlsSolver<AngularFlux>::Descriptor( void ) {

    return "petsc-peierls-gmres";
}


//--------------------------------------------------------------------------------------------------------
//! \brief  Returns the string descriptor for an objects type.
//--------------------------------------------------------------------------------------------------------
template< class AngularFlux >
const std::string PETScPeierlsSolver<AngularFlux>::GetDescriptor( void ) const {

    return Descriptor();
}


//------------------------------------------------------------------------------------------------------------
//! \brief  Initializes an PETScPeierlsSolver object.
//!
//! \param[in]  domain_decomposition    Contains parameters of spatial discretization.
//! \param[in]  ordinate_set            Contains parameters of angular discretization.
//! \param[in]  input_list              Contains additional parameters.
//! \param[in]  dt                      Initial timestep size.
//------------------------------------------------------------------------------------------------------------
template< class AngularFlux >
PETScPeierlsSolver<AngularFlux>::PETScPeierlsSolver (

    const DomainDecomposition & domain_decomposition,
    const Quadrule::OrdinateSet & ordinate_set,
    const ParameterList & input_list,
    const double dt                 // = std::numeric_limits<double>::infinity()
) :
    ImplicitSolver<AngularFlux>( domain_decomposition, ordinate_set, input_list, dt ),

    dt_save{ std::nan("0") },
    sigma_t_save{ nullptr },
    sigma_s_save{ nullptr },
    use_dsa{ false }
{
    PRINT_STATUS( "Executing PETScPeierlsSolver<%s>::%s.\n", typeid(AngularFlux).name(), __func__ )

    // Setup sweep.
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

    // Setup PETSc variables.
    const PetscInt N = this->phi.at(0)->GetLocalArrayDim();

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
                          (void(*)(void)) PETScPeierlsSolver<AngularFlux>::IPLS );

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

# if SPACE_DIMS == 1

    // Test for DSA preconditioner options.
    PetscBool petsc_dsa = PETSC_FALSE;

    PetscOptionsGetBool( nullptr, nullptr, "-dsa", &petsc_dsa, nullptr );
    try {  this->use_dsa = input_list.GetValue<bool>( "use_dsa" );  } catch (...) {/* emtpy */}

    this->use_dsa |= petsc_dsa;

    // Setup variables for DSA preconditioner.
    if ( this->use_dsa ) {

    # if defined (ENABLE_MPI)
        if ( Global::MPI_num_ranks > 1 ) {

            std::string error_message = "DSA preconditioner enabled with more than one MPI rank, "
                                        "but DSA preconditioner only supports one MPI rank.\n";

            PRINT_ERROR( error_message.c_str() )
            throw std::invalid_argument( error_message );
        }
    # endif

        this->phi.push_back( std::shared_ptr<ScalarFlux>(
                ScalarFlux::Create( domain_decomposition, input_list )
            ) );

        // Create vectors for diffusion solve.
        VecDuplicate( this->transp_vec_sol, &this->diff_vec_sol );
        VecDuplicate( this->transp_vec_sol, &this->diff_vec_rhs );

        // Obtain and setup preconditioner context for transport solve.
        KSPGetPC( this->transp_solver, &this->transp_precond );

        PCSetType( this->transp_precond, PCSHELL );
        PCShellSetApply( this->transp_precond, DSA_Apply );
        PCShellSetName( this->transp_precond, "DSA" );

        // Setup diffusion matrix.
        MatCreate(
            # if defined (ENABLE_MPI)
                Global::MPI_cart_comm,
            # else
                PETSC_COMM_WORLD,
            # endif
                &this->diff_mat );

        MatSetType( this->diff_mat, MATSEQAIJ );
        MatCreateSeqAIJ(
            # if defined (ENABLE_MPI)
                Global::MPI_cart_comm,
            # else
                PETSC_COMM_WORLD,
            # endif
                N, N, 3 * this->phi.at(0)->DOFPerCell(), nullptr, &this->diff_mat );

        MatSetOption( this->diff_mat, MAT_IGNORE_ZERO_ENTRIES, PETSC_TRUE );

        // Setup context for diffusion solve.
        KSPCreate(
            # if defined (ENABLE_MPI)
                Global::MPI_cart_comm,
            # else
                PETSC_COMM_WORLD,
            # endif
                &this->diff_solver );

        KSPSetOperators( this->diff_solver, this->diff_mat, this->diff_mat );
        KSPSetType( this->diff_solver, KSPCG );
//         KSPSetInitialGuessNonzero( this->diff_solver, PETSC_TRUE );
        KSPSetFromOptions( this->diff_solver );

        // Set parameter list with default values.
        this->diff_params = std::shared_ptr<ParameterList>( new ParameterList() );

        this->diff_params->SetValue( "dsa_rtol",    "1.0e-60"   );
        this->diff_params->SetValue( "dsa_abs_tol", "1.0e-8"    );
        this->diff_params->SetValue( "dsa_dtol",    "1.0e+8"    );
        this->diff_params->SetValue( "dsa_max_its", "2000"      );

    } // end if ( this->use_dsa )

# endif // if SPACE_DIMS == 1

    PETScPeierlsSolver<AngularFlux>::SetParameters( input_list );
}


//------------------------------------------------------------------------------------------------------------
//! \brief  Virtual destructor.
//------------------------------------------------------------------------------------------------------------
template< class AngularFlux >
PETScPeierlsSolver<AngularFlux>::~PETScPeierlsSolver ( void ) {

    PRINT_STATUS( "Executing PETScPeierlsSolver<%s>::%s.\n", typeid(AngularFlux).name(), __func__ )

    KSPDestroy( &this->transp_solver );
    MatDestroy( &this->transp_mat );
    VecDestroy( &this->transp_vec_sol );
    VecDestroy( &this->transp_vec_rhs );

    if ( this->use_dsa ) {

        KSPDestroy( &this->diff_solver );
        VecDestroy( &this->diff_vec_sol );
        VecDestroy( &this->diff_vec_rhs );
    }
}


//------------------------------------------------------------------------------------------------------------
//! \brief  Prints a summary of the object's configuration to the logging interface.
//------------------------------------------------------------------------------------------------------------
template< class AngularFlux >
void PETScPeierlsSolver<AngularFlux>::Print (

    const std::string & prefix // = "  "

) const {

    PRINT_STATUS( "Executing PETScPeierlsSolver<%s>::%s.\n", typeid(AngularFlux).name(), __func__ )

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

    // Print preconditioner information, tolerances, etc.
    if ( this->use_dsa ) {

        PRINT_LOG( "%s%-*s  %s\n", prefix.c_str(), Global::col_width, "Preconditioner:", "dsa" )

        PRINT_LOG( "%s%-*s % .4e\n", prefix.c_str(), Global::col_width, "DSA rtol:",
                   this->diff_params->template GetValue<double>( "dsa_rtol" ) )

        PRINT_LOG( "%s%-*s % .4e\n", prefix.c_str(), Global::col_width, "DSA abstol:",
                   this->diff_params->template GetValue<double>( "dsa_abs_tol" ) )

        PRINT_LOG( "%s%-*s % .4e\n", prefix.c_str(), Global::col_width, "DSA dtol:",
                   this->diff_params->template GetValue<double>( "dsa_dtol" ) )

        PRINT_LOG( "%s%-*s % d\n", prefix.c_str(), Global::col_width, "DSA maxits:",
                   this->diff_params->template GetValue<int>( "dsa_max_its" ) )

    } else {

        PRINT_LOG( "%s%-*s  %s\n", prefix.c_str(), Global::col_width, "Preconditioner:", "none" )
    }

    this->sweep_operator->Print( prefix );
}


//------------------------------------------------------------------------------------------------------------
//! \brief  Sets parameters given a ParameterList.
//!
//! \param[in]  input_list  Contains parameters to set in the object.
//------------------------------------------------------------------------------------------------------------
template< class AngularFlux >
void PETScPeierlsSolver<AngularFlux>::SetParameters (

    const ParameterList & input_list
) {
    PRINT_STATUS( "Executing PETScPeierlsSolver<%s>::%s.\n", typeid(AngularFlux).name(), __func__ )

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

    // DSA options.
    if ( this->use_dsa ) {

        for ( const std::string & key : { "dsa_rtol", "dsa_abs_tol", "dsa_dtol", "dsa_max_its" } ) {

            try {
                this->diff_params->SetValue( key,
                                             input_list.template GetValue<std::string>( key ) );
            } catch (...) {/* empty */}
        }

        PetscReal dsa_rtol      = this->diff_params->template GetValue<PetscReal>( "dsa_rtol" ),
                  dsa_abs_tol   = this->diff_params->template GetValue<PetscReal>( "dsa_abs_tol" ),
                  dsa_dtol      = this->diff_params->template GetValue<PetscReal>( "dsa_dtol" );

        PetscInt dsa_max_its    = this->diff_params->template GetValue<PetscInt>( "dsa_max_its" );

        KSPSetTolerances( this->diff_solver, dsa_rtol, dsa_abs_tol, dsa_dtol, dsa_max_its );
    }

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
void PETScPeierlsSolver<AngularFlux>::SetInitialGuess (

    const AngularFlux * const guess    // = nullptr
) {
    PRINT_STATUS( "Executing PETScPeierlsSolver<%s>::%s.\n", typeid(AngularFlux).name(), __func__ )

    if ( guess == nullptr ) {

        PRINT_STATUS( "Setting initial guess to zero.\n" )

        VecSet( this->transp_vec_sol, 0.0 );

    // Otherwise use guess.
    } else {

        PRINT_STATUS( "Setting initial guess from given object.\n" )

        ScalarFlux & phi_ref = *this->phi.at(0);

        TransportOperator::Pmv( 1.0, 0.0, *guess, phi_ref );
        phi_ref.PackExternalVector( this->transp_vec_sol, 1.0, 0.0 );
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
void PETScPeierlsSolver<AngularFlux>::Solve (

    AngularFlux & source,
    const CrossSection & sigma_t,
    const CrossSection & sigma_s,
    const double dt,                                    // = std::numeric_limits<double>::infinity(),
    const AngularFluxInit * const initial_condition     // = nullptr
) {

    PRINT_STATUS( "Executing PETScPeierlsSolver<%s>::%s.\n", typeid(AngularFlux).name(), __func__ )

    PetscInt its;
    PetscReal rnorm;
    KSPConvergedReason reason;

    // --- Initialization. ---------------------------------------------------------------------------- //

    // Alias solve parameters for simplicity.
    AngularFlux & psi_ref = *this->psi;
    ScalarFlux & phi_ref = *this->phi.at(0);

    // Save parameters for solve.
    this->dt_save = dt;
    this->sigma_t_save = &sigma_t;
    this->sigma_s_save = &sigma_s;
    instance = this;

# if SPACE_DIMS == 1

    // Setup preconditioner if necessary.
    if ( this->use_dsa )
        this->DSA_ComputeMatrix( sigma_t, sigma_s, dt );

# endif

    // --- Compute right hand side of system. --------------------------------------------------------- //

    // L⁻¹Q
    psi_ref.Copy( source );
    this->sweep_operator->Apply( psi_ref, sigma_t, dt, initial_condition );

    // PL⁻¹Q
    TransportOperator::Pmv( 1.0, 0.0, psi_ref, phi_ref );
    phi_ref.PackExternalVector( this->transp_vec_rhs, 1.0, 0.0 );

    // --- Solve for the scalar flux. ----------------------------------------------------------------- //

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

    // --- Compute angular flux from scalar flux. ----------------------------------------------------- //

    phi_ref.UnpackExternalVector( this->transp_vec_sol );

    // L⁻¹(Q + SPΨ).
    TransportOperator::Smv( 1.0, 1.0, sigma_s, phi_ref, source );
    this->sweep_operator->Apply( source, sigma_t, dt, initial_condition );

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
PetscErrorCode PETScPeierlsSolver<AngularFlux>::IPLS (

    PETScMat,   // A,
    PETScVec x,
    PETScVec y
) {

    PRINT_STATUS( "Executing PETScPeierlsSolver<%s>::%s.\n", typeid(AngularFlux).name(), __func__ )

    Global::TMR_PETSc.Stop();

    // Alias solve parameters for simplicity.
    AngularFlux & psi = *instance->psi;
    ScalarFlux & phi = *instance->phi.at(0);

    const double dt = instance->dt_save;
    const CrossSection & sigma_t = *instance->sigma_t_save;
    const CrossSection & sigma_s = *instance->sigma_s_save;

    phi.UnpackExternalVector( x );

    // (I - PL⁻¹S)(PΨ)
    TransportOperator::Smv( 1.0, 0.0, sigma_s, phi, psi );
    instance->sweep_operator->Apply( psi, sigma_t, dt, nullptr );
    TransportOperator::Pmv( -1.0, 1.0, psi, phi );

    phi.PackExternalVector( y, 1.0, 0.0 );

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


# if SPACE_DIMS == 1


//------------------------------------------------------------------------------------------------------------
//! \brief  Applies the diffusion preconditioner \f$ \mathcal{I} - \mathcal{D}^{-1} \sigma_{\mathrm{s}} \f$ to
//!         a given vector.
//!
//! \param[in]  x       PETSc vector containing the coefficients to which the preconditioner matrix is
//!                     applied.
//! \param[out] y       Upon return contains the coefficients resulting from applying the DSA preconditioner
//!                     to the vector \pp{x}.
//!
//! \return     Returns a PetscErrorCode of \f$ 0 \f$ upon successful execution.
//------------------------------------------------------------------------------------------------------------
template< class AngularFlux >
PetscErrorCode PETScPeierlsSolver<AngularFlux>::DSA_Apply (

    PETScPC,
    PETScVec x,
    PETScVec y
) {

    PRINT_STATUS( "Executing PETScPeierlsSolver<%s>::%s.\n", typeid(AngularFlux).name(), __func__ )

    Global::TMR_PETSc.Stop();

    ScalarFlux & phi1 = *instance->phi.at(0);
    ScalarFlux & phi2 = *instance->phi.at(1);

    KSPConvergedReason reason;

    // Compute right hand side of diffusion system.
    phi1.UnpackExternalVector( x );
    instance->DSA_ComputeRHS( phi1, phi2 );

    // Solve diffusion system.
    Global::TMR_DSA_Solve.Start();

    KSPSolve( instance->diff_solver, instance->diff_vec_rhs, instance->diff_vec_sol );

    KSPGetConvergedReason( instance->diff_solver, &reason );

    Global::TMR_DSA_Solve.Stop();

# if LOGLEVEL >= 2

    PetscInt its;
    PetscReal rnorm;

    KSPGetIterationNumber( instance->diff_solver, &its );
    KSPGetResidualNorm( instance->diff_solver, &rnorm );

    PRINT_LOG( "DSA Iterations = %5d \tFinal Residual = % 10.4e \tStop Reason = %s\n",
               its, rnorm, KSPConvergedReasons[reason] )

# endif // if LOGLEVEL >= 2

    if ( reason < 0 ) {

        std::string error_message =   "DSA diffusion solve failed to converge with reason: "
                                    + std::string( KSPConvergedReasons[reason] )
                                    + ".\n";

        PRINT_ERROR( error_message.c_str() )
        throw std::runtime_error( error_message );
    }

    // Compute final result using diffusion solution.
    phi2.UnpackExternalVector( instance->diff_vec_sol );
    ScalarFlux::AXPY( 1.0, 1.0, phi2, phi1 );
    phi1.PackExternalVector( y, 1.0, 0.0 );

    Global::TMR_PETSc.Start();

    return 0;
}


//------------------------------------------------------------------------------------------------------------
//! \brief  Computes the right side of the DSA diffusion system.
//!
//! This default implementation throws an exception of type \c std::runtime_error, and should be overridden
//! through specialization based on the flux object type.
//------------------------------------------------------------------------------------------------------------
template< class AngularFlux >
void PETScPeierlsSolver<AngularFlux>::DSA_ComputeRHS (

    const ScalarFlux &,
    ScalarFlux &
) {

    PRINT_STATUS( "Executing PETScPeierlsSolver<%s>::%s.\n", typeid(AngularFlux).name(), __func__ )

    throw std::runtime_error(   "Implementation of PETScPeierlsSolver<"
                              + std::string( typeid(AngularFlux).name() )
                              + ">::"
                              + std::string(__func__)
                              + " incomplete." );
}


//------------------------------------------------------------------------------------------------------------
//! \brief  Computes the coefficients of the matrix for the diffusion solve in the DSA preconditioner.
//!
//! This default implementation throws an exception of type \c std::runtime_error, and should be overridden
//! through specialization based on the flux object type.
//------------------------------------------------------------------------------------------------------------
template< class AngularFlux >
void PETScPeierlsSolver<AngularFlux>::DSA_ComputeMatrix (

    const CrossSection &,
    const CrossSection &,
    const double
) {

    PRINT_STATUS( "Executing PETScPeierlsSolver<%s>::%s.\n", typeid(AngularFlux).name(), __func__ )

    throw std::runtime_error(   "Implementation of PETScPeierlsSolver<"
                              + std::string( typeid(AngularFlux).name() )
                              + ">::"
                              + std::string(__func__)
                              + " incomplete." );
}


//------------------------------------------------------------------------------------------------------------
//! \brief  Computes penalty parameters for the SIP-DG diffusion matrix.
//!
//! This default implementation throws an exception of type \c std::runtime_error, and should be overridden
//! through specialization based on the flux object type.
//------------------------------------------------------------------------------------------------------------
template< class AngularFlux >
double PETScPeierlsSolver<AngularFlux>::Kappa (

    const int64_t,
    const int64_t,
    const CrossSection &,
    const double

) const {

    PRINT_STATUS( "Executing PETScPeierlsSolver<%s>::%s.\n", typeid(AngularFlux).name(), __func__ )

    throw std::runtime_error(   "Implementation of PETScPeierlsSolver<"
                              + std::string( typeid(AngularFlux).name() )
                              + ">::"
                              + std::string(__func__)
                              + " incomplete." );
}


# endif // if SPACE_DIMS == 1


//============================================================================================================
//=== EXPLICIT INSTANTIATIONS ================================================================================
//============================================================================================================

# include "linear_solvers/Abstract/ImplicitSolver/ImplicitSolverDerivedFactory.hpp"
# include "objects/RKDG/OrdinateFlux.hpp"
# include "objects/STDG/OrdinateFlux.hpp"


//!
//! \brief  Instantiation of PETScPeierlsSolver<RKDG::OrdinateFlux> class.
//!
template class PETScPeierlsSolver<RKDG::OrdinateFlux>;

//!
//! \brief  Instantiation of factory for PETScPeierlsSolver<RKDG::OrdinateFlux>.
//!
template class
ImplicitSolver<RKDG::OrdinateFlux>::ImplicitSolverDerivedFactory<
        PETScPeierlsSolver<RKDG::OrdinateFlux>
    >;


//!
//! \brief  Instantiation of PETScPeierlsSolver<STDG::OrdinateFlux> class.
//!
template class PETScPeierlsSolver<STDG::OrdinateFlux>;

//!
//! \brief  Instantiation of factory for PETScPeierlsSolver<STDG::OrdinateFlux>.
//!
template class
ImplicitSolver<STDG::OrdinateFlux>::ImplicitSolverDerivedFactory<
        PETScPeierlsSolver<STDG::OrdinateFlux>
    >;


# endif // if defined (ENABLE_PETSC)
