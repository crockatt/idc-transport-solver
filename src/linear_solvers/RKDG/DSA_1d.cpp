//------------------------------------------------------------------------------------------------------------
//! \file   linear_solvers/RKDG/DSA_1d.cpp
//! \brief  Implements a diffusion synthetic acceleration (DSA) preconditioner for the transport system using
//!         a modified interior penalty (MIP) construction.
//!
//! \todo   Citations for MIP DSA method.
//!
//! \authors    Michael Crockatt
//! \date       July 2017
//------------------------------------------------------------------------------------------------------------

# if defined (DONT_COMPILE_THIS)

# include <cmath>
# include <cstdint>
# include <cstdlib>
# include <limits>

# if defined (ENABLE_PETSC) || defined (DOXYCOMPILE)
    # include <petscksp.h>
# endif

# if ( SPACE_DIMS == 1 && defined (ENABLE_PETSC) && !defined (ENABLE_MPI) ) || defined (DOXYCOMPILE)

# include "linear_solvers/RKDG/ImplicitSolver_old.hpp"
# include "utils/CLog.hpp"
# include "utils/global.hpp"


using namespace RKDG;


namespace RKDG {

    //!
    //! \brief  Global variable used to store the solver class during an implicit solve.
    //!
    //! Necessary for accessing the solver class from RKDG::ImplicitSolverOLD::DSA_Apply(), which is called by
    //! PETSc.
    //!
    //! Variable is defined in \ref linear_solvers/RKDG/ImplicitSolverOLD.cpp
    //!
    extern ImplicitSolverOLD * s_this;
}


//============================================================================================================
//=== DIFFUSION SOLVER =======================================================================================
//============================================================================================================


//------------------------------------------------------------------------------------------------------------
//! \brief  Computes penalty parameters for the SIP-DG diffusion matrix.
//!
//! \param[in]  i
//! \param[in]  j
//! \param[in]  sigma_t     Total cross section \f$ \sigma_{\mathrm{t}} \f$.
//! \param[in]  dt          Timestep size. For steady-state problems, use \pp{dt} = \c inf.
//!
//! \return     Returns the desired penalty parameter.
//------------------------------------------------------------------------------------------------------------
double ImplicitSolverOLD::Kappa (

    const int64_t i,
    const int64_t j,
    const RKDG::CrossSection & sigma_t,
    const double dt

) const {

    double temp_val = 4.0 * DG_degree * (DG_degree + 1)
                      / ( 3.0 * dx(0) * ( sigma_t(i,0) + sigma_t(j,0) + 2.0 / dt ) );

//     return 200.0;
//     return temp;

    return std::max( temp_val, 0.25 );
}


//------------------------------------------------------------------------------------------------------------
//! \brief  Evaluates the \f$ L^2 \f$ inner product of derivatives of Legendre polynomials.
//!
//! \param[in]  a       Degree of first Legendre polynomial.
//! \param[in]  d       Degree of second Legendre polynomial.
//!
//! \return     Returns the value of \f[ \int_{-1} ^1 P'_a (x) P'_d (x) \, dx \f] which is equal to
//!             \f$ \beta ( \beta + 1 ) \f$ if \pp{a} and \pp{d} have the same parity and \f$ 0 \f$ otherwise,
//!             where \f$ \beta = \min ( a,d ) \f$.
//------------------------------------------------------------------------------------------------------------
template< typename T >
T B (

    const T a,
    const T d
) {

    T beta = std::min( a, d );

    return ( (beta + 1) * beta ) * ((T) ((a % 2) == (d % 2)));
}


//------------------------------------------------------------------------------------------------------------
//! \brief  Computes the coefficients of the SIP-DG matrix for the DSA diffusion solve.
//!
//! The material cross sections here use a piecewise constant approximation irregardless of the type of
//! approximation used elsewhere.
//!
//! The computed diffusion matrix is stored in ImplicitSolverOLD::D.
//!
//! \attention  This routine uses the indexing function associated with the first RKDG::ScalarFlux object
//!             stored at ImplicitSolverOLD::phi. Hence the ScalarFlux objects must be allocated to this pointer
//!             _before_ calling this routine.
//!
//! \param[in]  sigma_t     Total cross section \f$ \sigma_{\mathrm{t}} \f$.
//! \param[in]  sigma_s     Scattering cross section \f$ \sigma_{\mathrm{s}} \f$.
//! \param[in]  dt          Timestep size. For steady-state problems, use \pp{dt} = \c inf.
//------------------------------------------------------------------------------------------------------------
void ImplicitSolverOLD::ComputeDiffusionMatrix (

    const CrossSection & sigma_t,
    const CrossSection & sigma_s,
    const double dt
) {

    Global::TMR_DSA_Assemble.Start();

    MatZeroEntries( this->D );

    double coeff, val;
    PetscInt m_row, m_col;

    // Set values for the matrix.
    for ( PetscInt i = 0; i <  this->nx(0);     ++i ) {
    for ( PetscInt d = 0; d <= this->DG_degree; ++d ) {

        m_row = this->phi->IndexNoGC(i,d);

        // --- Cell integral terms. ------------------------------------------------------------------- //

        m_col = this->phi->IndexNoGC(i,d);
        val = this->dx(0) * ( sigma_t(i+1,0) - sigma_s(i+1,0) + 1.0 / dt ) / (2*d + 1);
        MatSetValues( this->D, 1, &m_row, 1, &m_col, &val, ADD_VALUES );

        coeff = 2.0 / ( 3.0 * this->dx(0) * ( sigma_t(i+1,0) + 1.0 / dt ) );

        for ( PetscInt a = 0; a <= this->DG_degree; ++a ) {

            m_col = this->phi->IndexNoGC(i,a);
            val = coeff * B(d,a);
            MatSetValues( this->D, 1, &m_row, 1, &m_col, &val, ADD_VALUES );
        }

        // --- Cell interface (inter-cell) flux terms. ------------------------------------------------ //

        // Left boundary cell.
        if ( i == 0 ) {

            // Penalty flux.
            for ( PetscInt a = 0; a <= this->DG_degree; ++a ) {

                m_col = this->phi->IndexNoGC(i+1,a);
                val = neg1pow(a+1) * Kappa( i+1, i+2, sigma_t, dt );
                MatSetValues( this->D, 1, &m_row, 1, &m_col, &val, ADD_VALUES );

                m_col = this->phi->IndexNoGC(i,a);
                val = Kappa( i+1, i+1, sigma_t, dt );
                MatSetValues( this->D, 1, &m_row, 1, &m_col, &val, ADD_VALUES );
            }

            // Test normal derivative flux.
            coeff = d * (d + 1) / ( 6.0 * this->dx(0) );

            for ( PetscInt a = 0; a <= this->DG_degree; ++a ) {

                m_col = this->phi->IndexNoGC(i,a);
                val = - coeff / ( sigma_t(i+1,0) + 1.0 / dt );
                MatSetValues( this->D, 1, &m_row, 1, &m_col, &val, ADD_VALUES );

                m_col = this->phi->IndexNoGC(i+1,a);
                val = 2.0 * coeff * neg1pow(a) / ( sigma_t(i+1,0) + sigma_t(i+2,0) + 2.0 / dt );
                MatSetValues( this->D, 1, &m_row, 1, &m_col, &val, ADD_VALUES );
            }

            // Unknown normal derivative flux.
            coeff = -1.0 / ( 6.0 * this->dx(0) );

            for ( PetscInt a = 0; a <= this->DG_degree; ++a ) {

                m_col = this->phi->IndexNoGC(i+1,a);
                val = 2.0 * coeff * neg1pow(a) * a * (a + 1) / ( sigma_t(i+1,0) + sigma_t(i+2,0) + 2.0 / dt );
                MatSetValues( this->D, 1, &m_row, 1, &m_col, &val, ADD_VALUES );

                m_col = this->phi->IndexNoGC(i,a);
                val = coeff * a * (a + 1) / ( sigma_t(i+1,0) + 1.0 / dt );
                MatSetValues( this->D, 1, &m_row, 1, &m_col, &val, ADD_VALUES );
            }

        // Right boundary cell
        } else if ( i == this->nx(0) -1 ) {

            // Penalty flux.
            for ( PetscInt a = 0; a <= this->DG_degree; ++a ) {

                m_col = this->phi->IndexNoGC(i-1,a);
                val = neg1pow(d+1) * Kappa( i, i+1, sigma_t, dt );
                MatSetValues( this->D, 1, &m_row, 1, &m_col, &val, ADD_VALUES );

                m_col = this->phi->IndexNoGC(i,a);
                val = neg1pow(d+a) * Kappa( i+1, i+1, sigma_t, dt );
                MatSetValues( this->D, 1, &m_row, 1, &m_col, &val, ADD_VALUES );
            }

            // Test normal derivative flux.
            coeff = neg1pow(d) * d * (d + 1) / ( 6.0 * this->dx(0) );

            for ( PetscInt a = 0; a <= this->DG_degree; ++a ) {

                m_col = this->phi->IndexNoGC(i-1,a);
                val = - 2.0 * coeff / ( sigma_t(i,0) + sigma_t(i+1,0) + 2.0 / dt );
                MatSetValues( this->D, 1, &m_row, 1, &m_col, &val, ADD_VALUES );

                m_col = this->phi->IndexNoGC(i,a);
                val = coeff * neg1pow(a) / ( sigma_t(i+1,0) + 1.0 / dt );
                MatSetValues( this->D, 1, &m_row, 1, &m_col, &val, ADD_VALUES );
            }

            // Unknown normal derivative flux.
            coeff = 1.0 / ( 6.0 * this->dx(0) );

            for ( PetscInt a = 0; a <= this->DG_degree; ++a ) {

                m_col = this->phi->IndexNoGC(i,a);
                val = coeff * neg1pow(d+a) * a * (a + 1) / ( sigma_t(i+1,0) + 1.0 / dt );
                MatSetValues( this->D, 1, &m_row, 1, &m_col, &val, ADD_VALUES );

                m_col = this->phi->IndexNoGC(i-1,a);
                val = 2.0 * coeff * neg1pow(d) * a * (a + 1) / ( sigma_t(i,0) + sigma_t(i+1,0) + 2.0 / dt );
                MatSetValues( this->D, 1, &m_row, 1, &m_col, &val, ADD_VALUES );
            }

        // Interior cells.
        } else {

            // Penalty flux.
            for ( PetscInt a = 0; a <= this->DG_degree; ++a ) {

                m_col = this->phi->IndexNoGC(i-1,a);
                val = Kappa( i, i+1, sigma_t, dt ) * neg1pow(d+1);
                MatSetValues( this->D, 1, &m_row, 1, &m_col, &val, ADD_VALUES );

                m_col = this->phi->IndexNoGC(i,a);
                val = Kappa( i+1, i+1, sigma_t, dt ) * ( 1 + neg1pow(d+a) );
                MatSetValues( this->D, 1, &m_row, 1, &m_col, &val, ADD_VALUES );

                m_col = this->phi->IndexNoGC(i+1,a);
                val = Kappa( i+1, i+2, sigma_t, dt ) * neg1pow(a+1);
                MatSetValues( this->D, 1, &m_row, 1, &m_col, &val, ADD_VALUES );
            }

            // Test normal derivative flux.
            coeff = d * (d + 1) / ( 6.0 * this->dx(0) );

            for ( PetscInt a = 0; a <= this->DG_degree; ++a ) {

                m_col = this->phi->IndexNoGC(i-1,a);
                val = 2.0 * coeff * neg1pow(d+1) / ( sigma_t(i,0) + sigma_t(i+1,0) + 2.0 / dt );
                MatSetValues( this->D, 1, &m_row, 1, &m_col, &val, ADD_VALUES );

                m_col = this->phi->IndexNoGC(i,a);
                val = coeff * ( neg1pow(d+a) - 1 ) / ( sigma_t(i+1,0) + 1.0 / dt );
                MatSetValues( this->D, 1, &m_row, 1, &m_col, &val, ADD_VALUES );

                m_col = this->phi->IndexNoGC(i+1,a);
                val = 2.0 * coeff * neg1pow(a) / ( sigma_t(i+1,0) + sigma_t(i+2,0) + 2.0 / dt );
                MatSetValues( this->D, 1, &m_row, 1, &m_col, &val, ADD_VALUES );
            }

            // Unknown normal derivative flux.
            coeff = 1.0 / ( 6.0 * this->dx(0) );

            for ( PetscInt a = 0; a <= this->DG_degree; ++a ) {

                m_col = this->phi->IndexNoGC(i-1,a);
                val = 2.0 * coeff * neg1pow(d) * a * (a + 1) / ( sigma_t(i,0) + sigma_t(i+1,0) + 2.0 / dt );
                MatSetValues( this->D, 1, &m_row, 1, &m_col, &val, ADD_VALUES );

                m_col = this->phi->IndexNoGC(i,a);
                val = coeff * ( neg1pow(d+a) - 1 ) * a * (a + 1) / ( sigma_t(i+1,0) + 1.0 / dt );
                MatSetValues( this->D, 1, &m_row, 1, &m_col, &val, ADD_VALUES );

                m_col = this->phi->IndexNoGC(i+1,a);
                val = 2.0 * coeff * neg1pow(a+1) * a * (a + 1) / ( sigma_t(i+1,0) +  sigma_t(i+2,0) + 2.0 / dt );
                MatSetValues( this->D, 1, &m_row, 1, &m_col, &val, ADD_VALUES );
            }
        }

        // --- Boundary flux terms. ------------------------------------------------------------------- //

        // Left boundary.
        if ( i == 0 ) {

            // Penalty flux.
            for ( PetscInt a = 0; a <= this->DG_degree; ++a ) {

                m_col = this->phi->IndexNoGC(i,a);
                val = Kappa( i+1, i+1, sigma_t, dt ) * neg1pow(d+a);
                MatSetValues( this->D, 1, &m_row, 1, &m_col, &val, ADD_VALUES );
            }

            // Test normal derivative flux.
            coeff = neg1pow(d+1) * d * (d + 1) / ( 6.0 * this->dx(0) * ( sigma_t(i+1,0) + 1.0 / dt ) );

            for ( PetscInt a = 0; a <= this->DG_degree; ++a ) {

                m_col = this->phi->IndexNoGC(i,a);
                val = coeff * neg1pow(a);
                MatSetValues( this->D, 1, &m_row, 1, &m_col, &val, ADD_VALUES );
            }

            // Unknown normal derivative flux.
            coeff = neg1pow(d+1) / ( 6.0 * this->dx(0) * ( sigma_t(i+1,0) + 1.0 / dt ) );

            for ( PetscInt a = 0; a <= this->DG_degree; ++a ) {

                m_col = this->phi->IndexNoGC(i,a);
                val = coeff * neg1pow(a) * a * (a + 1);
                MatSetValues( this->D, 1, &m_row, 1, &m_col, &val, ADD_VALUES );
            }

        // Right boundary
        } else if ( i == this->nx(0) - 1 ) {

            // Penalty flux.
            for ( PetscInt a = 0; a <= this->DG_degree; ++a ) {

                m_col = this->phi->IndexNoGC(i,a);
                val = Kappa( i+1, i+1, sigma_t, dt );
                MatSetValues( this->D, 1, &m_row, 1, &m_col, &val, ADD_VALUES );
            }

            // Test normal derivative flux.
            coeff = -d * (d + 1) / ( 6.0 * this->dx(0) * ( sigma_t(i+1,0) + 1.0 / dt ) );

            for ( PetscInt a = 0; a <= this->DG_degree; ++a ) {

                m_col = this->phi->IndexNoGC(i,a);
                val = coeff;
                MatSetValues( this->D, 1, &m_row, 1, &m_col, &val, ADD_VALUES );
            }

            // Unknown normal derivative flux.
            coeff = -1.0 / ( 6.0 * this->dx(0) * ( sigma_t(i+1,0) + 1.0 / dt ) );

            for ( PetscInt a = 0; a <= this->DG_degree; ++a ) {

                m_col = this->phi->IndexNoGC(i,a);
                val = coeff * a * (a + 1);
                MatSetValues( this->D, 1, &m_row, 1, &m_col, &val, ADD_VALUES );
            }
        }
    }}

    MatAssemblyBegin( this->D, MAT_FINAL_ASSEMBLY );
    MatAssemblyEnd( this->D, MAT_FINAL_ASSEMBLY );

    Global::TMR_DSA_Assemble.Stop();
}


//------------------------------------------------------------------------------------------------------------
//! \brief   Computes the RHS vector for the DSA diffusion system.
//!
//! The RHS vector is stored in ImplicitSolverOLD::v.
//------------------------------------------------------------------------------------------------------------
void ImplicitSolverOLD::DSA_ComputeRHS (

    const RKDG::ScalarFlux & phi,
    RKDG::ScalarFlux & Q
) {

    Q.ZeroDensity();

    Global::TMR_DSA_Solve.Start();

    # pragma omp parallel for
    for ( int64_t i = 0; i <= Q.nx(0) + 1; ++i ) {
    for ( int64_t d = 0; d <= Q.DG_degree; ++d ) {

        Q(i,d) = Q.dx(0) * phi(i,d) * (*sigma_s)(i,0) / (2*d + 1);
    }}

    Global::TMR_DSA_Solve.Stop();

    Q.PackPETScVec( this->v );
}


//------------------------------------------------------------------------------------------------------------
//! \brief  Applies the DSA preconditioner \f$ \mathcal{I} - \mathcal{D}^{-1} \sigma_{\mathrm{s}} \f$ to the
//!         given vector.
//!
//! \attention  This function is passed to PETSc as a function pointer to create a shell preconditioner for
//!             the GMRES transport solver.
//!
//! \param[in]  pc      PETSc shell precondition object which is to be applied.
//! \param[in]  x       PETSc vector containing the coefficients to which the preconditioner matrix is applied.
//! \param[out] y       Upon return contains the coefficients resulting from applying the DSA preconditioner
//!                     to the vector \pp{x}.
//!
//! \return     Returns a PetscErrorCode of \f$ 0 \f$ upon successful execution.
//!
//! \see    RKDG::ImplicitSolverOLD::ComputeDiffusionMatrix()
//! \see    RKDG::ImplicitSolverOLD::DSA_ComputeRHS()
//------------------------------------------------------------------------------------------------------------
PetscErrorCode ImplicitSolverOLD::DSA_Apply (

    PC pc,
    Vec x,
    Vec y
) {

    Global::TMR_PETSc.Stop();

    ScalarFlux & phi1 = s_this->phi[0];
    ScalarFlux & phi2 = s_this->phi[1];

    PetscInt its;
    PetscReal rnorm;
    KSPConvergedReason reason;

    // --- First solve the diffusion equation. -------------------------------------------------------- //

    // Compute right hand vector for the diffusion system.
    phi1.UnpackPETScVec( x );
    s_this->DSA_ComputeRHS( phi1, phi2 );

    // Solve the diffusion system.
    Global::TMR_DSA_Solve.Start();

    KSPSolve( s_this->diffusion_solver, s_this->v, s_this->u );

    // Get information about the solve.
    KSPGetIterationNumber( s_this->diffusion_solver, &its );
    KSPGetResidualNorm( s_this->diffusion_solver, &rnorm );
    KSPGetConvergedReason( s_this->diffusion_solver, &reason );

    Global::TMR_DSA_Solve.Stop();

# if LOGLEVEL >= 2
    PRINT_LOG( "DSA Iterations = %5d \tFinal Residual = % 10.4e \tStop Reason = %s\n",
               its, rnorm, KSPConvergedReasons[reason] )
# endif

    if ( reason < 0 ) {

        std::string error_message = "DSA diffusion solve failed to converge with reason '"
                                    + std::string( KSPConvergedReasons[reason] )
                                    + "'.\n";

        PRINT_ERROR( error_message.c_str() )
        throw std::runtime_error( error_message );
    }

    // --- Then compute final result using solution of diffusion equation. ---------------------------- //

    phi2.UnpackPETScVec( s_this->u );
    DensityFunction::AXPY( 1.0, phi2, phi1 );
    phi1.PackPETScVec( y );

    Global::TMR_PETSc.Start();

    return 0;
}


//============================================================================================================
//=== SOURCE ITERATION TRANSPORT SOLVER ======================================================================
//============================================================================================================

# if 0

//------------------------------------------------------------------------------------------------------------
// SI_solver
//      Solves the linear transport system (L - SP)Ψ = Q using source iteration:
//
//              Ψⁿ⁺¹ = L⁻¹ ( SP Ψⁿ + Q ).
//
// INPUT:
//      Q           - Pointer to array containing the coefficients of the source term of the transport system.
//
// OUTPUT:
//      Q           - Upon return contains the vector of coefficients Ψ which solves the transport system.
//------------------------------------------------------------------------------------------------------------
static
void SI_solver (

    double * const Q
) {

    double *psi_n1   = s_1tempPsi;
    double *psi_n    = s_2tempPsi;
    double *phi_one  = s_1tempPhi;
    double *phi_half = s_2tempPhi;
    double *swapr    = NULL;

    double *rhs2;
    double *phi2;

    int64_t iter = 0;
    double error    = std::numeric_limits<double>::max();
    double dPhiNorm = std::numeric_limits<double>::max();
    double tempVal;
    bool useDSA = s_useDSA;

    PetscInt its;
    PetscReal rnorm;
    KSPConvergedReason reason;


    if ( useDSA ) {

        memset( phi_one, 0, g_sizeofPhi );

        #pragma omp parallel for
        for ( int64_t i = 1 ; i <= g_nx      ; ++i ) {
        for ( int64_t d = 0 ; d <= g_dgOrder ; ++d ) {
        for ( int64_t q = 0 ; q <  g_nq      ; ++q ) {

            phi_one[ IPHI(i,d) ] += g_w[q] * psi_n1[ IPSI(i,q,d) ];

        }}}
    }

    while (
        ( iter  < s_SImaxIts ) &&
        ( error > s_SItol    )
    ) {

        ++iter;

        // Swap iterate pointers.
        swapr  = psi_n1;
        psi_n1 = psi_n;
        psi_n  = swapr;


        // Turn off DSA if norm of correction is below threshold.
        if (
            ( useDSA                 ) &&
            ( dPhiNorm < s_DSAabsTol )
        ) {

            printf( "Turning off DSA (iter = %" PRId64", dPhiNorm = %.4e)\n", iter, dPhiNorm );
            fprintf( g_outfp, "Turning off DSA (iter = %" PRId64", dPhiNorm = %.4e)\n", iter, dPhiNorm );

            useDSA = false;

        }


        // --- Compute next iterate. ------------------------------------------------------------------ //

        // Q.
        memcpy( psi_n1, Q, g_sizeofPsi );

        if ( !useDSA ) {

            memset( phi_one, 0, g_sizeofPhi );

            // PΨⁿ
            #pragma omp parallel for
            for ( int64_t i = 1 ; i <= g_nx      ; ++i ) {
            for ( int64_t d = 0 ; d <= g_dgOrder ; ++d ) {
            for ( int64_t q = 0 ; q <  g_nq      ; ++q ) {

                phi_one[ IPHI(i,d) ] += g_w[q] * psi_n[ IPSI(i,q,d) ];

            }}}

        }

        // SPΨⁿ + Q
        #pragma omp parallel for private( tempVal )
        for ( int64_t i = 1 ; i <= g_nx      ; ++i ) {
        for ( int64_t d = 0 ; d <= g_dgOrder ; ++d ) {

            tempVal = 0.0;

            for ( int64_t a = 0 ; a <= g_dgOrder ; ++a ) {
            for ( int64_t u = 0 ; u <= g_dgOrder ; ++u ) {

                tempVal += (2*d + 1) / 4.0 * g_TPI[ ITPI(u,a,d) ]
                           * s_sigmaS[ IPHI(i,u) ] * phi_one[ IPHI(i,a) ];

            }}

            for ( int64_t q = 0 ; q < g_nq ; ++q ) { psi_n1[ IPSI(i,q,d) ] += tempVal; }

        }}

        // L⁻¹ ( SPΨⁿ + Q ).
        g_solverTime.stop();
        sweep( psi_n1, s_sigmaT );
        g_solverTime.start();


        // --- Measure error between successive iterates. --------------------------------------------- //

        error = 0.0;

        for ( int64_t i = 1 ; i <= g_nx      ; ++i ) {
        for ( int64_t d = 0 ; d <= g_dgOrder ; ++d ) {
        for ( int64_t q = 0 ; q <  g_nq      ; ++q ) {

            tempVal = psi_n[ IPSI(i,q,d) ] - psi_n1[ IPSI(i,q,d) ];

            error += tempVal * tempVal;

        }}}

        error = sqrt( error );

        if ( !(iter % 10) ) {

            printf( "Iteration:  %5" PRId64"       Error: %9.4e\t", iter, error );
            fprintf( g_outfp, "Iteration:  %5" PRId64"       Error: %9.4e\t", iter, error );

            if ( useDSA ) {

                printf( "dPhiNorm = %.4e", dPhiNorm );
                fprintf( g_outfp, "dPhiNorm = %.4e", dPhiNorm );

            }

            printf( "\n" );
            fprintf( g_outfp, "\n" );

        }


        // --- If using DSA, solve for and apply the diffusion correction. ---------------------------- //

        if ( useDSA ) {

            // Save PΨⁿ⁺¹ to phi_half to be corrected later.
            memset( phi_half, 0, g_sizeofPhi );

            #pragma omp parallel for
            for ( int64_t i = 1 ; i <= g_nx      ; ++i ) {
            for ( int64_t d = 0 ; d <= g_dgOrder ; ++d ) {
            for ( int64_t q = 0 ; q <  g_nq      ; ++q ) {

                phi_half[ IPHI(i,d) ] += g_w[q] * psi_n1[ IPSI(i,q,d) ];

            }}}

            // Compute ( phi_half - phi_one ) for right hand side of diffusion solve.
            #pragma omp parallel for
            for ( int64_t i = 1 ; i <= g_nx      ; ++i ) {
            for ( int64_t d = 0 ; d <= g_dgOrder ; ++d ) {

                phi_one[ IPHI(i,d) ] = phi_half[ IPHI(i,d) ] - phi_one[ IPHI(i,d) ];

            }}

            // Compute right hand vector for the diffusion system.
            VecGetArray( s_v2, &rhs2 );
            D_computeRHS( phi_one, rhs2 );
            VecRestoreArray( s_v2, &rhs2 );

            // Compute the matrix for the diffusion system.
            setupD_matrix( &s_D2 );

            // Solve the diffusion system.
            KSPSolve( s_ksp2, s_v2, s_u2 );

            // Get information about the diffusion solve.
            KSPGetIterationNumber( s_ksp2, &its );
            KSPGetResidualNorm( s_ksp2, &rnorm );
            KSPGetConvergedReason( s_ksp2, &reason );

            // Print information about diffusion solve.
            printf(
                "D Iterations: %5d       Final Residual: %9.4e       Stop Reason: %d\n",
                its, rnorm, reason
            );
            fprintf(
                g_outfp,
                "D Iterations: %5d       Final Residual: %9.4e       Stop Reason: %d\n",
                its, rnorm, reason
            );

            if ( reason < 0 ) {

                fprintf( stderr, __RED "ERROR" __RESET " at %s:%d -- Diffusion solve failed to converge.\n",
                         __FILE__, __LINE__ );
                fprintf( g_outfp, "ERROR at %s:%d -- Diffusion solve failed to converge.\n",
                         __FILE__, __LINE__ );

                exit( EXIT_FAILURE );

            }

            // Use result from the diffusion solve to correct phi_half (stored in phi_one).
            memset( phi_one, 0, g_sizeofPhi );
            VecGetArray( s_u2, &phi2 );

            for ( int64_t i = 1 ; i <= g_nx      ; ++i ) {
            for ( int64_t d = 0 ; d <= g_dgOrder ; ++d ) {

                phi_one[ IPHI(i,d) ] = phi_half[ IPHI(i,d) ] + phi2[ IPHI2(i-1,d) ];

            }}

            dPhiNorm = 0.0;

            for ( int64_t i = 1 ; i <= g_nx      ; ++i ) {
            for ( int64_t d = 0 ; d <= g_dgOrder ; ++d ) {

                tempVal = phi2[ IPHI2(i-1,d) ];

                dPhiNorm += tempVal * tempVal;

            }}

            dPhiNorm = sqrt( dPhiNorm );

            VecRestoreArray( s_u2, &phi2 );

        }

    }


    // --- Source iterations finished. ---------------------------------------------------------------- //

    memcpy( Q, psi_n1, g_sizeofPsi );


    // Update cumulative iteration counters.
    s_totalIts += iter;
    ++s_totalSolves;


    // Print information about solve.
    if ( s_useDSA ) {

        printf( "\n" );
        fprintf( g_outfp, "\n" );
    }

    printf( "Iterations: %5" PRId64"       Final Error: %9.4e\n", iter, error );
    fprintf( g_outfp, "Iterations: %5" PRId64"       Final Error: %9.4e\n", iter, error );

    if ( s_useDSA ) {

        printf( "\n" );
        fprintf( g_outfp, "\n" );
    }

}

# endif // if 0



//============================================================================================================
//=== INITIALIZATION AND CLEANUP ROUTINES ====================================================================
//============================================================================================================

# if 0
//------------------------------------------------------------------------------------------------------------
// initSolver
//      Initialization routine for the collisional solvers.
//
// INPUT:
//      argc        - argc from int main (used to activate PETSc command line switches).
//      argv        - argv from int main (used to activate PETSc command line switches).
//      KSP_rtol    - PETSc tolerances for the GMRES transport solve. Also used to set the tolerances for the
//      KSP_abstol  : SI transport solve.
//      KSP_dtol    :
//      KSP_maxits  :
//      useDSA      - Boolean value determining whether DSA should be applied.
//      DSA_rtol    - PETSc tolerances for the SIP-DG diffusion solve.
//      DSA_abstol  : NOTE: This parameter is also used as a tolerance to turn off DSA for the source
//                  :       iteration transport solver.
//      DSA_dtol    :
//      DSA_maxits  :
//
// RETURN:
//      Returns 0 upon successful execution.
//------------------------------------------------------------------------------------------------------------
int initSolver (

    int argc,
    char **argv,
    const double KSP_rtol,
    const double KSP_abstol,
    const double KSP_dtol,
    const int KSP_maxits,
    const bool useDSA,
    const double DSA_rtol,
    const double DSA_abstol,
    const double DSA_dtol,
    const int DSA_maxits
) {

    PetscMPIInt size;
    PetscScalar one;
    void *ctx;
    PetscErrorCode ierr;
    PetscInt n;
    PetscBool dsa = PETSC_FALSE;

    KSPType   t_solver,  d_solver;
    PetscReal t_rtol,    d_rtol;
    PetscReal t_abstol,  d_abstol;
    PetscReal t_dtol,    d_dtol;
    PetscInt  t_maxits,  d_maxits;
    PetscInt  t_restart;


    // Initialize iteration counters.
    s_totalIts    = 0;
    s_totalSolves = 0;

    // Transfer value of DSA boolean flag.
    s_useDSA = useDSA;

    // Allocate memory for auxiliary variables.
    s_sizeofPhi2 =        (g_nx)     * (g_dgOrder + 1) * sizeof(double);
    s_sizeofPsi2 = g_nq * (g_nx)     * (g_dgOrder + 1) * sizeof(double);

    s_1tempPsi = (double*) malloc( g_sizeofPsi );

    memset( s_1tempPsi, 0, g_sizeofPsi );

    if ( g_solverType == SOLVER_SI ) {

        s_2tempPsi = (double*) malloc( g_sizeofPsi );
        s_1tempPhi = (double*) malloc( g_sizeofPhi );
        s_2tempPhi = (double*) malloc( g_sizeofPhi );

        memset( s_2tempPsi, 0, g_sizeofPsi );
        memset( s_1tempPhi, 0, g_sizeofPhi );
        memset( s_2tempPhi, 0, g_sizeofPhi );

    }


    // PETSc initialization
    ctx  = NULL;
    ierr = 0;
    one  = 1.0;
    n    = (g_dgOrder + 1) * (g_nx);

    PetscInitialize( &argc, &argv, NULL, NULL );
    ierr = MPI_Comm_size( PETSC_COMM_WORLD, &size );    CHKERRQ(ierr);

    if ( size != 1 ) { SETERRQ( PETSC_COMM_WORLD, 1, "This is a uniprocessor example only!" ); }


    // --- Setup GMRES Krylov solver. ----------------------------------------------------------------- //

    // Setup vectors.
    ierr = VecCreate( PETSC_COMM_WORLD, &s_x1 );    CHKERRQ(ierr);
    ierr = VecSetSizes( s_x1, PETSC_DECIDE, n );    CHKERRQ(ierr);
    ierr = VecSetFromOptions( s_x1 );               CHKERRQ(ierr);
    ierr = VecDuplicate( s_x1, &s_b1 );             CHKERRQ(ierr);

    // Set initial value for vectors.
    ierr = VecSet( s_b1, one );     CHKERRQ(ierr);
    ierr = VecSet( s_x1, one );     CHKERRQ(ierr);

    // Setup shell matrix.
    ierr = MatCreateShell( PETSC_COMM_WORLD, n, n, PETSC_DETERMINE, PETSC_DETERMINE, ctx, &s_A1 );
           CHKERRQ(ierr);
    ierr = MatShellSetOperation( s_A1, MATOP_MULT, (void(*)(void)) GMRES_computeLHS );
           CHKERRQ(ierr);

    // Create and setup KSP context.
    ierr = KSPCreate( PETSC_COMM_WORLD, &s_ksp1 );      CHKERRQ(ierr);
    ierr = KSPSetOperators( s_ksp1, s_A1, s_A1 );       CHKERRQ(ierr);

    // Function KSPSetOperators changed in version 3.5.x. Additional argument below no longer needed.
    //, DIFFERENT_NONZERO_PATTERN );    CHKERRQ(ierr);

    // Set solver type and basic options.
    ierr = KSPSetType( s_ksp1, KSPLGMRES );                                         CHKERRQ(ierr);
    ierr = KSPSetInitialGuessNonzero( s_ksp1, PETSC_TRUE );                         CHKERRQ(ierr);
    ierr = KSPSetTolerances( s_ksp1, KSP_rtol, KSP_abstol, KSP_dtol, KSP_maxits );     CHKERRQ(ierr);

    // Get preconditioner context to setup DSA.
    ierr = KSPGetPC( s_ksp1, &s_pc1 );      CHKERRQ(ierr);

    // Check for -dsa switch to turn on DSA preconditioner for GMRES transport solve.
    PetscOptionsGetBool( NULL, "-dsa", &dsa, NULL );

    // Use dsa method if -dsa switch found or useDSA variable set.
    if ( dsa || useDSA ) {

        PCSetType( s_pc1, PCSHELL );
        PCShellSetApply( s_pc1, computeDSA );
        PCShellSetName( s_pc1, "DSA" );

    // Otherwise do not use preconditioner.
    } else {

        PCSetType( s_pc1, PCNONE );

    }

    // Set solver options from command line arguments.
    KSPSetFromOptions( s_ksp1 );


    // --- Setup diffusion solver for DSA. ------------------------------------------------------------ //

    // Setup vectors.
    ierr = VecCreate( PETSC_COMM_WORLD, &s_u2 );    CHKERRQ(ierr);
    ierr = VecSetSizes( s_u2, PETSC_DECIDE, n );    CHKERRQ(ierr);
    ierr = VecSetFromOptions( s_u2 );               CHKERRQ(ierr);
    ierr = VecDuplicate( s_u2, &s_v2 );             CHKERRQ(ierr);

    // Setup matrix.
    MatCreate( PETSC_COMM_WORLD, &s_D2 );
    MatSetType( s_D2, MATSEQAIJ );
    MatCreateSeqAIJ( PETSC_COMM_WORLD, n, n, 3 * (g_dgOrder + 1), NULL, &s_D2 );
    MatSetOption( s_D2, MAT_IGNORE_ZERO_ENTRIES, PETSC_TRUE );

    // Set initial value for vectors.
    ierr = VecSet( s_u2, one );                     CHKERRQ(ierr);
    ierr = VecSet( s_x1, one );                     CHKERRQ(ierr);

    // Create and setup KSP context.
    ierr = KSPCreate( PETSC_COMM_WORLD, &s_ksp2 );  CHKERRQ(ierr);
    ierr = KSPSetOperators( s_ksp2, s_D2, s_D2 );   CHKERRQ(ierr);
    ierr = KSPSetType( s_ksp2, KSPCG );             CHKERRQ(ierr);

    ierr = KSPSetInitialGuessNonzero( s_ksp2, PETSC_TRUE );     CHKERRQ(ierr);
    ierr = KSPSetTolerances( s_ksp2, DSA_rtol, DSA_abstol, DSA_dtol, DSA_maxits );     CHKERRQ(ierr);


    // --- Seup source iteration (SI) solver. --------------------------------------------------------- //

    s_SImaxIts  = KSP_maxits;
    s_SItol     = KSP_abstol;
    s_DSAabsTol = DSA_abstol;


    // --- Print solver information. ------------------------------------------------------------------ //

    KSPGetTolerances( s_ksp1, &t_rtol, &t_abstol, &t_dtol, &t_maxits );
    KSPGetTolerances( s_ksp2, &d_rtol, &d_abstol, &d_dtol, &d_maxits );

    KSPGetType( s_ksp1, &t_solver );
    KSPGetType( s_ksp2, &d_solver );

    KSPGMRESGetRestart( s_ksp1, &t_restart );

    switch ( g_solverType ) {

        case SOLVER_SWEEP:

            printf( "  %-20s %s\n", "Solver Type:", " sweep only" );
            fprintf( g_outfp, "  %-20s %s\n", "Solver Type:", " sweep only" );

            break;

        case SOLVER_GMRES:

            // Solver type for transport solve.
            printf( "  %-20s  %s\n", "Solver Type:", t_solver );
            fprintf( g_outfp, "  %-20s  %s\n", "Solver Type:", t_solver );

            // Tolerances for transport solve.
            printf( "  %-20s % .4e\n", "Transp. rtol:",    t_rtol    );
            printf( "  %-20s % .4e\n", "Transp. abstol:",  t_abstol  );
            printf( "  %-20s % .4e\n", "Transp. dtol:",    t_dtol    );
            printf( "  %-20s % d\n",   "Transp. maxits:",  t_maxits  );
            printf( "  %-20s % d\n",   "Transp. restart:", t_restart );

            fprintf( g_outfp, "  %-20s % .4e\n", "Transp. rtol:",    t_rtol    );
            fprintf( g_outfp, "  %-20s % .4e\n", "Transp. abstol:",  t_abstol  );
            fprintf( g_outfp, "  %-20s % .4e\n", "Transp. dtol:",    t_dtol    );
            fprintf( g_outfp, "  %-20s % d\n",   "Transp. maxits:",  t_maxits  );
            fprintf( g_outfp, "  %-20s % d\n",   "Transp. restart:", t_restart );

            printf( "\n" );
            fprintf( g_outfp, "\n" );

            if ( s_useDSA ) {

                // Use DSA preconditioner.
                printf( "  %-20s %s\n", "Preconditioner:", " DSA" );
                fprintf( g_outfp, "  %-20s %s\n", "Preconditioner:", " DSA" );

                // Solver type for diffusion solve.
                printf( "  %-20s  %s\n", "Solver Type:", d_solver );
                fprintf( g_outfp, "  %-20s  %s\n", "Solver Type:", d_solver );

                // Tolerances for diffusion solve.
                printf( "  %-20s % .4e\n", "Diff. rtol:",   d_rtol   );
                printf( "  %-20s % .4e\n", "Diff. abstol:", d_abstol );
                printf( "  %-20s % .4e\n", "Diff. dtol:",   d_dtol   );
                printf( "  %-20s % d\n",   "Diff. maxits:", d_maxits );

                fprintf( g_outfp, "  %-20s % .4e\n", "Diff. rtol:",   d_rtol   );
                fprintf( g_outfp, "  %-20s % .4e\n", "Diff. abstol:", d_abstol );
                fprintf( g_outfp, "  %-20s % .4e\n", "Diff. dtol:",   d_dtol   );
                fprintf( g_outfp, "  %-20s % d\n",   "Diff. maxits:", d_maxits );


            } else {

                // Do not use any preconditioner.
                printf( "  %-20s %s\n", "Preconditioner:", " none" );
                fprintf( g_outfp, "  %-20s %s\n", "Preconditioner:", " none" );

            }

            break;

        case SOLVER_SI:

            // Solver type for transport solve.
            printf( "  %-20s %s\n", "Solver Type:", " source iteration" );
            fprintf( g_outfp, "  %-20s %s\n", "Solver Type:", " source iteration" );

            // Tolerances for transport solve.
            printf( "  %-20s % .4e\n", "Transp. tol:", s_SItol );
            printf( "  %-20s % d\n", "Transp. maxits:", s_SImaxIts );

            fprintf( g_outfp, "  %-20s % .4e\n", "Transp. tol:", s_SItol );
            fprintf( g_outfp, "  %-20s % d\n", "Transp. maxits:", s_SImaxIts );

            printf( "\n" );
            fprintf( g_outfp, "\n" );

            if ( s_useDSA ) {

                // Use DSA preconditioner.
                printf( "  %-20s %s\n", "Preconditioner:", " DSA" );
                fprintf( g_outfp, "  %-20s %s\n", "Preconditioner:", " DSA" );

                // Solver type for diffusion solve.
                printf( "  %-20s  %s\n", "Solver Type:", d_solver );
                fprintf( g_outfp, "  %-20s  %s\n", "Solver Type:", d_solver );

                // Tolerances for diffusion solve.
                printf( "  %-20s % .4e\n", "Diff. rtol:",   d_rtol   );
                printf( "  %-20s % .4e\n", "Diff. abstol:", d_abstol );
                printf( "  %-20s % .4e\n", "Diff. dtol:",   d_dtol   );
                printf( "  %-20s % d\n",   "Diff. maxits:", d_maxits );

                fprintf( g_outfp, "  %-20s % .4e\n", "Diff. rtol:",   d_rtol   );
                fprintf( g_outfp, "  %-20s % .4e\n", "Diff. abstol:", d_abstol );
                fprintf( g_outfp, "  %-20s % .4e\n", "Diff. dtol:",   d_dtol   );
                fprintf( g_outfp, "  %-20s % d\n",   "Diff. maxits:", d_maxits );


            } else {

                // Do not use any preconditioner.
                printf( "  %-20s %s\n", "Preconditioner:", " none" );
                fprintf( g_outfp, "  %-20s %s\n", "Preconditioner:", " none" );

            }

            break;

        default:

            fprintf( stderr, __RED "ERROR" __RESET " at %s:%d -- Invalid value for solver type.\n",
                     __FILE__, __LINE__ );
            fprintf( g_outfp, "ERROR at %s:%d -- Invalid value for solver type.\n",
                     __FILE__, __LINE__ );

            exit( EXIT_FAILURE );

    }


    return 0;

}
# endif // if 0


# endif // if SPACE_DIMS == 1 && defined (ENABLE_PETSC) && !defined (ENABLE_MPI)

# endif // if defined (DONT_COMPILE_THIS)
