//------------------------------------------------------------------------------------------------------------
//! \file   linear_solvers/STDG/DSA_1d.cpp
//! \brief  Implements a diffusion synthetic acceleration (DSA) preconditioner for the transport system using
//!         a modified interior penalty (MIP) construction.
//!
//! \todo   Citations for MIP DSA method.
//!
//! \authors Michael Crockatt
//! \date    July 2017
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

# include "linear_solvers/STDG/ImplicitSolver_old.hpp"
# include "utils/CLog.hpp"
# include "utils/global.hpp"


namespace STDG {

    //!
    //! \brief  Global variable used to store the solver class during an implicit solve.
    //!
    //! Necessary for accessing the solver class from STDG::ImplicitSolverOLD::DSA_Apply(), which is called by
    //! PETSc.
    //!
    //! Variable is defined in \ref linear_solvers/STDG/ImplicitSolverOLD.cpp
    //!
    extern ImplicitSolverOLD * s_this;
}


static const double MIN_CROSS = 1.0e-3;


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
double STDG::ImplicitSolverOLD::Kappa (

    const int64_t i,
    const int64_t j,
    const RKDG::CrossSection & sigma_t

) const {

    double temp_val = 2.0 * DG_degree_x * (DG_degree_x + 1)
                      / ( 3.0 * dx(0) * std::max( MIN_CROSS, ( sigma_t(i,0) + sigma_t(j,0) ) / 2.0 ) );

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
//! \attention  This routine uses the indexing function associated with the first STDG::ScalarFlux object
//!             stored at ImplicitSolverOLD::phi. Hence the ScalarFlux objects must be allocated to this pointer
//!             _before_ calling this routine.
//!
//! \param[in]  sigma_t     Total cross section \f$ \sigma_{\mathrm{t}} \f$.
//! \param[in]  sigma_s     Scattering cross section \f$ \sigma_{\mathrm{s}} \f$.
//! \param[in]  dt          Timestep size. For steady-state problems, use \pp{dt} = \c inf.
//------------------------------------------------------------------------------------------------------------
void STDG::ImplicitSolverOLD::ComputeDiffusionMatrix (

    const RKDG::CrossSection & sigma_t,
    const RKDG::CrossSection & sigma_s,
    const double dt
) {

    Global::TMR_DSA_Assemble.Start();

    MatZeroEntries( this->D );

    double coeff, val;
    PetscInt m_row, m_col;

    // Set values for the matrix.
    for ( PetscInt i = 0; i <  this->nx(0);       ++i ) {
    for ( PetscInt d = 0; d <= this->DG_degree_x; ++d ) {
    for ( PetscInt s = 0; s <= this->DG_degree_t; ++s ) {

        m_row = this->phi->IndexNoGC(i,d,s);

        // --- Cell integral terms. ------------------------------------------------------------------- //

        m_col = this->phi->IndexNoGC(i,d,s);
        val = this->dx(0) * ( sigma_t(i+1,0) - sigma_s(i+1,0) ) / ( (2*d + 1) * (2*s + 1) );
        MatSetValues( this->D, 1, &m_row, 1, &m_col, &val, ADD_VALUES );

        coeff = 2.0 / ( 3.0 * this->dx(0) * std::max( MIN_CROSS, sigma_t(i+1,0) ) );

        for ( PetscInt a = 0; a <= this->DG_degree_x; ++a ) {

            m_col = this->phi->IndexNoGC(i,a,s);
            val = coeff * B(d,a) / (2*s + 1);
            MatSetValues( this->D, 1, &m_row, 1, &m_col, &val, ADD_VALUES );
        }

        // --- Temporal derivative terms. ------------------------------------------------------------- //

        coeff = this->dx(0) / ( this->dt * (2*d + 1) );

        for ( PetscInt r = 0; r <= this->DG_degree_t; ++r ) {

            m_col = this->phi->IndexNoGC(i,d,r);
            MatSetValues( this->D, 1, &m_row, 1, &m_col, &coeff, ADD_VALUES );
        }

        coeff = -2.0 * this->dx(0) / ( this->dt * (2*d + 1) );

        for ( PetscInt r = s-1; r >= 0; r -= 2 ) {

            m_col = this->phi->IndexNoGC(i,d,r);
            MatSetValues( this->D, 1, &m_row, 1, &m_col, &coeff, ADD_VALUES );
        }

        // --- Cell interface (inter-cell) flux terms. ------------------------------------------------ //

        // Left boundary cell.
        if ( i == 0 ) {

            // Penalty flux.
            for ( PetscInt a = 0; a <= this->DG_degree_x; ++a ) {

                m_col = this->phi->IndexNoGC(i+1,a,s);
                val = neg1pow(a+1) * Kappa( i+1, i+2, sigma_t ) / (2*s + 1);
                MatSetValues( this->D, 1, &m_row, 1, &m_col, &val, ADD_VALUES );

                m_col = this->phi->IndexNoGC(i,a,s);
                val = Kappa( i+1, i+1, sigma_t ) / (2*s + 1);
                MatSetValues( this->D, 1, &m_row, 1, &m_col, &val, ADD_VALUES );
            }

            // Test normal derivative flux.
            coeff = d * (d + 1) / ( 6.0 * (2*s + 1) * this->dx(0) );

            for ( PetscInt a = 0; a <= this->DG_degree_x; ++a ) {

                m_col = this->phi->IndexNoGC(i,a,s);
                val = - coeff / std::max( MIN_CROSS, sigma_t(i+1,0) );
                MatSetValues( this->D, 1, &m_row, 1, &m_col, &val, ADD_VALUES );

                m_col = this->phi->IndexNoGC(i+1,a,s);
                val = coeff * neg1pow(a) / std::max( MIN_CROSS, ( sigma_t(i+1,0) + sigma_t(i+2,0) ) / 2.0 );
                MatSetValues( this->D, 1, &m_row, 1, &m_col, &val, ADD_VALUES );
            }

            // Unknown normal derivative flux.
            coeff = - 1.0 / ( 6.0 * (2*s + 1) * this->dx(0) );

            for ( PetscInt a = 0; a <= this->DG_degree_x; ++a ) {

                m_col = this->phi->IndexNoGC(i+1,a,s);
                val = coeff * neg1pow(a) * a * (a + 1)
                      / std::max( MIN_CROSS, ( sigma_t(i+1,0) + sigma_t(i+2,0) ) / 2.0 );
                MatSetValues( this->D, 1, &m_row, 1, &m_col, &val, ADD_VALUES );

                m_col = this->phi->IndexNoGC(i,a,s);
                val = coeff * a * (a + 1) / std::max( MIN_CROSS, sigma_t(i+1,0) );
                MatSetValues( this->D, 1, &m_row, 1, &m_col, &val, ADD_VALUES );
            }

        // Right boundary cell
        } else if ( i == this->nx(0) -1 ) {

            // Penalty flux.
            for ( PetscInt a = 0; a <= this->DG_degree_x; ++a ) {

                m_col = this->phi->IndexNoGC(i-1,a,s);
                val = neg1pow(d+1) * Kappa( i, i+1, sigma_t ) / (2*s + 1);
                MatSetValues( this->D, 1, &m_row, 1, &m_col, &val, ADD_VALUES );

                m_col = this->phi->IndexNoGC(i,a,s);
                val = neg1pow(d+a) * Kappa( i+1, i+1, sigma_t ) / (2*s + 1);
                MatSetValues( this->D, 1, &m_row, 1, &m_col, &val, ADD_VALUES );
            }

            // Test normal derivative flux.
            coeff = neg1pow(d) * d * (d + 1) / ( 6.0 * (2*s + 1) * this->dx(0) );

            for ( PetscInt a = 0; a <= this->DG_degree_x; ++a ) {

                m_col = this->phi->IndexNoGC(i-1,a,s);
                val = -coeff / std::max( MIN_CROSS, ( sigma_t(i,0) + sigma_t(i+1,0) ) / 2.0 );
                MatSetValues( this->D, 1, &m_row, 1, &m_col, &val, ADD_VALUES );

                m_col = this->phi->IndexNoGC(i,a,s);
                val = coeff * neg1pow(a) / std::max( MIN_CROSS, sigma_t(i+1,0) );
                MatSetValues( this->D, 1, &m_row, 1, &m_col, &val, ADD_VALUES );
            }

            // Unknown normal derivative flux.
            coeff = 1.0 / ( 6.0 * (2*s + 1) * this->dx(0) );

            for ( PetscInt a = 0; a <= this->DG_degree_x; ++a ) {

                m_col = this->phi->IndexNoGC(i,a,s);
                val = coeff * neg1pow(d+a) * a * (a + 1) / std::max( MIN_CROSS, sigma_t(i+1,0) );
                MatSetValues( this->D, 1, &m_row, 1, &m_col, &val, ADD_VALUES );

                m_col = this->phi->IndexNoGC(i-1,a,s);
                val = coeff * neg1pow(d) * a * (a + 1)
                      / std::max( MIN_CROSS, ( sigma_t(i,0) + sigma_t(i+1,0) ) / 2.0 );
                MatSetValues( this->D, 1, &m_row, 1, &m_col, &val, ADD_VALUES );
            }

        // Interior cells.
        } else {

            // Penalty flux.
            for ( PetscInt a = 0; a <= this->DG_degree_x; ++a ) {

                m_col = this->phi->IndexNoGC(i-1,a,s);
                val = Kappa( i, i+1, sigma_t ) * neg1pow(d+1) / (2*s + 1);
                MatSetValues( this->D, 1, &m_row, 1, &m_col, &val, ADD_VALUES );

                m_col = this->phi->IndexNoGC(i,a,s);
                val = Kappa( i+1, i+1, sigma_t ) * ( 1 + neg1pow(d+a) ) / (2*s + 1);
                MatSetValues( this->D, 1, &m_row, 1, &m_col, &val, ADD_VALUES );

                m_col = this->phi->IndexNoGC(i+1,a,s);
                val = Kappa( i+1, i+2, sigma_t ) * neg1pow(a+1) / (2*s + 1);
                MatSetValues( this->D, 1, &m_row, 1, &m_col, &val, ADD_VALUES );
            }

            // Test normal derivative flux.
            coeff = d * (d + 1) / ( 6.0 * (2*s + 1) * this->dx(0) );

            for ( PetscInt a = 0; a <= this->DG_degree_x; ++a ) {

                m_col = this->phi->IndexNoGC(i-1,a,s);
                val = coeff * neg1pow(d+1) / std::max( MIN_CROSS, ( sigma_t(i,0) + sigma_t(i+1,0) ) / 2.0 );
                MatSetValues( this->D, 1, &m_row, 1, &m_col, &val, ADD_VALUES );

                m_col = this->phi->IndexNoGC(i,a,s);
                val = coeff * ( neg1pow(d+a) - 1 ) / std::max( MIN_CROSS, sigma_t(i+1,0) );
                MatSetValues( this->D, 1, &m_row, 1, &m_col, &val, ADD_VALUES );

                m_col = this->phi->IndexNoGC(i+1,a,s);
                val = coeff * neg1pow(a) / std::max( MIN_CROSS, ( sigma_t(i+1,0) + sigma_t(i+2,0) ) / 2.0 );
                MatSetValues( this->D, 1, &m_row, 1, &m_col, &val, ADD_VALUES );
            }

            // Unknown normal derivative flux.
            coeff = 1.0 / ( 6.0 * (2*s + 1) * this->dx(0) );

            for ( PetscInt a = 0; a <= this->DG_degree_x; ++a ) {

                m_col = this->phi->IndexNoGC(i-1,a,s);
                val = coeff * neg1pow(d) * a * (a + 1)
                      / std::max( MIN_CROSS, ( sigma_t(i,0) + sigma_t(i+1,0) ) / 2.0 );
                MatSetValues( this->D, 1, &m_row, 1, &m_col, &val, ADD_VALUES );

                m_col = this->phi->IndexNoGC(i,a,s);
                val = coeff * ( neg1pow(d+a) - 1 ) * a * (a + 1) / std::max( MIN_CROSS, sigma_t(i+1,0) );
                MatSetValues( this->D, 1, &m_row, 1, &m_col, &val, ADD_VALUES );

                m_col = this->phi->IndexNoGC(i+1,a,s);
                val = coeff * neg1pow(a+1) * a * (a + 1)
                      / std::max( MIN_CROSS, ( sigma_t(i+1,0) + sigma_t(i+2,0) ) / 2.0 );
                MatSetValues( this->D, 1, &m_row, 1, &m_col, &val, ADD_VALUES );
            }
        }

        // --- Boundary flux terms. ------------------------------------------------------------------- //

        // Left boundary.
        if ( i == 0 ) {

            // Penalty flux.
            for ( PetscInt a = 0; a <= this->DG_degree_x; ++a ) {

                m_col = this->phi->IndexNoGC(i,a,s);
                val = Kappa( i+1, i+1, sigma_t ) * neg1pow(d+a) / (2*s + 1);
                MatSetValues( this->D, 1, &m_row, 1, &m_col, &val, ADD_VALUES );
            }

            // Test normal derivative flux.
            coeff = neg1pow(d+1) * d * (d + 1)
                    / ( 6.0 * (2*s + 1) * this->dx(0) * std::max( MIN_CROSS, sigma_t(i+1,0) ) );

            for ( PetscInt a = 0; a <= this->DG_degree_x; ++a ) {

                m_col = this->phi->IndexNoGC(i,a,s);
                val = coeff * neg1pow(a);
                MatSetValues( this->D, 1, &m_row, 1, &m_col, &val, ADD_VALUES );
            }

            // Unknown normal derivative flux.
            coeff = neg1pow(d+1) / ( 6.0 * (2*s + 1) * this->dx(0) * std::max( MIN_CROSS, sigma_t(i+1,0) ) );

            for ( PetscInt a = 0; a <= this->DG_degree_x; ++a ) {

                m_col = this->phi->IndexNoGC(i,a,s);
                val = coeff * neg1pow(a) * a * (a + 1);
                MatSetValues( this->D, 1, &m_row, 1, &m_col, &val, ADD_VALUES );
            }

        // Right boundary
        } else if ( i == this->nx(0) -1 ) {

            // Penalty flux.
            for ( PetscInt a = 0; a <= this->DG_degree_x; ++a ) {

                m_col = this->phi->IndexNoGC(i,a,s);
                val = Kappa( i+1, i+1, sigma_t ) / (2*s + 1);
                MatSetValues( this->D, 1, &m_row, 1, &m_col, &val, ADD_VALUES );
            }

            // Test normal derivative flux.
            coeff = -d * (d + 1) / ( 6.0 * (2*s + 1) * this->dx(0) * std::max( MIN_CROSS, sigma_t(i+1,0) ) );

            for ( PetscInt a = 0; a <= this->DG_degree_x; ++a ) {

                m_col = this->phi->IndexNoGC(i,a,s);
                val = coeff;
                MatSetValues( this->D, 1, &m_row, 1, &m_col, &val, ADD_VALUES );
            }

            // Unknown normal derivative flux.
            coeff = -1.0 / ( 6.0 * (2*s + 1) * this->dx(0) * std::max( MIN_CROSS, sigma_t(i+1,0) ) );

            for ( PetscInt a = 0; a <= this->DG_degree_x; ++a ) {

                m_col = this->phi->IndexNoGC(i,a,s);
                val = coeff * a * (a + 1);
                MatSetValues( this->D, 1, &m_row, 1, &m_col, &val, ADD_VALUES );
            }
        }
    }}}

    MatAssemblyBegin( this->D, MAT_FINAL_ASSEMBLY );
    MatAssemblyEnd( this->D, MAT_FINAL_ASSEMBLY );

    Global::TMR_DSA_Assemble.Stop();
}


//------------------------------------------------------------------------------------------------------------
//! \brief   Computes the RHS vector for the DSA diffusion system.
//!
//! The RHS vector is stored in ImplicitSolverOLD::v.
//------------------------------------------------------------------------------------------------------------
void STDG::ImplicitSolverOLD::DSA_ComputeRHS (

    const STDG::ScalarFlux & phi,
    STDG::ScalarFlux & Q
) {

    Q.ZeroDensity();

    Global::TMR_DSA_Solve.Start();

    # pragma omp parallel for
    for ( int64_t i = 0; i <= Q.nx(0) + 1;   ++i ) {
    for ( int64_t d = 0; d <= Q.DG_degree_x; ++d ) {
    for ( int64_t s = 0; s <= Q.DG_degree_t; ++s ) {

        Q(i,d,s) = Q.dx(0) * phi(i,d,s) * (*sigma_s)(i,0) / ( (2*d + 1) * (2*s + 1) );
    }}}

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
//! \see    STDG::ImplicitSolverOLD::ComputeDiffusionMatrix()
//! \see    STDG::ImplicitSolverOLD::DSA_ComputeRHS()
//------------------------------------------------------------------------------------------------------------
PetscErrorCode STDG::ImplicitSolverOLD::DSA_Apply (

    PC pc,
    Vec x,
    Vec y
) {

    Global::TMR_PETSc.Stop();

    STDG::ScalarFlux & phi1 = s_this->phi[0];
    STDG::ScalarFlux & phi2 = s_this->phi[1];

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


# endif // if SPACE_DIMS == 1 && defined (ENABLE_PETSC) && !defined (ENABLE_MPI)

# endif // if defined (DONT_COMPILE_THIS)
