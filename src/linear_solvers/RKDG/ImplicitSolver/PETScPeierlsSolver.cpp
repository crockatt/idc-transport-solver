//------------------------------------------------------------------------------------------------------------
//! \file   linear_solvers/RKDG/ImplicitSolver/PETScPeierlsSolver.hpp
//! \brief  Contains member function specializations for
//!         SweepOperator<RKDG::OrdinateFlux>::PETScPeierlsSolver.
//!
//! \author Michael M. Crockatt
//! \date   February 2018
//------------------------------------------------------------------------------------------------------------

# include "linear_solvers/Abstract/ImplicitSolver/PETScPeierlsSolver.hpp"


# if SPACE_DIMS == 1  && defined (ENABLE_PETSC)


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
template<>
double PETScPeierlsSolver<RKDG::OrdinateFlux>::Kappa (

    const int64_t i,
    const int64_t j,
    const CrossSection & sigma_t,
    const double dt

) const {

    const double DG_degree = this->phi.at(0)->DG_degree;

    double temp_val = 4.0 * DG_degree * (DG_degree + 1)
                      / ( 3.0 * this->dx(0) * ( sigma_t(i,0) + sigma_t(j,0) + 2.0 / dt ) );

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
template< typename IntegerType >
IntegerType B (

    const IntegerType a,
    const IntegerType d
) {

    IntegerType beta = std::min( a, d );

    return ( (beta + 1) * beta ) * ((IntegerType) ((a % 2) == (d % 2)));
}


//------------------------------------------------------------------------------------------------------------
//! \brief  Computes the coefficients of the matrix for the diffusion solve in the DSA preconditioner.
//!
//! The material cross sections here use a piecewise constant approximation irregardless of the type of
//! approximation used elsewhere.
//!
//! The computed diffusion matrix is stored in PETScPeierlsSolver::diff_mat.
//!
//! \attention  This routine uses the indexing function associated with the first ScalarFlux object stored
//!             at PETScPeierlsSolver::phi. This means that the ScalarFlux objects must be allocated _before_
//!             calling this routine in order to obtain the correct diffusion matrix.
//!
//! \attention  The indexing function used here is Abstract::DensityFunction::IndexNoGC. This is  sufficient
//!             for the 1D RKDG implementation. Implementations for STDG and multi-dimensional cases will
//!             require defining an appropriate IndexNoGC function in the corresponding concrete class.
//!
//! \param[in]  sigma_t     Total cross section \f$ \sigma_{\mathrm{t}} \f$.
//! \param[in]  sigma_s     Scattering cross section \f$ \sigma_{\mathrm{s}} \f$.
//! \param[in]  dt          Timestep size. For steady-state problems, use \pp{dt} = \c inf.
//------------------------------------------------------------------------------------------------------------
template<>
void PETScPeierlsSolver<RKDG::OrdinateFlux>::DSA_ComputeMatrix (

    const CrossSection & sigma_t,
    const CrossSection & sigma_s,
    const double dt
) {

    Global::TMR_DSA_Assemble.Start();

    MatZeroEntries( this->diff_mat );

    ScalarFlux & phi_ref = *this->phi.at(0);

    PetscReal coeff = 0.0, val = 0.0;
    PetscInt m_row = 0, m_col = 0;

    // Set values for the matrix.
    for ( PetscInt i = 1; i <= phi_ref.nx(0);     ++i ) {
    for ( PetscInt d = 0; d <= phi_ref.DG_degree; ++d ) {

        m_row = phi_ref.IndexNoGC(i,d);

        // --- Cell integral terms. ------------------------------------------------------------------- //

        m_col = phi_ref.IndexNoGC(i,d);
        val = this->dx(0) * ( sigma_t(i,0) - sigma_s(i,0) + 1.0 / dt ) / (2*d + 1);
        MatSetValues( this->diff_mat, 1, &m_row, 1, &m_col, &val, ADD_VALUES );

        coeff = 2.0 / ( 3.0 * this->dx(0) * ( sigma_t(i,0) + 1.0 / dt ) );

        for ( PetscInt a = 0; a <= phi_ref.DG_degree; ++a ) {

            m_col = phi_ref.IndexNoGC(i,a);
            val = coeff * B(d,a);
            MatSetValues( this->diff_mat, 1, &m_row, 1, &m_col, &val, ADD_VALUES );
        }

        // --- Cell interface (inter-cell) flux terms. ------------------------------------------------ //

        // Left boundary cell.
        if ( i == 1 ) {

            // Penalty flux.
            for ( PetscInt a = 0; a <= phi_ref.DG_degree; ++a ) {

                m_col = phi_ref.IndexNoGC(i+1,a);
                val = neg1pow(a+1) * Kappa( i, i+1, sigma_t, dt );
                MatSetValues( this->diff_mat, 1, &m_row, 1, &m_col, &val, ADD_VALUES );

                m_col = phi_ref.IndexNoGC(i,a);
                val = Kappa( i, i, sigma_t, dt );
                MatSetValues( this->diff_mat, 1, &m_row, 1, &m_col, &val, ADD_VALUES );
            }

            // Test normal derivative flux.
            coeff = d * (d + 1) / ( 6.0 * this->dx(0) );

            for ( PetscInt a = 0; a <= phi_ref.DG_degree; ++a ) {

                m_col = phi_ref.IndexNoGC(i,a);
                val = - coeff / ( sigma_t(i,0) + 1.0 / dt );
                MatSetValues( this->diff_mat, 1, &m_row, 1, &m_col, &val, ADD_VALUES );

                m_col = phi_ref.IndexNoGC(i+1,a);
                val = 2.0 * coeff * neg1pow(a) / ( sigma_t(i,0) + sigma_t(i+1,0) + 2.0 / dt );
                MatSetValues( this->diff_mat, 1, &m_row, 1, &m_col, &val, ADD_VALUES );
            }

            // Unknown normal derivative flux.
            coeff = -1.0 / ( 6.0 * this->dx(0) );

            for ( PetscInt a = 0; a <= phi_ref.DG_degree; ++a ) {

                m_col = phi_ref.IndexNoGC(i+1,a);
                val = 2.0 * coeff * neg1pow(a) * a * (a + 1) / ( sigma_t(i,0) + sigma_t(i+1,0) + 2.0 / dt );
                MatSetValues( this->diff_mat, 1, &m_row, 1, &m_col, &val, ADD_VALUES );

                m_col = phi_ref.IndexNoGC(i,a);
                val = coeff * a * (a + 1) / ( sigma_t(i,0) + 1.0 / dt );
                MatSetValues( this->diff_mat, 1, &m_row, 1, &m_col, &val, ADD_VALUES );
            }

        // Right boundary cell
        } else if ( i == phi_ref.nx(0) ) {

            // Penalty flux.
            for ( PetscInt a = 0; a <= phi_ref.DG_degree; ++a ) {

                m_col = phi_ref.IndexNoGC(i-1,a);
                val = neg1pow(d+1) * Kappa( i-1, i, sigma_t, dt );
                MatSetValues( this->diff_mat, 1, &m_row, 1, &m_col, &val, ADD_VALUES );

                m_col = phi_ref.IndexNoGC(i,a);
                val = neg1pow(d+a) * Kappa( i, i, sigma_t, dt );
                MatSetValues( this->diff_mat, 1, &m_row, 1, &m_col, &val, ADD_VALUES );
            }

            // Test normal derivative flux.
            coeff = neg1pow(d) * d * (d + 1) / ( 6.0 * this->dx(0) );

            for ( PetscInt a = 0; a <= phi_ref.DG_degree; ++a ) {

                m_col = phi_ref.IndexNoGC(i-1,a);
                val = - 2.0 * coeff / ( sigma_t(i-1,0) + sigma_t(i,0) + 2.0 / dt );
                MatSetValues( this->diff_mat, 1, &m_row, 1, &m_col, &val, ADD_VALUES );

                m_col = phi_ref.IndexNoGC(i,a);
                val = coeff * neg1pow(a) / ( sigma_t(i,0) + 1.0 / dt );
                MatSetValues( this->diff_mat, 1, &m_row, 1, &m_col, &val, ADD_VALUES );
            }

            // Unknown normal derivative flux.
            coeff = 1.0 / ( 6.0 * this->dx(0) );

            for ( PetscInt a = 0; a <= phi_ref.DG_degree; ++a ) {

                m_col = phi_ref.IndexNoGC(i,a);
                val = coeff * neg1pow(d+a) * a * (a + 1) / ( sigma_t(i,0) + 1.0 / dt );
                MatSetValues( this->diff_mat, 1, &m_row, 1, &m_col, &val, ADD_VALUES );

                m_col = phi_ref.IndexNoGC(i-1,a);
                val = 2.0 * coeff * neg1pow(d) * a * (a + 1) / ( sigma_t(i-1,0) + sigma_t(i,0) + 2.0 / dt );
                MatSetValues( this->diff_mat, 1, &m_row, 1, &m_col, &val, ADD_VALUES );
            }

        // Interior cells.
        } else {

            // Penalty flux.
            for ( PetscInt a = 0; a <= phi_ref.DG_degree; ++a ) {

                m_col = phi_ref.IndexNoGC(i-1,a);
                val = Kappa( i-1, i, sigma_t, dt ) * neg1pow(d+1);
                MatSetValues( this->diff_mat, 1, &m_row, 1, &m_col, &val, ADD_VALUES );

                m_col = phi_ref.IndexNoGC(i,a);
                val = Kappa( i, i, sigma_t, dt ) * ( 1 + neg1pow(d+a) );
                MatSetValues( this->diff_mat, 1, &m_row, 1, &m_col, &val, ADD_VALUES );

                m_col = phi_ref.IndexNoGC(i+1,a);
                val = Kappa( i, i+1, sigma_t, dt ) * neg1pow(a+1);
                MatSetValues( this->diff_mat, 1, &m_row, 1, &m_col, &val, ADD_VALUES );
            }

            // Test normal derivative flux.
            coeff = d * (d + 1) / ( 6.0 * this->dx(0) );

            for ( PetscInt a = 0; a <= phi_ref.DG_degree; ++a ) {

                m_col = phi_ref.IndexNoGC(i-1,a);
                val = 2.0 * coeff * neg1pow(d+1) / ( sigma_t(i-1,0) + sigma_t(i,0) + 2.0 / dt );
                MatSetValues( this->diff_mat, 1, &m_row, 1, &m_col, &val, ADD_VALUES );

                m_col = phi_ref.IndexNoGC(i,a);
                val = coeff * ( neg1pow(d+a) - 1 ) / ( sigma_t(i,0) + 1.0 / dt );
                MatSetValues( this->diff_mat, 1, &m_row, 1, &m_col, &val, ADD_VALUES );

                m_col = phi_ref.IndexNoGC(i+1,a);
                val = 2.0 * coeff * neg1pow(a) / ( sigma_t(i,0) + sigma_t(i+1,0) + 2.0 / dt );
                MatSetValues( this->diff_mat, 1, &m_row, 1, &m_col, &val, ADD_VALUES );
            }

            // Unknown normal derivative flux.
            coeff = 1.0 / ( 6.0 * this->dx(0) );

            for ( PetscInt a = 0; a <= phi_ref.DG_degree; ++a ) {

                m_col = phi_ref.IndexNoGC(i-1,a);
                val = 2.0 * coeff * neg1pow(d) * a * (a + 1) / ( sigma_t(i-1,0) + sigma_t(i,0) + 2.0 / dt );
                MatSetValues( this->diff_mat, 1, &m_row, 1, &m_col, &val, ADD_VALUES );

                m_col = phi_ref.IndexNoGC(i,a);
                val = coeff * ( neg1pow(d+a) - 1 ) * a * (a + 1) / ( sigma_t(i,0) + 1.0 / dt );
                MatSetValues( this->diff_mat, 1, &m_row, 1, &m_col, &val, ADD_VALUES );

                m_col = phi_ref.IndexNoGC(i+1,a);
                val = 2.0 * coeff * neg1pow(a+1) * a * (a + 1) / ( sigma_t(i,0) +  sigma_t(i+1,0) + 2.0 / dt );
                MatSetValues( this->diff_mat, 1, &m_row, 1, &m_col, &val, ADD_VALUES );
            }
        }

        // --- Boundary flux terms. ------------------------------------------------------------------- //

        // Left boundary.
        if ( i == 1 ) {

            // Penalty flux.
            for ( PetscInt a = 0; a <= phi_ref.DG_degree; ++a ) {

                m_col = phi_ref.IndexNoGC(i,a);
                val = Kappa( i, i, sigma_t, dt ) * neg1pow(d+a);
                MatSetValues( this->diff_mat, 1, &m_row, 1, &m_col, &val, ADD_VALUES );
            }

            // Test normal derivative flux.
            coeff = neg1pow(d+1) * d * (d + 1) / ( 6.0 * this->dx(0) * ( sigma_t(i,0) + 1.0 / dt ) );

            for ( PetscInt a = 0; a <= phi_ref.DG_degree; ++a ) {

                m_col = phi_ref.IndexNoGC(i,a);
                val = coeff * neg1pow(a);
                MatSetValues( this->diff_mat, 1, &m_row, 1, &m_col, &val, ADD_VALUES );
            }

            // Unknown normal derivative flux.
            coeff = neg1pow(d+1) / ( 6.0 * this->dx(0) * ( sigma_t(i,0) + 1.0 / dt ) );

            for ( PetscInt a = 0; a <= phi_ref.DG_degree; ++a ) {

                m_col = phi_ref.IndexNoGC(i,a);
                val = coeff * neg1pow(a) * a * (a + 1);
                MatSetValues( this->diff_mat, 1, &m_row, 1, &m_col, &val, ADD_VALUES );
            }

        // Right boundary
        } else if ( i == phi_ref.nx(0) ) {

            // Penalty flux.
            for ( PetscInt a = 0; a <= phi_ref.DG_degree; ++a ) {

                m_col = phi_ref.IndexNoGC(i,a);
                val = Kappa( i, i, sigma_t, dt );
                MatSetValues( this->diff_mat, 1, &m_row, 1, &m_col, &val, ADD_VALUES );
            }

            // Test normal derivative flux.
            coeff = -d * (d + 1) / ( 6.0 * this->dx(0) * ( sigma_t(i,0) + 1.0 / dt ) );

            for ( PetscInt a = 0; a <= phi_ref.DG_degree; ++a ) {

                m_col = phi_ref.IndexNoGC(i,a);
                val = coeff;
                MatSetValues( this->diff_mat, 1, &m_row, 1, &m_col, &val, ADD_VALUES );
            }

            // Unknown normal derivative flux.
            coeff = -1.0 / ( 6.0 * this->dx(0) * ( sigma_t(i,0) + 1.0 / dt ) );

            for ( PetscInt a = 0; a <= phi_ref.DG_degree; ++a ) {

                m_col = phi_ref.IndexNoGC(i,a);
                val = coeff * a * (a + 1);
                MatSetValues( this->diff_mat, 1, &m_row, 1, &m_col, &val, ADD_VALUES );
            }
        }
    }}

    MatAssemblyBegin( this->diff_mat, MAT_FINAL_ASSEMBLY );
    MatAssemblyEnd( this->diff_mat, MAT_FINAL_ASSEMBLY );

    Global::TMR_DSA_Assemble.Stop();
}


//------------------------------------------------------------------------------------------------------------
//! \brief  Computes the right side of the DSA diffusion system.
//!
//! \param[in]  phi_in  Scalar flux distribution to invert diffusion operator onto.
//! \param[in]  Q       Temporary object used as workspace to construct right hand vector.
//------------------------------------------------------------------------------------------------------------
template<>
void PETScPeierlsSolver<RKDG::OrdinateFlux>::DSA_ComputeRHS (

    const ScalarFlux & phi_in,
    ScalarFlux & Q
) {

    Q.ZeroDensity();

    Global::TMR_DSA_Solve.Start();

    # pragma omp parallel for
    for ( int64_t i = 0; i <= Q.nx(0) + 1; ++i ) {
    for ( int64_t d = 0; d <= Q.DG_degree; ++d ) {

        Q(i,d) = Q.dx(0) * phi_in(i,d) * (*sigma_s_save)(i,0) / (2*d + 1);
    }}

    Global::TMR_DSA_Solve.Stop();

    Q.PackExternalVector( this->diff_vec_rhs, 1.0, 0.0 );
}


# endif // if SPACE_DIMS == 1
