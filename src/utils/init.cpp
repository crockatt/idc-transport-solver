//------------------------------------------------------------------------------------------------------------
//! \file   utils/init_1d.cpp
//! \brief  Implements routines for initializing 1D test problems.
//!
//! \author Michael M. Crockatt
//! \date   June 2018
//------------------------------------------------------------------------------------------------------------

# include <algorithm>
# include <cinttypes>
# include <cmath>
# include <cstdint>
# include <cstdio>
# include <cstdlib>
# include <cstring>
# include <limits>
# include <memory>
# include <sys/types.h>
# include <unistd.h>

# include "objects/DomainDecomposition.hpp"
# include "objects/RKDG/CrossSection.hpp"
# include "objects/RKDG/OrdinateFlux.hpp"
# include "objects/RKDG/ScalarFlux.hpp"
# include "time_integrators/Abstract/TimeIntegrator.hpp"
# include "utils/CLog.hpp"
# include "utils/Factory/Factory.hpp"
# include "utils/global.hpp"
# include "utils/init.hpp"
# include "utils/Problems/Problem.hpp"
# include "utils/Quadrule/Quadrule.hpp"


using namespace RKDG;
using namespace Quadrule;


//------------------------------------------------------------------------------------------------------------
//! \brief  Computes the \f$ L^2 \f$ projection of the specified problem component into the given object.
//!
//! \param[in]  prob            Problem object containing parameters to compute projection for.
//! \param[in]  dest            Object to store projected coefficients in.
//! \param[in]  eval_fcn        Pointer to function specifying which problem parameter to compute projection
//!                             for.
//! \param[in]  apply_limiter   (optional) <br>
//!                             Specifies whether a positivity-preserving limiter should be applied. Defaults
//!                             to true.
//! \param[in]  proj_order      Number of Gauss quadrature nodes to compute projections with.
//------------------------------------------------------------------------------------------------------------
static void ComputeProjection (

    const Problem & prob,
    RKDG::DensityFunction & dest,
    const ProblemEvalFcn eval_fcn,
    const bool apply_limiter = true,
    const int64_t proj_order = 14
) {

    const double TOL = 1.0e-14;
    const int64_t sample_pts = 20;

    // Zero destination object.
    dest.ZeroDensity();

    // Generate Gaussian quadrature.
    double * const GL_nodes   = new double[ proj_order ];
    double * const GL_weights = new double[ proj_order ];

    ComputeQuadrature( proj_order, GL_nodes, GL_weights, NodesType::GaussLegendre );


# if SPACE_DIMS == 1

    # pragma omp parallel for
    for ( int64_t i = 1; i <= dest.nx(0); ++i ) {

        // Compute coefficients of projection.
        for ( int64_t l = 0; l < proj_order; ++l ) {

            // Scale quadrature node to current spatial cell.
            const double x_l = dest.ax(0) + (i-1) * dest.dx(0) + (GL_nodes[l] + 1.0) * dest.dx(0) / 2.0;

            // Evaluate function at quadrature node.
            const double f_val = (prob.*eval_fcn)( x_l );

            // Accumulate result for each coefficient.
            for ( int64_t d = 0; d <= dest.DG_degree; ++d )
                dest(i,d) += GL_weights[l] * Legendre( d, GL_nodes[l] ) * (2*d + 1) * f_val / 2.0;
        }

        // Do not apply limiter unless flag is set.
        if ( !apply_limiter ) {  continue;  }

        // If cell average is negative, don't bother limiting.
        if ( dest(i,0) < 0.0 ) {

            PRINT_WARNING( "Negative total density in cell %" PRId64 "\n", i )
            continue;
        }

        // Only apply limiter if negative values are detected.
        bool is_negative = false;

        // Sample projected solution and check for negative values.
        for ( int64_t l = 0; l <= sample_pts; ++l ) {

            const double x = -1.0 + (2.0 / sample_pts) * l;

            double pt_val = 0.0;

            for ( int64_t d = 0; d <= dest.DG_degree; ++d )
                pt_val += dest(i,d) * Legendre(d,x);

            is_negative |= ( pt_val < 0.0 );
        }

        if ( !is_negative ) {  continue;  }

    # if LOGLEVEL >= 2
        PRINT_WARNING( "Limiting necessary in cell %" PRId64 "\n", i )
    # endif

        // If limiting is necessary, find suitable limiting parameter. Optimal limiting parameter is computed
        // using a bisection algorithm.
        double theta_min = 0.0,
               theta_max = 1.0;

        while ( (theta_max - theta_min) > TOL ) {

            // Determine if midpoint yields positive solution.
            const double theta_mid = 0.5 * (theta_max + theta_min);

            is_negative = false;

            for ( int64_t l = 0; l <= sample_pts; ++l ) {

                const double x = -1.0 + (2.0 / sample_pts) * l;

                double pt_val = 0.0;

                for ( int64_t d = 0; d <= dest.DG_degree; ++d )
                    pt_val += ( d == 0 ? 1.0 : theta_mid ) * dest(i,d) * Legendre(d,x);

                is_negative |= ( pt_val < 0.0 );
            }

            // Bisect.
            if ( is_negative ) {  theta_max = theta_mid;  }
            else               {  theta_min = theta_mid;  }
        }

        // Apply limiter to high-order moments.
        for ( int64_t d = 1; d <= dest.DG_degree; ++d )
            dest(i,d) *= theta_min;
    }

# elif SPACE_DIMS == 2

    const int64_t (& nx) [SPACE_DIMS] = dest.nx();

    # pragma omp parallel for schedule(dynamic,1) collapse(2)
    for ( int64_t i = 1; i <= nx[0]; ++i ) {
    for ( int64_t j = 1; j <= nx[1]; ++j ) {

//         PRINT_ERROR( "Projecting cell (%" PRId64 ",%" PRId64 ")\n", i,j )

        // Compute coefficients of projection.
        for ( int64_t l = 0; l < proj_order; ++l ) {
        for ( int64_t m = 0; m < proj_order; ++m ) {

            // Scale quadrature node to current spatial cell.
            const double x_l = dest.ax(0) + (i-1) * dest.dx(0) + (GL_nodes[l] + 1.0) * dest.dx(0) / 2.0;
            const double y_m = dest.ax(1) + (j-1) * dest.dx(1) + (GL_nodes[m] + 1.0) * dest.dx(1) / 2.0;

            // Evaluate function at quadrature node.
            const double f_val = (prob.*eval_fcn)( x_l, y_m );

            // Accumulate result for each coefficient.
            for ( int64_t d = 0; d <= dest.DG_degree; ++d ) {
            for ( int64_t e = 0; e <= dest.DG_degree; ++e ) {

                dest(i,j,d,e) += GL_weights[l] * GL_weights[m]
                                 * Legendre( d, GL_nodes[l] ) * Legendre( e, GL_nodes[m] )
                                 * (2*d + 1) * (2*e + 1) / 4.0
                                 * f_val;
            }}
        }}

        // Do not apply limiter unless flag is set.
        if ( !apply_limiter ) {  continue;  }

        // If cell average is negative, don't bother limiting.
        if ( dest(i,j,0,0) < 0.0 ) {

            PRINT_WARNING( "Negative total density in cell %" PRId64 "\n", i )
            continue;
        }

        // Only apply limiter if negative values are detected.
        bool is_negative = false;

        // Sample projected solution and check for negative values.
        for ( int64_t l = 0; l <= sample_pts; ++l ) {
        for ( int64_t m = 0; m <= sample_pts; ++m ) {

            const double x = -1.0 + (2.0 / sample_pts) * l;
            const double y = -1.0 + (2.0 / sample_pts) * m;

            double pt_val = 0.0;

            for ( int64_t d = 0; d <= dest.DG_degree; ++d ) {
            for ( int64_t e = 0; e <= dest.DG_degree; ++e ) {

                pt_val += dest(i,j,d,e) * Legendre(d,x) * Legendre(e,y);
            }}

            is_negative |= ( pt_val < 0.0 );
        }}

        if ( !is_negative ) {  continue;  }

    # if LOGLEVEL >= 2
        PRINT_WARNING( "Limiting necessary in cell (%" PRId64 ",%" PRId64 ")\n", i,j )
    # endif

        // If limiting is necessary, find suitable limiting parameter. Optimal limiting parameter is computed
        // using a bisection algorithm.
        double theta_min = 0.0,
               theta_max = 1.0;

        while ( (theta_max - theta_min) > TOL ) {

            // Determine if midpoint yields positive solution.
            const double theta_mid = 0.5 * (theta_max + theta_min);

            is_negative = false;

            for ( int64_t l = 0; l <= sample_pts; ++l ) {
            for ( int64_t m = 0; m <= sample_pts; ++m ) {

                const double x = -1.0 + (2.0 / sample_pts) * l;
                const double y = -1.0 + (2.0 / sample_pts) * m;

                double pt_val = 0.0;

                for ( int64_t d = 0; d <= dest.DG_degree; ++d ) {
                for ( int64_t e = 0; e <= dest.DG_degree; ++e ) {

                    pt_val += ( d == 0 && e == 0 ? 1.0 : theta_mid )
                              * dest(i,j,d,e) * Legendre(d,x) * Legendre(e,y);
                }}

                is_negative |= ( pt_val < 0.0 );
            }}

            // Bisect.
            if ( is_negative ) {  theta_max = theta_mid;  }
            else               {  theta_min = theta_mid;  }
        }

        // Apply limiter to high-order moments.
        for ( int64_t d = 0; d <= dest.DG_degree; ++d ) {
        for ( int64_t e = 0; e <= dest.DG_degree; ++e ) {

            if ( d == 0 && e == 0 ) {  continue;  }

            dest(i,j,d,e) *= theta_min;
        }}
    }}

# elif SPACE_DIMS == 3

    const int64_t (& nx) [SPACE_DIMS] = dest.nx();

    # pragma omp parallel for schedule(dynamic,1) collapse(3)
    for ( int64_t i = 1; i <= nx[0]; ++i ) {
    for ( int64_t j = 1; j <= nx[1]; ++j ) {
    for ( int64_t k = 1; k <= nx[2]; ++k ) {

//         PRINT_ERROR( "Projecting cell (%" PRId64 ",%" PRId64 ",%" PRId64 ")\n", i,j,k )

        // Compute coefficients of projection.
        for ( int64_t l = 0; l < proj_order; ++l ) {
        for ( int64_t m = 0; m < proj_order; ++m ) {
        for ( int64_t n = 0; n < proj_order; ++n ) {

            // Scale quadrature node to current spatial cell.
            const double x_l = dest.ax(0) + (i-1) * dest.dx(0) + (GL_nodes[l] + 1.0) * dest.dx(0) / 2.0;
            const double y_m = dest.ax(1) + (j-1) * dest.dx(1) + (GL_nodes[m] + 1.0) * dest.dx(1) / 2.0;
            const double z_n = dest.ax(2) + (k-1) * dest.dx(2) + (GL_nodes[n] + 1.0) * dest.dx(2) / 2.0;

            // Evaluate function at quadrature node.
            const double f_val = (prob.*eval_fcn)( x_l, y_m, z_n );

            // Accumulate result for each coefficient.
            for ( int64_t d = 0; d <= dest.DG_degree; ++d ) {
            for ( int64_t e = 0; e <= dest.DG_degree; ++e ) {
            for ( int64_t f = 0; f <= dest.DG_degree; ++f ) {

                dest(i,j,d,e) +=   GL_weights[l] * Legendre( d, GL_nodes[l] ) * (2*d + 1)
                                 * GL_weights[m] * Legendre( e, GL_nodes[m] ) * (2*e + 1)
                                 * GL_weights[n] * Legendre( f, GL_nodes[n] ) * (2*f + 1)
                                 * f_val / 8.0;
            }}}
        }}}

        // Do not apply limiter unless flag is set.
        if ( !apply_limiter ) {  continue;  }

        // If cell average is negative, don't bother limiting.
        if ( dest(i,j,k,0,0,0) < 0.0 ) {

            PRINT_WARNING( "Negative total density in cell %" PRId64 "\n", i )
            continue;
        }

        // Only apply limiter if negative values are detected.
        bool is_negative = false;

        // Sample projected solution and check for negative values.
        for ( int64_t l = 0; l <= sample_pts; ++l ) {
        for ( int64_t m = 0; m <= sample_pts; ++m ) {
        for ( int64_t n = 0; n <= sample_pts; ++n ) {

            const double x = -1.0 + (2.0 / sample_pts) * l;
            const double y = -1.0 + (2.0 / sample_pts) * m;
            const double z = -1.0 + (2.0 / sample_pts) * n;

            double pt_val = 0.0;

            for ( int64_t d = 0; d <= dest.DG_degree; ++d ) {
            for ( int64_t e = 0; e <= dest.DG_degree; ++e ) {
            for ( int64_t f = 0; f <= dest.DG_degree; ++f ) {

                pt_val += dest(i,j,k,d,e,f) * Legendre(d,x) * Legendre(e,y) * Legendre(f,z);
            }}}

            is_negative |= ( pt_val < 0.0 );
        }}}

        if ( !is_negative ) {  continue;  }

    # if LOGLEVEL >= 2
        PRINT_WARNING( "Limiting necessary in cell (%" PRId64 ",%" PRId64 ",%" PRId64 ")\n", i,j,k )
    # endif

        // If limiting is necessary, find suitable limiting parameter. Optimal limiting parameter is computed
        // using a bisection algorithm.
        double theta_min = 0.0,
               theta_max = 1.0;

        while ( (theta_max - theta_min) > TOL ) {

            // Determine if midpoint yields positive solution.
            const double theta_mid = 0.5 * (theta_max + theta_min);

            bool is_negative = false;

            for ( int64_t l = 0; l <= sample_pts; ++l ) {
            for ( int64_t m = 0; m <= sample_pts; ++m ) {
            for ( int64_t n = 0; n <= sample_pts; ++n ) {

                const double x = -1.0 + (2.0 / sample_pts) * l;
                const double y = -1.0 + (2.0 / sample_pts) * m;
                const double z = -1.0 + (2.0 / sample_pts) * n;

                double pt_val = 0.0;

                for ( int64_t d = 0; d <= dest.DG_degree; ++d ) {
                for ( int64_t e = 0; e <= dest.DG_degree; ++e ) {
                for ( int64_t f = 0; f <= dest.DG_degree; ++f ) {

                    pt_val += ( d == 0 && e == 0 && f == 0 ? 1.0 : theta_mid )
                              * dest(i,j,k,d,e,f) * Legendre(d,x) * Legendre(e,y) * Legendre(f,z);
                }}}

                is_negative |= ( pt_val < 0.0 );
            }}}

            // Bisect.
            if ( is_negative ) {  theta_max = theta_mid;  }
            else               {  theta_min = theta_mid;  }
        }

        // Apply limiter to high-order moments.
        for ( int64_t d = 0; d <= dest.DG_degree; ++d ) {
        for ( int64_t e = 0; e <= dest.DG_degree; ++e ) {
        for ( int64_t f = 0; f <= dest.DG_degree; ++f ) {

            if ( d == 0 && e == 0 && f == 0 ) {  continue;  }

            dest(i,j,k,d,e,f) *= theta_min;
        }}}
    }}}

# endif // if SPACE_DIMS == ?


    // Cleanup.
    delete [] GL_nodes;
    delete [] GL_weights;
}


//============================================================================================================
//=== MAIN INITIALIZATION ROUTINE ============================================================================
//============================================================================================================


//------------------------------------------------------------------------------------------------------------
//! \brief  Performs initialization of material cross sections, initial conditions, and sources.
//!
//! Creates required objects and sets initial conditions. Prints a problem summary to the global logging
//! interface.
//!
//! \param[out]     psi             Upon return contains the initial condition for the angular flux.
//! \param[out]     source          Upon return contains the angular flux corresponding to the source term.
//! \param[out]     sigma_t         Upon return contains the total cross section \f$ \sigma_{\mathrm{t}} \f$.
//! \param[out]     sigma_s         Upon return contains the scattering cross section \f$ \sigma_{\mathrm{s}} \f$.
//! \param[out]     spatial_params  Upon return contains the parameters of the spatial mesh.
//------------------------------------------------------------------------------------------------------------
void Init (

    OrdinateFlux & psi,
    OrdinateFlux & source,
    CrossSection & sigma_t,
    CrossSection & sigma_s,
    DomainDecomposition & spatial_params
) {

    // Determine filtered input list to use for initialization.
    ParameterList param_list = Global::input_list;

    try {

        if ( Global::input_list.GetValue<std::string>( "hybrid_method" ) != "none" )
            param_list = Abstract::TimeIntegrator::MakeUncollidedList( Global::input_list );

    } catch (...) {/* empty */}

    // Initialize problem.
    auto prob = std::shared_ptr<Problem>( Factory<Problem>::CreateObject( param_list ) );
    prob->Print();
    prob->OverrideOptions( Global::input_list );
    prob->OverrideOptions( param_list );

    // Setup spatial mesh.
    {
        try {  Global::periodic = param_list.GetValue<bool>( "periodic" );  } catch (...) {/* empty */}

        /*  Determine reflecting boundary conditions.
         *
         *  Note: 1D implementation does not support reflecting boundary conditions.
         */
    # if SPACE_DIMS >= 2

        try {
            if ( param_list.GetValue<bool>( "reflect_x_min" ) )
                Global::reflecting_boundaries |= BoundaryEdge::X_Min;

        } catch (...) {/* empty */}

        try {
            if ( param_list.GetValue<bool>( "reflect_x_max" ) )
                Global::reflecting_boundaries |= BoundaryEdge::X_Max;

        } catch (...) {/* empty */}

        try {
            if ( param_list.GetValue<bool>( "reflect_y_min" ) )
                Global::reflecting_boundaries |= BoundaryEdge::Y_Min;

        } catch (...) {/* empty */}

        try {
            if ( param_list.GetValue<bool>( "reflect_y_max" ) )
                Global::reflecting_boundaries |= BoundaryEdge::Y_Max;

        } catch (...) {/* empty */}

    # endif // if SPACE_DIMS >= 2

    # if SPACE_DIMS == 3

        try {
            if ( param_list.GetValue<bool>( "reflect_z_min" ) )
                Global::reflecting_boundaries |= BoundaryEdge::Z_Min;

        } catch (...) {/* empty */}

        try {
            if ( param_list.GetValue<bool>( "reflect_z_max" ) )
                Global::reflecting_boundaries |= BoundaryEdge::Z_Max;

        } catch (...) {/* empty */}

    # endif // if SPACE_DIMS == 3

    # if SPACE_DIMS >= 2

        // Implementation is not designed to handle combinations of periodic and reflecting boundary
        // conditions.
        if (    Global::periodic
             && BitmaskHasAny( Global::reflecting_boundaries, BoundaryEdge::All )
        ) {
            std::string error_message = "Combined reflecting and periodic boundary conditions not "
                                        "supported.\n";

            PRINT_ERROR( error_message.c_str() )
            throw std::runtime_error( error_message );
        }

    # endif // if SPACE_DIMS >= 2

        // Read in MPI domain decomposition parameters and construct Cartesian communicator.
    # if defined (ENABLE_MPI)

        int MPI_num_blocks [SPACE_DIMS] =
            # if SPACE_DIMS == 1
                {1}
            # elif SPACE_DIMS == 2
                {1,1}
            # elif SPACE_DIMS == 3
                {1,1,1}
            # endif
                ;

        int MPI_block_coords [SPACE_DIMS] =
            # if SPACE_DIMS == 1
                {0}
            # elif SPACE_DIMS == 2
                {0,0}
            # elif SPACE_DIMS == 3
                {0,0,0}
            # endif
                ;

        int periods [SPACE_DIMS] =
            # if SPACE_DIMS == 1
                {0}
            # elif SPACE_DIMS == 2
                {0,0}
            # elif SPACE_DIMS == 3
                {0,0,0}
            # endif
                ;

        param_list.GetValue( "num_mpi_blocks_x", MPI_num_blocks[0] );
    # if SPACE_DIMS >= 2
        param_list.GetValue( "num_mpi_blocks_y", MPI_num_blocks[1] );
    # endif
    # if SPACE_DIMS == 3
        param_list.GetValue( "num_mpi_blocks_z", MPI_num_blocks[2] );
    # endif

        char err_str [MPI_MAX_ERROR_STRING];
        int err_len;

        if ( Global::periodic ) {

            for ( int64_t i = 0; i < SPACE_DIMS; ++i )
                periods[i] = 1;
        }

        int MPI_err = MPI_Cart_create( MPI_COMM_WORLD, SPACE_DIMS, MPI_num_blocks, periods, 1,
                                       &Global::MPI_cart_comm );

        if ( MPI_err ) {

            MPI_Error_string( MPI_err, err_str, &err_len );
            std::string error_string = std::string( err_str );
            std::replace( error_string.begin(), error_string.end(), '\n', ' ' );

            PRINT_ERROR( "MPI_Cart_create returned error '%s'.\nAborting run.\n", error_string.c_str() )
            MPI_Abort( MPI_COMM_WORLD, EXIT_FAILURE );
        }

        MPI_Cart_coords( Global::MPI_cart_comm, Global::MPI_rank, SPACE_DIMS, MPI_block_coords );

    # endif // if defined (ENABLE_MPI)

        spatial_params = DomainDecomposition( Global::input_list );
    }

    // Setup additional parameters.
    bool pwc_cross_sections = false,
         symmetric_reduce   = false,
         output_init        = false;

    int64_t ang_order;
    OrdinateSet::OrdinateType ordinate_type;

# if SPACE_DIMS == 2

    try         {  param_list.GetValue( "ordinate_sym_reduce", symmetric_reduce );  }
    catch (...) {  symmetric_reduce = false;                                        }

# endif // if SPACE_DIMS == 2

    try         {  param_list.GetValue( "piecewise_cross", pwc_cross_sections );  }
    catch (...) {  pwc_cross_sections = false;                                    }

    try         {  param_list.GetValue( "output_init", output_init );  }
    catch (...) {  output_init = false;                                }

    param_list.GetValue( "ang_order", ang_order );

    try {
        param_list.GetValue( "ordinate_type", ordinate_type, OrdinateSet::String_to_OrdinateType );
    } catch (...) {
        ordinate_type = OrdinateSet::OrdinateType::
            # if SPACE_DIMS == 1
                GaussLegendre
            # else
                ChebyshevLegendre
            # endif
                ;
    }

    // Output parameters to logging interface.
# if defined (ENABLE_MPI)
    if ( Global::MPI_rank == 0 )
# endif
    {
        PRINT_LOG( "\n" )

        PRINT_LOG( "  %-*s  [ % .2e, % .2e ]\n", Global::col_width, "x Spatial Domain:",
                   spatial_params.global_ax(0), spatial_params.global_bx(0) )

    # if SPACE_DIMS >= 2

        PRINT_LOG( "  %-*s  [ % .2e, % .2e ]\n", Global::col_width, "y Spatial Domain:",
                   spatial_params.global_ax(1), spatial_params.global_bx(1) )

    # endif // if SPACE_DIMS >= 2

    # if SPACE_DIMS == 3

        PRINT_LOG( "  %-*s  [ % .2e, % .2e ]\n", Global::col_width, "z Spatial Domain:",
                   spatial_params.global_ax(2), spatial_params.global_bx(2) )

    # endif // if SPACE_DIMS == 3

        PRINT_LOG( "  %-*s  %s\n", Global::col_width, "Periodic:", ( Global::periodic ? "yes" : "no" ) )

    # if SPACE_DIMS >= 2

        PRINT_LOG( "  %-*s  %s\n", Global::col_width, "X_Min Reflect:",
                   ( BitmaskHasAny( Global::reflecting_boundaries, BoundaryEdge::X_Min ) ? "yes" : "no" ) )

        PRINT_LOG( "  %-*s  %s\n", Global::col_width, "X_Max Reflect:",
                   ( BitmaskHasAny( Global::reflecting_boundaries, BoundaryEdge::X_Max ) ? "yes" : "no" ) )

        PRINT_LOG( "  %-*s  %s\n", Global::col_width, "Y_Min Reflect:",
                   ( BitmaskHasAny( Global::reflecting_boundaries, BoundaryEdge::Y_Min ) ? "yes" : "no" ) )

        PRINT_LOG( "  %-*s  %s\n", Global::col_width, "Y_Max Reflect:",
                   ( BitmaskHasAny( Global::reflecting_boundaries, BoundaryEdge::Y_Max ) ? "yes" : "no" ) )

    # endif // if SPACE_DIMS >= 2

    # if SPACE_DIMS == 3

        PRINT_LOG( "  %-*s  %s\n", Global::col_width, "Z_Min Reflect:",
                   ( BitmaskHasAny( Global::reflecting_boundaries, BoundaryEdge::Z_Min ) ? "yes" : "no" ) )

        PRINT_LOG( "  %-*s  %s\n", Global::col_width, "Z_Max Reflect:",
                   ( BitmaskHasAny( Global::reflecting_boundaries, BoundaryEdge::Z_Max ) ? "yes" : "no" ) )

    # endif // if SPACE_DIMS == 3

    # if defined (ENABLE_MPI)

        PRINT_LOG( "  %-*s % " PRId64"\n", Global::col_width,
                   "Num x MPI Blocks:", spatial_params.MPI_num_blocks(0) )

    # if SPACE_DIMS >= 2

        PRINT_LOG( "  %-*s % " PRId64"\n", Global::col_width,
                   "Num y MPI Blocks:", spatial_params.MPI_num_blocks(1) )

    # endif // if SPACE_DIMS >= 2

    # if SPACE_DIMS == 3

        PRINT_LOG( "  %-*s % " PRId64"\n", Global::col_width,
                   "Num z MPI Blocks:", spatial_params.MPI_num_blocks(2) )

    # endif // if SPACE_DIMS == 3
    # endif // if defined (ENABLE_MPI)

        PRINT_LOG( "  %-*s % " PRId64"\n", Global::col_width,
                   "Num x Spatial Cells:", spatial_params.global_nx(0) )

    # if SPACE_DIMS >= 2

        PRINT_LOG( "  %-*s % " PRId64"\n", Global::col_width,
                   "Num y Spatial Cells:", spatial_params.global_nx(1) )

    # endif // if SPACE_DIMS >= 2

    # if SPACE_DIMS == 3

        PRINT_LOG( "  %-*s % " PRId64"\n", Global::col_width,
                   "Num z Spatial Cells:", spatial_params.global_nx(2) )

    # endif // if SPACE_DIMS == 3

        PRINT_LOG( "  %-*s % .4e\n", Global::col_width, "Delta x:", spatial_params.dx(0) )

    # if SPACE_DIMS >= 2

        PRINT_LOG( "  %-*s % .4e\n", Global::col_width, "Delta y:", spatial_params.dx(1) )

    # endif // if SPACE_DIMS >= 2

    # if SPACE_DIMS == 3

        PRINT_LOG( "  %-*s % .4e\n", Global::col_width, "Delta z:", spatial_params.dx(2) )

    # endif // if SPACE_DIMS == 3

        PRINT_LOG( "  %-*s % " PRId64"\n", Global::col_width, "Max DG Degree (x):", Global::DG_degree  )

        PRINT_LOG( "  %-*s  %-20s\n", Global::col_width, "Cross-Sections:",
                   ( pwc_cross_sections ? "piecewise constant" : "high order" ) )

        PRINT_LOG( "\n" )
    }

    // --- Compute problem parameters: cross section, initial condition, and source. ------------------ //

    std::string filename;
    ScalarFlux temp_phi {};
    ComputeTPI();

    char nx_buff[32] = "";
    std::sprintf( nx_buff, "%07" PRId64 , spatial_params.global_nx(0) );

    // Determine degree of spatial approximation of cross sections.
    const int64_t sigma_degree = ( pwc_cross_sections ? 0 : Global::DG_degree );

    // Configure objects with given parameters.
    sigma_t.Reconfigure( spatial_params, sigma_degree );
    sigma_s.Reconfigure( spatial_params, sigma_degree );

    temp_phi.Reconfigure( spatial_params, Global::DG_degree );

    psi.Reconfigure( spatial_params, Global::DG_degree, ang_order, symmetric_reduce, ordinate_type );
    source.Reconfigure( spatial_params, Global::DG_degree, ang_order, symmetric_reduce, ordinate_type );

    // Compute total cross section.
    filename = "sigma_t_dg"
               + std::to_string( sigma_t.DG_degree )
               + "_nx"
               + std::string( nx_buff )
               + ".cs"
               + std::to_string( SPACE_DIMS );

    // Try reading from disk.
    try {

        if ( !output_init ) {  throw std::range_error( "Aborting reading from disk.\n" );  }

        sigma_t.ReadFromDisk( filename );

    // If reading fails, compute projection and (optionally) write result to disk.
    } catch (...) {

        PRINT_NOTE( "Computing total cross section.\n" )

        ComputeProjection( *prob, sigma_t, &Problem::EvalTotalCross );

        if ( output_init ) {  sigma_t.WriteToDisk( filename );  }
    }

    sigma_t.ComputeTensor();

    // Compute scattering cross section.
    filename = "sigma_s_dg"
               + std::to_string( sigma_s.DG_degree )
               + "_nx"
               + std::string( nx_buff )
               + ".cs"
               + std::to_string( SPACE_DIMS );

    // Try reading from disk.
    try {

        if ( !output_init ) {  throw std::range_error( "Aborting reading from disk.\n" );  }

        sigma_s.ReadFromDisk( filename );

    // If reading fails, compute projection and (optionally) write result to disk.
    } catch (...) {

        PRINT_NOTE( "Computing scattering cross section.\n" )

        ComputeProjection( *prob, sigma_s, &Problem::EvalScatterCross );

        if ( output_init ) {  sigma_s.WriteToDisk( filename );  }
    }

    sigma_s.ComputeTensor();

    // Compute isotropic initial condition.
    filename = "initial_dg"
               + std::to_string( temp_phi.DG_degree )
               + "_nx"
               + std::string( nx_buff )
               + ".sd"
               + std::to_string( SPACE_DIMS );

    // Try reading from disk.
    try {

        if ( !output_init ) {  throw std::range_error( "Aborting reading from disk.\n" );  }

        temp_phi.ReadFromDisk( filename );

    // If reading fails, compute projection and (optionally) write result to disk.
    } catch (...) {

        PRINT_NOTE( "Computing initial condition.\n" )

        ComputeProjection( *prob, temp_phi, &Problem::EvalInitialCondition );

        if ( output_init ) {  temp_phi.WriteToDisk( filename );  }
    }

    // Distribute initial condition to all ordinates.
    {   const int64_t (& nx) [SPACE_DIMS] = psi.nx();

    # if SPACE_DIMS == 1

        # pragma omp parallel for
        for ( int64_t i = 1; i <= nx[0];         ++i ) {
        for ( int64_t q = 0; q <  psi.nq();      ++q ) {
        for ( int64_t d = 0; d <= psi.DG_degree; ++d ) {

            psi(q,i,d) = temp_phi(i,d) / 2.0;
        }}}

    # elif SPACE_DIMS == 2

        # pragma omp parallel for
        for ( int64_t i = 1; i <= nx[0];         ++i ) {
        for ( int64_t j = 1; j <= nx[1];         ++j ) {
        for ( int64_t q = 0; q <  psi.nq();      ++q ) {
        for ( int64_t d = 0; d <= psi.DG_degree; ++d ) {
        for ( int64_t e = 0; e <= psi.DG_degree; ++e ) {

            psi(q,i,j,d,e) = temp_phi(i,j,d,e) / (4.0 * M_PI);
        }}}}}

    # elif SPACE_DIMS == 3

        # pragma omp parallel for
        for ( int64_t i = 1; i <= nx[0];         ++i ) {
        for ( int64_t j = 1; j <= nx[1];         ++j ) {
        for ( int64_t k = 1; k <= nx[2];         ++k ) {
        for ( int64_t q = 0; q <  psi.nq();      ++q ) {
        for ( int64_t d = 0; d <= psi.DG_degree; ++d ) {
        for ( int64_t e = 0; e <= psi.DG_degree; ++e ) {
        for ( int64_t f = 0; f <= psi.DG_degree; ++f ) {

            psi(q,i,j,k,d,e,f) = temp_phi(i,j,k,d,e,f) / (4.0 * M_PI);
        }}}}}}}

    # endif // if SPACE_DIMS == ?
    }

    psi.SynchronizeHalos();

    // Compute isotropic source.
    filename = "source_dg"
               + std::to_string( temp_phi.DG_degree )
               + "_nx"
               + std::string( nx_buff )
               + ".sd"
               + std::to_string( SPACE_DIMS );

    // Try reading from disk.
    try {

        if ( !output_init ) {  throw std::range_error( "Aborting reading from disk.\n" );  }

        temp_phi.ReadFromDisk( filename );

    // If reading fails, compute projection and (optionally) write result to disk.
    } catch (...) {

        PRINT_NOTE( "Computing source term.\n" )

        ComputeProjection( *prob, temp_phi, &Problem::EvalSource );

        if ( output_init ) {  temp_phi.WriteToDisk( filename );  }
    }

    prob->SetInflowBoundaries( temp_phi );

    // Distribute initial condition to all ordinates.
    {   const int64_t (& nx) [SPACE_DIMS] = source.nx();

    # if SPACE_DIMS == 1

        # pragma omp parallel for
        for ( int64_t i = 0; i <= nx[0] + 1;        ++i ) {
        for ( int64_t q = 0; q <  source.nq();      ++q ) {
        for ( int64_t d = 0; d <= source.DG_degree; ++d ) {

            source(q,i,d) = temp_phi(i,d) / 2.0;
        }}}

    # elif SPACE_DIMS == 2

        # pragma omp parallel for
        for ( int64_t i = 0; i <= nx[0] + 1;        ++i ) {
        for ( int64_t j = 0; j <= nx[1] + 1;        ++j ) {
        for ( int64_t q = 0; q <  source.nq();      ++q ) {
        for ( int64_t d = 0; d <= source.DG_degree; ++d ) {
        for ( int64_t e = 0; e <= source.DG_degree; ++e ) {

            source(q,i,j,d,e) = temp_phi(i,j,d,e) / (4.0 * M_PI);
        }}}}}

    # elif SPACE_DIMS == 3

        # pragma omp parallel for
        for ( int64_t i = 0; i <= nx[0] + 1;        ++i ) {
        for ( int64_t j = 0; j <= nx[1] + 1;        ++j ) {
        for ( int64_t k = 0; k <= nx[2] + 1;        ++k ) {
        for ( int64_t q = 0; q <  source.nq();      ++q ) {
        for ( int64_t d = 0; d <= source.DG_degree; ++d ) {
        for ( int64_t e = 0; e <= source.DG_degree; ++e ) {
        for ( int64_t f = 0; f <= source.DG_degree; ++f ) {

            source(q,i,j,k,d,e,f) = temp_phi(i,j,k,d,e,f) / (4.0 * M_PI);
        }}}}}}}

    # endif // if SPACE_DIMS == ?
    }

    source.SynchronizeHalos();

}
