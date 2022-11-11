//------------------------------------------------------------------------------------------------------------
//! \file   operators/RKDG/TransportOperator.cpp
//! \brief  Implementation file for RKDG transport operator routines.
//!
//! \author Michael M. Crockatt
//! \date   February 2018
//------------------------------------------------------------------------------------------------------------

# include <algorithm>

# include "operators/RKDG/TransportOperator.hpp"
# include "utils/CLog.hpp"
# include "utils/SIMD.hpp"


namespace TransportOperator {


//============================================================================================================
//=== MULTI-DIMENSIONAL OPERATOR ROUTINES ====================================================================
//============================================================================================================


//------------------------------------------------------------------------------------------------------------
//! \brief  Computes \f$ \vec{y} \gets \beta \vec{y} + \alpha \mathcal{P} \vec{x} \f$ where
//!         \f$ \mathcal{P} \f$ is the operator that integrates over the angular dimension.
//!
//! Computes the generalized matrix-vector product using the operator \f$ \mathcal{P} \f$.
//!
//! \param[in]      alpha       Scalar to augment \f$ \mathcal{P} \f$.
//! \param[in]      beta        Scalar to augment \f$ y \f$.
//! \param[in]      x           Angular flux to which \f$ \mathcal{P} \f$ is applied.
//! \param[in,out]  y           Scalar flux that is scaled and into which the result is stored.
//! \param[in]      op_domain   (optional) <br>
//!                             Domain over which to apply copy operation. Defaults to OpDomain::All.
//------------------------------------------------------------------------------------------------------------
void Pmv (

    const double alpha,
    const double beta,
    const RKDG::OrdinateFlux & x,
    RKDG::ScalarFlux & y,
    const OpDomain op_domain        // = OpDomain::All
) {

    PRINT_STATUS( "Applying RKDG transport operator P.\n" )

# if defined (STRICT_CHECK)

    if ( !BitmaskHasAll( op_domain, OpDomain::Interior ) ) {

        std::string error_message = "OpDomain passed to '" + std::string(__func__)
                                    + "' does not contain flag 'OpDomain::Interior'\n";

        PRINT_ERROR( error_message.c_str() )
        throw std::out_of_range( error_message );
    }

    if (    !DomainDecomposition::AreMatching( x, y )
         || x.DG_degree != y.DG_degree
    ) {
        std::string error_message
            = "Parameter mismatch between RKDG::OrdinateFlux and RKDG::ScalarFlux objects in '"
            + std::string(__func__) + "'.\n";

        PRINT_ERROR( error_message.c_str() )
        throw std::runtime_error( error_message );
    }

# endif // if defined (STRICT_CHECK)

    Global::TMR_Pmv.Start();

    /* Determine offsets for handling boundary and halo cells.
     *
     * Offsets indicate whether the ghost cells along the specified local domain edge are included (1) or not
     * (0) in the operation.
     *
     * First index specifies spatial dimension.
     * Second index specifies either:
     *      index 0:    offset for lower bound of loop
     *      index 1:    offset for upper bound of loop
     * in the given spatial dimension.
     */
    int64_t offset[][2] = {     { 0, 0 }
                        # if SPACE_DIMS >= 2
                            ,   { 0, 0 }
                        # endif
                        # if SPACE_DIMS == 3
                            ,   { 0, 0 }
                        # endif
                          };

    x.DetermineOffsets( op_domain, offset );

    // Update status of halo cells in output object.
    if (/*
         * Set halos dirty if:
         *
         *  1.  Interior halo cells ARE NOT included in the operation; AND
         *  2.  (a) boundaries ARE periodic; OR
         *      (b) there is at least one reflecting boundary condition; OR
         *      (c) there is more than one MPI rank.
         *
         *  Otherwise preserve original values of input objects.
         */
            !BitmaskHasAny( op_domain, OpDomain::Halo )
         && (
                 Global::periodic
              || BitmaskHasAny( Global::reflecting_boundaries, BoundaryEdge::All )
        # if defined (ENABLE_MPI)
              || Global::MPI_num_ranks > 1
        # endif
        )
    ) {
        y.MarkAllHalosDirty();

    } else {

        y.halo_cells_dirty |= x.halo_cells_dirty;
    }

    const int64_t (& nx) [SPACE_DIMS] = y.nx();

# if SPACE_DIMS == 1

    # pragma omp parallel for schedule(static)
    for ( int64_t i = 1 - offset[0][0]; i <= nx[0] + offset[0][1]; ++i ) {

        for ( int64_t d = 0; d <= x.DG_degree; ++d )
            y(i,d) *= beta;

        for ( int64_t q = 0; q <  x.nq();      ++q ) {
        for ( int64_t d = 0; d <= x.DG_degree; ++d ) {

            y(i,d) += alpha * x.w(q) * x(q,i,d);
        }}
    }

# elif SPACE_DIMS == 2

    # pragma omp parallel for collapse(2) schedule(static)
    for ( int64_t i = 1 - offset[0][0]; i <= nx[0] + offset[0][1]; ++i ) {
    for ( int64_t j = 1 - offset[1][0]; j <= nx[1] + offset[1][1]; ++j ) {

        for ( int64_t d = 0; d <= x.DG_degree; ++d ) {
        for ( int64_t e = 0; e <= x.DG_degree; ++e ) {

            y(i,j,d,e) *= beta;
        }}

        for ( int64_t q = 0; q <  x.nq();      ++q ) {
        for ( int64_t d = 0; d <= x.DG_degree; ++d ) {
        for ( int64_t e = 0; e <= x.DG_degree; ++e ) {

            y(i,j,d,e) += alpha * x.w(q) * x(q,i,j,d,e);
        }}}
    }}

# elif SPACE_DIMS == 3

    # pragma omp parallel for collapse(3) schedule(static)
    for ( int64_t i = 1 - offset[0][0]; i <= nx[0] + offset[0][1]; ++i ) {
    for ( int64_t j = 1 - offset[1][0]; j <= nx[1] + offset[1][1]; ++j ) {
    for ( int64_t k = 1 - offset[2][0]; k <= nx[2] + offset[2][1]; ++k ) {

        for ( int64_t d = 0; d <= x.DG_degree; ++d ) {
        for ( int64_t e = 0; e <= x.DG_degree; ++e ) {
        for ( int64_t f = 0; f <= x.DG_degree; ++f ) {

            y(i,j,k,d,e,f) *= beta;
        }}}

        for ( int64_t q = 0; q <  x.nq();      ++q ) {
        for ( int64_t d = 0; d <= x.DG_degree; ++d ) {
        for ( int64_t e = 0; e <= x.DG_degree; ++e ) {
        for ( int64_t f = 0; f <= x.DG_degree; ++f ) {

            y(i,j,k,d,e,f) += alpha * x.w(q) * x(q,i,j,k,d,e,f);
        }}}}
    }}}

# endif // if SPACE_DIMS == ?

    Global::TMR_Pmv.Stop();
}


//------------------------------------------------------------------------------------------------------------
//! \brief  Computes \f$ \vec{y} \gets \beta \vec{y} + \alpha \mathcal{S} \vec{x} \f$ where
//!         \f$ \mathcal{S} \f$ is the angular redistribution operator.
//!
//! Computes the generalized matrix-vector product using the operator \f$ \mathcal{S} \f$.
//!
//! \param[in]      alpha       Scalar to augment \f$ \mathcal{S} \f$.
//! \param[in]      beta        Scalar to augment \f$ y \f$.
//! \param[in]      sigma       Cross section that defines \f$ \mathcal{S} \f$.
//! \param[in]      x           Scalar flux to which \f$ \mathcal{S} \f$ is applied.
//! \param[in,out]  y           Angular flux that is scaled and into which the result is stored.
//! \param[in]      op_domain   (optional) <br>
//!                             Domain over which to apply copy operation. Defaults to OpDomain::All.
//------------------------------------------------------------------------------------------------------------
void Smv (

    const double alpha,
    const double beta,
    const RKDG::CrossSection & sigma,
    const RKDG::ScalarFlux & x,
    RKDG::OrdinateFlux & y,
    const OpDomain op_domain            // = OpDomain::Full
) {

    PRINT_STATUS( "Applying RKDG transport operator S.\n" )

# if defined (STRICT_CHECK)

    if ( !BitmaskHasAll( op_domain, OpDomain::Interior ) ) {

        std::string error_message = "OpDomain passed to '" + std::string(__func__)
                                    + "' does not contain flag 'OpDomain::Interior'\n";

        PRINT_ERROR( error_message.c_str() )
        throw std::out_of_range( error_message );
    }

    if (    !DomainDecomposition::AreMatching( x, y )
         || x.DG_degree != y.DG_degree
    ) {
        std::string error_message
            = "Parameter mismatch between RKDG::OrdinateFlux and RKDG::ScalarFlux objects in '"
            + std::string(__func__) + "'.\n";

        PRINT_ERROR( error_message.c_str() )
        throw std::runtime_error( error_message );
    }

# endif // if defined (STRICT_CHECK)

    Global::TMR_AF_Smv.Start();

    /* Determine offsets for handling boundary and halo cells.
     *
     * Offsets indicate whether the ghost cells along the specified local domain edge are included (1) or not
     * (0) in the operation.
     *
     * First index specifies spatial dimension.
     * Second index specifies either:
     *      index 0:    offset for lower bound of loop
     *      index 1:    offset for upper bound of loop
     * in the given spatial dimension.
     */
    int64_t offset[][2] = {     { 0, 0 }
                        # if SPACE_DIMS >= 2
                            ,   { 0, 0 }
                        # endif
                        # if SPACE_DIMS == 3
                            ,   { 0, 0 }
                        # endif
                          };

    x.DetermineOffsets( op_domain, offset );

    // Update status of halo cells in output object.
    if (/*
         * Set halos dirty if:
         *
         *  1.  Interior halo cells ARE NOT included in the operation; AND
         *  2.  (a) boundaries ARE periodic; OR
         *      (b) there is at least one reflecting boundary condition; OR
         *      (c) there is more than one MPI rank.
         *
         *  Otherwise preserve original values of input objects.
         */
            !BitmaskHasAny( op_domain, OpDomain::Halo )
         && (
                 Global::periodic
              || BitmaskHasAny( Global::reflecting_boundaries, BoundaryEdge::All )
        # if defined (ENABLE_MPI)
              || Global::MPI_num_ranks > 1
        # endif
        )
    ) {
        y.MarkAllHalosDirty();

    } else {

        y.halo_cells_dirty |= x.halo_cells_dirty;
        y.upwind_halo_cells_dirty |= x.halo_cells_dirty;
    }

# if SPACE_DIMS == 1

    # pragma omp parallel for schedule(static)
    for ( int64_t i = 1 - offset[0][0]; i <= y.nx(0) + offset[0][1]; ++i ) {
    for ( int64_t d = 0; d <= x.DG_degree; ++d ) {

        double psi_val = 0.0;

        if ( sigma.DG_degree == 0 ) {

            psi_val = sigma(i,0) * x(i,d) / 2.0;

        } else {

            for ( int64_t u = 0; u <= x.DG_degree; ++u )
                psi_val += sigma(i,d,u) * x(i,u);

            psi_val *= (2*d + 1) / 4.0;
        }

        psi_val *= alpha;

        for ( int64_t q = 0; q < y.nq(); ++q )
            y(q,i,d) = beta * y(q,i,d) + psi_val;
    }}

# elif SPACE_DIMS == 2

//
// \brief  Indexing function for temporary array psi_val.
//
// \param[in]  d   Degree of test function with respect to \f$ x_1 \f$.
// \param[in]  e   Degree of test function with respect to \f$ x_2 \f$.
//
# define IPSIVAL(d,e) (psi_val[ (d) + (DG_degree + 1)*(e) ])

    const int64_t DG_degree = x.DG_degree;

    // Number of bytes of temporary memory to allocate per thread.
    const size_t bytes_allocd = (DG_degree + 1)*(DG_degree +1) * sizeof(double);

    # pragma omp parallel
    {
    # if defined (ENABLE_HWLOC) && defined (USE_HWLOC_ALLOC)

        int tid
            # if defined (_OPENMP)
                = omp_get_thread_num();
            # else
                = 0;
            # endif

        double * const psi_val = (double *)
            hwloc_alloc_membind( Global::machine_topology, bytes_allocd,
                                 Global::thread_masks[tid], HWLOC_MEMBIND_BIND, HWLOC_MEMBIND_THREAD );

    # else // if defined (ENABLE_HWLOC) && defined (USE_HWLOC_ALLOC)

        double * const psi_val = (double *)
            # if defined (USE_ALIGNED_ALLOC)
                aligned_alloc( ALIGNED_ALLOC_ALIGNMENT, bytes_allocd );
            # else
                std::malloc( bytes_allocd );
            # endif

    # endif // if defined (ENABLE_HWLOC) && defined (USE_HWLOC_ALLOC)

        const int64_t (& nx) [SPACE_DIMS] = y.nx();

        # pragma omp for collapse(2) schedule(static)
        for ( int64_t i = 1 - offset[0][0]; i <= nx[0] + offset[0][1]; ++i ) {
        for ( int64_t j = 1 - offset[1][0]; j <= nx[1] + offset[1][1]; ++j ) {

            for ( int64_t d = 0; d <= x.DG_degree; ++d ) {
            for ( int64_t e = 0; e <= x.DG_degree; ++e ) {

                IPSIVAL(d,e) = 0.0;

                if ( sigma.DG_degree == 0 ) {

                    IPSIVAL(d,e) = sigma(i,j,0,0) * x(i,j,d,e) / (4.0 * M_PI);

                } else {

                    for ( int64_t u = 0; u <= x.DG_degree; ++u ) {
                    for ( int64_t v = 0; v <= x.DG_degree; ++v ) {

                        IPSIVAL(d,e) += sigma(i,j,d,e,u,v) * x(i,j,u,v);
                    }}

                    IPSIVAL(d,e) *= (2*d + 1)*(2*e + 1) / (16.0 * M_PI);
                }

                IPSIVAL(d,e) *= alpha;
            }}

            for ( int64_t q = 0; q < y.nq(); ++q ) {
            for ( int64_t d = 0; d <= x.DG_degree; ++d ) {
            for ( int64_t e = 0; e <= x.DG_degree; ++e ) {

                y(q,i,j,d,e) = beta * y(q,i,j,d,e) + IPSIVAL(d,e);
            }}}
        }}

    # if defined (ENABLE_HWLOC) && defined (USE_HWLOC_ALLOC)
        hwloc_free( Global::machine_topology, psi_val, bytes_allocd );
    # else
        std::free( psi_val );
    # endif
    }

    # undef IXCOEFF

# elif SPACE_DIMS == 3

    const int64_t (& nx) [SPACE_DIMS] = y.nx();

    # pragma omp parallel for collapse(3) schedule(static)
    for ( int64_t i = 1 - offset[0][0]; i <= nx[0] + offset[0][1]; ++i ) {
    for ( int64_t j = 1 - offset[1][0]; j <= nx[1] + offset[1][1]; ++j ) {
    for ( int64_t k = 1 - offset[2][0]; k <= nx[2] + offset[2][1]; ++k ) {
    for ( int64_t d = 0; d <= x.DG_degree; ++d ) {
    for ( int64_t e = 0; e <= x.DG_degree; ++e ) {
    for ( int64_t f = 0; f <= x.DG_degree; ++f ) {

        double psi_val = 0.0;

        if ( sigma.DG_degree == 0 ) {

            psi_val = sigma(i,j,k,0,0,0) * x(i,j,k,d,e,f) / (4.0 * M_PI);

        } else {

            for ( int64_t u = 0; u <= x.DG_degree; ++u ) {
            for ( int64_t v = 0; v <= x.DG_degree; ++v ) {
            for ( int64_t w = 0; w <= x.DG_degree; ++w ) {

                psi_val += sigma(i,j,k,d,e,f,u,v,w) * x(i,j,k,u,v,w);
            }}}

            psi_val *= (2*d + 1)*(2*e + 1)*(2*f + 1) / (32.0 * M_PI);
        }

        psi_val *= alpha;

        for ( int64_t q = 0; q < y.nq(); ++q )
            y(q,i,j,k,d,e,f) = beta * y(q,i,j,k,d,e,f) + psi_val;
    }}}}}}

# endif // if SPACE_DIMS == ?

    Global::TMR_AF_Smv.Stop();
}


//============================================================================================================
//=== ONE-DIMENSIONAL OPERATOR ROUTINES ======================================================================
//============================================================================================================

# if SPACE_DIMS == 1


//------------------------------------------------------------------------------------------------------------
//! \brief  Computes \f$ \vec{y} \gets \beta \vec{y} + \alpha \mathcal{L} \vec{x} \f$ where
//!         \f$ \mathcal{L} \f$ is the advection-absorption operator with total cross-section given by
//!         \pp{sigma}.
//!         Computes the 1D case.
//!
//! Generalized matrix-vector product using the operator \f$ \mathcal{L} \f$.
//!
//! \param[in]      alpha       Scalar to augment \f$ \mathcal{L} \f$.
//! \param[in]      beta        Scalar to augment \f$ \vec{y} \f$.
//! \param[in]      sigma       Cross section which defines \f$ \mathcal{L} \f$.
//! \param[in]      x           Angular flux to which \f$ \mathcal{L} \f$ is applied.
//! \param[in,out]  y           Angular flux that is scaled and into which the result is stored.
//------------------------------------------------------------------------------------------------------------
void Lmv (

    const double alpha,
    const double beta,
    const RKDG::CrossSection & sigma,
    RKDG::OrdinateFlux & x,
    RKDG::OrdinateFlux & y
) {

    PRINT_STATUS( "Applying RKDG transport operator L.\n" )

    if ( x.UpwindHalosAreDirty() ) {  x.SynchronizeHalos();  }

    Global::TMR_AF_Lmv.Start();

    # pragma omp parallel
    {
        double coeff;
        double * const x_coeff = new double[ x.DG_degree + 1 ];

        # define IXCOEFF(d) (x_coeff[(d)])

        # pragma omp for
        for ( int64_t q = 0; q <  y.nq();  ++q ) {
        for ( int64_t i = 1; i <= y.nx(0); ++i ) {

            int64_t ii = i;

            if ( y.xi(q) > 0 ) {  ii = y.nx(0) - i + 1;  }

            // Save local DG coefficients into x_coeff array.
            for ( int64_t d = 0; d <= y.DG_degree; ++d ) {

                IXCOEFF(d) = x(q,ii,d);
                y(q,ii,d) *= beta;
            }

            // Set value from absorption term.
            if ( sigma.DG_degree == 0 ) {

                for ( int64_t d = 0; d <= y.DG_degree; ++d )
                    y(q,ii,d) += alpha * sigma(ii,0) * IXCOEFF(d);

            } else {

                for ( int64_t d = 0; d <= y.DG_degree; ++d ) {
                for ( int64_t u = 0; u <= y.DG_degree; ++u ) {

                    y(q,ii,d) += alpha * (2*d + 1) * sigma(ii,d,u) * IXCOEFF(u) / 2.0;
                }}
            }

            // Include components of advection term for ξ angular direction.
            coeff = alpha * y.xi(q) / y.dx(0);

            if ( y.xi(q) < 0 ) {

                for ( int64_t d = 0; d <= y.DG_degree; ++d ) {
                for ( int64_t a = 0; a <= y.DG_degree; ++a ) {

                    y(q,ii,d) += coeff * (2*d + 1) * ( neg1pow(a) * x(q,ii+1,a)
                                                       + neg1pow(d+a+1) * IXCOEFF(a) );
                }}

            } else {

                for ( int64_t d = 0; d <= y.DG_degree; ++d ) {
                for ( int64_t a = 0; a <= y.DG_degree; ++a ) {

                    y(q,ii,d) += coeff * (2*d + 1) * ( IXCOEFF(a) + neg1pow(d+1) * x(q,ii-1,a) );
                }}
            }

            coeff *= -2.0;

            for ( int64_t d = 0; d <= y.DG_degree; ++d ) {
            for ( int64_t a = d-1; a >= 0; a -= 2 ) {

                y(q,ii,d) += coeff * (2*d + 1) * IXCOEFF(a);
            }}
        }}

        delete [] x_coeff;
        # undef IXCOEFF
    }

    if (/*
         * Set halos dirty if:
         *
         *  1.  Boundaries ARE periodic; OR
         *  2.  There is at least one reflecting boundary condition; OR
         *  3.  There is more than one MPI rank.
         */
            Global::periodic
         || BitmaskHasAny( Global::reflecting_boundaries, BoundaryEdge::All )
    # if defined (ENABLE_MPI)
         || Global::MPI_num_ranks > 1
    # endif
    ) {
        y.MarkAllHalosDirty();
    }

    Global::TMR_AF_Lmv.Stop();
}


# endif // if SPACE_DIMS == 1


//============================================================================================================
//=== TWO-DIMENSIONAL OPERATOR ROUTINES ======================================================================
//============================================================================================================

# if SPACE_DIMS == 2 || defined (DOXYCOMPILE)


//------------------------------------------------------------------------------------------------------------
//! \brief  Computes \f$ \vec{y} \gets \beta \vec{y} + \alpha \mathcal{L} \vec{x} \f$ where
//!         \f$ \mathcal{L} \f$ is the advection-absorption operator with total cross-section given by
//!         \pp{sigma}.
//!         Computes the 2D case.
//!
//! Generalized matrix-vector product using the operator \f$ \mathcal{L} \f$. Computation is performed using a
//! wavefront sweep technique over the spatial cells that permits \pp{x} and \pp{y} to be the same
//! OrdinateFlux object.
//!
//! \attention  This routine is automatically called by TransportOperator::Lmv if \pp{x} and \pp{y} are
//!             determined to be the same object.
//!
//! \param[in]      alpha       Scalar to augment \f$ \mathcal{L} \f$.
//! \param[in]      beta        Scalar to augment \f$ \vec{y} \f$.
//! \param[in]      sigma       Cross section which defines \f$ \mathcal{L} \f$.
//! \param[in]      x           Angular flux to which \f$ \mathcal{L} \f$ is applied.
//! \param[in,out]  y           Angular flux that is scaled and into which the result is stored.
//------------------------------------------------------------------------------------------------------------
static void Lmv_InPlace (

    const double alpha,
    const double beta,
    const RKDG::CrossSection & sigma,
    RKDG::OrdinateFlux & x,
    RKDG::OrdinateFlux & y
) {

    PRINT_STATUS( "Applying RKDG transport operator L using in-place routine.\n" )

    if ( x.UpwindHalosAreDirty() ) {  x.SynchronizeHalos();  }

    Global::TMR_AF_Lmv.Start();

    const int64_t num_diags = y.nx(0) + y.nx(1) - 1;

# if defined (ENABLE_SIMD_BLOCKING)

//!
//! \brief  Indexing function for temporary array for caching values from result.
//!
//! \param[in]  k   Index of temporary vector. <br/>
//!                 0: Input values from current cell,
//!                 1: Input values from x-upwind cell,
//!                 2: Input values from y-upwind cell,
//!                 3: Output values for current cell.
//! \param[in]  d   Degree of test function with respect to \f$ x_1 \f$.
//! \param[in]  e   Degree of test function with respect to \f$ x_2 \f$.
//! \param[in]  l   Index for system in SIMD block. In [ 0, LMV_SIMD_LEN ).
//!
# define IXCOEFF(k,d,e,l) (x_coeff[ (l) + (LMV_SIMD_LEN)*( \
                                    (e) + (DG_degree + 1)*( \
                                    (d) + (DG_degree + 1)*(k) \
                           ) ) ])

    const int64_t DG_degree = x.DG_degree;

# if defined (ENABLE_SIMD_BLOCKING)

    // Determine ordinate blocking structure.
    int64_t num_in_quad [4];                // Number of ordinates in each quadrant.
    int64_t leftover_angs [4*LMV_SIMD_LEN]; // Contains indices of leftover ordinates after blocking.
    int64_t num_leftover_angs = 0;          // Total number of leftover ordinates.

    for ( int64_t quad = 0; quad < 4; ++quad ) {

        num_in_quad[quad] = y.Quadrants(quad + 1) - y.Quadrants(quad);

        for ( int64_t q = y.Quadrants(quad + 1) - num_in_quad[quad] % ((int64_t)(LMV_SIMD_LEN));
                      q < y.Quadrants(quad + 1);
                      q++
        ) {  leftover_angs[ num_leftover_angs++ ] = q;  }
    }

# endif // if defined (ENABLE_SIMD_BLOCKING)

    // Number of bytes of temporary memory to allocate per thread.
    const size_t bytes_allocd = 4 * LMV_SIMD_LEN * (DG_degree + 1)*(DG_degree +1) * sizeof(double);

    # pragma omp parallel
    {
        double coeff;
        alignas( LMV_SIMD_LEN * sizeof(double) ) double coeff_v [LMV_SIMD_LEN];

    # if defined (ENABLE_HWLOC) && defined (USE_HWLOC_ALLOC)

        int tid
            # if defined (_OPENMP)
                = omp_get_thread_num();
            # else
                = 0;
            # endif

        double * const x_coeff = (double *)
            hwloc_alloc_membind( Global::machine_topology, bytes_allocd, Global::thread_masks[tid],
                                 HWLOC_MEMBIND_BIND, HWLOC_MEMBIND_THREAD );
    # else

        double * const x_coeff = (double *) aligned_alloc( ALIGNED_ALLOC_ALIGNMENT, bytes_allocd );

    # endif // if defined (ENABLE_HWLOC) && defined (USE_HWLOC_ALLOC)

        // Sweep using diagonal wavefront.
        for ( int64_t diag = 1; diag <= num_diags; ++diag ) {

            int64_t diag_length = std::min( diag, std::min( num_diags - diag + 1, std::min( y.nx(0), y.nx(1) ) ) );

            // First peel off blocks of type SameCell_ContiguousAngles.
            for ( int64_t quad = 0; quad < 4; ++quad ) {

                const int64_t quad_bounds [] = { y.Quadrants(quad), y.Quadrants(quad + 1) };

            # pragma omp for collapse(2) nowait
            for ( int64_t q = quad_bounds[0];
                          q < quad_bounds[1] - (num_in_quad[quad] % LMV_SIMD_LEN) - (LMV_SIMD_LEN - 1);
                          q += LMV_SIMD_LEN
            ) {
            for ( int64_t cell = 1; cell <= diag_length; ++cell ) {

                // Determine sign on angles in block.
                bool all_xi_negative = true;
                bool all_eta_negative = true;

                for ( int64_t l = 0; l < LMV_SIMD_LEN; ++l ) {

                    all_xi_negative  &= x.xi(q+l)  <= 0.0;
                    all_eta_negative &= x.eta(q+l) <= 0.0;
                }

                const int64_t xoffset = ( all_xi_negative  ? +1 : -1 );
                const int64_t yoffset = ( all_eta_negative ? +1 : -1 );

                // Compute spatial cell indices.
                int64_t i = std::min( diag, y.nx(0) ) - cell + 1;
                int64_t j = diag - i + 1;

                // Reverse direction of sweep for positive angles.
                if ( !all_xi_negative  ) {  i = y.nx(0) - i + 1;  }
                if ( !all_eta_negative ) {  j = y.nx(1) - j + 1;  }

                // Load values into temporary storage.
                for ( int64_t d = 0; d <= DG_degree; ++d ) {
                for ( int64_t e = 0; e <= DG_degree; ++e ) {
                # pragma omp simd aligned( x_coeff : LMV_SIMD_LEN * sizeof(double) )
                for ( int64_t l = 0; l < LMV_SIMD_LEN; ++l ) {

                    IXCOEFF(0,d,e,l) = x(q+l,i,j,d,e);

                    IXCOEFF(1,d,e,l) = x(q+l, i + xoffset, j, d,e);
                    IXCOEFF(2,d,e,l) = x(q+l, i, j + yoffset, d,e);

                    IXCOEFF(3,d,e,l) = 0.0;
                }}}

                // Set value from absorption term.
                if ( sigma.DG_degree == 0 ) {

                    const double sigma_val = sigma(i,j,0,0);

                    for ( int64_t d = 0; d <= DG_degree; ++d ) {
                    for ( int64_t e = 0; e <= DG_degree; ++e ) {
                    # pragma omp simd aligned( x_coeff : LMV_SIMD_LEN * sizeof(double) )
                    for ( int64_t l = 0; l < LMV_SIMD_LEN; ++l ) {

                        IXCOEFF(3,d,e,l) += sigma_val * IXCOEFF(0,d,e,l);
                    }}}

                } else {

                    for ( int64_t d = 0; d <= DG_degree; ++d ) {
                    for ( int64_t e = 0; e <= DG_degree; ++e ) {
                    for ( int64_t u = 0; u <= DG_degree; ++u ) {
                    for ( int64_t v = 0; v <= DG_degree; ++v ) {

                        const double sigma_val = 0.25 * (2*d + 1)*(2*e + 1) * sigma(i,j,d,e,u,v);

                        # pragma omp simd aligned( x_coeff : LMV_SIMD_LEN * sizeof(double) )
                        for ( int64_t l = 0; l < LMV_SIMD_LEN; ++l )
                            IXCOEFF(3,d,e,l) += sigma_val * IXCOEFF(0,u,v,l);
                    }}}}
                }

                // Include components of advection term for ξ angular direction.
                # pragma omp simd aligned( x_coeff, coeff_v : LMV_SIMD_LEN * sizeof(double) )
                for ( int64_t l = 0; l < LMV_SIMD_LEN; ++l )
                    coeff_v[l] = x.xi(q+l) / x.dx(0);

                if ( all_xi_negative ) {

                    for ( int64_t d = 0; d <= DG_degree; ++d ) {
                    for ( int64_t a = 0; a <= DG_degree; ++a ) {
                    for ( int64_t e = 0; e <= DG_degree; ++e ) {
                    # pragma omp simd aligned( x_coeff, coeff_v : LMV_SIMD_LEN * sizeof(double) )
                    for ( int64_t l = 0; l < LMV_SIMD_LEN; ++l ) {

                        IXCOEFF(3,d,e,l) += coeff_v[l] * (2*d + 1) * ( neg1pow(a) * IXCOEFF(1,a,e,l)
                                                                       + neg1pow(d+a+1) * IXCOEFF(0,a,e,l) );
                    }}}}

                } else {

                    for ( int64_t d = 0; d <= DG_degree; ++d ) {
                    for ( int64_t a = 0; a <= DG_degree; ++a ) {
                    for ( int64_t e = 0; e <= DG_degree; ++e ) {
                    # pragma omp simd aligned( x_coeff, coeff_v : LMV_SIMD_LEN * sizeof(double) )
                    for ( int64_t l = 0; l < LMV_SIMD_LEN; ++l ) {

                        IXCOEFF(3,d,e,l) += coeff_v[l] * (2*d + 1) * ( IXCOEFF(0,a,e,l)
                                                                       + neg1pow(d+1) * IXCOEFF(1,a,e,l) );
                    }}}}
                }

                for ( int64_t d = 0; d <= DG_degree; ++d ) {
                for ( int64_t a = d-1; a >= 0; a -= 2 ) {
                for ( int64_t e = 0; e <= DG_degree; ++e ) {
                # pragma omp simd aligned( x_coeff, coeff_v : LMV_SIMD_LEN * sizeof(double) )
                for ( int64_t l = 0; l < LMV_SIMD_LEN; ++l ) {

                    IXCOEFF(3,d,e,l) -= ( 2 * (2*d + 1) ) * coeff_v[l] * IXCOEFF(0,a,e,l);
                }}} }

                // Include components of advection term for η angular direction.
                # pragma omp simd aligned( x_coeff, coeff_v : LMV_SIMD_LEN * sizeof(double) )
                for ( int64_t l = 0; l < LMV_SIMD_LEN; ++l )
                    coeff_v[l] = x.eta(q+l) / x.dx(1);

                if ( all_eta_negative ) {

                    for ( int64_t d = 0; d <= DG_degree; ++d ) {
                    for ( int64_t e = 0; e <= DG_degree; ++e ) {
                    for ( int64_t b = 0; b <= DG_degree; ++b ) {
                    # pragma omp simd aligned( x_coeff, coeff_v : LMV_SIMD_LEN * sizeof(double) )
                    for ( int64_t l = 0; l < LMV_SIMD_LEN; ++l ) {

                        IXCOEFF(3,d,e,l) += coeff_v[l] * (2*e + 1) * ( neg1pow(b) * IXCOEFF(2,d,b,l)
                                                                       + neg1pow(e+b+1) * IXCOEFF(0,d,b,l) );
                    }}}}

                } else {

                    for ( int64_t d = 0; d <= DG_degree; ++d ) {
                    for ( int64_t e = 0; e <= DG_degree; ++e ) {
                    for ( int64_t b = 0; b <= DG_degree; ++b ) {
                    # pragma omp simd aligned( x_coeff, coeff_v : LMV_SIMD_LEN * sizeof(double) )
                    for ( int64_t l = 0; l < LMV_SIMD_LEN; ++l ) {

                        IXCOEFF(3,d,e,l) += coeff_v[l] * (2*e + 1) * ( IXCOEFF(0,d,b,l)
                                                                       + neg1pow(e+1) * IXCOEFF(2,d,b,l) );
                    }}}}
                }

                for ( int64_t d = 0; d <= DG_degree; ++d ) {
                for ( int64_t e = 0; e <= DG_degree; ++e ) {
                for ( int64_t b = e-1; b >= 0; b -= 2 ) {
                # pragma omp simd aligned( x_coeff, coeff_v : LMV_SIMD_LEN * sizeof(double) )
                for ( int64_t l = 0; l < LMV_SIMD_LEN; ++l ) {

                    IXCOEFF(3,d,e,l) -= ( 2 * (2*e + 1) ) * coeff_v[l] * IXCOEFF(0,d,b,l);
                }}}}

                // Store output values.
                for ( int64_t d = 0; d <= DG_degree; ++d ) {
                for ( int64_t e = 0; e <= DG_degree; ++e ) {
                # pragma omp simd aligned( x_coeff : LMV_SIMD_LEN * sizeof(double) )
                for ( int64_t l = 0; l < LMV_SIMD_LEN; ++l ) {

                    y(q+l,i,j,d,e) = beta * y(q+l,i,j,d,e) + alpha * IXCOEFF(3,d,e,l);
                }}}
            }}}

            // Then compute remaining blocks of type SameAngle.
            # pragma omp for collapse(2) nowait
            for ( int64_t cell = 1; cell <= diag_length; cell += LMV_SIMD_LEN ) {
            for ( int64_t q_leftover = 0; q_leftover < num_leftover_angs; ++q_leftover ) {

                const int64_t q = leftover_angs[ q_leftover ];

                SIMD_BlkIdx<LMV_SIMD_LEN> idx {};
                idx.type = BlockType::SameAngle;

                // Truncate last block.
                idx.len = std::min( ((int64_t)(LMV_SIMD_LEN)), diag_length - cell + 1 );

                for ( int64_t l = 0; l < idx.len; ++l ) {

                    // Compute spatial cell indices.
                    int64_t i = std::min( diag, y.nx(0) ) - (cell + l) + 1;
                    int64_t j = diag - i + 1;

                    // Reverse direction of sweep for positive angles.
                    if ( y.xi(q)  > 0 ) {  i = y.nx(0) - i + 1;  }
                    if ( y.eta(q) > 0 ) {  j = y.nx(1) - j + 1;  }

                    idx.i[l] = i;
                    idx.j[l] = j;
                    idx.q[l] = q;
                }

                // Pad extra spaces in truncated SIMD blocks with duplicate values.
                for ( int64_t l = idx.len; l < LMV_SIMD_LEN; ++l ) {

                    idx.i[l] = idx.i[idx.len - 1];
                    idx.j[l] = idx.j[idx.len - 1];
                    idx.q[l] = idx.q[idx.len - 1];
                }

                const int64_t xoffset = ( y.xi(q)  < 0 ? +1 : -1 );
                const int64_t yoffset = ( y.eta(q) < 0 ? +1 : -1 );

                // Load values into temporary storage.
                for ( int64_t d = 0; d <= DG_degree; ++d ) {
                for ( int64_t e = 0; e <= DG_degree; ++e ) {
                for ( int64_t l = 0; l < LMV_SIMD_LEN; ++l ) {

                    IXCOEFF(0,d,e,l) = x(q, idx.i[l], idx.j[l], d,e);

                    IXCOEFF(1,d,e,l) = x(q, idx.i[l] + xoffset, idx.j[l], d,e);
                    IXCOEFF(2,d,e,l) = x(q, idx.i[l], idx.j[l] + yoffset, d,e);

                    IXCOEFF(3,d,e,l) = 0.0;
                }}}

                // Set value from absorption term.
                if ( sigma.DG_degree == 0 ) {

                    for ( int64_t d = 0; d <= DG_degree; ++d ) {
                    for ( int64_t e = 0; e <= DG_degree; ++e ) {
                    for ( int64_t l = 0; l < LMV_SIMD_LEN; ++l ) {

                        IXCOEFF(3,d,e,l) += sigma( idx.i[l], idx.j[l], 0,0) * IXCOEFF(0,d,e,l);
                    }}}

                } else {

                    for ( int64_t u = 0; u <= DG_degree; ++u ) {
                    for ( int64_t v = 0; v <= DG_degree; ++v ) {
                    for ( int64_t d = 0; d <= DG_degree; ++d ) {
                    for ( int64_t e = 0; e <= DG_degree; ++e ) {
                    for ( int64_t l = 0; l < LMV_SIMD_LEN; ++l ) {

                        IXCOEFF(3,d,e,l) += 0.25 * (2*d + 1)*(2*e + 1) * sigma( idx.i[l], idx.j[l], d,e,u,v)
                                            * IXCOEFF(0,u,v,l);
                    }}}}}
                }

                // Include components of advection term for ξ angular direction.
                coeff = y.xi(q) / y.dx(0);

                if ( y.xi(q) < 0 ) {

                    for ( int64_t d = 0; d <= DG_degree; ++d ) {
                    for ( int64_t a = 0; a <= DG_degree; ++a ) {
                    for ( int64_t e = 0; e <= DG_degree; ++e ) {
                    # pragma omp simd aligned( x_coeff : LMV_SIMD_LEN * sizeof(double) )
                    for ( int64_t l = 0; l < LMV_SIMD_LEN; ++l ) {

                        IXCOEFF(3,d,e,l) += coeff * (2*d + 1) * ( neg1pow(a) * IXCOEFF(1,a,e,l)
                                                                  + neg1pow(d+a+1) * IXCOEFF(0,a,e,l) );
                    }}}}

                } else {

                    for ( int64_t d = 0; d <= DG_degree; ++d ) {
                    for ( int64_t a = 0; a <= DG_degree; ++a ) {
                    for ( int64_t e = 0; e <= DG_degree; ++e ) {
                    # pragma omp simd aligned( x_coeff : LMV_SIMD_LEN * sizeof(double) )
                    for ( int64_t l = 0; l < LMV_SIMD_LEN; ++l ) {

                        IXCOEFF(3,d,e,l) += coeff * (2*d + 1) * ( IXCOEFF(0,a,e,l)
                                                                  + neg1pow(d+1) * IXCOEFF(1,a,e,l) );
                    }}}}
                }

                for ( int64_t d = 0; d <= DG_degree; ++d ) {
                for ( int64_t a = d-1; a >= 0; a -= 2 ) {
                for ( int64_t e = 0; e <= DG_degree; ++e ) {
                # pragma omp simd aligned( x_coeff : LMV_SIMD_LEN * sizeof(double) )
                for ( int64_t l = 0; l < LMV_SIMD_LEN; ++l ) {

                    IXCOEFF(3,d,e,l) -= ( 2 * (2*d + 1) ) * coeff * IXCOEFF(0,a,e,l);
                }}} }

                // Include components of advection term for η angular direction.
                coeff = y.eta(q) / y.dx(1);

                if ( y.eta(q) < 0 ) {

                    for ( int64_t d = 0; d <= DG_degree; ++d ) {
                    for ( int64_t e = 0; e <= DG_degree; ++e ) {
                    for ( int64_t b = 0; b <= DG_degree; ++b ) {
                    # pragma omp simd aligned( x_coeff : LMV_SIMD_LEN * sizeof(double) )
                    for ( int64_t l = 0; l < LMV_SIMD_LEN; ++l ) {

                        IXCOEFF(3,d,e,l) += coeff * (2*e + 1) * ( neg1pow(b) * IXCOEFF(2,d,b,l)
                                                                  + neg1pow(e+b+1) * IXCOEFF(0,d,b,l) );
                    }}}}

                } else {

                    for ( int64_t d = 0; d <= DG_degree; ++d ) {
                    for ( int64_t e = 0; e <= DG_degree; ++e ) {
                    for ( int64_t b = 0; b <= DG_degree; ++b ) {
                    # pragma omp simd aligned( x_coeff : LMV_SIMD_LEN * sizeof(double) )
                    for ( int64_t l = 0; l < LMV_SIMD_LEN; ++l ) {

                        IXCOEFF(3,d,e,l) += coeff * (2*e + 1) * ( IXCOEFF(0,d,b,l)
                                                                  + neg1pow(e+1) * IXCOEFF(2,d,b,l) );
                    }}}}
                }

                for ( int64_t d = 0; d <= DG_degree; ++d ) {
                for ( int64_t e = 0; e <= DG_degree; ++e ) {
                for ( int64_t b = e-1; b >= 0; b -= 2 ) {
                # pragma omp simd aligned( x_coeff : LMV_SIMD_LEN * sizeof(double) )
                for ( int64_t l = 0; l < LMV_SIMD_LEN; ++l ) {

                    IXCOEFF(3,d,e,l) -= ( 2 * (2*e + 1) ) * coeff * IXCOEFF(0,d,b,l);
                }}}}

                // Store output values.
                for ( int64_t d = 0; d <= DG_degree; ++d ) {
                for ( int64_t e = 0; e <= DG_degree; ++e ) {
                for ( int64_t l = 0; l < idx.len; ++l ) {

                    y(q, idx.i[l], idx.j[l], d,e) = beta * y(q, idx.i[l], idx.j[l], d,e)
                                                    + alpha * IXCOEFF(3,d,e,l);
                }}}
            }}

            # pragma omp barrier
        }

    # if defined (ENABLE_HWLOC) && defined (USE_HWLOC_ALLOC)
        hwloc_free( Global::machine_topology, x_coeff, bytes_allocd );
    # else
        std::free( x_coeff );
    # endif

        # undef IXCOEFF

    } // end OMP parallel region.

# else // if defined (ENABLE_SIMD_BLOCKING)

    // Number of bytes of temporary memory to allocate per thread.
    const size_t bytes_allocd = (x.DG_degree + 1)*(x.DG_degree +1) * sizeof(double);

    # pragma omp parallel
    {
        double coeff;

    # if defined (ENABLE_HWLOC) && defined (USE_HWLOC_ALLOC)

        int tid
            # if defined (_OPENMP)
                = omp_get_thread_num();
            # else
                = 0;
            # endif

        double * const x_coeff = (double *)
            hwloc_alloc_membind( Global::machine_topology, bytes_allocd, Global::thread_masks[tid],
                                 HWLOC_MEMBIND_BIND, HWLOC_MEMBIND_THREAD );
    # else

        double * const x_coeff = (double *)
            # if defined (USE_ALIGNED_ALLOC)
                aligned_alloc( ALIGNED_ALLOC_ALIGNMENT, bytes_allocd );
            # else
                std::malloc( bytes_allocd );
            # endif

    # endif // if defined (ENABLE_HWLOC) && defined (USE_HWLOC_ALLOC)

        # define IXCOEFF(d,e) (x_coeff[ (e) + (x.DG_degree + 1)*(d) ])

        // Sweep using diagonal wavefront.
        for ( int64_t diag = 1; diag <= num_diags; ++diag ) {

            int64_t diag_length = std::min( diag, std::min( num_diags - diag + 1, std::min( y.nx(0), y.nx(1) ) ) );

            const int64_t nq = y.nq();

            # pragma omp for collapse(2)
            for ( int64_t cell = 1; cell <= diag_length; ++cell ) {
            for ( int64_t q = 0; q < nq; ++q ) {

                int64_t i = std::min( diag, y.nx(0) ) - cell + 1;
                int64_t j = diag - i + 1;

                // Reverse direction of sweep for negative angles.
                if ( y.xi(q)  > 0 ) {  i = y.nx(0) - i + 1;  }
                if ( y.eta(q) > 0 ) {  j = y.nx(1) - j + 1;  }

                // Save local DG coefficients into x_coeff array.
                for ( int64_t d = 0; d <= x.DG_degree; ++d ) {
                for ( int64_t e = 0; e <= x.DG_degree; ++e ) {

                    IXCOEFF(d,e) = x(q,i,j,d,e);
                    y(q,i,j,d,e) *= beta;
                }}

                // Set value from absorption term.
                if ( sigma.DG_degree == 0 ) {

                    for ( int64_t d = 0; d <= x.DG_degree; ++d ) {
                    for ( int64_t e = 0; e <= x.DG_degree; ++e ) {

                        y(q,i,j,d,e) += alpha * sigma(i,j,0,0) * IXCOEFF(d,e);
                    }}

                } else {

                    for ( int64_t d = 0; d <= x.DG_degree; ++d ) {
                    for ( int64_t e = 0; e <= x.DG_degree; ++e ) {
                    for ( int64_t u = 0; u <= x.DG_degree; ++u ) {
                    for ( int64_t v = 0; v <= x.DG_degree; ++v ) {

                        y(q,i,j,d,e) += alpha * (2*d + 1)*(2*e + 1) * sigma(i,j,d,e,u,v) * IXCOEFF(u,v) / 4.0;
                    }}}}
                }

                // Include components of advection term for ξ angular direction.
                coeff = alpha * x.xi(q) / x.dx(0);

                if ( x.xi(q) < 0 ) {

                    for ( int64_t d = 0; d <= x.DG_degree; ++d ) {
                    for ( int64_t a = 0; a <= x.DG_degree; ++a ) {
                    for ( int64_t e = 0; e <= x.DG_degree; ++e ) {

                        y(q,i,j,d,e) += coeff * (2*d + 1) * ( neg1pow(a) * x(q,i+1,j,a,e)
                                                              + neg1pow(d+a+1) * IXCOEFF(a,e) );
                    }}}

                } else {

                    for ( int64_t d = 0; d <= x.DG_degree; ++d ) {
                    for ( int64_t a = 0; a <= x.DG_degree; ++a ) {
                    for ( int64_t e = 0; e <= x.DG_degree; ++e ) {

                        y(q,i,j,d,e) += coeff * (2*d + 1) * ( IXCOEFF(a,e) + neg1pow(d+1) * x(q,i-1,j,a,e) );
                    }}}
                }

                coeff *= -2.0;

                for ( int64_t d = 0; d <= x.DG_degree; ++d ) {
                for ( int64_t a = d-1; a >= 0; a -= 2 ) {
                for ( int64_t e = 0; e <= x.DG_degree; ++e ) {

                    y(q,i,j,d,e) += coeff * (2*d + 1) * IXCOEFF(a,e);
                }}}

                // Include components of advection term for η angular direction.
                coeff = alpha * x.eta(q) / x.dx(1);

                if ( x.eta(q) < 0 ) {

                    for ( int64_t d = 0; d <= x.DG_degree; ++d ) {
                    for ( int64_t e = 0; e <= x.DG_degree; ++e ) {
                    for ( int64_t b = 0; b <= x.DG_degree; ++b ) {

                        y(q,i,j,d,e) += coeff * (2*e + 1) * ( neg1pow(b) * x(q,i,j+1,d,b)
                                                              + neg1pow(e+b+1) * IXCOEFF(d,b) );
                    }}}

                } else {

                    for ( int64_t d = 0; d <= x.DG_degree; ++d ) {
                    for ( int64_t e = 0; e <= x.DG_degree; ++e ) {
                    for ( int64_t b = 0; b <= x.DG_degree; ++b ) {

                        y(q,i,j,d,e) += coeff * (2*e + 1) * ( IXCOEFF(d,b) + neg1pow(e+1) * x(q,i,j-1,d,b) );
                    }}}
                }

                coeff *= -2.0;

                for ( int64_t d = 0; d <= x.DG_degree; ++d ) {
                for ( int64_t e = 0; e <= x.DG_degree; ++e ) {
                for ( int64_t b = e-1; b >= 0; b -= 2 ) {

                    y(q,i,j,d,e) += coeff * (2*e + 1) * IXCOEFF(d,b);
                }}}
            }}
        }

    # if defined (ENABLE_HWLOC) && defined (USE_HWLOC_ALLOC)
        hwloc_free( Global::machine_topology, x_coeff, bytes_allocd );
    # else
        std::free( x_coeff );
    # endif

        # undef IXCOEFF

    } // end OMP parallel region.

# endif // if defined (ENABLE_SIMD_BLOCKING)

    if (/*
         * Set halos dirty if:
         *
         *  1.  Boundaries ARE periodic; OR
         *  2.  There is at least one reflecting boundary condition; OR
         *  3.  There is more than one MPI rank.
         */
            Global::periodic
         || BitmaskHasAny( Global::reflecting_boundaries, BoundaryEdge::All )
    # if defined (ENABLE_MPI)
         || Global::MPI_num_ranks > 1
    # endif
    ) {
        y.MarkAllHalosDirty();
    }

    Global::TMR_AF_Lmv.Stop();
}


//------------------------------------------------------------------------------------------------------------
//! \brief  Computes \f$ \vec{y} \gets \beta \vec{y} + \alpha \mathcal{L} \vec{x} \f$ where
//!         \f$ \mathcal{L} \f$ is the advection-absorption operator with total cross section given by
//!         \pp{sigma}.
//!         Computes the 2D case.
//!
//! \param[in]      alpha       Scalar to augment \f$ \mathcal{L} \f$.
//! \param[in]      beta        Scalar to augment \f$ \vec{y} \f$.
//! \param[in]      sigma       Cross section that defines \f$ \mathcal{L} \f$.
//! \param[in]      x           Angular flux to which \f$ \mathcal{L} \f$ is applied.
//! \param[in,out]  y           Angular flux that is scaled and into which the result is stored.
//------------------------------------------------------------------------------------------------------------
void Lmv (

    const double alpha,
    const double beta,
    const RKDG::CrossSection & sigma,
    RKDG::OrdinateFlux & x,
    RKDG::OrdinateFlux & y
) {

    if ( x == y ) {

        Lmv_InPlace( alpha, beta, sigma, x, y );
        return;
    }

    PRINT_STATUS( "Applying RKDG transport operator L.\n" )

    if ( x.UpwindHalosAreDirty() ) {  x.SynchronizeHalos();  }

    Global::TMR_AF_Lmv.Start();

# if defined (ENABLE_SIMD_BLOCKING)

//!
//! \brief  Indexing function for temporary array for caching values from result.
//!
//! \param[in]  k   Index of temporary vector. <br/>
//!                 0: Input values from current cell,
//!                 1: Input values from x-upwind cell,
//!                 2: Input values from y-upwind cell,
//!                 3: Output values for current cell.
//! \param[in]  d   Degree of test function with respect to \f$ x_1 \f$.
//! \param[in]  e   Degree of test function with respect to \f$ x_2 \f$.
//! \param[in]  l   Index for system in SIMD block. In [ 0, LMV_SIMD_LEN ).
//!
# define IXCOEFF(k,d,e,l) (x_coeff[ (l) + (LMV_SIMD_LEN)*( \
                                    (e) + (DG_degree + 1)*( \
                                    (d) + (DG_degree + 1)*(k) \
                           ) ) ])

    const int64_t DG_degree = x.DG_degree;

    // Determine how much temporary memory should be allocated.
    bool use_temp_mem = false;
    size_t bytes_mallocd = 0;

    for ( int64_t quad = 0; quad < 4; ++quad )
        use_temp_mem |= 0 != ( x.Quadrants(quad +1) - x.Quadrants(quad) ) / ((int64_t)(LMV_SIMD_LEN));

    if ( use_temp_mem )
        bytes_mallocd = 4 * LMV_SIMD_LEN * (DG_degree +1)*(DG_degree +1) * sizeof(double);

    # pragma omp parallel
    {
        double coeff;
        alignas( LMV_SIMD_LEN * sizeof(double) ) double coeff_v [LMV_SIMD_LEN];
        double * x_coeff = nullptr;

        if ( use_temp_mem ) {
        # if defined (ENABLE_HWLOC) && defined (USE_HWLOC_ALLOC)

            int tid
                # if defined (_OPENMP)
                    = omp_get_thread_num();
                # else
                    = 0;
                # endif

            x_coeff = (double *)
                hwloc_alloc_membind( Global::machine_topology, bytes_mallocd, Global::thread_masks[tid],
                                     HWLOC_MEMBIND_BIND, HWLOC_MEMBIND_THREAD );
        # else

            x_coeff = (double *) aligned_alloc( ALIGNED_ALLOC_ALIGNMENT, bytes_mallocd );

        # endif // if defined (ENABLE_HWLOC) && defined (USE_HWLOC_ALLOC)
        }

        const int64_t (& nx) [SPACE_DIMS] = y.nx();

        # pragma omp for collapse(2) schedule(static)
        for ( int64_t i = 1; i <= nx[0]; ++i ) {
        for ( int64_t j = 1; j <= nx[1]; ++j ) {
        for ( int64_t quad = 0; quad < 4; ++quad ) {
        for ( int64_t q = y.Quadrants(quad); q < y.Quadrants(quad + 1); q += LMV_SIMD_LEN ) {

            const int64_t len = std::min( ((int64_t)(LMV_SIMD_LEN)), y.Quadrants(quad + 1) - q );

            // --- Full blocks. ----------------------------------------------------------------------- //

            if ( len == LMV_SIMD_LEN ) {

                // Determine sign on angles in block.
                bool all_xi_negative = true;
                bool all_eta_negative = true;

                for ( int64_t l = 0; l < LMV_SIMD_LEN; ++l ) {

                    all_xi_negative  &= x.xi(q+l)  <= 0.0;
                    all_eta_negative &= x.eta(q+l) <= 0.0;
                }

                const int64_t xoffset = ( all_xi_negative  ? +1 : -1 );
                const int64_t yoffset = ( all_eta_negative ? +1 : -1 );

                // Load values into temporary storage.
                for ( int64_t d = 0; d <= DG_degree; ++d ) {
                for ( int64_t e = 0; e <= DG_degree; ++e ) {
                # pragma omp simd aligned( x_coeff : LMV_SIMD_LEN * sizeof(double) )
                for ( int64_t l = 0; l < LMV_SIMD_LEN; ++l ) {

                    IXCOEFF(0,d,e,l) = x(q+l,i,j,d,e);

                    IXCOEFF(1,d,e,l) = x(q+l, i + xoffset, j, d,e);
                    IXCOEFF(2,d,e,l) = x(q+l, i, j + yoffset, d,e);

                    IXCOEFF(3,d,e,l) = 0.0;
                }}}

                // Set value from absorption term.
                if ( sigma.DG_degree == 0 ) {

                    const double sigma_val = sigma(i,j,0,0);

                    for ( int64_t d = 0; d <= DG_degree; ++d ) {
                    for ( int64_t e = 0; e <= DG_degree; ++e ) {
                    # pragma omp simd aligned( x_coeff : LMV_SIMD_LEN * sizeof(double) )
                    for ( int64_t l = 0; l < LMV_SIMD_LEN; ++l ) {

                        IXCOEFF(3,d,e,l) += sigma_val * IXCOEFF(0,d,e,l);
                    }}}

                } else {

                    for ( int64_t d = 0; d <= DG_degree; ++d ) {
                    for ( int64_t e = 0; e <= DG_degree; ++e ) {
                    for ( int64_t u = 0; u <= DG_degree; ++u ) {
                    for ( int64_t v = 0; v <= DG_degree; ++v ) {

                        const double sigma_val = 0.25 * (2*d + 1)*(2*e + 1) * sigma(i,j,d,e,u,v);

                        # pragma omp simd aligned( x_coeff : LMV_SIMD_LEN * sizeof(double) )
                        for ( int64_t l = 0; l < LMV_SIMD_LEN; ++l )
                            IXCOEFF(3,d,e,l) += sigma_val * IXCOEFF(0,u,v,l);
                    }}}}
                }

                // Include components of advection term for ξ angular direction.
                # pragma omp simd aligned( x_coeff, coeff_v : LMV_SIMD_LEN * sizeof(double) )
                for ( int64_t l = 0; l < LMV_SIMD_LEN; ++l )
                    coeff_v[l] = x.xi(q+l) / x.dx(0);

                if ( all_xi_negative ) {

                    for ( int64_t d = 0; d <= DG_degree; ++d ) {
                    for ( int64_t a = 0; a <= DG_degree; ++a ) {
                    for ( int64_t e = 0; e <= DG_degree; ++e ) {
                    # pragma omp simd aligned( x_coeff, coeff_v : LMV_SIMD_LEN * sizeof(double) )
                    for ( int64_t l = 0; l < LMV_SIMD_LEN; ++l ) {

                        IXCOEFF(3,d,e,l) += coeff_v[l] * (2*d + 1) * ( neg1pow(a) * IXCOEFF(1,a,e,l)
                                                                       + neg1pow(d+a+1) * IXCOEFF(0,a,e,l) );
                    }}}}

                } else {

                    for ( int64_t d = 0; d <= DG_degree; ++d ) {
                    for ( int64_t a = 0; a <= DG_degree; ++a ) {
                    for ( int64_t e = 0; e <= DG_degree; ++e ) {
                    # pragma omp simd aligned( x_coeff, coeff_v : LMV_SIMD_LEN * sizeof(double) )
                    for ( int64_t l = 0; l < LMV_SIMD_LEN; ++l ) {

                        IXCOEFF(3,d,e,l) += coeff_v[l] * (2*d + 1) * ( IXCOEFF(0,a,e,l)
                                                                       + neg1pow(d+1) * IXCOEFF(1,a,e,l) );
                    }}}}
                }

                for ( int64_t d = 0; d <= DG_degree; ++d ) {
                for ( int64_t a = d-1; a >= 0; a -= 2 ) {
                for ( int64_t e = 0; e <= DG_degree; ++e ) {
                # pragma omp simd aligned( x_coeff, coeff_v : LMV_SIMD_LEN * sizeof(double) )
                for ( int64_t l = 0; l < LMV_SIMD_LEN; ++l ) {

                    IXCOEFF(3,d,e,l) -= ( 2 * (2*d + 1) ) * coeff_v[l] * IXCOEFF(0,a,e,l);
                }}}}

                // Include components of advection term for η angular direction.
                # pragma omp simd aligned( x_coeff, coeff_v : LMV_SIMD_LEN * sizeof(double) )
                for ( int64_t l = 0; l < LMV_SIMD_LEN; ++l )
                    coeff_v[l] = x.eta(q+l) / x.dx(1);

                if ( all_eta_negative ) {

                    for ( int64_t d = 0; d <= DG_degree; ++d ) {
                    for ( int64_t e = 0; e <= DG_degree; ++e ) {
                    for ( int64_t b = 0; b <= DG_degree; ++b ) {
                    # pragma omp simd aligned( x_coeff, coeff_v : LMV_SIMD_LEN * sizeof(double) )
                    for ( int64_t l = 0; l < LMV_SIMD_LEN; ++l ) {

                        IXCOEFF(3,d,e,l) += coeff_v[l] * (2*e + 1) * ( neg1pow(b) * IXCOEFF(2,d,b,l)
                                                                       + neg1pow(e+b+1) * IXCOEFF(0,d,b,l) );
                    }}}}

                } else {

                    for ( int64_t d = 0; d <= DG_degree; ++d ) {
                    for ( int64_t e = 0; e <= DG_degree; ++e ) {
                    for ( int64_t b = 0; b <= DG_degree; ++b ) {
                    # pragma omp simd aligned( x_coeff, coeff_v : LMV_SIMD_LEN * sizeof(double) )
                    for ( int64_t l = 0; l < LMV_SIMD_LEN; ++l ) {

                        IXCOEFF(3,d,e,l) += coeff_v[l] * (2*e + 1) * ( IXCOEFF(0,d,b,l)
                                                                       + neg1pow(e+1) * IXCOEFF(2,d,b,l) );
                    }}}}
                }

                for ( int64_t d = 0; d <= DG_degree; ++d ) {
                for ( int64_t e = 0; e <= DG_degree; ++e ) {
                for ( int64_t b = e-1; b >= 0; b -= 2 ) {
                # pragma omp simd aligned( x_coeff, coeff_v : LMV_SIMD_LEN * sizeof(double) )
                for ( int64_t l = 0; l < LMV_SIMD_LEN; ++l ) {

                    IXCOEFF(3,d,e,l) -= ( 2 * (2*e + 1) ) * coeff_v[l] * IXCOEFF(0,d,b,l);
                }}}}

                // Store output values.
                for ( int64_t d = 0; d <= DG_degree; ++d ) {
                for ( int64_t e = 0; e <= DG_degree; ++e ) {
                # pragma omp simd aligned( x_coeff : LMV_SIMD_LEN * sizeof(double) )
                for ( int64_t l = 0; l < LMV_SIMD_LEN; ++l ) {

                    y(q+l,i,j,d,e) = beta * y(q+l,i,j,d,e) + alpha * IXCOEFF(3,d,e,l);
                }}}

            // --- Partial blocks. -------------------------------------------------------------------- //

            } else {
            for ( int64_t l = 0; l < len; ++l ) {

                // Scale coefficients of output vector.
                for ( int64_t d = 0; d <= x.DG_degree; ++d ) {
                for ( int64_t e = 0; e <= x.DG_degree; ++e ) {

                    y(q+l,i,j,d,e) *= beta;
                }}

                // Set value from absorption term.
                if ( sigma.DG_degree == 0 ) {

                    for ( int64_t d = 0; d <= x.DG_degree; ++d ) {
                    for ( int64_t e = 0; e <= x.DG_degree; ++e ) {

                        y(q+l,i,j,d,e) += alpha * sigma(i,j,0,0) * x(q+l,i,j,d,e);
                    }}

                } else {

                    for ( int64_t d = 0; d <= x.DG_degree; ++d ) {
                    for ( int64_t e = 0; e <= x.DG_degree; ++e ) {
                    for ( int64_t u = 0; u <= x.DG_degree; ++u ) {
                    for ( int64_t v = 0; v <= x.DG_degree; ++v ) {

                        y(q+l,i,j,d,e) += 0.25 * alpha * (2*d + 1)*(2*e + 1) * sigma(i,j,d,e,u,v)
                                          * x(q+l,i,j,u,v);
                    }}}}
                }

                // Include components of advection term for ξ angular direction.
                coeff = alpha * x.xi(q+l) / x.dx(0);

                if ( x.xi(q+l) < 0 ) {

                    for ( int64_t d = 0; d <= x.DG_degree; ++d ) {
                    for ( int64_t a = 0; a <= x.DG_degree; ++a ) {
                    for ( int64_t e = 0; e <= x.DG_degree; ++e ) {

                        y(q+l,i,j,d,e) += coeff * (2*d + 1) * ( neg1pow(a) * x(q+l,i+1,j,a,e)
                                                                + neg1pow(d+a+1) * x(q+l,i,j,a,e) );
                    }}}

                } else {

                    for ( int64_t d = 0; d <= x.DG_degree; ++d ) {
                    for ( int64_t a = 0; a <= x.DG_degree; ++a ) {
                    for ( int64_t e = 0; e <= x.DG_degree; ++e ) {

                        y(q+l,i,j,d,e) += coeff * (2*d + 1) * ( x(q+l,i,j,a,e)
                                                                + neg1pow(d+1) * x(q+l,i-1,j,a,e) );
                    }}}
                }

                coeff *= -2.0;

                for ( int64_t d = 0; d <= x.DG_degree; ++d ) {
                for ( int64_t a = d-1; a >= 0; a -= 2 ) {
                for ( int64_t e = 0; e <= x.DG_degree; ++e ) {

                    y(q+l,i,j,d,e) += coeff * (2*d + 1) * x(q+l,i,j,a,e);
                }}}

                // Include components of advection term for η angular direction.
                coeff = alpha * x.eta(q+l) / x.dx(1);

                if ( x.eta(q+l) < 0 ) {

                    for ( int64_t d = 0; d <= x.DG_degree; ++d ) {
                    for ( int64_t e = 0; e <= x.DG_degree; ++e ) {
                    for ( int64_t b = 0; b <= x.DG_degree; ++b ) {

                        y(q+l,i,j,d,e) += coeff * (2*e + 1) * ( neg1pow(b) * x(q+l,i,j+1,d,b)
                                                                + neg1pow(e+b+1) * x(q+l,i,j,d,b) );
                    }}}

                } else {

                    for ( int64_t d = 0; d <= x.DG_degree; ++d ) {
                    for ( int64_t e = 0; e <= x.DG_degree; ++e ) {
                    for ( int64_t b = 0; b <= x.DG_degree; ++b ) {

                        y(q+l,i,j,d,e) += coeff * (2*e + 1) * ( x(q+l,i,j,d,b)
                                                                + neg1pow(e+1) * x(q+l,i,j-1,d,b) );
                    }}}
                }

                coeff *= -2.0;

                for ( int64_t d = 0; d <= x.DG_degree; ++d ) {
                for ( int64_t e = 0; e <= x.DG_degree; ++e ) {
                for ( int64_t b = e-1; b >= 0; b -= 2 ) {

                    y(q+l,i,j,d,e) += coeff * (2*e + 1) * x(q+l,i,j,d,b);
                }}}
            }}
        }}}}

        if ( use_temp_mem ) {
        # if defined (ENABLE_HWLOC) && defined (USE_HWLOC_ALLOC)
            hwloc_free( Global::machine_topology, x_coeff, bytes_mallocd );
        # else
            std::free( x_coeff );
        # endif
        }
    }

# undef IXCOEFF

# else // if defined (ENABLE_SIMD_BLOCKING)

    double coeff;

    const int64_t (& nx) [SPACE_DIMS] = y.nx();

    # pragma omp parallel for private(coeff) collapse(2) schedule(static)
    for ( int64_t i = 1; i <= nx[0]; ++i ) {
    for ( int64_t j = 1; j <= nx[1]; ++j ) {
    for ( int64_t q = 0; q <  y.nq();  ++q ) {

        // Scale coefficients of output vector.
        for ( int64_t d = 0; d <= x.DG_degree; ++d ) {
        for ( int64_t e = 0; e <= x.DG_degree; ++e ) {

            y(q,i,j,d,e) *= beta;
        }}

        // Set value from absorption term.
        if ( sigma.DG_degree == 0 ) {

            for ( int64_t d = 0; d <= x.DG_degree; ++d ) {
            for ( int64_t e = 0; e <= x.DG_degree; ++e ) {

                y(q,i,j,d,e) += alpha * sigma(i,j,0,0) * x(q,i,j,d,e);
            }}

        } else {

            for ( int64_t d = 0; d <= x.DG_degree; ++d ) {
            for ( int64_t e = 0; e <= x.DG_degree; ++e ) {
            for ( int64_t u = 0; u <= x.DG_degree; ++u ) {
            for ( int64_t v = 0; v <= x.DG_degree; ++v ) {

                y(q,i,j,d,e) += alpha * (2*d + 1)*(2*e + 1) * sigma(i,j,d,e,u,v) * x(q,i,j,u,v) / 4.0;
            }}}}
        }

        // Include components of advection term for ξ angular direction.
        coeff = alpha * x.xi(q) / x.dx(0);

        if ( x.xi(q) < 0 ) {

            for ( int64_t d = 0; d <= x.DG_degree; ++d ) {
            for ( int64_t a = 0; a <= x.DG_degree; ++a ) {
            for ( int64_t e = 0; e <= x.DG_degree; ++e ) {

                y(q,i,j,d,e) += coeff * (2*d + 1) * ( neg1pow(a) * x(q,i+1,j,a,e)
                                                      + neg1pow(d+a+1) * x(q,i,j,a,e) );
            }}}

        } else {

            for ( int64_t d = 0; d <= x.DG_degree; ++d ) {
            for ( int64_t a = 0; a <= x.DG_degree; ++a ) {
            for ( int64_t e = 0; e <= x.DG_degree; ++e ) {

                y(q,i,j,d,e) += coeff * (2*d + 1) * ( x(q,i,j,a,e) + neg1pow(d+1) * x(q,i-1,j,a,e) );
            }}}
        }

        coeff *= -2.0;

        for ( int64_t d = 0; d <= x.DG_degree; ++d ) {
        for ( int64_t a = d-1; a >= 0; a -= 2 ) {
        for ( int64_t e = 0; e <= x.DG_degree; ++e ) {

            y(q,i,j,d,e) += coeff * (2*d + 1) * x(q,i,j,a,e);
        }}}

        // Include components of advection term for η angular direction.
        coeff = alpha * x.eta(q) / x.dx(1);

        if ( x.eta(q) < 0 ) {

            for ( int64_t d = 0; d <= x.DG_degree; ++d ) {
            for ( int64_t e = 0; e <= x.DG_degree; ++e ) {
            for ( int64_t b = 0; b <= x.DG_degree; ++b ) {

                y(q,i,j,d,e) += coeff * (2*e + 1) * ( neg1pow(b) * x(q,i,j+1,d,b)
                                                      + neg1pow(e+b+1) * x(q,i,j,d,b) );
            }}}

        } else {

            for ( int64_t d = 0; d <= x.DG_degree; ++d ) {
            for ( int64_t e = 0; e <= x.DG_degree; ++e ) {
            for ( int64_t b = 0; b <= x.DG_degree; ++b ) {

                y(q,i,j,d,e) += coeff * (2*e + 1) * ( x(q,i,j,d,b) + neg1pow(e+1) * x(q,i,j-1,d,b) );
            }}}
        }

        coeff *= -2.0;

        for ( int64_t d = 0; d <= x.DG_degree; ++d ) {
        for ( int64_t e = 0; e <= x.DG_degree; ++e ) {
        for ( int64_t b = e-1; b >= 0; b -= 2 ) {

            y(q,i,j,d,e) += coeff * (2*e + 1) * x(q,i,j,d,b);
        }}}
    }}}

# endif // if defined (ENABLE_SIMD_BLOCKING)

    if (/*
         * Set halos dirty if:
         *
         *  1.  Boundaries ARE periodic; OR
         *  2.  There is at least one reflecting boundary condition; OR
         *  3.  There is more than one MPI rank.
         */
            Global::periodic
         || BitmaskHasAny( Global::reflecting_boundaries, BoundaryEdge::All )
    # if defined (ENABLE_MPI)
         || Global::MPI_num_ranks > 1
    # endif
    ) {
        y.MarkAllHalosDirty();
    }

    Global::TMR_AF_Lmv.Stop();
}


# endif // if SPACE_DIMS == 2


} // namespace TransportOperator
