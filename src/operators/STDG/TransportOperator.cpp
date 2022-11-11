//------------------------------------------------------------------------------------------------------------
//! \file   operators/STDG/TransportOperator.cpp
//! \brief  Implementation file for RKDG transport operator routines.
//!
//! \author Michael M. Crockatt
//! \date   October 2017
//------------------------------------------------------------------------------------------------------------

# include <algorithm>

# include "operators/STDG/TransportOperator.hpp"
# include "utils/CLog.hpp"


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
    const STDG::OrdinateFlux & x,
    STDG::ScalarFlux & y,
    const OpDomain op_domain        // = OpDomain::All
) {

    PRINT_STATUS( "Applying STDG transport operator P.\n" )

# if defined (STRICT_CHECK)

    if ( !BitmaskHasAll( op_domain, OpDomain::Interior ) ) {

        std::string error_message = "OpDomain passed to '" + std::string(__func__)
                                    + "' does not contain flag 'OpDomain::Interior'\n";

        PRINT_ERROR( error_message.c_str() )
        throw std::out_of_range( error_message );
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

        for ( int64_t d = 0; d <= x.DG_degree_x; ++d ) {
        for ( int64_t s = 0; s <= x.DG_degree_t; ++s ) {

            y(i,d,s) *= beta;
        }}

        for ( int64_t q = 0; q <  x.nq();        ++q ) {
        for ( int64_t d = 0; d <= x.DG_degree_x; ++d ) {
        for ( int64_t s = 0; s <= x.DG_degree_t; ++s ) {

            y(i,d,s) += alpha * x.w(q) * x(q,i,d,s);
        }}}
    }

# elif SPACE_DIMS == 2

    # pragma omp parallel for collapse(2) schedule(static)
    for ( int64_t i = 1 - offset[0][0]; i <= nx[0] + offset[0][1]; ++i ) {
    for ( int64_t j = 1 - offset[1][0]; j <= nx[1] + offset[1][1]; ++j ) {

        for ( int64_t d = 0; d <= x.DG_degree_x; ++d ) {
        for ( int64_t e = 0; e <= x.DG_degree_x; ++e ) {
        for ( int64_t s = 0; s <= x.DG_degree_t; ++s ) {

            y(i,j,d,e,s) *= beta;
        }}}

        for ( int64_t q = 0; q <  x.nq();        ++q ) {
        for ( int64_t d = 0; d <= x.DG_degree_x; ++d ) {
        for ( int64_t e = 0; e <= x.DG_degree_x; ++e ) {
        for ( int64_t s = 0; s <= x.DG_degree_t; ++s ) {

            y(i,j,d,e,s) += alpha * x.w(q) * x(q,i,j,d,e,s);
        }}}}
    }}

# elif SPACE_DIMS == 3

    # pragma omp parallel for collapse(3) schedule(static)
    for ( int64_t i = 1 - offset[0][0]; i <= nx[0] + offset[0][1]; ++i ) {
    for ( int64_t j = 1 - offset[1][0]; j <= nx[1] + offset[1][1]; ++j ) {
    for ( int64_t k = 1 - offset[2][0]; k <= nx[2] + offset[2][1]; ++k ) {

        for ( int64_t d = 0; d <= x.DG_degree_x; ++d ) {
        for ( int64_t e = 0; e <= x.DG_degree_x; ++e ) {
        for ( int64_t f = 0; f <= x.DG_degree_x; ++f ) {
        for ( int64_t s = 0; s <= x.DG_degree_t; ++s ) {

            y(i,j,k,d,e,f,s) *= beta;
        }}}}

        for ( int64_t q = 0; q <  x.nq();        ++q ) {
        for ( int64_t d = 0; d <= x.DG_degree_x; ++d ) {
        for ( int64_t e = 0; e <= x.DG_degree_x; ++e ) {
        for ( int64_t f = 0; f <= x.DG_degree_x; ++f ) {
        for ( int64_t s = 0; s <= x.DG_degree_t; ++s ) {

            y(i,j,k,d,e,f,s) += alpha * x->w(q) * x(q,i,j,k,d,e,f,s);
        }}}}}
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
    const STDG::ScalarFlux & x,
    STDG::OrdinateFlux & y,
    const OpDomain op_domain            // = OpDomain::Full
) {

    PRINT_STATUS( "Applying STDG transport operator S.\n" )

# if defined (STRICT_CHECK)

    if ( !BitmaskHasAll( op_domain, OpDomain::Interior ) ) {

        std::string error_message = "OpDomain passed to '" + std::string(__func__)
                                    + "' does not contain flag 'OpDomain::Interior'\n";

        PRINT_ERROR( error_message.c_str() )
        throw std::out_of_range( error_message );
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

    const int64_t (& nx) [SPACE_DIMS] = y.nx();

# if SPACE_DIMS == 1

    # pragma omp parallel for schedule(static)
    for ( int64_t i = 1 - offset[0][0]; i <= nx[0] + offset[0][1]; ++i ) {
    for ( int64_t d = 0; d <= x.DG_degree_x; ++d ) {
    for ( int64_t s = 0; s <= x.DG_degree_t; ++s ) {

        double psi_val = 0.0;

        if ( sigma.DG_degree == 0 ) {

            psi_val = sigma(i,0) * x(i,d,s) / 2.0;

        } else {

            for ( int64_t u = 0; u <= x.DG_degree_x; ++u )
                psi_val += sigma(i,d,u) * x(i,u,s);

            psi_val *= (2*d + 1) / 4.0;
        }

        psi_val *= alpha;

        for ( int64_t q = 0; q < y.nq(); ++q )
            y(q,i,d,s) = beta * y(q,i,d,s) + psi_val;
    }}}

# elif SPACE_DIMS == 2

    # pragma omp parallel for collapse(2) schedule(static)
    for ( int64_t i = 1 - offset[0][0]; i <= nx[0] + offset[0][1]; ++i ) {
    for ( int64_t j = 1 - offset[1][0]; j <= nx[1] + offset[1][1]; ++j ) {
    for ( int64_t d = 0; d <= x.DG_degree_x; ++d ) {
    for ( int64_t e = 0; e <= x.DG_degree_x; ++e ) {
    for ( int64_t s = 0; s <= x.DG_degree_t; ++s ) {

        double psi_val = 0.0;

        if ( sigma.DG_degree == 0 ) {

            psi_val = sigma(i,j,0,0) * x(i,j,d,e,s) / (4.0 * M_PI);

        } else {

            for ( int64_t u = 0; u <= x.DG_degree_x; ++u ) {
            for ( int64_t v = 0; v <= x.DG_degree_x; ++v ) {

                psi_val += sigma(i,j,d,e,u,v) * x(i,j,u,v,s);
            }}

            psi_val *= (2*d + 1)*(2*e + 1) / (16.0 * M_PI);
        }

        psi_val *= alpha;

        for ( int64_t q = 0; q < y.nq(); ++q )
            y(q,i,j,d,e,s) = beta * y(q,i,j,d,e,s) + psi_val;
    }}}}}

# elif SPACE_DIMS == 3

    # pragma omp parallel for collapse(3) schedule(static)
    for ( int64_t i = 1 - offset[0][0]; i <= nx[0] + offset[0][1]; ++i ) {
    for ( int64_t j = 1 - offset[1][0]; j <= nx[1] + offset[1][1]; ++j ) {
    for ( int64_t k = 1 - offset[2][0]; k <= nx[2] + offset[2][1]; ++k ) {
    for ( int64_t d = 0; d <= x.DG_degree_x; ++d ) {
    for ( int64_t e = 0; e <= x.DG_degree_x; ++e ) {
    for ( int64_t f = 0; f <= x.DG_degree_x; ++f ) {
    for ( int64_t s = 0; s <= x.DG_degree_t; ++s ) {

        double psi_val = 0.0;

        if ( sigma.DG_degree == 0 {

            psi_val = sigma(i,j,k,0,0,0) * x(i,j,k,d,e,f,s) / (4.0 * M_PI);

        } else {

            for ( int64_t u = 0; u <= x.DG_degree_x; ++u ) {
            for ( int64_t v = 0; v <= x.DG_degree_x; ++v ) {
            for ( int64_t w = 0; w <= x.DG_degree_x; ++w ) {

                psi_val += sigma(i,j,k,d,e,f,u,v,w) * x(i,j,k,u,v,w,s);
            }}}

            psi_val *= (2*d + 1)*(2*e + 1)*(2*f + 1) / (32.0 * M_PI);
        }

        psi_val *= alpha;

        for ( int64_t q = 0; q < y.nq(); ++q )
            y(q,i,j,k,d,e,f,s) = beta * y(q,i,j,k,d,e,f,s) + psi_val;
    }}}}}}}

# endif // if SPACE_DIMS == ?

    Global::TMR_AF_Smv.Stop();
}


} // namespace TransportOperator
