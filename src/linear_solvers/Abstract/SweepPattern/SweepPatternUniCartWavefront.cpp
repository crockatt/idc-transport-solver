//------------------------------------------------------------------------------------------------------------
//! \file   linear_solvers/Abstract/SweepPattern/SweepPatternUniCartWavefront.cpp
//! \brief  Contains implementations and instantiations of methods from SweepPatternUniCartWavefront class.
//!
//! \author Michael M. Crockatt
//! \date   March 2018
//------------------------------------------------------------------------------------------------------------

# if SPACE_DIMS == 2


# include <typeinfo>

# if defined (ENABLE_HWLOC) || defined (DOXYCOMPILE)
    # include <hwloc.h>
# endif

# if defined (_OPENMP)
    # include <omp.h>
# endif

# include "linear_solvers/Abstract/SweepOperator.hpp"
# include "linear_solvers/Abstract/SweepPattern/SweepPatternUniCartWavefront.hpp"
# include "linear_solvers/Abstract/SweepCellSolve/SweepCellSolveManager.hpp"
# include "utils/SIMD.hpp"


//============================================================================================================
//=== MEMBER DEFINITIONS FOR SweepPatternUniCartWavefront CLASS ==============================================
//============================================================================================================


//--------------------------------------------------------------------------------------------------------
//! \brief  Returns the string descriptor of the class type.
//--------------------------------------------------------------------------------------------------------
template< class OrdinateFlux >
const std::string SweepOperator<OrdinateFlux>::SweepPatternUniCartWavefront::Descriptor( void ) {

    return "wavefront";
}


//--------------------------------------------------------------------------------------------------------
//! \brief  Returns the string descriptor for an objects type.
//--------------------------------------------------------------------------------------------------------
template< class OrdinateFlux >
const std::string SweepOperator<OrdinateFlux>::SweepPatternUniCartWavefront::GetDescriptor( void ) const {

    return Descriptor();
}


//------------------------------------------------------------------------------------------------------------
//! \brief  Virtual destructor.
//------------------------------------------------------------------------------------------------------------
template< class OrdinateFlux >
SweepOperator<OrdinateFlux>::SweepPatternUniCartWavefront::~SweepPatternUniCartWavefront( void ) {

    PRINT_STATUS( "Executing SweepOperator<%s>::SweepPatternUniCartWavefront::%s.\n",
                  typeid(OrdinateFlux).name(), __func__ )

    /* empty */
}


//------------------------------------------------------------------------------------------------------------
//! \brief  Prints a summary of the object's configuration to the logging interface.
//!
//! \param[in]  prefix  String prepended to each line of output, intended to be used for controlling
//!                     indentation of output.
//!
//! \todo   Update Print implementation to use map for string descriptors.
//------------------------------------------------------------------------------------------------------------
template< class OrdinateFlux >
void SweepOperator<OrdinateFlux>::SweepPatternUniCartWavefront::Print (

    const std::string & prefix // = "  "

) const {

    PRINT_STATUS( "Executing SweepOperator<%s>::SweepPatternUniCartWavefront::%s.\n",
                  typeid(OrdinateFlux).name(), __func__ )

    PRINT_LOG( "%s%-*s  %s\n", prefix.c_str(),
               Global::col_width, "Sweep Pattern:", this->Descriptor().c_str() )

    this->SweepPatternUniCart::Print( prefix );
}


//------------------------------------------------------------------------------------------------------------
//! \brief  Sets parameters given in a ParameterList.
//------------------------------------------------------------------------------------------------------------
template< class OrdinateFlux >
void SweepOperator<OrdinateFlux>::SweepPatternUniCartWavefront::SetParameters (

    const ParameterList & input_list
) {

    PRINT_STATUS( "Executing SweepOperator<%s>::SweepPatternUniCartWavefront::%s.\n",
                  typeid(OrdinateFlux).name(), __func__ )

    this->SweepPatternUniCart::SetParameters( input_list );
}


//------------------------------------------------------------------------------------------------------------
//! \brief  Performs a transport sweep.
//!
//! \param[in]      source      Contains the coefficients of the source term \f$ Q \f$ for the transport
//!                             system.
//! \param[out]     result      Upon return, contains the coefficients \f$ \Psi \f$ that satisfy the transport
//!                             equation.
//! \param[in]      sigma       The total cross section \f$ \sigma_{\mathrm{t}} \f$.
//! \param[in]      initial     (optional) <br>
//!                             Pointer to object containing the initial condition for the timestep.
//!                             If this value is null (default) then the initial condition is assumed to be
//!                             zero and this pointer is not dereferenced.
//------------------------------------------------------------------------------------------------------------
template< class OrdinateFlux >
void SweepOperator<OrdinateFlux>::SweepPatternUniCartWavefront::Sweep (

    const OrdinateFlux & source,
    OrdinateFlux & result,
    const RKDG::CrossSection & sigma,
    const RKDG::OrdinateFlux * const initial // = nullptr

) const {

    PRINT_STATUS( "Executing SweepOperator<%s>::SweepPatternUniCartWavefront::%s.\n",
                  typeid(OrdinateFlux).name(), __func__ )

    // Reference to SweepCellSolve object to call for local solves.
    const SweepCellSolveManager & cell_solve = *this->sw_op.sweep_cell_solve;

    // Number of diagonal wavefronts in spatial mesh.
    const int64_t num_diags = this->sw_op.nx(0) + this->sw_op.nx(1) - 1;

    // Number of bytes of temporary memory to allocate per thread.
    const size_t bytes_allocd = SWEEP_SIMD_LEN * sizeof(double) * (   cell_solve.GetDimofA()
                                                                    + cell_solve.GetDimofB()
                                                                    + cell_solve.GetDimofWork() );
    // Construct list of angular blocks to loop over.
    auto angle_blocks = this->ConstructSIMDAngleBlocks();

    # pragma omp parallel
    {
        // Allocate temporary memory for linear solves.
    # if defined (ENABLE_HWLOC) && defined (USE_HWLOC_ALLOC)

        int tid
            # if defined (_OPENMP)
                = omp_get_thread_num();
            # else
                = 0;
            # endif

        double * const temporary_memory = (double *)
            hwloc_alloc_membind( Global::machine_topology, bytes_allocd,
                                 Global::thread_masks[tid], HWLOC_MEMBIND_BIND, HWLOC_MEMBIND_THREAD );

    # else // if defined (ENABLE_HWLOC) && defined (USE_HWLOC_ALLOC)

        double * const temporary_memory = (double *)
            # if defined (USE_ALIGNED_ALLOC)
                aligned_alloc( ALIGNED_ALLOC_ALIGNMENT, bytes_allocd );
            # else
                std::malloc( bytes_allocd );
            # endif

    # endif // if defined (ENABLE_HWLOC) && defined (USE_HWLOC_ALLOC)

        double * const A    = temporary_memory;
        double * const B    = A + SWEEP_SIMD_LEN * cell_solve.GetDimofA();
        void   * const WORK = B + SWEEP_SIMD_LEN * cell_solve.GetDimofB();

    # if SPACE_DIMS == 2

        /*
         *  Spatial dimensions: 2
         */

        // Used to specify the signs associated with the coordinates of the ordinates in a given quadrant.
        bool signs[SPACE_DIMS] = { 0, 0 };

        if ( this->sweep_octants_in_sequence ) {
        /*
         *  Sweep quadrants sequentially, reflecting boundary values between quadrants.
         */
            for ( int64_t quad_index : { 0, 1, 2, 3 } ) {

                // Determine quadrant to sweep based on quadrant ordering.
                const int64_t quad = this->sweep_order[quad_index];

                // Determine signs on ordinates in current quadrant.
                Quadrule::OrdinateSet::SetSigns( quad, signs );

                // Get list of SIMD blocks for current quadrant.
                const auto & quad_angle_blocks = angle_blocks.at( quad );

                // Number of SIMD blocks for current quadrant.
                const int64_t num_angle_blocks = quad_angle_blocks.size();

                // Sweep using diagonal wavefront.
                for ( int64_t diag = 1; diag <= num_diags; ++diag ) {

                    const int64_t diag_length = this->ComputeDiagLength( diag );

                    # pragma omp for collapse(2)
                    for ( int64_t cell = 1; cell <= diag_length; ++cell ) {
                    for ( int64_t q_blk = 0; q_blk < num_angle_blocks; ++q_blk ) {

                        // Copy block containing angular component and length of block from list.
                        SIMD_BlkIdx<0> idx = quad_angle_blocks[ q_blk ];

                        // Compute indices of spatial cell to update.
                        int64_t i = std::min( diag, this->sw_op.nx(0) ) - cell + 1;
                        int64_t j = diag - i + 1;

                        // Reverse direction of sweep for negative angles.
                        if ( signs[0] /* ξ < 0 */ ) {  i = this->sw_op.nx(0) - i + 1;  }
                        if ( signs[1] /* η < 0 */ ) {  j = this->sw_op.nx(1) - j + 1;  }

                        // Set spatial index in SIMD block object.
                        idx.i = i;
                        idx.j = j;

                        // Do solve.
                        cell_solve.ConstructAndSolve( A, B, idx, sigma, source, result, WORK, initial );
                    }}

                } // end for diag

                if ( BitmaskHasAny( Global::reflecting_boundaries, BoundaryEdge::All ) ) {

                    # pragma omp barrier

                    # pragma omp master
                    Global::TMR_AF_Lsv.Stop();

                    result.ReflectBoundaries();
                    # pragma omp barrier

                    # pragma omp master
                    Global::TMR_AF_Lsv.Start();
                }

                # pragma omp barrier

            } // end for quad_index.

        } else {
        /*
         *  Sweep all quadrants simultaneously.
         */

            // Assemble list of SIMD blocks across all quadrants.
            int64_t num_total_blocks = 0;
            for ( auto & v : angle_blocks ) {  num_total_blocks += v.size();  }

            std::vector< SIMD_BlkIdx<0> > all_angle_blocks;
            all_angle_blocks.reserve( num_total_blocks );

            for ( auto & v : angle_blocks )
                all_angle_blocks.insert( all_angle_blocks.end(), v.begin(), v.end() );

            // Sweep using diagonal wavefront.
            for ( int64_t diag = 1; diag <= num_diags; ++diag ) {

                const int64_t diag_length = this->ComputeDiagLength( diag );

                # pragma omp for collapse(2)
                for ( int64_t cell = 1; cell <= diag_length; ++cell ) {
                for ( int64_t q_blk = 0; q_blk < num_total_blocks; ++q_blk ) {

                    // Copy block containing angular component and length of block from list.
                    SIMD_BlkIdx<0> idx = all_angle_blocks[ q_blk ];

                    // Determine signs on ordinates in block.
                    Quadrule::OrdinateSet::SetSigns( idx.quadrant, signs );

                    // Compute indices of spatial cell to update.
                    int64_t i = std::min( diag, this->sw_op.nx(0) ) - cell + 1;
                    int64_t j = diag - i + 1;

                    // Reverse direction of sweep for negative angles.
                    if ( signs[0] /* ξ < 0 */ ) {  i = this->sw_op.nx(0) - i + 1;  }
                    if ( signs[1] /* η < 0 */ ) {  j = this->sw_op.nx(1) - j + 1;  }

                    // Set spatial index in SIMD block object.
                    idx.i = i;
                    idx.j = j;

                    // Do solve.
                    cell_solve.ConstructAndSolve( A, B, idx, sigma, source, result, WORK, initial );
                }}

            } // end for diag

            /*
             *  NOTE:  This verison of the sweep is only used when there are no reflecting boundary conditions
             *          along the edges of the current MPI block. In this case, OrdinateFlux::ReflectBoundaries
             *          does nothing, and we therefore don't need to bother calling it here.
             */
        }

    # elif SPACE_DIMS == 3

        /*
         *  Spatial dimensions: 3
         */

        # warning "Implementation of SweepOperator<OrdinateFlux>::SweepPatternUniCartWavefront::Sweep incomplete."

    # endif // if SPACE_DIMS == ?

    # if defined (ENABLE_HWLOC) && defined (USE_HWLOC_ALLOC)
        hwloc_free( Global::machine_topology, temporary_memory, bytes_allocd );
    # else
        std::free( temporary_memory );
    # endif
    } // end omp parallel

    if (    Global::periodic
    # if defined (ENABLE_MPI)
         || Global::MPI_num_ranks > 1
    # endif
    ) {
        result.halo_cells_dirty = true;
        result.upwind_halo_cells_dirty = true;
    }
}


//============================================================================================================
//=== EXPLICIT INSTANTIATIONS ================================================================================
//============================================================================================================

# include "linear_solvers/Abstract/SweepPattern/SweepPatternDerivedFactory.hpp"
# include "objects/RKDG/OrdinateFlux.hpp"
# include "objects/STDG/OrdinateFlux.hpp"


//!
//! \brief  Instantiation of SweepOperator<RKDG::OrdinateFlux>::SweepPatternUniCartWavefront class.
//!
template class SweepOperator<RKDG::OrdinateFlux>::SweepPatternUniCartWavefront;

//!
//! \brief  Instantiation of factory for SweepOperator<RKDG::OrdinateFlux>::SweepPatternUniCartWavefront.
//!
template class
SweepOperator<RKDG::OrdinateFlux>::SweepPattern::SweepPatternDerivedFactory<
        SweepOperator<RKDG::OrdinateFlux>::SweepPatternUniCartWavefront
    >;


//!
//! \brief  Instantiation of SweepOperator<STDG::OrdinateFlux>::SweepPatternUniCartWavefront class.
//!
template class SweepOperator<STDG::OrdinateFlux>::SweepPatternUniCartWavefront;

//!
//! \brief  Instantiation of factory for SweepOperator<STDG::OrdinateFlux>::SweepPatternUniCartWavefront.
//!
template class
SweepOperator<STDG::OrdinateFlux>::SweepPattern::SweepPatternDerivedFactory<
        SweepOperator<STDG::OrdinateFlux>::SweepPatternUniCartWavefront
    >;


# endif // if SPACE_DIMS == 2
