//------------------------------------------------------------------------------------------------------------
//! \file   linear_solvers/Abstract/SweepPattern/SweepPatternUniCartWavefrontMPI.cpp
//! \brief  Contains implementations and instantiations of methods from SweepPatternUniCartWavefrontMPI class.
//!
//! \author Michael M. Crockatt
//! \date   March 2018
//------------------------------------------------------------------------------------------------------------

# if defined (ENABLE_MPI)
# if SPACE_DIMS == 2


# include <typeinfo>
# include <time.h>

# if defined (ENABLE_HWLOC) || defined (DOXYCOMPILE)
    # include <hwloc.h>
# endif

# if defined (_OPENMP)
    # include <omp.h>
# endif

# include "linear_solvers/Abstract/SweepOperator.hpp"
# include "linear_solvers/Abstract/SweepPattern/SweepPatternUniCartWavefrontMPI.hpp"
# include "linear_solvers/Abstract/SweepCellSolve/SweepCellSolveManager.hpp"
# include "utils/SIMD.hpp"


//============================================================================================================
//=== MEMBER DEFINITIONS FOR SweepPatternUniCartWavefrontMPI CLASS ===========================================
//============================================================================================================


//--------------------------------------------------------------------------------------------------------
//! \brief  Returns the string descriptor of the class type.
//--------------------------------------------------------------------------------------------------------
template< class OrdinateFlux >
const std::string SweepOperator<OrdinateFlux>::SweepPatternUniCartWavefrontMPI::Descriptor( void ) {

    return "wavefront-mpi";
}


//--------------------------------------------------------------------------------------------------------
//! \brief  Returns the string descriptor for an objects type.
//--------------------------------------------------------------------------------------------------------
template< class OrdinateFlux >
const std::string SweepOperator<OrdinateFlux>::SweepPatternUniCartWavefrontMPI::GetDescriptor( void ) const {

    return Descriptor();
}


//------------------------------------------------------------------------------------------------------------
//! \brief  Virtual destructor.
//------------------------------------------------------------------------------------------------------------
template< class OrdinateFlux >
SweepOperator<OrdinateFlux>::SweepPatternUniCartWavefrontMPI::~SweepPatternUniCartWavefrontMPI( void ) {

    PRINT_STATUS( "Executing SweepOperator<%s>::SweepPatternUniCartWavefrontMPI::%s.\n",
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
void SweepOperator<OrdinateFlux>::SweepPatternUniCartWavefrontMPI::Print (

    const std::string & prefix // = "  "

) const {

    PRINT_STATUS( "Executing SweepOperator<%s>::SweepPatternUniCartWavefrontMPI::%s.\n",
                  typeid(OrdinateFlux).name(), __func__ )

    PRINT_LOG( "%s%-*s  %s\n", prefix.c_str(),
               Global::col_width, "Sweep Pattern:",  this->Descriptor().c_str() )

    this->SweepPatternUniCart::Print( prefix );
}


//------------------------------------------------------------------------------------------------------------
//! \brief  Sets parameters given in a ParameterList.
//------------------------------------------------------------------------------------------------------------
template< class OrdinateFlux >
void SweepOperator<OrdinateFlux>::SweepPatternUniCartWavefrontMPI::SetParameters (

    const ParameterList & input_list
) {

    PRINT_STATUS( "Executing SweepOperator<%s>::SweepPatternUniCartWavefrontMPI::%s.\n",
                  typeid(OrdinateFlux).name(), __func__ )

    this->SweepPatternUniCart::SetParameters( input_list );
}


//------------------------------------------------------------------------------------------------------------
//! \brief  Performs a transport sweep.
//!
//! Implements a very basic sequential (non-pipelined) distributed memory sweep algorithm.
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
void SweepOperator<OrdinateFlux>::SweepPatternUniCartWavefrontMPI::Sweep (

    const OrdinateFlux & source,
    OrdinateFlux & result,
    const RKDG::CrossSection & sigma,
    const RKDG::OrdinateFlux * const initial // = nullptr

) const {

    PRINT_STATUS( "Executing SweepOperator<%s>::SweepPatternUniCartWavefrontMPI::%s.\n",
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
        // Allocate temporary memory for solving linear systems.
    # if defined (ENABLE_HWLOC) && defined (USE_HWLOC_ALLOC)

        hwloc_cpuset_t thread_mask = hwloc_bitmap_alloc();
        hwloc_get_cpubind( Global::machine_topology, thread_mask, HWLOC_CPUBIND_THREAD );

        double * const temporary_memory = (double *)
            hwloc_alloc_membind( Global::machine_topology, bytes_allocd, thread_mask, HWLOC_MEMBIND_BIND,
                                    HWLOC_MEMBIND_THREAD );

        hwloc_bitmap_free( thread_mask );

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

        // Used to specify the signs associated with the coordinates of the ordinates in a given quadrant.
        bool signs[SPACE_DIMS] = { 0, 0 };

        for ( int64_t quad_index : { 0, 1, 2, 3 } ) {

            // Determine quadrant to sweep based on quadrant ordering.
            const int64_t quad = this->sweep_order[quad_index];

            // Determine signs on ordinates in current quadrant.
            Quadrule::OrdinateSet::SetSigns( quad, signs );

            // Get list of SIMD blocks for current quadrant.
            const auto & quad_angle_blocks = angle_blocks.at( quad );

            // Number of SIMD blocks for current quadrant.
            const int64_t num_angle_blocks = quad_angle_blocks.size();

            // Construct resized MPI_Datatype objects in each dimension for the current quadrant.
            MPI_Datatype cell_X, cell_Y;

            # pragma omp master
            {
                MPI_Datatype cell_contig;

                const int64_t count = result.SpatialDOFPerCell()
                                      * ( this->sw_op.Quadrants(quad + 1) - this->sw_op.Quadrants(quad) );

                const int64_t stride_X = result.CellStride(1);
                const int64_t stride_Y = result.CellStride(0);

                MPI_Type_contiguous( count, MPI_DOUBLE, &cell_contig );
                MPI_Type_commit( &cell_contig );

                MPI_Type_create_resized( cell_contig, 0, stride_X * sizeof(double), &cell_X );
                MPI_Type_create_resized( cell_contig, 0, stride_Y * sizeof(double), &cell_Y );

                MPI_Type_commit( &cell_X );
                MPI_Type_commit( &cell_Y );

                MPI_Type_free( &cell_contig );
            }

            MPI_Request comm_reqs [] = { MPI_REQUEST_NULL, MPI_REQUEST_NULL };

            // Compute shifts for quadrant.
            int rank_source_X, rank_dest_X, direction_X = 0, disp_X = 1,
                rank_source_Y, rank_dest_Y, direction_Y = 1, disp_Y = 1;

            if ( signs[0] /* ξ < 0 */ ) {  disp_X = -1;  }
            if ( signs[1] /* η < 0 */ ) {  disp_Y = -1;  }

            # pragma omp master
            {
                MPI_Cart_shift( Global::MPI_cart_comm, direction_X, disp_X, &rank_source_X, &rank_dest_X );
                MPI_Cart_shift( Global::MPI_cart_comm, direction_Y, disp_Y, &rank_source_Y, &rank_dest_Y );
            }

            // Wait for incoming data.
            # pragma omp master
            if ( rank_source_X != MPI_PROC_NULL ) {

                // Coordinates of halo cells to receive data into for current quadrant.
                const int64_t ii = ( signs[0] /* ξ < 0 */ ? this->sw_op.nx(0) + 1 : 0 );
                const int64_t jj = 1;

                void * buff_X  = result.PointerAtQuadrantCell( quad, ii,jj );
                int    count_X = this->sw_op.nx(1);
                int    tag_X   = quad;

                MPI_Irecv( buff_X, count_X, cell_X, rank_source_X, tag_X, Global::MPI_cart_comm,
                           &comm_reqs[0] );
            }

            # pragma omp master
            if ( rank_source_Y != MPI_PROC_NULL ) {

                // Coordinates of halo cells to receive data into for current quadrant.
                const int64_t ii = 1;
                const int64_t jj = ( signs[1] /* η < 0 */ ? this->sw_op.nx(1) + 1 : 0 );

                void * buff_Y  = result.PointerAtQuadrantCell( quad, ii,jj );
                int    count_Y = this->sw_op.nx(0);
                int    tag_Y   = quad;

                MPI_Irecv( buff_Y, count_Y, cell_Y, rank_source_Y, tag_Y, Global::MPI_cart_comm,
                           &comm_reqs[1] );
            }

            # pragma omp master
            MPI_Waitall( 2, comm_reqs, MPI_STATUSES_IGNORE );

            # pragma omp barrier

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

            // Send new boundary data.
            # pragma omp master
            if ( rank_dest_X != MPI_PROC_NULL ) {

                // Coordinate of cell to send data from for current quadrant.
                const int64_t ii = ( signs[0] /* ξ < 0 */ ? 1 : this->sw_op.nx(0) );
                const int64_t jj = 1;

                void * buff_X  = result.PointerAtQuadrantCell( quad, ii,jj );
                int    count_X = this->sw_op.nx(1);
                int    tag_X   = quad;

                MPI_Isend( buff_X, count_X, cell_X, rank_dest_X, tag_X, Global::MPI_cart_comm,
                           &comm_reqs[0] );
            }

            # pragma omp master
            if ( rank_dest_Y != MPI_PROC_NULL ) {

                // Coordinates of cells to send data from for current quadrant.
                const int64_t ii = 1;
                const int64_t jj = ( signs[1] /* η < 0 */ ? 1 : this->sw_op.nx(1) );

                void * buff_Y  = result.PointerAtQuadrantCell( quad, ii,jj);
                int    count_Y = this->sw_op.nx(0);
                int    tag_Y   = quad;

                MPI_Isend( buff_Y, count_Y, cell_Y, rank_dest_Y, tag_Y, Global::MPI_cart_comm,
                           &comm_reqs[1] );
            }

            # pragma omp master
            MPI_Waitall( 2, comm_reqs, MPI_STATUSES_IGNORE );

            # pragma omp barrier
        }

    # if defined (ENABLE_HWLOC) && defined (USE_HWLOC_ALLOC)
        hwloc_free( Global::machine_topology, temporary_memory, bytes_allocd );
    # else
        std::free( temporary_memory );
    # endif
    }

    if ( Global::MPI_num_ranks > 1 )
        result.halo_cells_dirty = true;

    result.upwind_halo_cells_dirty = false;
}


//============================================================================================================
//=== EXPLICIT INSTANTIATIONS ================================================================================
//============================================================================================================

# include "linear_solvers/Abstract/SweepPattern/SweepPatternDerivedFactory.hpp"
# include "objects/RKDG/OrdinateFlux.hpp"
# include "objects/STDG/OrdinateFlux.hpp"


//!
//! \brief  Instantiation of SweepOperator<RKDG::OrdinateFlux>::SweepPatternUniCartWavefrontMPI class.
//!
template class SweepOperator<RKDG::OrdinateFlux>::SweepPatternUniCartWavefrontMPI;

//!
//! \brief  Instantiation of factory for SweepOperator<RKDG::OrdinateFlux>::SweepPatternUniCartWavefrontMPI.
//!
template class
SweepOperator<RKDG::OrdinateFlux>::SweepPattern::SweepPatternDerivedFactory<
        SweepOperator<RKDG::OrdinateFlux>::SweepPatternUniCartWavefrontMPI
    >;


//!
//! \brief  Instantiation of SweepOperator<STDG::OrdinateFlux>::SweepPatternUniCartWavefrontMPI class.
//!
template class SweepOperator<STDG::OrdinateFlux>::SweepPatternUniCartWavefrontMPI;

//!
//! \brief  Instantiation of factory for SweepOperator<STDG::OrdinateFlux>::SweepPatternUniCartWavefrontMPI.
//!
template class
SweepOperator<STDG::OrdinateFlux>::SweepPattern::SweepPatternDerivedFactory<
        SweepOperator<STDG::OrdinateFlux>::SweepPatternUniCartWavefrontMPI
    >;


# endif // if SPACE_DIMS == 2
# endif // if defined (ENABLE_MPI)
