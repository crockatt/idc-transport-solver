//------------------------------------------------------------------------------------------------------------
//! \file   linear_solvers/Abstract/SweepPattern/SweepPatternUniCartKBA2.cpp
//! \brief  Contains implementations and instantiations of methods from SweepPatternUniCartKBA2 class.
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
# include "linear_solvers/Abstract/SweepPattern/SweepPatternUniCartKBA2.hpp"
# include "linear_solvers/Abstract/SweepCellSolve/SweepCellSolveManager.hpp"
# include "utils/SIMD.hpp"


//!
//! \brief  Sleep time used for spinlocks to prevent hogging of the synchronization variables.
//!
const struct timespec s_sleep_ts {

    0,      // tv_sec
    1'000   // tv_nsec
};


//============================================================================================================
//=== MEMBER DEFINITIONS FOR SweepPatternUniCartKBA2 CLASS ====================================================
//============================================================================================================


//--------------------------------------------------------------------------------------------------------
//! \brief  Returns the string descriptor of the class type.
//--------------------------------------------------------------------------------------------------------
template< class OrdinateFlux >
const std::string SweepOperator<OrdinateFlux>::SweepPatternUniCartKBA2::Descriptor( void ) {

    return "kba2";
}


//--------------------------------------------------------------------------------------------------------
//! \brief  Returns the string descriptor for an objects type.
//--------------------------------------------------------------------------------------------------------
template< class OrdinateFlux >
const std::string SweepOperator<OrdinateFlux>::SweepPatternUniCartKBA2::GetDescriptor( void ) const {

    return Descriptor();
}


//------------------------------------------------------------------------------------------------------------
//! \brief  Virtual destructor.
//------------------------------------------------------------------------------------------------------------
template< class OrdinateFlux >
SweepOperator<OrdinateFlux>::SweepPatternUniCartKBA2::~SweepPatternUniCartKBA2( void ) {

    PRINT_STATUS( "Executing SweepOperator<%s>::SweepPatternUniCartKBA2::%s.\n",
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
void SweepOperator<OrdinateFlux>::SweepPatternUniCartKBA2::Print (

    const std::string & prefix // = "  "

) const {

    PRINT_STATUS( "Executing SweepOperator<%s>::SweepPatternUniCartKBA2::%s.\n",
                  typeid(OrdinateFlux).name(), __func__ )

    PRINT_LOG( "%s%-*s  %s\n", prefix.c_str(),
               Global::col_width, "Sweep Pattern:",  this->Descriptor().c_str() )

    PRINT_LOG( "%s%-*s % d\n", prefix.c_str(),
               Global::col_width, "KBA block size:", this->KBA_block_size     )

    this->SweepPatternUniCart::Print( prefix );
}


//------------------------------------------------------------------------------------------------------------
//! \brief  Sets parameters given in a ParameterList.
//------------------------------------------------------------------------------------------------------------
template< class OrdinateFlux >
void SweepOperator<OrdinateFlux>::SweepPatternUniCartKBA2::SetParameters (

    const ParameterList & input_list
) {

    PRINT_STATUS( "Executing SweepOperator<%s>::SweepPatternUniCartKBA2::%s.\n",
                  typeid(OrdinateFlux).name(), __func__ )

    try         {  this->KBA_block_size = input_list.GetValue<int64_t>( "kba_block_size" );  }
    catch (...) {/* empty */}

    this->SweepPatternUniCart::SetParameters( input_list );
}


//------------------------------------------------------------------------------------------------------------
//! \brief  Performs a transport sweep.
//!
//! ### Rationale and Standards Compliance
//!
//! This implementation is motivated by the observation that most MPI implementations are unable to
//! progress background communications well enough to enable effective overlapping of communication and
//! computation \cite Denis2016. This implementation seeks to alleviate this shortcoming through the use of
//! a dedicated communication threads.
//!
//! Separation of communication and computation threads is achieved using nested OpenMP parallel regions. The
//! outermost parallel region spawns only two threads. These threads will become the master threads of the
//! communication and computation thread teams. Thread IDs are obtained in the outermost parallel region from
//! the OpenMP library routine \c omp_get_thread_num and stored in \c tid_outer. The innermost parallel region
//! is used to spawn multiple work threads in the computation team. The communication thread creates a team
//! bound to the innermost parallel region consisting of only one thread (itself). Thread IDs obtained in the
//! innermost parallel region are stored in \c tid_inner. In the innermost parallel region, threads are
//! uniquely identified by the _pair_ (\c tid_outer, \c tid_inner).
//!
//! References to sections of the OpenMP standard (version 4.5, November 2015) that are relevant to the
//! compliance of this implementation with said standard are listed below.
//! Reference locations are given in the format:
//!
//! > §\<section\> p.\<page\>:\<line\>
//!
//! #### §2.5 \c parallel Construct
//!
//! - p. 48:13-14
//! > The binding thread set for a \c parallel region is the encountering thread. The encountering thread
//! > becomes the master thread of the new team.
//! - p. 48:16-17
//! > When a thread encounters a \c parallel construct, a team of threads is created to execute the
//! > \c parallel region.
//! - p. 48:19-22
//! > The thread that encountered the \c parallel construct becomes the master thread of the new team, with a
//! > thread number of zero for the duration of the new \c parallel region. All threads in the new team,
//! > including the master thread, execute the region.
//! - p. 48:26-28
//! > Within a \c parallel region, thread numbers uniquely identify each thread. Thread numbers are
//! > consecutive whole numbers ranging from zero for the master thread up to one less than the number of
//! > threads in the team.
//! - p. 49:11-13
//! > If a thread in a team executing a \c parallel region encounters another \c parallel directive, it
//! > creates a new team ... and it becomes the master of that new team.
//!
//! #### §2.7.1 Loop Construct
//!
//! - p. 58:17-20
//! > The binding thread set for a loop region is the current team. A loop region binds to the innermost
//! > enclosing \c parallel region. Only the threads of the team executing the binding \c parallel region
//! > participate in the execution of the loop iterations and the implied barrier of the loop region...
//!
//! #### §2.13.3 \c barrier Construct
//!
//! - p. 152:5-6
//! > The binding thread set for a \c barrier region is the current team. A \c barrier region binds to the
//! > innermost enclosing parallel region.
//!
//! #### §3.2.2 <code>omp_get_num_threads</code>
//!
//! - p. 233:2-3
//! > The binding region for an \c omp_get_num_threads region is the innermost enclosing \c parallel region.
//! - p. 233:5-7
//! > The \c omp_get_num_threads routine returns the number of threads in the team executing the \c parallel
//! > region to which the routine region binds.
//!
//! #### §3.2.4 <code>omp_get_thread_num</code>
//!
//! - p. 235:7-8
//! > The binding thread set for an \c omp_get_thread_num region is the current team. The binding region for
//! > an \c omp_get_thread_num region is the innermost enclosing \c parallel region.
//! - p. 235:10-12
//! > The \c omp_get_thread_num routine returns the thread number of the calling thread, within the team
//! > executing the \c parallel region to which the routine region binds. The thread number is an integer
//! > between 0 and one less than the value returned by \c omp_get_num_threads, inclusive.
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
void SweepOperator<OrdinateFlux>::SweepPatternUniCartKBA2::Sweep (

    const OrdinateFlux & source,
    OrdinateFlux & result,
    const RKDG::CrossSection & sigma,
    const RKDG::OrdinateFlux * const initial // = nullptr

) const {

    PRINT_STATUS( "Executing SweepOperator<%s>::SweepPatternUniCartKBA2::%s.\n",
                  typeid(OrdinateFlux).name(), __func__ )

    // Reference to SweepCellSolve object to call for local solves.
    const SweepCellSolveManager & cell_solve = *this->sw_op.sweep_cell_solve;

    /* Counters used to synchronize communication and computation thread teams.
     *
     * NOTE: Cells are counted across all angular quadrants.
     */
    int64_t cells_received = 0;
    int64_t cells_computed = 0;

    // Number of bytes of temporary memory to allocate per thread.
    const size_t bytes_allocd = SWEEP_SIMD_LEN * sizeof(double) * (   cell_solve.GetDimofA()
                                                                    + cell_solve.GetDimofB()
                                                                    + cell_solve.GetDimofWork() );
    // Construct list of angular blocks to loop over.
    auto angle_blocks = this->ConstructSIMDAngleBlocks();

    # pragma omp parallel num_threads(2)
    /************************************************************************************************
     |  This outermost parallel region creates exactly two threads: one to handle communication     |
     |  tasks and one that will spawn a thread team within a nested parallel region to handle       |
     |  computation.                                                                                |
     |                                                                                              |
     |  The thread IDs for the outermost parallel region are organized as                           |
     |                                                                                              |
     |  tid_outer = 0: Communication team                                                           |
     |              1: Computation team                                                             |
     |                                                                                              |
     ************************************************************************************************/
    {
        const int tid_outer = omp_get_thread_num();

    // "if" clause forces communication team to include only one thread.
    # pragma omp parallel if( tid_outer != 0 )
    {
        // const int tid_inner = omp_get_thread_num();

        switch ( tid_outer ) {

        case 0:
        /********************************************************************************************
         |  COMMUNICATION TEAM                                                                      |
         |                                                                                          |
         |  Variables in this section use the following conventions:                                |
         |                                                                                          |
         |  1. Variables ending in X or Y indicate that the variable is associated with messages    |
         |     based on a shift in the X or Y dimension of the Cartesian topology, respectively.    |
         |                                                                                          |
         |     Note that for a shift in the X (Y) dimension, the cells along the corresponding      |
         |     edges of the spatial domain have different j (i) indices; i.e., the index used is    |
         |     the opposite of that typically used for the dimension indicated by the shift.        |
         ********************************************************************************************/
        {

            // Used to specify the signs associated with the coordinates of the ordinates in a given
            // quadrant.
            bool signs[SPACE_DIMS] = { 0, 0 };

            // Array of MPI_Request objects used for pipeline messages.
            MPI_Request comm_reqs [] = { MPI_REQUEST_NULL,   // 0: Used for receive operations, Y-shift send.
                                         MPI_REQUEST_NULL }; // 1: Used for X-shift send operations.

            for ( int64_t quad_index : { 0, 1, 2, 3 } ) {

                // Determine next quadrant to sweep based on quadrant ordering.
                const int64_t quad = this->sweep_order[quad_index];

                // Determine signs on ordinates in current quadrant.
                Quadrule::OrdinateSet::SetSigns( quad, signs );

                // Construct resized MPI_Datatype objects in each dimension for the current quadrant.
                MPI_Datatype cell_X, cell_Y;

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

                // Compute shifts for quadrant.
                int rank_source_X, rank_dest_X, direction_X = 0, disp_X = 1,
                    rank_source_Y, rank_dest_Y, direction_Y = 1, disp_Y = 1;

                if ( signs[0] /* ξ < 0 */ ) {  disp_X = -1;  }
                if ( signs[1] /* η < 0 */ ) {  disp_Y = -1;  }

                MPI_Cart_shift( Global::MPI_cart_comm, direction_X, disp_X, &rank_source_X, &rank_dest_X );
                MPI_Cart_shift( Global::MPI_cart_comm, direction_Y, disp_Y, &rank_source_Y, &rank_dest_Y );

                // Wait for incoming boundary data with shift in Y dimension.
                if ( rank_source_Y != MPI_PROC_NULL ) {

                    // Coordinates of halo cells to receive data into for current quadrant.
                    const int64_t ii = 1;
                    const int64_t jj = ( signs[1] /* η < 0 */ ? this->sw_op.nx(1) + 1 : 0 );

                    void * buff_Y  = result.PointerAtQuadrantCell( quad, ii,jj );
                    int    count_Y = this->sw_op.nx(0);
                    int    tag_Y   = quad;

                    MPI_Recv( buff_Y, count_Y, cell_Y, rank_source_Y, tag_Y, Global::MPI_cart_comm,
                              MPI_STATUS_IGNORE );
                }

                cells_received += 1;
                # pragma omp flush( cells_received )

                // Short-circuit: no receive operations are required along pipeline for this rank and quadrant.
                if ( rank_source_X == MPI_PROC_NULL ) {

                    cells_received += this->sw_op.nx(1);
                    # pragma omp flush( cells_received )
                }

                // Traverse pipeline
                for ( int64_t blk_j = 1; blk_j <= this->sw_op.nx(1); blk_j += KBA_block_size ) {

                    // Height of current block.
                    int count_X = std::min( this->KBA_block_size, this->sw_op.nx(1) - blk_j + 1 );

                    // Coordinate of halo cells for send/receive for current quadrant.
                    const int64_t jj = ( signs[1] /* η < 0 */ ? this->sw_op.nx(1) - blk_j - count_X + 2
                                                              : blk_j );

                    // Tag for send/receive operations for current block of pipeline.
                    int tag_X = (jj - 1) + this->sw_op.nx(1) * quad;

                    // Post receive for incoming boundary data with shift in X dimension.
                    if ( rank_source_X != MPI_PROC_NULL ) {

                        // Coordinate of halo cell for receive for current quadrant.
                        const int64_t ii = ( signs[0] /* ξ < 0 */ ? this->sw_op.nx(0) + 1 : 0 );

                        void * buff_X = result.PointerAtQuadrantCell( quad, ii,jj );

                        MPI_Recv( buff_X, count_X, cell_X, rank_source_X, tag_X, Global::MPI_cart_comm,
                                  MPI_STATUS_IGNORE );

                        // Increment receive count.
                        cells_received += count_X;
                        # pragma omp flush( cells_received )
                    }

                    // Wait for pending send operations to complete.
                    MPI_Waitall( 2, comm_reqs, MPI_STATUSES_IGNORE );

                    // Post send for boundary data with shift in X dimension.
                    if ( rank_dest_X != MPI_PROC_NULL ) {

                        // Wait for computation to finish.
                    # if defined (ENABLE_KBA_WAIT_TIMING)
                        Global::TMR_KBA_comm_wait.Start();
                    # endif

                        while ( cells_computed - quad_index * this->sw_op.nx(1) < blk_j + count_X - 1 ) {

                            nanosleep( &s_sleep_ts, nullptr );

                            # pragma omp flush( cells_computed )
                        }

                    # if defined (ENABLE_KBA_WAIT_TIMING)
                        Global::TMR_KBA_comm_wait.Stop();
                    # endif

                        // Coordinate of cell to send data from for current quadrant.
                        const int64_t ii = ( signs[0] /* ξ < 0 */ ? 1 : this->sw_op.nx(0) );

                        void * buff_X = result.PointerAtQuadrantCell( quad, ii,jj );

                        MPI_Isend( buff_X, count_X, cell_X, rank_dest_X, tag_X, Global::MPI_cart_comm,
                                   &comm_reqs[0] );
                    }
                } // end for blk_j.

                // Send boundary data with shift in Y dimension.
                if ( rank_dest_Y != MPI_PROC_NULL ) {

                    // Wait for computation to finish.
                # if defined (ENABLE_KBA_WAIT_TIMING)
                    Global::TMR_KBA_comm_wait.Start();
                # endif

                    while ( cells_computed - quad_index * this->sw_op.nx(1) < this->sw_op.nx(1) ) {

                        nanosleep( &s_sleep_ts, nullptr );

                        # pragma omp flush( cells_computed )
                    }

                # if defined (ENABLE_KBA_WAIT_TIMING)
                    Global::TMR_KBA_comm_wait.Stop();
                # endif

                    // Coordinates of cells to send data from for current quadrant.
                    const int64_t ii = 1;
                    const int64_t jj = ( signs[1] /* η < 0 */ ? 1 : this->sw_op.nx(1) );

                    void * buff_Y  = result.PointerAtQuadrantCell( quad, ii,jj);
                    int    count_Y = this->sw_op.nx(0);
                    int    tag_Y   = quad;

                    MPI_Isend( buff_Y, count_Y, cell_Y, rank_dest_Y, tag_Y, Global::MPI_cart_comm,
                               &comm_reqs[0] );
                }

                // Wait for remaining send operations to complete.
                MPI_Waitall( 2, comm_reqs, MPI_STATUSES_IGNORE );

                // Cleanup MPI_Datatype objects.
                MPI_Type_free( &cell_X );
                MPI_Type_free( &cell_Y );

            } // end for quad index.
        } break; // end case 0 (communication team).


        case 1:
        /********************************************************************************************
         *     COMPUTATION TEAM                                                                     *
         ********************************************************************************************/
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

            // Used to specify the signs associated with the coordinates of the ordinates in a given
            // quadrant.
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

                // Wait for incoming boundary data with shift in Y dimension.
                # pragma omp master
                {
                    while ( cells_received - quad_index * ( this->sw_op.nx(1) + 1 )
                            < 1
                    ) {
                        nanosleep( &s_sleep_ts, nullptr );

                        # pragma omp flush( cells_received )
                    }
                }
                # pragma omp barrier

            // Traverse pipeline.
            for ( int64_t blk_j = 1; blk_j <= this->sw_op.nx(1); blk_j += this->KBA_block_size ) {

                // Height of current block.
                const int64_t blk_height = std::min( this->KBA_block_size, this->sw_op.nx(1) - blk_j + 1 );

                // Number of diagonals in current block
                const int64_t num_diags = this->sw_op.nx(0) + blk_height - 1;

                // Wait for receive.
                # pragma omp master
                {
                /*
                 * Because the condition blk_j > 1 is used, the timer TMR_KBA_comp_wait records
                 * only the portion of the wait time between blocks; i.e., the time spent waiting
                 * for the first messages to arrive from neighboring ranks is not included.
                 */
                # if defined (ENABLE_KBA_WAIT_TIMING)
                    if ( blk_j > 1 ) {  Global::TMR_KBA_comp_wait.Start();  }
                # endif

                    while ( cells_received - quad_index * ( this->sw_op.nx(1) + 1 ) - 1
                            < blk_j + blk_height - 1
                    ) {
                        nanosleep( &s_sleep_ts, nullptr );

                        # pragma omp flush( cells_received )
                    }

                # if defined (ENABLE_KBA_WAIT_TIMING)
                    if ( blk_j > 1 ) {  Global::TMR_KBA_comp_wait.Stop();  }
                # endif
                }
                # pragma omp barrier

                // Sweep using diagonal wavefront.
                for ( int64_t diag = 1; diag <= num_diags; ++diag ) {

                    const int64_t diag_length = std::min( diag,
                                                          std::min( num_diags - diag + 1,
                                                                    std::min( this->sw_op.nx(0), blk_height )
                                                        ) );

                    # pragma omp for collapse(2)
                    for ( int64_t cell = 1; cell <= diag_length; ++cell ) {
                    for ( int64_t q_blk = 0; q_blk < num_angle_blocks; ++q_blk ) {

                        // Copy block containing angular component and length of block from list.
                        SIMD_BlkIdx<0> idx = quad_angle_blocks[ q_blk ];

                        // Compute indices of spatial cell to update.
                        int64_t i = std::min( diag, this->sw_op.nx(0) ) - cell + 1;
                        int64_t j = blk_j + diag - i;

                        // Reverse direction of sweep for negative angles.
                        if ( signs[0] /* ξ < 0 */ ) {  i = this->sw_op.nx(0) - i + 1;  }
                        if ( signs[1] /* η < 0 */ ) {  j = this->sw_op.nx(1) - j + 1;  }

                        // Set spatial index in SIMD block object.
                        idx.i = i;
                        idx.j = j;

                        // Do solve.
                        cell_solve.ConstructAndSolve( A, B, idx, sigma, source, result, WORK, initial );
                    }}

                } // end for diag.


                // Increment count.
                # pragma omp master
                {
                    cells_computed += blk_height;
                    # pragma omp flush( cells_computed )
                }

            } // end for blk_j.

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

        # if defined (ENABLE_HWLOC) && defined (USE_HWLOC_ALLOC)
            hwloc_free( Global::machine_topology, temporary_memory, bytes_allocd );
        # else
            std::free( temporary_memory );
        # endif
        } break; // end case 1 (computation team).

        /********************************************************************************************
         *  Invalid outer thread ID.                                                                *
         ********************************************************************************************/
        default:
        {   std::string error_message =   "Variable tid_outer is '"
                                        + std::to_string(tid_outer)
                                        + "' in '"
                                        + std::string(__func__)
                                        + ": should be 0 or 1.\n";

            PRINT_ERROR( error_message.c_str() )
        }
        } // end switch tid_outer.
    }} // end parallel regions.

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
//! \brief  Instantiation of SweepOperator<RKDG::OrdinateFlux>::SweepPatternUniCartKBA2 class.
//!
template class SweepOperator<RKDG::OrdinateFlux>::SweepPatternUniCartKBA2;

//!
//! \brief  Instantiation of factory for SweepOperator<RKDG::OrdinateFlux>::SweepPatternUniCartKBA2.
//!
template class
SweepOperator<RKDG::OrdinateFlux>::SweepPattern::SweepPatternDerivedFactory<
        SweepOperator<RKDG::OrdinateFlux>::SweepPatternUniCartKBA2
    >;


//!
//! \brief  Instantiation of SweepOperator<STDG::OrdinateFlux>::SweepPatternUniCartKBA2 class.
//!
template class SweepOperator<STDG::OrdinateFlux>::SweepPatternUniCartKBA2;

//!
//! \brief  Instantiation of factory for SweepOperator<STDG::OrdinateFlux>::SweepPatternUniCartKBA2.
//!
template class
SweepOperator<STDG::OrdinateFlux>::SweepPattern::SweepPatternDerivedFactory<
        SweepOperator<STDG::OrdinateFlux>::SweepPatternUniCartKBA2
    >;


# endif // if SPACE_DIMS == 2
# endif // if defined (ENABLE_MPI)
