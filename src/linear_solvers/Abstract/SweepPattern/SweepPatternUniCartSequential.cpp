//------------------------------------------------------------------------------------------------------------
//! \file   linear_solvers/Abstract/SweepPattern/SweepPatternUniCartSequential.cpp
//! \brief  Contains implementations and instantiations of methods from SweepPatternUniCartSequential class.
//!
//! \author Michael M. Crockatt
//! \date   March 2018
//------------------------------------------------------------------------------------------------------------

# include <typeinfo>
# include <vector>

# if defined (ENABLE_HWLOC) || defined (DOXYCOMPILE)
    # include <hwloc.h>
# endif

# if defined (_OPENMP)
    # include <omp.h>
# endif

# include "linear_solvers/Abstract/SweepOperator.hpp"
# include "linear_solvers/Abstract/SweepPattern/SweepPatternUniCartSequential.hpp"
# include "linear_solvers/Abstract/SweepCellSolve/SweepCellSolveManager.hpp"
# include "utils/SIMD.hpp"


//============================================================================================================
//=== MEMBER DEFINITIONS FOR SweepPatternUniCartSequential CLASS ==============================================
//============================================================================================================


//--------------------------------------------------------------------------------------------------------
//! \brief  Returns the string descriptor of the class type.
//--------------------------------------------------------------------------------------------------------
template< class OrdinateFlux >
const std::string SweepOperator<OrdinateFlux>::SweepPatternUniCartSequential::Descriptor( void ) {

    return "sequential";
}


//--------------------------------------------------------------------------------------------------------
//! \brief  Returns the string descriptor for an objects type.
//--------------------------------------------------------------------------------------------------------
template< class OrdinateFlux >
const std::string SweepOperator<OrdinateFlux>::SweepPatternUniCartSequential::GetDescriptor( void ) const {

    return Descriptor();
}


//------------------------------------------------------------------------------------------------------------
//! \brief  Virtual destructor.
//------------------------------------------------------------------------------------------------------------
template< class OrdinateFlux >
SweepOperator<OrdinateFlux>::SweepPatternUniCartSequential::~SweepPatternUniCartSequential( void ) {

    PRINT_STATUS( "Executing SweepOperator<%s>::SweepPatternUniCartSequential::%s.\n",
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
void SweepOperator<OrdinateFlux>::SweepPatternUniCartSequential::Print (

    const std::string & prefix // = "  "

) const {

    PRINT_STATUS( "Executing SweepOperator<%s>::SweepPatternUniCartSequential::%s.\n",
                  typeid(OrdinateFlux).name(), __func__ )

    PRINT_LOG( "%s%-*s  %s\n", prefix.c_str(),
               Global::col_width, "Sweep Pattern:", this->Descriptor().c_str() )

    this->SweepPatternUniCart::Print( prefix );
}


//------------------------------------------------------------------------------------------------------------
//! \brief  Sets parameters given in a ParameterList.
//------------------------------------------------------------------------------------------------------------
template< class OrdinateFlux >
void SweepOperator<OrdinateFlux>::SweepPatternUniCartSequential::SetParameters (

    const ParameterList & input_list
) {

    PRINT_STATUS( "Executing SweepOperator<%s>::SweepPatternUniCartSequential::%s.\n",
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
void SweepOperator<OrdinateFlux>::SweepPatternUniCartSequential::Sweep (

    const OrdinateFlux & source,
    OrdinateFlux & result,
    const RKDG::CrossSection & sigma,
    const RKDG::OrdinateFlux * const initial // = nullptr

) const {

    PRINT_STATUS( "Executing SweepOperator<%s>::SweepPatternUniCartSequential::%s.\n",
                  typeid(OrdinateFlux).name(), __func__ )

    // Reference to SweepCellSolve object to call for local solves.
    const SweepCellSolveManager & cell_solve = *this->sw_op.sweep_cell_solve;

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

    # if SPACE_DIMS == 1

        /*
         *  Spatial dimensions: 1
         */

        const int64_t num_angle_blocks = angle_blocks.size();

        // Sweep cell-by-cell.
        for ( int64_t i = 1; i <= this->sw_op.nx(0); ++i ) {

            # pragma omp for
            for ( int64_t q_blk = 0; q_blk < num_angle_blocks; ++q_blk ) {

                // Copy block containing angular component and length of block from list.
                SIMD_BlkIdx<0> idx = angle_blocks[ q_blk ];

                int64_t ii = i;

                // Reverse direction of sweep for negative angles.
                if ( this->sw_op.xi( idx.q ) < 0 ) {  ii = this->sw_op.nx(0) - i + 1;  }

                // Set spatial index in SIMD block object.
                idx.i = ii;

                // Do solve.
                cell_solve.ConstructAndSolve( A, B, idx, sigma, source, result, WORK, initial );

            } // end for q_blk.
        } // end for i.

    # elif SPACE_DIMS == 2

        /*
         *  Spatial dimensions: 2
         */

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

            // Sweep cell-by-cell
            for ( int64_t i = 1; i <= this->sw_op.nx(0); ++i ) {
            for ( int64_t j = 1; j <= this->sw_op.nx(1); ++j ) {

                # pragma omp for
                for ( int64_t q_blk = 0; q_blk < num_angle_blocks; ++q_blk ) {

                    // Copy block containing angular component and length of block from list.
                    SIMD_BlkIdx<0> idx = quad_angle_blocks[ q_blk ];

                    int64_t ii = i;
                    int64_t jj = j;

                    // Reverse direction of sweep for negative angles.
                    if ( signs[0] /* ξ < 0 */ ) {  ii = this->sw_op.nx(0) - i + 1;  }
                    if ( signs[1] /* η < 0 */ ) {  jj = this->sw_op.nx(1) - j + 1;  }

                    // Set spatial index in SIMD block object.
                    idx.i = ii;
                    idx.j = jj;

                    // Do solve.
                    cell_solve.ConstructAndSolve( A, B, idx, sigma, source, result, WORK, initial );
                }
            }} // end for i,j

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

    # elif SPACE_DIMS == 3

        /*
         *  Spatial dimensions: 3
         */

        # warning "Implementation of SweepOperator<OrdinateFlux>::SweepPatternUniCartSequential::Sweep incomplete."

        throw std::runtime_error( "Implementation of SweepOperator<"
                                  + std::string( typeid(OrdinateFlux).name() )
                                  + ">::SweepPatternUniCartSequential::" + std::string(__func__)
                                  + " incomplete." );

    # endif // if SPACE_DIMS == ?

    /* --- Cleanup. ----------------------------------------------------------------------------------- */

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
//! \brief  Instantiation of SweepOperator<RKDG::OrdinateFlux>::SweepPatternUniCartSequential class.
//!
template class SweepOperator<RKDG::OrdinateFlux>::SweepPatternUniCartSequential;

//!
//! \brief  Instantiation of factory for SweepOperator<RKDG::OrdinateFlux>::SweepPatternUniCartSequential.
//!
template class
SweepOperator<RKDG::OrdinateFlux>::SweepPattern::SweepPatternDerivedFactory<
        SweepOperator<RKDG::OrdinateFlux>::SweepPatternUniCartSequential
    >;


//!
//! \brief  Instantiation of SweepOperator<STDG::OrdinateFlux>::SweepPatternUniCartSequential class.
//!
template class SweepOperator<STDG::OrdinateFlux>::SweepPatternUniCartSequential;

//!
//! \brief  Instantiation of factory for SweepOperator<STDG::OrdinateFlux>::SweepPatternUniCartSequential.
//!
template class
SweepOperator<STDG::OrdinateFlux>::SweepPattern::SweepPatternDerivedFactory<
        SweepOperator<STDG::OrdinateFlux>::SweepPatternUniCartSequential
    >;
