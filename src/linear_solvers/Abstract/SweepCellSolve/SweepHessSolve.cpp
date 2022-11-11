//------------------------------------------------------------------------------------------------------------
//! \file   linear_solvers/Abstract/SweepCellSolve/SweepHessSolve.cpp
//! \brief  Contains implementation and instantiations of methods from SweepHessSolve class template.
//!
//! \author Michael M. Crockatt
//! \date   March 2018
//------------------------------------------------------------------------------------------------------------

# include <typeinfo>

# include "linear_solvers/Abstract/SweepCellSolve/SweepHessSolve.hpp"
# include "linear_solvers/LinearAlgebra.h"
# include "linear_solvers/SIMD_LinearAlgebra.hpp"


//============================================================================================================
//=== MEMBER DEFINITIONS FOR SweepHessSolve CLASS ============================================================
//============================================================================================================


//--------------------------------------------------------------------------------------------------------
//! \brief  Returns the string descriptor of the class type.
//--------------------------------------------------------------------------------------------------------
template< class OrdinateFlux, int64_t SIMD_length >
const std::string SweepHessSolve< OrdinateFlux, SIMD_length >::Descriptor( void ) {

    return "hessenberg";
}


//--------------------------------------------------------------------------------------------------------
//! \brief  Returns the string descriptor for an objects type.
//--------------------------------------------------------------------------------------------------------
template< class OrdinateFlux, int64_t SIMD_length >
const std::string SweepHessSolve< OrdinateFlux, SIMD_length >::GetDescriptor( void ) const {

    return Descriptor();
}


//------------------------------------------------------------------------------------------------------------
//! \brief  Construct solver object from given parameters.
//!
//! \param[in]  enclosing       Reference to instance of enclosing SweepOperator class.
//! \param[in]  input_list      Contains parameters for initialization of object.
//! \param[in]  dt_in           Initial timestep size.
//------------------------------------------------------------------------------------------------------------
template< class OrdinateFlux, int64_t SIMD_length >
SweepHessSolve< OrdinateFlux, SIMD_length >::SweepHessSolve (

    const SweepOperator<OrdinateFlux> & enclosing,
    const ParameterList & input_list,
    const double dt_in              // = std::numeric_limits<double>::infinity()
) :
    SweepLinearSolveSpec< OrdinateFlux, SIMD_length >( enclosing, input_list, dt_in )
{

    PRINT_STATUS( "Executing SweepHessSolve<%s,%" PRId64 ">::%s.\n",
                  typeid(OrdinateFlux).name(), SIMD_length, __func__ )

    // Force zone-wise constant approximation of cross sections.
    this->force_zwc = true;

    // Allocate memory for precomputed Hessenberg reductions.
    this->dimof_tau = enclosing.nq();
    this->tau = new double*[ this->dimof_tau ];

    for ( int64_t q = 0; q < this->dimof_tau; ++q )
        this->tau[q] = nullptr;

# if SPACE_DIMS == 1

    /*
     *  Spatial dimensions: 1
     *  SIMD blocking:      Yes
     */

    for ( int64_t sign = 0; sign <= 1; ++sign ) {
    /*
     *  sign = 0:   negative angles
     *       = 1:   positive angles
     */

        const int64_t q_min = this->sw_op.nq() * sign / 2;
        const int64_t q_max = this->sw_op.nq() * (sign + 1) / 2;

        # pragma omp parallel for
        for ( int64_t q = q_min; q < q_max; q += SIMD_length ) {

            this->tau[q] = (double *)
                # if defined (USE_ALIGNED_ALLOC)
                    aligned_alloc( ALIGNED_ALLOC_ALIGNMENT, SIMD_length * this->sizeof_A );
                # else
                    std::malloc( SIMD_length * this->sizeof_A );
                # endif
        }
    }

# elif SPACE_DIMS == 2

    /*
     *  Spatial dimensions: 2
     *  SIMD blocking:      Yes
     */

    for ( int64_t quad = 0; quad < 4; ++quad ) {

        const int64_t q_min = enclosing.Quadrants(quad);
        const int64_t q_max = enclosing.Quadrants(quad + 1);

        # pragma omp parallel for
        for ( int64_t q = q_min; q < q_max; q += SIMD_length ) {

            this->tau[q] = (double *)
                # if defined (USE_ALIGNED_ALLOC)
                    aligned_alloc( ALIGNED_ALLOC_ALIGNMENT, SIMD_length * this->sizeof_B );
                # else
                    std::malloc( SIMD_length * this->sizeof_B );
                # endif
        }
    }

# elif SPACE_DIMS == 3

    /*
     *  Spatial dimensions: 3
     *  SIMD blocking:      Yes
     */

    # warning "Implementation of SweepOperator<OrdinateFlux,SIMD_length>::SweepHessSolve incomplete."

    throw std::runtime_error( "Implementation of SweepHessSolve<"
                              + std::string( typeid(OrdinateFlux).name() ) + ","
                              + std::to_string(SIMD_length) + ">::" + std::string(__func__)
                              + " incomplete." );

# endif // if SPACE_DIMS == ?

    // Compute Hessenberg reductions of angular components of matrices.
    SweepHessSolve< OrdinateFlux, SIMD_length >::ComputeMatrices();
}


//--------------------------------------------------------------------------------------------------------
//! \brief  Prints a summary of the object's configuration to the logging interface.
//!
//! \param[in]  prefix  String prepended to each line of output, intended to be used for controlling
//!                     indentation of output.
//--------------------------------------------------------------------------------------------------------
template< class OrdinateFlux, int64_t SIMD_length >
void SweepHessSolve< OrdinateFlux, SIMD_length >::Print (

    const std::string & prefix // = "  "

) const {

    PRINT_STATUS( "Executing SweepHessSolve<%s,%" PRId64 ">::%s.\n",
                  typeid(OrdinateFlux).name(), SIMD_length, __func__ )

    PRINT_LOG( "%s%-*s  %s\n", prefix.c_str(),
               Global::col_width, "Sweep Solve:", this->Descriptor().c_str() )
}


//--------------------------------------------------------------------------------------------------------
//! \brief  Virtual destructor.
//--------------------------------------------------------------------------------------------------------
template< class OrdinateFlux, int64_t SIMD_length >
SweepHessSolve< OrdinateFlux, SIMD_length >::~SweepHessSolve ( void ) {

    PRINT_STATUS( "Executing SweepHessSolve<%s,%" PRId64 ">::%s.\n",
                  typeid(OrdinateFlux).name(), SIMD_length, __func__ )

    if ( this->tau != nullptr )
    {
        for ( int64_t q = 0; q < this->dimof_tau; ++q )
            if ( this->tau[q] != nullptr )
                std::free( this->tau[q] );

        delete [] this->tau;
    }
}


//--------------------------------------------------------------------------------------------------------
//! \brief  Solves the system for the updated coefficients.
//!
//! \param[in]      A           Pointer to memory location containing matrix for linear system to be solved.
//! \param[in,out]  B           Initially contains the RHS vector for the linear system to be solved.
//!                             Upon return, contains the solution of the system.
//! \param[in]      work_ptr    Pointer to temporary workspace used by solver.
//--------------------------------------------------------------------------------------------------------
template< class OrdinateFlux, int64_t SIMD_length >
void SweepHessSolve< OrdinateFlux, SIMD_length >::SolveSystem (

    double * const A,
    double * const B,
    void * const work_ptr

) const {

    PRINT_STATUS( "Executing SweepHessSolve<%s,%" PRId64 ">::%s.\n",
                  typeid(OrdinateFlux).name(), SIMD_length, __func__ )

# if defined (DO_SWEEP_SUBTIMING)

    Global::TMR_linearSolve[
        # if defined (_OPENMP)
            omp_get_thread_num()
        # else
            0
        # endif
        ].Start();

# endif // if defined (DO_SWEEP_SUBTIMING)

# if defined (ENABLE_SIMD_BLOCKING)
    SIMD_UpperHessenbergSolve_CM<SIMD_length>
# else
    UpperHessenbergSolve_CM
# endif
    ( this->system_dim, A, B, static_cast<double *>( work_ptr ), nullptr /* no pivoting */ );
    //!
    //! \todo   Update call to Hessenberg solver.
    //!

# if defined (DO_SWEEP_SUBTIMING)

    Global::TMR_linearSolve[
        # if defined (_OPENMP)
            omp_get_thread_num()
        # else
            0
        # endif
        ].Stop();

# endif // if defined (DO_SWEEP_SUBTIMING)

    return;
}


//------------------------------------------------------------------------------------------------------------
//! \brief  Zeros the arrays for storing precomputed values (e.g., angular component of matrices).
//------------------------------------------------------------------------------------------------------------
template< class OrdinateFlux, int64_t SIMD_length >
void SweepHessSolve< OrdinateFlux, SIMD_length >::ZeroMatrices ( void ) const {

    PRINT_STATUS( "Executing SweepHessSolve<%s,%" PRId64 ">::%s.\n",
                  typeid(OrdinateFlux).name(), SIMD_length, __func__ )

# if SPACE_DIMS == 1

    /*
     *  Spatial dimensions: 1
     *  SIMD blocking:      Yes
     */

    for ( int64_t sign = 0; sign <= 1; ++sign ) {
    /*
     *  sign = 0:   negative angles
     *       = 1:   positive angles
     */

        const int64_t q_min = this->sw_op.nq() * sign / 2;
        const int64_t q_max = this->sw_op.nq() * (sign + 1) / 2;

        // SIMD blocks of contiguous angles.
        # pragma omp parallel for
        for ( int64_t q = q_min; q < q_max; q += SIMD_length )
            std::memset( this->tau[q], 0, SIMD_length * this->sizeof_B );
    }

# elif SPACE_DIMS == 2

    /*
     *  Spatial dimensions: 2
     *  SIMD blocking:      Yes
     */

    for ( int64_t quad = 0; quad < 4; ++quad ) {

        const int64_t q_min = this->sw_op.Quadrants(quad);
        const int64_t q_max = this->sw_op.Quadrants(quad + 1);

        // SIMD blocks of contiguous angles.
        # pragma omp parallel for
        for ( int64_t q = q_min; q < q_max; q += SIMD_length )
            std::memset( this->tau[q], 0, SIMD_length * this->sizeof_B );
    }

# elif SPACE_DIMS == 3

    /*
     *  Spatial dimensions: 3
     *  SIMD blocking:      Yes
     *
     *  ---------------------------
     *  ----- NOT IMPLEMENTED -----
     *  ---------------------------
     */

    # warning "Implementation of SweepOperator<OrdinateFlux,SIMD_length>::ZeroMatrices incomplete."

    throw std::runtime_error( "Implementation of SweepHessSolve<"
                              + std::string( typeid(OrdinateFlux).name() ) + ","
                              + std::to_string(SIMD_length) + ">::" + std::string(__func__)
                              + " incomplete." );

# endif // if SPACE_DIMS == ?

    // Pass up hierarchy to zero Aq arrays.
    SweepLinearSolveBase< OrdinateFlux, SIMD_length >::ZeroMatrices();
}


//-----------------------------------------------------------------------------------------------------------
//! \brief  Computes Hessenberg reductions of the angular components of the matrices that define the linear
//!         systems.
//!
//! Assumes:
//!     - Memory at SweepHessSolve::tau has already been allocated (this is done in
//!       SweepHessSolve::SweepHessSolve.
//!     - The matrices at SweepLinearSolve::Aq have already been computed (this function only modifies them).
//!
//-----------------------------------------------------------------------------------------------------------
template< class OrdinateFlux, int64_t SIMD_length >
void SweepHessSolve< OrdinateFlux, SIMD_length >::ComputeMatrices ( void ) const {

    PRINT_STATUS( "Executing SweepHessSolve<%s,%" PRId64 ">::%s.\n",
                  typeid(OrdinateFlux).name(), SIMD_length, __func__ )

    // Needed so that Aq matrices are recomputed when this function is called by SetDt.
    SweepLinearSolveSpec< OrdinateFlux, SIMD_length >::ComputeMatrices();

# if SPACE_DIMS == 1

    /*
     *  Spatial dimensions: 1
     *  SIMD blocking:      Yes
     */

     # pragma omp parallel
    {
        int N     = this->system_dim;
        int ILO   = 1;
        int IHI   = N;
        int LDA   = N;
        int LWORK = N*N;
        int INFO;

        double * const WORK = new double[ LWORK ];

        double * const Aq_temp  = new double[ this->dimof_A ];
        double * const tau_temp = new double[ this->dimof_B ];

        for ( int64_t sign = 0; sign <= 1; ++sign ) {
        /*
         *  sign = 0:   negative angles
         *       = 1:   positive angles
         */

            const int64_t q_min = this->sw_op.nq() * sign / 2;
            const int64_t q_max = this->sw_op.nq() * (sign + 1) / 2;

            # pragma omp for collapse(2)
            for ( int64_t q = q_min; q < q_max; q += SIMD_length ) {
            for ( int64_t l = 0; l < SIMD_length; ++l ) {

                // Compute system for each angle.
                if ( q+l < q_max ) {

                    std::memset( Aq_temp,  0, this->sizeof_A );
                    std::memset( tau_temp, 0, this->sizeof_B );

                    // Load Aq matrix into temporary array.
                    for ( int64_t i = 0; i < this->system_dim; ++i ) {
                    for ( int64_t j = 0; j < this->system_dim; ++j ) {

                        Aq_temp[ MatIdxDI_CM<1>(i,j) ] = this->Aq[q][ MatIdxDI_CM<SIMD_length>(i,j,l) ];
                    }}

                    // Compute Hessenberg decomposition...
                    dgehrd_( &N, &ILO, &IHI, Aq_temp, &LDA, tau_temp, WORK, &LWORK, &INFO );

                    // ... and then store it.
                    for ( int64_t i = 0; i < this->system_dim; ++i ) {
                    for ( int64_t j = 0; j < this->system_dim; ++j ) {

                        this->Aq[q][ MatIdxDI_CM<SIMD_length>(i,j,l) ] = Aq_temp[ MatIdxDI_CM<1>(i,j) ];
                    }}

                    for ( int64_t i = 0; i < this->system_dim; ++i )
                        this->tau[q][ VecIdxDI<SIMD_length>(i,l) ] = tau_temp[ VecIdxDI<1>(i) ];

                // Pad remaining elements of system stack with identity matrices.
                } else {

                    // Zeros.
                    for ( int64_t i = 0; i < this->system_dim; ++i ) {
                    for ( int64_t j = 0; j < this->system_dim; ++j ) {

                        this->Aq[q][ MatIdxDI_CM<SIMD_length>(i,j,l) ] = 0.0;
                    }
                        this->tau[q][ VecIdxDI<SIMD_length>(i,l) ] = 0.0;
                    }

                    // Ones.
                    for ( int64_t i = 0; i < this->system_dim; ++i )
                        this->Aq[q][ MatIdxDI_CM<SIMD_length>(i,i,l) ] = 1.0;
                }
            }}
        }

        delete [] WORK;
        delete [] Aq_temp;
        delete [] tau_temp;

    } // end omp parallel.

# elif SPACE_DIMS == 2

    /*
     *  Spatial dimensions: 2
     *  SIMD blocking:      Yes
     */

     # pragma omp parallel
    {
        int N     = this->system_dim;
        int ILO   = 1;
        int IHI   = N;
        int LDA   = N;
        int LWORK = N*N;
        int INFO;

        double * const WORK = new double[ LWORK ];

        double * const Aq_temp  = new double[ this->dimof_A ];
        double * const tau_temp = new double[ this->dimof_B ];

        for ( int64_t quad = 0; quad < 4; ++quad ) {

            const int64_t q_min = this->sw_op.Quadrants(quad);
            const int64_t q_max = this->sw_op.Quadrants(quad + 1);

            // Full blocks of contiguous angles.
            # pragma omp for collapse(2)
            for ( int64_t q = q_min; q < q_max; q += SIMD_length ) {
            for ( int64_t l = 0; l < SIMD_length; ++l ) {

                // Compute system for each angle.
                if ( q+l < q_max ) {

                    std::memset( Aq_temp,  0, this->sizeof_A );
                    std::memset( tau_temp, 0, this->sizeof_B );

                    // Load Aq matrix into temporary array.
                    for ( int64_t i = 0; i < this->system_dim; ++i ) {
                    for ( int64_t j = 0; j < this->system_dim; ++j ) {

                        Aq_temp[ MatIdxDI_CM<1>(i,j) ] = this->Aq[q][ MatIdxDI_CM<SIMD_length>(i,j,l) ];
                    }}

                    // Compute Hessenberg decomposition...
                    dgehrd_( &N, &ILO, &IHI, Aq_temp, &LDA, tau_temp, WORK, &LWORK, &INFO );

                    // ... and then store it.
                    for ( int64_t i = 0; i < this->system_dim; ++i ) {
                    for ( int64_t j = 0; j < this->system_dim; ++j ) {

                        this->Aq[q][ MatIdxDI_CM<SIMD_length>(i,j,l) ] = Aq_temp[ MatIdxDI_CM<1>(i,j) ];
                    }}

                    for ( int64_t i = 0; i < this->system_dim; ++i )
                        this->tau[q][ VecIdxDI<SIMD_length>(i,l) ] = tau_temp[ VecIdxDI<1>(i) ];

                // Pad remaining elements of system stack with identity matrices.
                } else {

                    // Zeros.
                    for ( int64_t i = 0; i < this->system_dim; ++i ) {
                    for ( int64_t j = 0; j < this->system_dim; ++j ) {

                        this->Aq[q][ MatIdxDI_CM<SIMD_length>(i,j,l) ] = 0.0;
                    }
                        this->tau[q][ VecIdxDI<SIMD_length>(i,l) ] = 0.0;
                    }

                    // Ones.
                    for ( int64_t i = 0; i < this->system_dim; ++i )
                        this->Aq[q][ MatIdxDI_CM<SIMD_length>(i,i,l) ] = 1.0;
                }
            }}
        }

        delete [] WORK;
        delete [] Aq_temp;
        delete [] tau_temp;

    } // end omp parallel.

# elif SPACE_DIMS == 3

    /*
     *  Spatial dimensions: 3
     *  SIMD blocking:      Yes
     *
     *  ---------------------------
     *  ----- NOT IMPLEMENTED -----
     *  ---------------------------
     */

    # warning "Implementation of SweepOperator<OrdinateFlux,SIMD_length>::ComputeMatrices incomplete."

    throw std::runtime_error( "Implementation of SweepHessSolve<"
                              + std::string( typeid(OrdinateFlux).name() ) + ","
                              + std::to_string(SIMD_length) + ">::" + std::string(__func__)
                              + " incomplete." );

# endif // if SPACE_DIMS == ?
}


//-----------------------------------------------------------------------------------------------------------
//! \brief  Computes the elements of the matrix/matrices for the linear system(s) to be solved at the given
//!         step of a sweep.
//!
//! Computes the local matrix/matrices to be solved by the sweep algorithm for the ordinate(s) and spatial
//! cell(s) specified by \pp{idx}.
//!
//! For solving steady-state problems one should use <tt>\pp{dt} = inf</tt>.
//!
//! \param[out]     A           Pointer to memory location to store computed matrix/matrices.
//! \param[in]      idx         Contains indices for the system(s) to compute matrix/matrices for.
//! \param[in]      sigma       Total cross section \f$ \sigma_{\mathrm{t}} \f$.
//! \param[out]     work_ptr    Pointer for temporary workspace.
//------------------------------------------------------------------------------------------------------------
template< class OrdinateFlux, int64_t SIMD_length >
void SweepHessSolve< OrdinateFlux, SIMD_length >::CalcA (

    double * const A,
    const SIMD_BlkIdx<0> & idx,
    const RKDG::CrossSection & sigma,
    void * const work_ptr

) const {

    PRINT_STATUS( "Executing SweepHessSolve<%s,%" PRId64 ">::%s.\n",
                  typeid(OrdinateFlux).name(), SIMD_length, __func__ )

    SweepLinearSolveSpec< OrdinateFlux, SIMD_length >::CalcA( A, idx, sigma, work_ptr );

    /*
     *  Spatial dimensions: Any
     *  SIMD blocking:      Yes
     */

    const int64_t q = idx.q;

    double * const WORK = static_cast<double *>( work_ptr );

    // Full blocks.

    for ( int64_t i = 0; i < this->system_dim; ++i ) {
/*  # pragma omp simd aligned( WORK : SWEEP_SIMD_LEN * sizeof(double) )
 *
 *  This simd clause causes a compilation error when using Intel's compiler. Specifically, using
 *
 *      icpc (ICC) 18.0.1 20171018
 *
 *  yields the following error:
 *
 *      SweepHessSolve.cpp(638): error: syntax error in omp clause
 *              # pragma omp simd aligned( WORK : SWEEP_SIMD_LEN * sizeof(double) )
 *                                                                                ^
 *
 *  The source of the problem is unclear. Since I am not sure that this clause offers a significant
 *  performance benefit I have just removed it.
 *
 *      -- Michael M. Crockatt  2018-02-25
 */
    for ( int64_t l = 0; l < SIMD_length; ++l ) {

        WORK[ VecIdxDI<SIMD_length>(i,l) ] = this->tau[q][ VecIdxDI<SIMD_length>(i,l) ];
    }}
}


//============================================================================================================
//=== EXPLICIT INSTANTIATIONS ================================================================================
//============================================================================================================

# include "objects/RKDG/OrdinateFlux.hpp"
# include "objects/STDG/OrdinateFlux.hpp"


# if SWEEP_SIMD_LEN >= (1 << 0)
    template class SweepHessSolve< RKDG::OrdinateFlux, (1 << 0) >;
    template class SweepHessSolve< STDG::OrdinateFlux, (1 << 0) >;

    template class SweepCellSolveDerivedFactory< RKDG::OrdinateFlux, (1 << 0), SweepHessSolve >;
    template class SweepCellSolveDerivedFactory< STDG::OrdinateFlux, (1 << 0), SweepHessSolve >;
# endif
# if SWEEP_SIMD_LEN >= (1 << 1)
    template class SweepHessSolve< RKDG::OrdinateFlux, (1 << 1) >;
    template class SweepHessSolve< STDG::OrdinateFlux, (1 << 1) >;

    template class SweepCellSolveDerivedFactory< RKDG::OrdinateFlux, (1 << 1), SweepHessSolve >;
    template class SweepCellSolveDerivedFactory< STDG::OrdinateFlux, (1 << 1), SweepHessSolve >;
# endif
# if SWEEP_SIMD_LEN >= (1 << 2)
    template class SweepHessSolve< RKDG::OrdinateFlux, (1 << 2) >;
    template class SweepHessSolve< STDG::OrdinateFlux, (1 << 2) >;

    template class SweepCellSolveDerivedFactory< RKDG::OrdinateFlux, (1 << 2), SweepHessSolve >;
    template class SweepCellSolveDerivedFactory< STDG::OrdinateFlux, (1 << 2), SweepHessSolve >;
# endif
# if SWEEP_SIMD_LEN >= (1 << 3)
    template class SweepHessSolve< RKDG::OrdinateFlux, (1 << 3) >;
    template class SweepHessSolve< STDG::OrdinateFlux, (1 << 3) >;

    template class SweepCellSolveDerivedFactory< RKDG::OrdinateFlux, (1 << 3), SweepHessSolve >;
    template class SweepCellSolveDerivedFactory< STDG::OrdinateFlux, (1 << 3), SweepHessSolve >;
# endif
# if SWEEP_SIMD_LEN >= (1 << 4)
    template class SweepHessSolve< RKDG::OrdinateFlux, (1 << 4) >;
    template class SweepHessSolve< STDG::OrdinateFlux, (1 << 4) >;

    template class SweepCellSolveDerivedFactory< RKDG::OrdinateFlux, (1 << 4), SweepHessSolve >;
    template class SweepCellSolveDerivedFactory< STDG::OrdinateFlux, (1 << 4), SweepHessSolve >;
# endif
# if SWEEP_SIMD_LEN >= (1 << 5)
    template class SweepHessSolve< RKDG::OrdinateFlux, (1 << 5) >;
    template class SweepHessSolve< STDG::OrdinateFlux, (1 << 5) >;

    template class SweepCellSolveDerivedFactory< RKDG::OrdinateFlux, (1 << 5), SweepHessSolve >;
    template class SweepCellSolveDerivedFactory< STDG::OrdinateFlux, (1 << 5), SweepHessSolve >;
# endif
# if SWEEP_SIMD_LEN > (1 << 5)
    # error "Maximum SIMD vector length supported is (1 << 5)."
# endif
