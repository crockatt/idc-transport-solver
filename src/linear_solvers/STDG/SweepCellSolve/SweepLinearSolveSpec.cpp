//------------------------------------------------------------------------------------------------------------
//! \file   linear_solvers/STDG/SweepCellSolve/SweepLinearSolveSpec.hpp
//! \brief  Contains partial specialization of SweepLinearSolveSpec class for STDG::OrdinateFlux objects.
//!
//! \author Michael M. Crockatt
//! \date   March 2018
//------------------------------------------------------------------------------------------------------------

# include <typeinfo>

# include "linear_solvers/STDG/SweepCellSolve/SweepCellSolve.hpp"
# include "linear_solvers/STDG/SweepCellSolve/SweepLinearSolveSpec.hpp"


//-----------------------------------------------------------------------------------------------------------
//! \brief  Construct solver object from given parameters.
//!
//! \param[in]  enclosing       Reference to instance of enclosing SweepOperator class.
//! \param[in]  input_list      Contains parameters for initialization of object.
//! \param[in]  dt_in           Initial timestep size.
//-----------------------------------------------------------------------------------------------------------
template< int64_t SIMD_length >
SweepLinearSolveSpec< STDG::OrdinateFlux, SIMD_length >::SweepLinearSolveSpec (

    const SweepOperator<STDG::OrdinateFlux> & enclosing,
    const ParameterList & input_list,
    const double dt_in              // = std::numeric_limits<double>::infinity()
) :
    SweepLinearSolveBase< STDG::OrdinateFlux, SIMD_length >( enclosing, input_list, dt_in )
{
    PRINT_STATUS( "Executing SweepLinearSolveSpec<STDG::OrdinateFlux,%" PRId64 ">::%s.\n",
                  SIMD_length, __func__ )

    // Compute angular components of matrices.
    this->SweepLinearSolveSpec< STDG::OrdinateFlux, SIMD_length >::ComputeMatrices();
}


//-----------------------------------------------------------------------------------------------------------
//! \brief  Computes the angular component of the matrices that define the linear systems.
//!
//! Assumes that the memory at SweepLinearSolveBase::Aq has already been allocated (this is done in
//! SweepLinearSolveBase::SweepLinearSolveBase.
//!
//-----------------------------------------------------------------------------------------------------------
template< int64_t SIMD_length >
void SweepLinearSolveSpec< STDG::OrdinateFlux, SIMD_length >::ComputeMatrices ( void ) const {

    PRINT_STATUS( "Executing SweepLinearSolveSpec<STDG::OrdinateFlux,%" PRId64 ">::%s.\n",
                  SIMD_length, __func__ )

    this->ZeroMatrices();

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
        for ( int64_t q  = q_min; q < q_max; q += SIMD_length ) {
        for ( int64_t l = 0; l < SIMD_length; ++l ) {

            // Only do cases where (q+l) < q_max.
            if ( (q+l) >= q_max ) {  continue;  }

            for ( int64_t d = 0; d <= this->DG_degree_x; ++d ) {
            for ( int64_t s = 0; s <= this->DG_degree_t; ++s ) {

                // Components from time derivative term.
                for ( int64_t r = s-1; r >= 0; r -= 2 )
                    this->Aq[q][ this->template MatIdx_CM<SIMD_length>(d,s,d,r,l) ]
                        -= 2.0 * this->sw_op.dx(0) * (2*s + 1);

                for ( int64_t r = 0; r <= this->DG_degree_t; ++r )
                    this->Aq[q][ this->template MatIdx_CM<SIMD_length>(d,s,d,r,l) ]
                        += this->sw_op.dx(0) * (2*s + 1);

                // Components from advection in ξ angular direction.
                for ( int64_t a = d-1; a >= 0; a -= 2 )
                    this->Aq[q][ this->template MatIdx_CM<SIMD_length>(d,s,a,s,l) ]
                        -= 2.0 * this->dt * this->sw_op.xi(q+l) * (2*d + 1);

                if ( this->sw_op.xi(q+l) < 0 ) {

                    for ( int64_t a = 0; a <= this->DG_degree_x; ++a )
                        this->Aq[q][ this->template MatIdx_CM<SIMD_length>(d,s,a,s,l) ]
                            += neg1pow(d+a+1) * this->dt * this->sw_op.xi(q+l) * (2*d + 1);
                } else {

                    for ( int64_t a = 0; a <= this->DG_degree_x; ++a )
                        this->Aq[q][ this->template MatIdx_CM<SIMD_length>(d,s,a,s,l) ]
                            += this->dt * this->sw_op.xi(q+l) * (2*d + 1);
                }
            }}
        }}
    }

# elif SPACE_DIMS == 2

    /*
     *  Spatial dimensions: 2
     *  SIMD blocking:      Yes
     */

    for ( int64_t quad = 0; quad < 4; ++quad ) {

        const int64_t q_min = this->sw_op.Quadrants(quad);
        const int64_t q_max = this->sw_op.Quadrants(quad + 1);

        # pragma omp parallel for
        for ( int64_t q = q_min; q < q_max; q += SIMD_length ) {
        for ( int64_t l = 0; l < SIMD_length; ++l ) {

            // Only do cases where (q+l) < q_max.
            if ( (q+l) >= q_max ) {  continue;  }

            const double xi_coeff  = this->dt * this->sw_op.dx(1) * this->sw_op.xi(q+l);
            const double eta_coeff = this->dt * this->sw_op.dx(0) * this->sw_op.eta(q+l);
            const double dx_coeff  = this->sw_op.dx(0) * this->sw_op.dx(1);

            for ( int64_t d = 0; d <= this->DG_degree_x; ++d ) {
            for ( int64_t e = 0; e <= this->DG_degree_x; ++e ) {
            for ( int64_t s = 0; s <= this->DG_degree_t; ++s ) {

                // Components from time derivative term.
                for ( int64_t r = s-1; r >= 0; r -= 2 )
                    this->Aq[q][ this->template MatIdx_CM<SIMD_length>(d,e,s,d,e,r,l) ]
                        -= 2.0 * dx_coeff * (2*s + 1);

                for ( int64_t r = 0; r <= this->DG_degree_t; ++r )
                    this->Aq[q][ this->template MatIdx_CM<SIMD_length>(d,e,s,d,e,r,l) ]
                        += dx_coeff * (2*s + 1);

                // Components from spatial derivative in ξ angular direction.
                for ( int64_t a = d-1; a >= 0; a -= 2 )
                    this->Aq[q][ this->template MatIdx_CM<SIMD_length>(d,e,s,a,e,s,l) ]
                        -= 2.0 * xi_coeff * (2*d + 1);

                if ( this->sw_op.xi(q+l) < 0 ) {

                    for ( int64_t a = 0; a <= this->DG_degree_x; ++a )
                        this->Aq[q][ this->template MatIdx_CM<SIMD_length>(d,e,s,a,e,s,l) ]
                            += neg1pow(d+a+1) * xi_coeff * (2*d + 1);
                } else {

                    for ( int64_t a = 0; a <= this->DG_degree_x; ++a )
                        this->Aq[q][ this->template MatIdx_CM<SIMD_length>(d,e,s,a,e,s,l) ]
                            += xi_coeff * (2*d + 1);
                }

                // Components from spatial derivative in η angular direction.
                for ( int64_t b = e-1; b >= 0; b -= 2 )
                    this->Aq[q][ this->template MatIdx_CM<SIMD_length>(d,e,s,d,b,s,l) ]
                        -= 2.0 * eta_coeff * (2*e + 1);

                if ( this->sw_op.eta(q+l) < 0 ) {

                    for ( int64_t b = 0; b <= this->DG_degree_x; ++b )
                        this->Aq[q][ this->template MatIdx_CM<SIMD_length>(d,e,s,d,b,s,l) ]
                            += neg1pow(e+b+1) * eta_coeff * (2*e + 1);
                } else {

                    for ( int64_t b = 0; b <= this->DG_degree_x; ++b )
                        this->Aq[q][ this->template MatIdx_CM<SIMD_length>(d,e,s,d,b,s,l) ]
                            += eta_coeff * (2*e + 1);
                }
            }}}
        }}
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

    # warning "Implementation of SweepLinearSolveSpec<STDG::OrdinateFlux,SIMD_length>::ComputeMatrices incomplete."

    throw std::runtime_error( "Implementation of SweepLinearSolveSpec<STDG::OrdinateFlux,"
                              + std::to_string(SIMD_length) + ">::" + std::string(__func__)
                              + " incomplete." );

# endif // if SPACE_DIMS == ?
}


//------------------------------------------------------------------------------------------------------------
//! \brief  Sets the timestep size for the sweep, recomputing values if needed.
//------------------------------------------------------------------------------------------------------------
template< int64_t SIMD_length >
void SweepLinearSolveSpec< STDG::OrdinateFlux, SIMD_length >::SetDt (

    const double dt_in
) {
    PRINT_STATUS( "Executing SweepLinearSolveSpec<STDG::OrdinateFlux,%" PRId64 ">::%s.\n",
                  SIMD_length, __func__ )

    this->SweepCellSolveBase<STDG::OrdinateFlux>::SetDt( dt_in );

    // STDG implementation requires re-computing angular components of systems when timestep changes.
    this->ComputeMatrices();
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
// \param[out]     work_ptr    Pointer for temporary workspace.
//------------------------------------------------------------------------------------------------------------
template< int64_t SIMD_length >
void SweepLinearSolveSpec< STDG::OrdinateFlux, SIMD_length >::CalcA (

    double * const A,
    const SIMD_BlkIdx<0> & idx,
    const RKDG::CrossSection & sigma,
    void * const

) const {

    PRINT_STATUS( "Executing SweepLinearSolveSpec<STDG::OrdinateFlux,%" PRId64 ">::%s.\n",
                  SIMD_length, __func__ )

# if defined (DO_SWEEP_SUBTIMING)

    Global::TMR_calcA[
        # if defined (_OPENMP)
            omp_get_thread_num()
        # else
            0
        # endif
        ].Start();

# endif // if defined (DO_SWEEP_SUBTIMING)

# if SPACE_DIMS == 1

    /*
     *  Spatial dimensions: 1
     *  SIMD blocking:      Yes
     */

    const int64_t q = idx.q;
    const int64_t i = idx.i;

    std::memcpy( A, this->Aq[q], this->sizeof_A * SIMD_length );

    // Augment with cross section.
    if (    this->force_zwc
         || sigma.DG_degree == 0
    ) {

        const double coeff = this->dt * this->sw_op.dx(0) * sigma(i,0);

        for ( int64_t d = 0; d <= this->DG_degree_x; ++d ) {
        for ( int64_t s = 0; s <= this->DG_degree_t; ++s ) {

            # pragma omp simd aligned( A : SWEEP_SIMD_LEN * SIZEOF_DOUBLE )
            for ( int64_t l = 0; l < SIMD_length; ++l )
                A[ this->template MatIdx_CM<SIMD_length>(d,s,d,s,l) ] += coeff;
        }}

    } else {

        const double coeff = this->dt * this->sw_op.dx(0) / 2.0;

        for ( int64_t d = 0; d <= this->DG_degree_x; ++d ) {
        for ( int64_t u = 0; u <= this->DG_degree_x; ++u ) {
        for ( int64_t s = 0; s <= this->DG_degree_t; ++s ) {

            const double coeff2 = coeff * sigma(i,d,u) * (2*d + 1);

            # pragma omp simd aligned( A : SWEEP_SIMD_LEN * SIZEOF_DOUBLE )
            for ( int64_t l = 0; l < SIMD_length; ++l )
                A[ this->template MatIdx_CM<SIMD_length>(d,s,u,s,l) ] += coeff2;
        }}}
    }

# elif SPACE_DIMS == 2

    /*
     *  Spatial dimensions: 2
     *  SIMD blocking:      Yes
     */

    const int64_t q = idx.q;
    const int64_t i = idx.i;
    const int64_t j = idx.j;

    double * const Aq_ptr = this->Aq[q];

//     # pragma omp simd aligned( A, Aq_ptr : SIMD_length * SIZEOF_DOUBLE )
//     for ( uint64_t i = 0; i < this->dimof_A * SIMD_length; ++i )
//         A[i] = Aq_ptr[i];

    std::memcpy( A, Aq_ptr, this->sizeof_A * SIMD_length );

    if (    this->force_zwc
            || sigma.DG_degree == 0
    ) {

        const double coeff = this->dt * this->sw_op.dx(0) * this->sw_op.dx(1) * sigma(i,j,0,0);

        for ( int64_t d = 0; d <= this->DG_degree_x; ++d ) {
        for ( int64_t e = 0; e <= this->DG_degree_x; ++e ) {
        for ( int64_t s = 0; s <= this->DG_degree_t; ++s ) {

            # pragma omp simd aligned( A : SWEEP_SIMD_LEN * SIZEOF_DOUBLE )
            for ( int64_t l = 0; l < SIMD_length; ++l )
                A[ this->template MatIdx_CM<SIMD_length>(d,e,s,d,e,s,l) ] += coeff;
        }}}

    } else {

        const double dx_coeff = this->dt * this->sw_op.dx(0) * this->sw_op.dx(1) / 4.0;

        for ( int64_t d = 0; d <= this->DG_degree_x; ++d ) {
        for ( int64_t e = 0; e <= this->DG_degree_x; ++e ) {
        for ( int64_t u = 0; u <= this->DG_degree_x; ++u ) {
        for ( int64_t v = 0; v <= this->DG_degree_x; ++v ) {

            const double coeff = dx_coeff * (2*d + 1) * (2*e + 1) * sigma(i,j,d,e,u,v);

            for ( int64_t s = 0; s <= this->DG_degree_t; ++s ) {
            # pragma omp simd aligned( A : SWEEP_SIMD_LEN * SIZEOF_DOUBLE )
            for ( int64_t l = 0; l < SIMD_length; ++l ) {

                A[ this->template MatIdx_CM<SIMD_length>(d,e,s,u,v,s,l) ] += coeff;
            }}
        }}}}
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

    # warning "Implementation of SweepLinearSolveSpec<STDG::OrdinateFlux,SIMD_length>::CalcA incomplete."

    throw std::runtime_error( "Implementation of SweepLinearSolveSpec<STDG::OrdinateFlux,"
                              + std::to_string(SIMD_length) + ">::" + std::string(__func__)
                              + " incomplete." );

# endif // if SPACE_DIMS == ?

# if defined (DO_SWEEP_SUBTIMING)

    Global::TMR_calcA[
        # if defined (_OPENMP)
            omp_get_thread_num()
        # else
            0
        # endif
        ].Stop();

# endif // if defined (DO_SWEEP_SUBTIMING)
}


# if 0
//------------------------------------------------------------------------------------------------------------
//! \brief  Computes the elements of the right-hand vectors \f$ b \f$ to be used for solving the linear
//!         systems at the given step of a sweep.
//!
//! \param[out]     B           Pointer to memory location to store computed vector.
//! \param[in]      idx         Contains indices for the system(s) to construct vector(s) for.
//! \param[in]      source      Contains the source term \f$ Q \f$ of the transport equation.
//! \param[in]      result      Contains the solution of the transport equation.
// \param[in]      work_ptr    Pointer to temporary workspace.
//! \param[in]      initial     (optional) <br>
//!                             Pointer to initial condition for the step. If this is null then the initial
//!                             condition is assumed to be zero and this pointer is not dereferenced.
//------------------------------------------------------------------------------------------------------------
template<>
void SweepLinearSolveSpec< STDG::OrdinateFlux, 1 >
::CalcB (

    double * const B,
    const SIMD_BlkIdx<0> & idx,
    const STDG::OrdinateFlux & source,
    const STDG::OrdinateFlux & result,
    void * const,
    const RKDG::OrdinateFlux * const initial

) const {

    PRINT_STATUS( "Executing SweepLinearSolveSpec<STDG::OrdinateFlux,1>::%s.\n", __func__ )

# if defined (DO_SWEEP_SUBTIMING)

    Global::TMR_calcB[
        # if defined (_OPENMP)
            omp_get_thread_num()
        # else
            0
        # endif
        ].Start();

# endif // if defined (DO_SWEEP_SUBTIMING)

# if SPACE_DIMS == 1

    /*
     *  Spatial dimensions: 1
     *  SIMD blocking:      No
     */

    const int64_t q = idx.q;
    const int64_t i = idx.i;

    for ( int64_t d = 0; d <= this->DG_degree_x; ++d ) {
    for ( int64_t s = 0; s <= this->DG_degree_t; ++s ) {

        // Set value from known source term.
        B[ this->template VecIdx<1>(d,s) ] = this->dt * this->sw_op.dx(0) * source(q,i,d,s);

        // Include known components of time derivative term.
        if ( initial != nullptr )
            B[ this->template VecIdx<1>(d,s) ] += neg1pow(s) * this->sw_op.dx(0) * (2*s + 1)
                                                  * (*initial)(q,i,d);

        // Include known components of advection term.
        if ( this->sw_op.xi(q) < 0 ) {

            for ( int64_t a = 0; a <= this->DG_degree_x; ++a )
                B[ this->template VecIdx<1>(d,s) ] += neg1pow(a+1) * this->dt * this->sw_op.xi(q) * (2*d + 1)
                                                      * result(q,i+1,a,s);
        } else {

            for ( int64_t a = 0; a <= this->DG_degree_x; ++a )
                B[ this->template VecIdx<1>(d,s) ]+= neg1pow(d) * this->dt * this->sw_op.xi(q) * (2*d + 1)
                                                     * result(q,i-1,a,s);
        }
    }}

# elif SPACE_DIMS == 2

    /*
     *  Spatial dimensions: 2
     *  SIMD blocking:      No
     */

    const int64_t q = idx.q;
    const int64_t i = idx.i;
    const int64_t j = idx.j;

    const double xi_coeff  = this->dt * this->sw_op.dx(1) * this->sw_op.xi(q);
    const double eta_coeff = this->dt * this->sw_op.dx(0) * this->sw_op.eta(q);

    for ( int64_t d = 0; d <= this->DG_degree_x; ++d ) {
    for ( int64_t e = 0; e <= this->DG_degree_x; ++e ) {
    for ( int64_t s = 0; s <= this->DG_degree_t; ++s ) {

        // Set value from known source term.
        B[ this->template VecIdx<1>(d,e,s) ]
            = this->dt * this->sw_op.dx(0) * this->sw_op.dx(1) * source(q,i,j,d,e,s);

        // Include known components of time derivative term.
        if ( initial != nullptr )
            B[ this->template VecIdx<1>(d,e,s) ]
                += neg1pow(s) * this->sw_op.dx(0) * this->sw_op.dx(1) * (2*s + 1) * (*initial)(q,i,j,d,e);

        // Include known components of advection term for ξ angular direction.
        if ( this->sw_op.xi(q) < 0 ) {

            for ( int64_t a = 0; a <= this->DG_degree_x; ++a )
                B[ this->template VecIdx<1>(d,e,s) ]
                    += neg1pow(a+1) * xi_coeff * (2*d + 1) * result(q,i+1,j,a,e,s);

        } else {

            for ( int64_t a = 0; a <= this->DG_degree_x; ++a )
                B[ this->template VecIdx<1>(d,e,s) ]
                    += neg1pow(d) * xi_coeff * (2*d + 1) * result(q,i-1,j,a,e,s);
        }

        // Include known components of advection term for η angular direction.
        if ( this->sw_op.eta(q) < 0 ) {

            for ( int64_t b = 0; b <= this->DG_degree_x; ++b )
                B[ this->template VecIdx<1>(d,e,s) ]
                    += neg1pow(b+1) * eta_coeff * (2*e + 1) * result(q,i,j+1,d,b,s);

        } else {

            for ( int64_t b = 0; b <= this->DG_degree_x; ++b )
                B[ this->template VecIdx<1>(d,e,s) ]
                += neg1pow(e) * eta_coeff * (2*e + 1) * result(q,i,j-1,d,b,s);
        }
    }}}

# elif SPACE_DIMS == 3

    /*
     *  Spatial dimensions: 3
     *  SIMD blocking:      No
     *
     *  ---------------------------
     *  ----- NOT IMPLEMENTED -----
     *  ---------------------------
     */

    # warning "Implementation of SweepLinearSolveSpec<STDG::OrdinateFlux,SIMD_length>::CalcB incomplete."

    throw std::runtime_error( "Implementation of SweepLinearSolveSpec<STDG::OrdinateFlux,1>::"
                              + std::string(__func__) + " incomplete." );

# endif // if SPACE_DIMS == ?

# if defined (DO_SWEEP_SUBTIMING)

    Global::TMR_calcB[
        # if defined (_OPENMP)
            omp_get_thread_num()
        # else
            0
        # endif
        ].Stop();

# endif // if defined (DO_SWEEP_SUBTIMING)
}
# endif // if 0


//------------------------------------------------------------------------------------------------------------
//! \brief  Computes the elements of the right-hand vectors \f$ b \f$ to be used for solving the linear
//!         systems at the given step of a sweep.
//!
//! \param[out]     B           Pointer to memory location to store computed vector.
//! \param[in]      idx         Contains indices for the system(s) to construct vector(s) for.
//! \param[in]      source      Contains the source term \f$ Q \f$ of the transport equation.
//! \param[in]      result      Contains the solution of the transport equation.
//! \param[in]      work_ptr    Pointer to temporary workspace.
//! \param[in]      initial     (optional) <br>
//!                             Pointer to initial condition for the step. If this is null then the initial
//!                             condition is assumed to be zero and this pointer is not dereferenced.
//------------------------------------------------------------------------------------------------------------
template< int64_t SIMD_length >
void SweepLinearSolveSpec< STDG::OrdinateFlux, SIMD_length >::CalcB (

    double * const B,
    const SIMD_BlkIdx<0> & idx,
    const STDG::OrdinateFlux & source,
    const STDG::OrdinateFlux & result,
    void * const work_ptr,
    const RKDG::OrdinateFlux * const initial

) const {

    PRINT_STATUS( "Executing SweepLinearSolveSpec<STDG::OrdinateFlux,%" PRId64 ">::%s.\n",
                  SIMD_length, __func__ )

# if defined (DO_SWEEP_SUBTIMING)

    Global::TMR_calcB[
        # if defined (_OPENMP)
            omp_get_thread_num()
        # else
            0
        # endif
        ].Start();

# endif // if defined (DO_SWEEP_SUBTIMING)

# if SPACE_DIMS == 1

    /*
     *  Spatial dimensions: 1
     *  SIMD blocking:      Yes
     */

//!
//! \brief  Indexing function for temporary array for caching values from result.
//!
//! \param[in]  k   Dimension for shift in spatial index. In { 0, 1 }.
//! \param[in]  d   Degree of test function with respect to \f$ x_1 \f$.
//! \param[in]  s   Degree of test function with respect to time.
//! \param[in]  l   Index for system in SIMD block. In [ 0, SIMD_length ).
//!
//! \see    VecIdx()
//!
# define IWORK(k,d,s,l) ( WORK[ this->template VecIdx<SIMD_length>(d,s,l) + (SIMD_length * this->dimof_B)*(k) ] )

    double * const WORK = static_cast<double *>( work_ptr );

    const int64_t q = idx.q;
    const int64_t i = idx.i;

    // Determine sign on angles in block.
    bool all_ang_negative = true;

    for ( int64_t l = 0; l < SIMD_length; ++l )
        all_ang_negative &= this->sw_op.xi( q+l ) <= 0.0;

    const int64_t xoffset = ( all_ang_negative ? +1 : -1 );

    // Load upwind values into temporary storage.
    for ( int64_t d = 0; d <= this->DG_degree_x; ++d ) {
    for ( int64_t s = 0; s <= this->DG_degree_t; ++s ) {
    for ( int64_t l = 0; l < SIMD_length; ++l ) {

        IWORK(0,d,s,l) = result( q+l, i + xoffset, d,s);
    }}}

    // Scale temporary values.
    for ( int64_t d = 0; d <= this->DG_degree_x; ++d ) {
    for ( int64_t s = 0; s <= this->DG_degree_t; ++s ) {
    # pragma omp simd aligned( WORK : SWEEP_SIMD_LEN * SIZEOF_DOUBLE )
    for ( int64_t l = 0; l < SIMD_length; ++l ) {

        IWORK(0,d,s,l) *= this->dt * this->sw_op.xi(q+l);
    }}}

    // Set value from known source term.
    for ( int64_t d = 0; d <= this->DG_degree_x; ++d ) {
    for ( int64_t s = 0; s <= this->DG_degree_t; ++s ) {
    # pragma omp simd aligned( WORK : SWEEP_SIMD_LEN * SIZEOF_DOUBLE )
    for ( int64_t l = 0; l < SIMD_length; ++l ) {

        B[ this->template VecIdx<SIMD_length>(d,s,l) ] = this->dt * this->sw_op.dx(0) * source(q+l,i,d,s);
    }}}

    // Include known components of time derivative term.
    if ( initial != nullptr ) {

        for ( int64_t d = 0; d <= this->DG_degree_x; ++d ) {
        for ( int64_t s = 0; s <= this->DG_degree_t; ++s ) {

            const double dx_coeff = neg1pow(s) * this->sw_op.dx(0) * (2*s + 1);

            # pragma omp simd aligned( WORK : SWEEP_SIMD_LEN * SIZEOF_DOUBLE )
            for ( int64_t l = 0; l < SIMD_length; ++l )
                B[ this->template VecIdx<SIMD_length>(d,s,l) ] += dx_coeff * (*initial)(q+l,i,d);
        }}
    }

    // Include upwind values in ξ angular direction.
    if ( all_ang_negative ) {

        for ( int64_t d = 0; d <= this->DG_degree_x; ++d ) {
        for ( int64_t a = 0; a <= this->DG_degree_x; ++a ) {
        for ( int64_t s = 0; s <= this->DG_degree_t; ++s ) {

            const int64_t coeff = neg1pow(a+1) * (2*d + 1);

            # pragma omp simd aligned( WORK : SWEEP_SIMD_LEN * SIZEOF_DOUBLE )
            for ( int64_t l = 0; l < SIMD_length; ++l )
                B[ this->template VecIdx<SIMD_length>(d,s,l) ] += IWORK(0,a,s,l) * coeff;
        }}}

    } else {

        for ( int64_t d = 0; d <= this->DG_degree_x; ++d ) {
        for ( int64_t a = 0; a <= this->DG_degree_x; ++a ) {
        for ( int64_t s = 0; s <= this->DG_degree_t; ++s ) {

            const int64_t coeff = neg1pow(d) * (2*d + 1);

            # pragma omp simd aligned( WORK : SWEEP_SIMD_LEN * SIZEOF_DOUBLE )
            for ( int64_t l = 0; l < SIMD_length; ++l )
                B[ this->template VecIdx<SIMD_length>(d,s,l) ] += IWORK(0,a,s,l) * coeff;
        }}}
    }

# elif SPACE_DIMS == 2

    /*
     *  Spatial dimensions: 2
     *  SIMD blocking:      Yes
     */

//!
//! \brief  Indexing function for temporary array for caching values from result.
//!
//! \param[in]  k   Dimension for shift in spatial index. In { 0, 1 }.
//! \param[in]  d   Degree of test function with respect to \f$ x_1 \f$.
//! \param[in]  e   Degree of test function with respect to \f$ x_2 \f$.
//! \param[in]  s   Degree of test function with respect to time.
//! \param[in]  l   Index for system in SIMD block. In [ 0, SIMD_length ).
//!
//! \see    VecIdx()
//!
# define IWORK(k,d,e,s,l) ( WORK[ this->template VecIdx<SIMD_length>(d,e,s,l) + (SIMD_length * this->dimof_B)*(k) ] )

    double * const WORK = static_cast<double *>( work_ptr );

    const int64_t q = idx.q;
    const int64_t i = idx.i;
    const int64_t j = idx.j;

    // Determine sign on angles in block.
    bool all_xi_negative  = true;
    bool all_eta_negative = true;

    for ( int64_t l = 0; l < SIMD_length; ++l ) {

        all_xi_negative  &= this->sw_op.xi(  q+l ) <= 0.0;
        all_eta_negative &= this->sw_op.eta( q+l ) <= 0.0;
    }

    const int64_t xoffset = ( all_xi_negative  ? +1 : -1 );
    const int64_t yoffset = ( all_eta_negative ? +1 : -1 );

    // Load upwind values into temporary storage.
    for ( int64_t d = 0; d <= this->DG_degree_x; ++d ) {
    for ( int64_t e = 0; e <= this->DG_degree_x; ++e ) {
    for ( int64_t s = 0; s <= this->DG_degree_t; ++s ) {
    for ( int64_t l = 0; l < SIMD_length; ++l ) {

        IWORK(0,d,e,s,l) = result( q+l, i + xoffset, j, d,e,s);
        IWORK(1,d,e,s,l) = result( q+l, i, j + yoffset, d,e,s);
    }}}}

    // Scale temporary values.
    for ( int64_t d = 0; d <= this->DG_degree_x; ++d ) {
    for ( int64_t e = 0; e <= this->DG_degree_x; ++e ) {
    for ( int64_t s = 0; s <= this->DG_degree_t; ++s ) {

        # pragma omp simd aligned( WORK : SWEEP_SIMD_LEN * SIZEOF_DOUBLE )
        for ( int64_t l = 0; l < SIMD_length; ++l )
            IWORK(0,d,e,s,l) *= this->dt * this->sw_op.dx(1) * this->sw_op.xi(q+l);

        # pragma omp simd aligned( WORK : SWEEP_SIMD_LEN * SIZEOF_DOUBLE )
        for ( int64_t l = 0; l < SIMD_length; ++l )
            IWORK(1,d,e,s,l) *= this->dt * this->sw_op.dx(0) * this->sw_op.eta(q+l);
    }}}

    // Set value from known source term.
    const double dx_coeff = this->dt * this->sw_op.dx(0) * this->sw_op.dx(1);

    for ( int64_t d = 0; d <= this->DG_degree_x; ++d ) {
    for ( int64_t e = 0; e <= this->DG_degree_x; ++e ) {
    for ( int64_t s = 0; s <= this->DG_degree_t; ++s ) {
    # pragma omp simd aligned( B : SWEEP_SIMD_LEN * SIZEOF_DOUBLE )
    for ( int64_t l = 0; l < SIMD_length; ++l ) {

        B[ this->template VecIdx<SIMD_length>(d,e,s,l) ] = dx_coeff * source( q+l, i,j,d,e,s);
    }}}}

    // Include known components of time derivative term.
    if ( initial != nullptr ) {

        for ( int64_t d = 0; d <= this->DG_degree_x; ++d ) {
        for ( int64_t e = 0; e <= this->DG_degree_x; ++e ) {
        for ( int64_t s = 0; s <= this->DG_degree_t; ++s ) {
        # pragma omp simd aligned( B : SWEEP_SIMD_LEN * SIZEOF_DOUBLE )
        for ( int64_t l = 0; l < SIMD_length; ++l ) {

            B[ this->template VecIdx<SIMD_length>(d,e,s,l) ]
                += neg1pow(s) * this->sw_op.dx(0) * this->sw_op.dx(1)* (2*s + 1) * (*initial)( q+l, i,j,d,e);
        }}}}
    }

    // Include upwind values in ξ angular direction.
    if ( all_xi_negative ) {

        for ( int64_t d = 0; d <= this->DG_degree_x; ++d ) {
        for ( int64_t e = 0; e <= this->DG_degree_x; ++e ) {
        for ( int64_t s = 0; s <= this->DG_degree_t; ++s ) {
        for ( int64_t a = 0; a <= this->DG_degree_x; ++a ) {

            const int64_t coeff = neg1pow(a+1) * (2*d + 1);

            # pragma omp simd aligned( B, WORK : SWEEP_SIMD_LEN * SIZEOF_DOUBLE )
            for ( int64_t l = 0; l < SIMD_length; ++l )
                B[ this->template VecIdx<SIMD_length>(d,e,s,l) ] += IWORK(0,a,e,s,l) * coeff;
        }}}}

    } else {

        for ( int64_t d = 0; d <= this->DG_degree_x; ++d ) {
        for ( int64_t e = 0; e <= this->DG_degree_x; ++e ) {
        for ( int64_t s = 0; s <= this->DG_degree_t; ++s ) {
        for ( int64_t a = 0; a <= this->DG_degree_x; ++a ) {

            const int64_t coeff = neg1pow(d) * (2*d + 1);

            # pragma omp simd aligned( B, WORK : SWEEP_SIMD_LEN * SIZEOF_DOUBLE )
            for ( int64_t l = 0; l < SIMD_length; ++l )
                B[ this->template VecIdx<SIMD_length>(d,e,s,l) ] += IWORK(0,a,e,s,l) * coeff;
        }}}}
    }

    // Include upwind values in η angular direction.
    if ( all_eta_negative ) {

        for ( int64_t d = 0; d <= this->DG_degree_x; ++d ) {
        for ( int64_t e = 0; e <= this->DG_degree_x; ++e ) {
        for ( int64_t s = 0; s <= this->DG_degree_t; ++s ) {
        for ( int64_t b = 0; b <= this->DG_degree_x; ++b ) {

            const int64_t coeff = neg1pow(b+1) * (2*e + 1);

            # pragma omp simd aligned( B, WORK : SWEEP_SIMD_LEN * SIZEOF_DOUBLE )
            for ( int64_t l = 0; l < SIMD_length; ++l )
                B[ this->template VecIdx<SIMD_length>(d,e,s,l) ] += IWORK(1,d,b,s,l) * coeff;
        }}}}

    } else {

        for ( int64_t d = 0; d <= this->DG_degree_x; ++d ) {
        for ( int64_t e = 0; e <= this->DG_degree_x; ++e ) {
        for ( int64_t s = 0; s <= this->DG_degree_t; ++s ) {
        for ( int64_t b = 0; b <= this->DG_degree_x; ++b ) {

            const int64_t coeff = neg1pow(e) * (2*e + 1);

            # pragma omp simd aligned( B, WORK : SWEEP_SIMD_LEN * SIZEOF_DOUBLE )
            for ( int64_t l = 0; l < SIMD_length; ++l )
                B[ this->template VecIdx<SIMD_length>(d,e,s,l) ] += IWORK(1,d,b,s,l) * coeff;
        }}}}
    }

# undef IWORK

# elif SPACE_DIMS == 3

    /*
     *  Spatial dimensions: 3
     *  SIMD blocking:      Yes
     *
     *  ---------------------------
     *  ----- NOT IMPLEMENTED -----
     *  ---------------------------
     */

    # warning "Implementation of SweepLinearSolveSpec<STDG::OrdinateFlux,SIMD_length>::CalcB incomplete."

    throw std::runtime_error( "Implementation of SweepLinearSolveSpec<STDG::OrdinateFlux,"
                              + std::to_string(SIMD_length) + ">::" + std::string(__func__)
                              + " incomplete." );

# endif // if SPACE_DIMS == ?

# if defined (DO_SWEEP_SUBTIMING)

    Global::TMR_calcB[
        # if defined (_OPENMP)
            omp_get_thread_num()
        # else
            0
        # endif
        ].Stop();

# endif // if defined (DO_SWEEP_SUBTIMING)
}


//============================================================================================================
//=== EXPLICIT INSTANTIATIONS ================================================================================
//============================================================================================================


# if SWEEP_SIMD_LEN >= (1 << 0)
    template class SweepLinearSolveSpec< STDG::OrdinateFlux, (1 << 0) >;
# endif
# if SWEEP_SIMD_LEN >= (1 << 1)
    template class SweepLinearSolveSpec< STDG::OrdinateFlux, (1 << 1) >;
# endif
# if SWEEP_SIMD_LEN >= (1 << 2)
    template class SweepLinearSolveSpec< STDG::OrdinateFlux, (1 << 2) >;
# endif
# if SWEEP_SIMD_LEN >= (1 << 3)
    template class SweepLinearSolveSpec< STDG::OrdinateFlux, (1 << 3) >;
# endif
# if SWEEP_SIMD_LEN >= (1 << 4)
    template class SweepLinearSolveSpec< STDG::OrdinateFlux, (1 << 4) >;
# endif
# if SWEEP_SIMD_LEN >= (1 << 5)
    template class SweepLinearSolveSpec< STDG::OrdinateFlux, (1 << 5) >;
# endif
# if SWEEP_SIMD_LEN > (1 << 5)
    # error "Maximum SIMD vector length supported is (1 << 5)."
# endif
