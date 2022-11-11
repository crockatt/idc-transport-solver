//------------------------------------------------------------------------------------------------------------
//! \file   linear_solvers/Abstract/SweepCellSolve/SweepLinearSolveBase.cpp
//! \brief  Contains implementations and instantiations of methods from SweepLinearSolveBase class.
//!
//! \author Michael M. Crockatt
//! \date   March 2018
//------------------------------------------------------------------------------------------------------------

# include <typeinfo>

# include "linear_solvers/Abstract/SweepCellSolve/SweepLinearSolveBase.hpp"


//------------------------------------------------------------------------------------------------------------
//! \brief  Construct solver object from given parameters.
//!
//! \param[in]  enclosing       Reference to instance of enclosing SweepOperator class.
//! \param[in]  input_list      Contains parameters for initialization of object.
//! \param[in]  dt_in           Initial timestep size.
//------------------------------------------------------------------------------------------------------------
template< class OrdinateFlux, int64_t SIMD_length >
SweepLinearSolveBase< OrdinateFlux, SIMD_length >::SweepLinearSolveBase (

    const SweepOperator<OrdinateFlux> & enclosing,
    const ParameterList & input_list,
    const double dt_in              // = std::numeric_limits<double>::infinity()
) :
    SweepCellSolve< OrdinateFlux, SIMD_length >( enclosing, input_list, dt_in ),

    Aq{ nullptr },

    dimof_A{ this->system_dim * this->system_dim },
    dimof_B{ this->system_dim },

    sizeof_A{ this->dimof_A * sizeof(double) },
    sizeof_B{ this->dimof_B * sizeof(double) },

    force_zwc{ false }
{
    PRINT_STATUS( "Executing SweepLinearSolveBase<%s,%" PRId64 ">::%s.\n",
                  typeid(OrdinateFlux).name(), SIMD_length, __func__ )

    // Allocate memory for precomputed components of matrices.
    this->Aq = new double*[ this->sw_op.nq() ];

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

            this->Aq[q] = (double *)
                # if defined (USE_ALIGNED_ALLOC)
                    aligned_alloc( ALIGNED_ALLOC_ALIGNMENT, SIMD_length * this->sizeof_A );
                # else
                    std::malloc( SIMD_length * this->sizeof_A );
                # endif

            for ( int64_t l = 1; l < SIMD_length; ++l )
                if ( q+l < q_max )
                    this->Aq[q+l] = nullptr;
        }
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

            this->Aq[q] = (double *)
                # if defined (USE_ALIGNED_ALLOC)
                    aligned_alloc( ALIGNED_ALLOC_ALIGNMENT, SIMD_length * this->sizeof_A );
                # else
                    std::malloc( SIMD_length * this->sizeof_A );
                # endif

            for ( int64_t l = 1; l < SIMD_length; ++l )
                if ( q+l < q_max )
                    this->Aq[q+l] = nullptr;
        }
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

    # warning "Implementation of SweepLinearSolveBase<OrdinateFlux,SIMD_length>::SweepLinearSolveBase incomplete."

    throw std::runtime_error( "Implementation of SweepLinearSolveBase<"
                              + std::string( typeid(OrdinateFlux).name() ) + "," + std::to_string(SIMD_length)
                              + ">::" + std::string(__func__) + " incomplete." );

# endif // if SPACE_DIMS == ?
}


//------------------------------------------------------------------------------------------------------------
//! \brief  Virtual destructor.
//------------------------------------------------------------------------------------------------------------
template< class OrdinateFlux, int64_t SIMD_length >
SweepLinearSolveBase< OrdinateFlux, SIMD_length >::~SweepLinearSolveBase ( void ) {

    PRINT_STATUS( "Executing SweepLinearSolveBase<%s,%" PRId64 ">::%s.\n",
                  typeid(OrdinateFlux).name(), SIMD_length, __func__ )

    if ( this->Aq != nullptr ) {

        for ( int64_t q = 0; q < this->sw_op.nq(); ++q ) {

            if ( this->Aq[q] != nullptr )
                std::free( this->Aq[q] );
        }

        delete [] this->Aq;
        this->Aq = nullptr;
    }
}


//------------------------------------------------------------------------------------------------------------
//! \brief  Returns the dimension of the array to be allocated for A.
//------------------------------------------------------------------------------------------------------------
template< class OrdinateFlux, int64_t SIMD_length >
size_t SweepLinearSolveBase< OrdinateFlux, SIMD_length >::GetDimofA ( void ) const {

    PRINT_STATUS( "Executing SweepLinearSolveBase<%s,%" PRId64 ">::%s.\n",
                  typeid(OrdinateFlux).name(), SIMD_length, __func__ )

    return this->dimof_A;
}


//------------------------------------------------------------------------------------------------------------
//! \brief  Returns the dimension of the array to be allocated for B.
//------------------------------------------------------------------------------------------------------------
template< class OrdinateFlux, int64_t SIMD_length >
size_t SweepLinearSolveBase< OrdinateFlux, SIMD_length >::GetDimofB ( void ) const {

    PRINT_STATUS( "Executing SweepLinearSolveBase<%s,%" PRId64 ">::%s.\n",
                  typeid(OrdinateFlux).name(), SIMD_length, __func__ )

    return this->dimof_B;
}


//------------------------------------------------------------------------------------------------------------
//! \brief  Returns the dimension of the array to be allocated for extra work space.
//------------------------------------------------------------------------------------------------------------
template< class OrdinateFlux, int64_t SIMD_length >
size_t SweepLinearSolveBase< OrdinateFlux, SIMD_length >::GetDimofWork ( void ) const {

    PRINT_STATUS( "Executing SweepLinearSolveBase<%s,%" PRId64 ">::%s.\n",
                  typeid(OrdinateFlux).name(), SIMD_length, __func__ )

    return 2 * this->dimof_B;
}


//------------------------------------------------------------------------------------------------------------
//! \brief  Zeros the arrays for storing precomputed values (e.g., angular component of matrices).
//------------------------------------------------------------------------------------------------------------
template< class OrdinateFlux, int64_t SIMD_length >
void SweepLinearSolveBase< OrdinateFlux, SIMD_length >::ZeroMatrices ( void ) const {

    PRINT_STATUS( "Executing SweepLinearSolveBase<%s,%" PRId64 ">::%s.\n",
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

        # pragma omp parallel for
        for ( int64_t q  = q_min; q < q_max; q += SIMD_length )
            std::memset( this->Aq[q], 0, SIMD_length * this->sizeof_A );
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
            std::memset( this->Aq[q], 0, SIMD_length * this->sizeof_A );
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

    # warning "Implementation of SweepLinearSolveBase<RKDG::OrdinateFlux>::ComputeMatrices incomplete."

    throw std::runtime_error( "Implementation of SweepLinearSolveBase<"
                              + std::string( typeid(OrdinateFlux).name() ) + "," + std::to_string(SIMD_length)
                              + ">::" + std::string(__func__) + " incomplete." );

# endif // if SPACE_DIMS == ?
}


//============================================================================================================
//=== EXPLICIT INSTANTIATIONS ================================================================================
//============================================================================================================

# include "objects/RKDG/OrdinateFlux.hpp"
# include "objects/STDG/OrdinateFlux.hpp"


# if SWEEP_SIMD_LEN >= (1 << 0)
    template class SweepLinearSolveBase< RKDG::OrdinateFlux, (1 << 0) >;
    template class SweepLinearSolveBase< STDG::OrdinateFlux, (1 << 0) >;
# endif
# if SWEEP_SIMD_LEN >= (1 << 1)
    template class SweepLinearSolveBase< RKDG::OrdinateFlux, (1 << 1) >;
    template class SweepLinearSolveBase< STDG::OrdinateFlux, (1 << 1) >;
# endif
# if SWEEP_SIMD_LEN >= (1 << 2)
    template class SweepLinearSolveBase< RKDG::OrdinateFlux, (1 << 2) >;
    template class SweepLinearSolveBase< STDG::OrdinateFlux, (1 << 2) >;
# endif
# if SWEEP_SIMD_LEN >= (1 << 3)
    template class SweepLinearSolveBase< RKDG::OrdinateFlux, (1 << 3) >;
    template class SweepLinearSolveBase< STDG::OrdinateFlux, (1 << 3) >;
# endif
# if SWEEP_SIMD_LEN >= (1 << 4)
    template class SweepLinearSolveBase< RKDG::OrdinateFlux, (1 << 4) >;
    template class SweepLinearSolveBase< STDG::OrdinateFlux, (1 << 4) >;
# endif
# if SWEEP_SIMD_LEN >= (1 << 5)
    template class SweepLinearSolveBase< RKDG::OrdinateFlux, (1 << 5) >;
    template class SweepLinearSolveBase< STDG::OrdinateFlux, (1 << 5) >;
# endif
# if SWEEP_SIMD_LEN > (1 << 5)
    # error "Maximum SIMD vector length supported is (1 << 5)."
# endif
