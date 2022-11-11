//------------------------------------------------------------------------------------------------------------
//! \file   linear_solvers/RKDG/SweepCellSolve/SweepCellSolve.hpp
//! \brief  Contains implementations of methods from partial specialization of SweepCellSolve class for
//!         RKDG::OrdinateFlux objects.
//!
//! \author Michael M. Crockatt
//! \date   March 2018
//------------------------------------------------------------------------------------------------------------

# include <typeinfo>

# include "linear_solvers/RKDG/SweepCellSolve/SweepCellSolve.hpp"
# include "objects/RKDG/OrdinateFlux.hpp"
# include "utils/global.hpp"


//-----------------------------------------------------------------------------------------------------------
//! \brief  Construct solver object from given parameters.
//!
//! SweepCellSolve::DG_degree is set based on the value determined by GetDGDegreeX.
//!
//! \param[in]  enclosing       Reference to instance of enclosing SweepOperator class.
//! \param[in]  input_list      Contains parameters for initialization of object.
//! \param[in]  dt_in           Initial timestep size.
//-----------------------------------------------------------------------------------------------------------
template< int64_t SIMD_length >
SweepCellSolve< RKDG::OrdinateFlux, SIMD_length >::SweepCellSolve (

    const SweepOperator<RKDG::OrdinateFlux> & enclosing,
    const ParameterList & input_list,
    const double dt_in              // = std::numeric_limits<double>::infinity()
) :
    SweepCellSolveBase<RKDG::OrdinateFlux>( enclosing, dt_in ),
    DG_degree{ GetDGDegreeX( input_list ) }
{
    PRINT_STATUS( "Executing SweepCellSolve<RKDG::OrdinateFlux,%" PRId64 ">::%s.\n", SIMD_length, __func__ )

    this->system_dim = IntPow( this->DG_degree + 1, SPACE_DIMS );
}


//------------------------------------------------------------------------------------------------------------
//! \brief  Virtual destructor.
//------------------------------------------------------------------------------------------------------------
template< int64_t SIMD_length >
SweepCellSolve< RKDG::OrdinateFlux, SIMD_length >::~SweepCellSolve ( void ) {

    PRINT_STATUS( "Executing SweepCellSolve<RKDG::OrdinateFlux,%" PRId64 ">::%s.\n", SIMD_length, __func__ )

    /* empty */
}


//------------------------------------------------------------------------------------------------------------
//! \brief  Stores the coefficients obtained from the solve into the output object.
//!
//! \param[in]      B       Pointer to memory location containing coefficients to store.
//! \param[out]     result  Object to store coefficients into.
//! \param[in]      idx     Contains indices for where to store coefficients.
//------------------------------------------------------------------------------------------------------------
template< int64_t SIMD_length >
void SweepCellSolve< RKDG::OrdinateFlux, SIMD_length >::StoreResult (

    const double * const B,
    const SIMD_BlkIdx<0> & idx,
    const RKDG::OrdinateFlux & result

) const {

    PRINT_STATUS( "Executing SweepCellSolve<RKDG::OrdinateFlux,%" PRId64 ">::%s.\n", SIMD_length, __func__ )

# if SPACE_DIMS == 1

    /*
     *  Spatial dimensions: 1
     *  SIMD blocking:      Yes
     */

    for ( int64_t d = 0; d <= this->DG_degree; ++d ) {
    # pragma omp simd aligned( B : SWEEP_SIMD_LEN * sizeof(double) )
    for ( int64_t l = 0; l < SIMD_length; ++l ) {

        result( idx.q + l, idx.i, d) = B[ this->VecIdx<SIMD_length>(d,l) ];
    }}

# elif SPACE_DIMS == 2

    /*
     *  Spatial dimensions: 2
     *  SIMD blocking:      Yes
     */

    for ( int64_t d = 0; d <= this->DG_degree; ++d ) {
    for ( int64_t e = 0; e <= this->DG_degree; ++e ) {
    # pragma omp simd aligned( B : SWEEP_SIMD_LEN * SIZEOF_DOUBLE )
    for ( int64_t l = 0; l < SIMD_length; ++l ) {

        result( idx.q + l, idx.i, idx.j, d,e) = B[ this->VecIdx<SIMD_length>(d,e,l) ];
    }}}

# elif SPACE_DIMS == 3

    /*
     *  Spatial dimensions: 3
     *  SIMD blocking:      Yes
     */

    for ( int64_t d = 0; d <= this->DG_degree; ++d ) {
    for ( int64_t e = 0; e <= this->DG_degree; ++e ) {
    for ( int64_t f = 0; f <= this->DG_degree; ++f ) {
    # pragma omp simd aligned( B : SWEEP_SIMD_LEN * sizeof(double) )
    for ( int64_t l = 0; l < SIMD_length; ++l ) {

        result( idx.q + l, idx.i, idx.j, idx.k, d,e,f) = B[ this->VecIdx<SIMD_length>(d,e,f,l) ];
    }}}}

# endif // if SPACE_DIMS == ?
}


//============================================================================================================
//=== EXPLICIT INSTANTIATIONS ================================================================================
//============================================================================================================


# if SWEEP_SIMD_LEN >= (1 << 0)
    template class SweepCellSolve< RKDG::OrdinateFlux, (1 << 0) >;
# endif
# if SWEEP_SIMD_LEN >= (1 << 1)
    template class SweepCellSolve< RKDG::OrdinateFlux, (1 << 1) >;
# endif
# if SWEEP_SIMD_LEN >= (1 << 2)
    template class SweepCellSolve< RKDG::OrdinateFlux, (1 << 2) >;
# endif
# if SWEEP_SIMD_LEN >= (1 << 3)
    template class SweepCellSolve< RKDG::OrdinateFlux, (1 << 3) >;
# endif
# if SWEEP_SIMD_LEN >= (1 << 4)
    template class SweepCellSolve< RKDG::OrdinateFlux, (1 << 4) >;
# endif
# if SWEEP_SIMD_LEN >= (1 << 5)
    template class SweepCellSolve< RKDG::OrdinateFlux, (1 << 5) >;
# endif
# if SWEEP_SIMD_LEN > (1 << 5)
    # error "Maximum SIMD vector length supported is (1 << 5)."
# endif
