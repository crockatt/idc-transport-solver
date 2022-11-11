//------------------------------------------------------------------------------------------------------------
//! \file   linear_solvers/Abstract/SweepCellSolve/SweepLapackLUSolve.cpp
//! \brief  Contains implementation and instantiations of methods from SweepLapackLUSolve class template.
//!
//! \author Michael M. Crockatt
//! \date   March 2018
//------------------------------------------------------------------------------------------------------------

# include <typeinfo>

# include "linear_solvers/Abstract/SweepCellSolve/SweepLapackLUSolve.hpp"
# include "linear_solvers/LinearAlgebra.h"


//============================================================================================================
//=== MEMBER DEFINITIONS FOR SweepLapackLUSolve CLASS ========================================================
//============================================================================================================


//--------------------------------------------------------------------------------------------------------
//! \brief  Returns the string descriptor of the class type.
//--------------------------------------------------------------------------------------------------------
template< class OrdinateFlux, int64_t SIMD_length >
const std::string SweepLapackLUSolve< OrdinateFlux, SIMD_length >::Descriptor( void ) {

    return "lapack-lu";
}


//--------------------------------------------------------------------------------------------------------
//! \brief  Returns the string descriptor for an objects type.
//--------------------------------------------------------------------------------------------------------
template< class OrdinateFlux, int64_t SIMD_length >
const std::string SweepLapackLUSolve< OrdinateFlux, SIMD_length >::GetDescriptor( void ) const {

    return Descriptor();
}


//--------------------------------------------------------------------------------------------------------
//! \brief  Prints a summary of the object's configuration to the logging interface.
//!
//! \param[in]  prefix  String prepended to each line of output, intended to be used for controlling
//!                     indentation of output.
//--------------------------------------------------------------------------------------------------------
template< class OrdinateFlux, int64_t SIMD_length >
void SweepLapackLUSolve< OrdinateFlux, SIMD_length >::Print (

    const std::string & prefix // = "  "

) const {

    PRINT_STATUS( "Executing SweepLapackLUSolve<%s,%" PRId64 ">::%s.\n",
                  typeid(OrdinateFlux).name(), SIMD_length, __func__ )

    PRINT_LOG( "%s%-*s  %s\n", prefix.c_str(),
               Global::col_width, "Sweep Solve:", this->Descriptor().c_str() )
}


//--------------------------------------------------------------------------------------------------------
//! \brief  Virtual destructor.
//--------------------------------------------------------------------------------------------------------
template< class OrdinateFlux, int64_t SIMD_length >
SweepLapackLUSolve< OrdinateFlux, SIMD_length >::~SweepLapackLUSolve ( void ) {

    PRINT_STATUS( "Executing SweepLapackLUSolve<%s,%" PRId64 ">::%s.\n",
                  typeid(OrdinateFlux).name(), SIMD_length, __func__ )

    /* empty */
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
void SweepLapackLUSolve< OrdinateFlux, SIMD_length >::SolveSystem (

    double * const A,
    double * const B,
    void * const work_ptr

) const {

    PRINT_STATUS( "Executing SweepLapackLUSolve<%s,%" PRId64 ">::%s.\n",
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

    int N = this->system_dim;
    int NRHS = 1;
    int INFO;

    dgesv_( &N, &NRHS, A, &N, static_cast<int *>( work_ptr ), B, &N, &INFO );

# if defined (STRICT_CHECK)

    if ( INFO ) {

        std::string error_message = "dgesv_ returned INFO = " + std::to_string(INFO) + ".\n";

        PRINT_ERROR( error_message.c_str() )
        throw std::runtime_error( error_message );
    }

# endif // if defined (STRICT_CHECK)

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


//============================================================================================================
//=== EXPLICIT INSTANTIATIONS ================================================================================
//============================================================================================================

# include "objects/RKDG/OrdinateFlux.hpp"
# include "objects/STDG/OrdinateFlux.hpp"


# if SWEEP_SIMD_LEN >= (1 << 0)
    template class SweepLapackLUSolve< RKDG::OrdinateFlux, (1 << 0) >;
    template class SweepLapackLUSolve< STDG::OrdinateFlux, (1 << 0) >;

    template class SweepCellSolveDerivedFactory< RKDG::OrdinateFlux, (1 << 0), SweepLapackLUSolve >;
    template class SweepCellSolveDerivedFactory< STDG::OrdinateFlux, (1 << 0), SweepLapackLUSolve >;
# endif
# if SWEEP_SIMD_LEN >= (1 << 1)
    template class SweepLapackLUSolve< RKDG::OrdinateFlux, (1 << 1) >;
    template class SweepLapackLUSolve< STDG::OrdinateFlux, (1 << 1) >;

    template class SweepCellSolveDerivedFactory< RKDG::OrdinateFlux, (1 << 1), SweepLapackLUSolve >;
    template class SweepCellSolveDerivedFactory< STDG::OrdinateFlux, (1 << 1), SweepLapackLUSolve >;
# endif
# if SWEEP_SIMD_LEN >= (1 << 2)
    template class SweepLapackLUSolve< RKDG::OrdinateFlux, (1 << 2) >;
    template class SweepLapackLUSolve< STDG::OrdinateFlux, (1 << 2) >;

    template class SweepCellSolveDerivedFactory< RKDG::OrdinateFlux, (1 << 2), SweepLapackLUSolve >;
    template class SweepCellSolveDerivedFactory< STDG::OrdinateFlux, (1 << 2), SweepLapackLUSolve >;
# endif
# if SWEEP_SIMD_LEN >= (1 << 3)
    template class SweepLapackLUSolve< RKDG::OrdinateFlux, (1 << 3) >;
    template class SweepLapackLUSolve< STDG::OrdinateFlux, (1 << 3) >;

    template class SweepCellSolveDerivedFactory< RKDG::OrdinateFlux, (1 << 3), SweepLapackLUSolve >;
    template class SweepCellSolveDerivedFactory< STDG::OrdinateFlux, (1 << 3), SweepLapackLUSolve >;
# endif
# if SWEEP_SIMD_LEN >= (1 << 4)
    template class SweepLapackLUSolve< RKDG::OrdinateFlux, (1 << 4) >;
    template class SweepLapackLUSolve< STDG::OrdinateFlux, (1 << 4) >;

    template class SweepCellSolveDerivedFactory< RKDG::OrdinateFlux, (1 << 4), SweepLapackLUSolve >;
    template class SweepCellSolveDerivedFactory< STDG::OrdinateFlux, (1 << 4), SweepLapackLUSolve >;
# endif
# if SWEEP_SIMD_LEN >= (1 << 5)
    template class SweepLapackLUSolve< RKDG::OrdinateFlux, (1 << 5) >;
    template class SweepLapackLUSolve< STDG::OrdinateFlux, (1 << 5) >;

    template class SweepCellSolveDerivedFactory< RKDG::OrdinateFlux, (1 << 5), SweepLapackLUSolve >;
    template class SweepCellSolveDerivedFactory< STDG::OrdinateFlux, (1 << 5), SweepLapackLUSolve >;
# endif
# if SWEEP_SIMD_LEN > (1 << 5)
    # error "Maximum SIMD vector length supported is (1 << 5)."
# endif
