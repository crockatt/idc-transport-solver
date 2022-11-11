//------------------------------------------------------------------------------------------------------------
//! \file   sweep_3d.c
//! \brief  Contains the implementations of 3D sweep algorithms for inverting the discrete ordinate
//!         advection-absorption operator \f$ \mathcal{L} \f$.
//!
//! Implements implicit transport sweeps for inverting the discrete ordinate advection-absorption operator
//! \f[ \mathcal{L} = \vec{\Omega} \cdot \nabla_{\vec{x}} + \sigma + \frac{1}{ \Delta t } \f]
//! in full three-dimensional space.
//! Two methods are used for solving the required local systems within the transport sweep:
//!     - _LU Solve:_
//!         The local linear systems are solved using a standard LU decomposition.
//!         Cross-sections are assumed to be piecewise constant.
//!     - _Hessenberg Solve:_
//!         The matrix for each local linear system can be divided into a sum of two components: one depending
//!         only on the spatial cell and one only depending on the angular ordinate. When the material
//!         cross-sections are constant on a spatial cell, the spatially dependent component is diagonal.
//!         By pre-computing a Hessenberg reduction of the angularly-dependent components, the local linear
//!         systems can be re-written as an upper Hessenberg system. This reduces the computational cost of
//!         each local solve from \f$ O(N^3) \f$ to \f$ O(N^2) \f$.
//!
//! \author Michael Crockatt
//! \date   November 2016
//------------------------------------------------------------------------------------------------------------


# include <inttypes.h>
# include <stdint.h>
# include <stdio.h>
# include <stdlib.h>
# include <string.h>
# include <unistd.h>

# if SPACE_DIMS == 3 || defined (DOXYCOMPILE)

# include "linear_solvers/linearAlgebra.h"
# include "linear_solvers/sweep.h"
# include "misc/global.hpp"
# include "structs/angularFlux.h"


// PETSc defines these as macros. (Their macro expansion causes compiler warnings, which I don't like.)
# if defined (MPI_Isend)
    # undef MPI_Isend
# endif
# if defined (MPI_Irecv)
    # undef MPI_Irecv
# endif
# if defined (MPI_Wait)
    # undef MPI_Wait
# endif
# if defined (MPI_Waitall)
    # undef MPI_Waitall
# endif


# if defined (MIN)
    # undef MIN
# endif
# define MIN(a,b) ( (a) < (b) ? (a) : (b) )

# if defined (MAX)
    # undef MAX
# endif
# define MAX(a,b) ( (a) > (b) ? (a) : (b) )


# if defined (USE_ALIGNED_ALLOC) && ! defined (ALIGNED_ALLOC_ALIGNMENT)
    # error "Macro ALIGNED_ALLOC_ALIGNMENT not defined!"
# endif


//------------------------------------------------------------------------------------------------------------
//! \brief  Macro for indexing array of matrices used for 3D sweeps in column-major format.
//!
//! Constructs the appropriate indexing format for accessing elements of the matrices stored in
//! \pp{context->A} in column-major format.
//!
//! \param[in]      context     sweepContext_t object containing the context of the current sweep.
//! \param[in]      q           Index of angular ordinate.
//! \param[in]      a           Degree of test function with respect to \f$ x \f$.
//! \param[in]      b           Degree of test function with respect to \f$ y \f$.
//! \param[in]      c           Degree of test function with respect to \f$ z \f$.
//! \param[in]      d           Degree of basis function with respect to \f$ x \f$.
//! \param[in]      e           Degree of basis function with respect to \f$ y \f$.
//! \param[in]      f           Degree of basis function with respect to \f$ z \f$.
//!
//! \see    SWC_matrixIndexBlock_columnMajor()
//------------------------------------------------------------------------------------------------------------
# define IA_BLOCK_COL(context,q,a,b,c,d,e,f) \
    ((context)->A[ SWC_matrixIndexBlock_columnMajor(context,q,a,b,c,d,e,f) ])


//------------------------------------------------------------------------------------------------------------
//! \brief  Macro for indexing array of matrices used for 3D sweeps in row-major format.
//!
//! Constructs the appropriate indexing format for accessing elements of the matrices stored in
//! \pp{context->A} in row-major format.
//!
//! \param[in]      context     sweepContext_t object containing the context of the current sweep.
//! \param[in]      q           Index of angular ordinate.
//! \param[in]      a           Degree of test function with respect to \f$ x \f$.
//! \param[in]      b           Degree of test function with respect to \f$ y \f$.
//! \param[in]      c           Degree of test function with respect to \f$ z \f$.
//! \param[in]      d           Degree of basis function with respect to \f$ x \f$.
//! \param[in]      e           Degree of basis function with respect to \f$ y \f$.
//! \param[in]      f           Degree of basis function with respect to \f$ z \f$.
//!
//! \see    SWC_matrixIndexBlock_rowMajor()
//------------------------------------------------------------------------------------------------------------
# define IA_BLOCK_ROW(context,q,a,b,c,d,e,f) \
    ((context)->A[ SWC_matrixIndexBlock_rowMajor(context,q,a,b,c,d,e,f) ])


//------------------------------------------------------------------------------------------------------------
//! \brief  Macro for indexing array of angular component matrices used for 3D sweeps in column-major format.
//!
//! Constructs the appropriate indexing format for accessing elements of the matrices stored in
//! \pp{context->Aq} in column-major format.
//!
//! \param[in]      context     sweepContext_t object containing the context of the current sweep.
//! \param[in]      q           Index of angular ordinate.
//! \param[in]      a           Degree of test function with respect to \f$ x \f$.
//! \param[in]      b           Degree of test function with respect to \f$ y \f$.
//! \param[in]      c           Degree of test function with respect to \f$ z \f$.
//! \param[in]      d           Degree of basis function with respect to \f$ x \f$.
//! \param[in]      e           Degree of basis function with respect to \f$ y \f$.
//! \param[in]      f           Degree of basis function with respect to \f$ z \f$.
//!
//! \see    SWC_matrixIndexBlock_columnMajor()
//------------------------------------------------------------------------------------------------------------
# define IAQ_BLOCK_COL(context,q,a,b,c,d,e,f) \
    ((context)->Aq[ SWC_matrixIndexBlock_columnMajor(context,q,a,b,c,d,e,f) ])


//------------------------------------------------------------------------------------------------------------
//! \brief  Macro for indexing array of angular component matrices used for 3D sweeps in row-major format.
//!
//! Constructs the appropriate indexing format for accessing elements of the matrices stored in
//! \pp{context->Aq} in row-major format.
//!
//! \param[in]      context     sweepContext_t object containing the context of the current sweep.
//! \param[in]      q           Index of angular ordinate.
//! \param[in]      a           Degree of test function with respect to \f$ x \f$.
//! \param[in]      b           Degree of test function with respect to \f$ y \f$.
//! \param[in]      c           Degree of test function with respect to \f$ z \f$.
//! \param[in]      d           Degree of basis function with respect to \f$ x \f$.
//! \param[in]      e           Degree of basis function with respect to \f$ y \f$.
//! \param[in]      f           Degree of basis function with respect to \f$ z \f$.
//!
//! \see    SWC_matrixIndexBlock_rowMajor()
//------------------------------------------------------------------------------------------------------------
# define IAQ_BLOCK_ROW(context,q,a,b,c,d,e,f) \
    ((context)->Aq[ SWC_matrixIndexBlock_rowMajor(context,q,a,b,c,d,e,f) ])


//------------------------------------------------------------------------------------------------------------
//! \brief  Macro for indexing angular component matrices used for 3D wavefront sweeps in column-major format.
//!
//! Constructs the appropriate indexing format for accessing elements of the matrices stored in
//! \pp{context->Aq_wave}, which are stored non-contiguously in memory. Matrices are indexed using
//! column-major format.
//!
//! \param[in]      context     sweepContext_t object containing the context of the current sweep.
//! \param[in]      q           Index of angular ordinate.
//! \param[in]      a           Degree of test function with respect to \f$ x \f$.
//! \param[in]      b           Degree of test function with respect to \f$ y \f$.
//! \param[in]      c           Degree of test function with respect to \f$ z \f$.
//! \param[in]      d           Degree of basis function with respect to \f$ x \f$.
//! \param[in]      e           Degree of basis function with respect to \f$ y \f$.
//! \param[in]      f           Degree of basis function with respect to \f$ z \f$.
//!
//! \see    SWC_matrixIndexSingle_columnMajor()
//------------------------------------------------------------------------------------------------------------
# define IAQ_NONCONTIG_COL(context,q,a,b,c,d,e,f) \
    ((context)->Aq_wave[(q)][ SWC_matrixIndexSingle_columnMajor(context,a,b,c,d,e,f) ])


//------------------------------------------------------------------------------------------------------------
//! \brief  Macro for indexing angular component matrices used for 3D wavefront sweeps in row-major format.
//!
//! Constructs the appropriate indexing format for accessing elements of the matrices stored in
//! \pp{context->Aq_wave}, which are stored non-contiguously in memory. Matrices are indexed using row-major
//! format.
//!
//! \param[in]      context     sweepContext_t object containing the context of the current sweep.
//! \param[in]      q           Index of angular ordinate.
//! \param[in]      a           Degree of test function with respect to \f$ x \f$.
//! \param[in]      b           Degree of test function with respect to \f$ y \f$.
//! \param[in]      c           Degree of test function with respect to \f$ z \f$.
//! \param[in]      d           Degree of basis function with respect to \f$ x \f$.
//! \param[in]      e           Degree of basis function with respect to \f$ y \f$.
//! \param[in]      f           Degree of basis function with respect to \f$ z \f$.
//!
//! \see    SWC_matrixIndexSingle_rowMajor()
//------------------------------------------------------------------------------------------------------------
# define IAQ_NONCONTIG_ROW(context,q,a,b,c,d,e,f) \
    ((context)->Aq_wave[(q)][ SWC_matrixIndexSingle_rowMajor(context,a,b,c,d,e,f) ])


//------------------------------------------------------------------------------------------------------------
//! \brief  Macro for indexing array of vectors used for 3D sweeps.
//!
//! Constructs the appropriate indexing format for accessing elements of the vectors stored in \pp{context-B}.
//!
//! \param[in]      context     sweepContext_t object containing the context of the current sweep.
//! \param[in]      q           Index of angular ordinate.
//! \param[in]      d           Degree of test function with respect to \f$ x \f$.
//! \param[in]      e           Degree of test function with respect to \f$ y \f$.
//! \param[in]      f           Degree of test function with respect to \f$ z \f$.
//!
//! \see    SWC_vectorIndexBlock()
//------------------------------------------------------------------------------------------------------------
# define IB_BLOCK(context,q,d,e,f) ((context)->B[ SWC_vectorIndexBlock(context,q,d,e,f) ])


//------------------------------------------------------------------------------------------------------------
//! \brief  Column-major indexing function for array of matrices used for 3D sweeps.
//!
//! Column-major indexing function for matrices used by the sweep algorithms. Assumes that the matrix to be
//! indexed is stored within an array of matrices stored contiguously in memory, with one matrix per angular
//! ordinate arranged in order of increasing ordinate index.
//!
//! \param[in]      context     sweepContext_t object containing the context of the current sweep.
//! \param[in]      q           Integral value corresponding to index of angular ordinate.
//! \param[in]      a           Degree of test function with respect to \f$ x \f$.
//! \param[in]      b           Degree of test function with respect to \f$ y \f$.
//! \param[in]      c           Degree of test function with respect to \f$ z \f$.
//! \param[in]      d           Degree of basis function with respect to \f$ x \f$.
//! \param[in]      e           Degree of basis function with respect to \f$ y \f$.
//! \param[in]      f           Degree of basis function with respect to \f$ z \f$.
//!
//! \see    SWC_vectorIndexBlock()
//------------------------------------------------------------------------------------------------------------
static inline size_t SWC_matrixIndexBlock_columnMajor (

    const sweepContext_t context,
    const int64_t q,
    const int64_t a,
    const int64_t b,
    const int64_t c,
    const int64_t d,
    const int64_t e,
    const int64_t f
) {
    return (
        c + (context->dgOrder + 1)*(
            b + (context->dgOrder + 1)*(
                a + (context->dgOrder + 1)*(
                    f + (context->dgOrder + 1)*(
                        e + (context->dgOrder + 1)*(
                            d + (context->dgOrder + 1)*(q)
                        )
                    )
                )
            )
        )
    );
}


//------------------------------------------------------------------------------------------------------------
//! \brief  Row-major indexing function for array of matrices used for 3D sweeps.
//!
//! Row-major indexing function for matrices used by the sweep algorithms. Assumes that the matrix to be
//! indexed is stored within an array of matrices stored contiguously in memory, with one matrix per angular
//! ordinate arranged in order of increasing ordinate index.
//!
//! \param[in]      context     sweepContext_t object containing the context of the current sweep.
//! \param[in]      q           Integral value corresponding to index of angular ordinate.
//! \param[in]      a           Degree of test function with respect to \f$ x \f$.
//! \param[in]      b           Degree of test function with respect to \f$ y \f$.
//! \param[in]      c           Degree of test function with respect to \f$ z \f$.
//! \param[in]      d           Degree of basis function with respect to \f$ x \f$.
//! \param[in]      e           Degree of basis function with respect to \f$ y \f$.
//! \param[in]      f           Degree of basis function with respect to \f$ z \f$.
//!
//! \see    SWC_vectorIndexBlock()
//------------------------------------------------------------------------------------------------------------
static inline size_t SWC_matrixIndexBlock_rowMajor (

    const sweepContext_t context,
    const int64_t q,
    const int64_t a,
    const int64_t b,
    const int64_t c,
    const int64_t d,
    const int64_t e,
    const int64_t f
) {
    return (
        f + (context->dgOrder + 1)*(
            e + (context->dgOrder + 1)*(
                d + (context->dgOrder + 1)*(
                    c + (context->dgOrder + 1)*(
                        b + (context->dgOrder + 1)*(
                            a + (context->dgOrder + 1)*(q)
                        )
                    )
                )
            )
        )
    );
}


//------------------------------------------------------------------------------------------------------------
//! \brief  Column-major indexing function for a single matrix used in the 3D sweep algorithms.
//!
//! Column-major indexing function for matrices used by the sweep algorithms. Indexes into the matrix assuming
//! it is stored as a single matrix.
//!
//! \param[in]      context     sweepContext_t object containing the context of the current sweep.
//! \param[in]      a           Degree of test function with respect to \f$ x \f$.
//! \param[in]      b           Degree of test function with respect to \f$ y \f$.
//! \param[in]      c           Degree of test function with respect to \f$ z \f$.
//! \param[in]      d           Degree of basis function with respect to \f$ x \f$.
//! \param[in]      e           Degree of basis function with respect to \f$ y \f$.
//! \param[in]      f           Degree of basis function with respect to \f$ z \f$.
//!
//! \see    SWC_vectorIndexSingle()
//! \see    SWC_matrixIndexSingle_rowMajor()
//------------------------------------------------------------------------------------------------------------
static inline size_t SWC_matrixIndexSingle_columnMajor (

    const sweepContext_t context,
    const int64_t a,
    const int64_t b,
    const int64_t c,
    const int64_t d,
    const int64_t e,
    const int64_t f
) {
    return (
        c + (context->dgOrder + 1)*(
            b + (context->dgOrder + 1)*(
                a + (context->dgOrder + 1)*(
                    f + (context->dgOrder + 1)*(
                        e + (context->dgOrder + 1)*(
                            d
                        )
                    )
                )
            )
        )
    );
}


//------------------------------------------------------------------------------------------------------------
//! \brief  Row-major indexing function for a single matrix used in the 3D sweep algorithms.
//!
//! Row-major indexing function for matrices used by the sweep algorithms. Indexes into the matrix assuming
//! it is stored as a single matrix.
//!
//! \param[in]      context     sweepContext_t object containing the context of the current sweep.
//! \param[in]      a           Degree of test function with respect to \f$ x \f$.
//! \param[in]      b           Degree of test function with respect to \f$ y \f$.
//! \param[in]      c           Degree of test function with respect to \f$ z \f$.
//! \param[in]      d           Degree of basis function with respect to \f$ x \f$.
//! \param[in]      e           Degree of basis function with respect to \f$ y \f$.
//! \param[in]      f           Degree of basis function with respect to \f$ z \f$.
//!
//! \see    SWC_vectorIndexSingle()
//! \see    SWC_matrixIndexSingle_columnMajor()
//------------------------------------------------------------------------------------------------------------
static inline size_t SWC_matrixIndexSingle_rowMajor (

    const sweepContext_t context,
    const int64_t a,
    const int64_t b,
    const int64_t c,
    const int64_t d,
    const int64_t e,
    const int64_t f
) {
    return (
        f + (context->dgOrder + 1)*(
            e + (context->dgOrder + 1)*(
                d + (context->dgOrder + 1)*(
                    c + (context->dgOrder + 1)*(
                        b + (context->dgOrder + 1)*(
                            a
                        )
                    )
                )
            )
        )
    );
}


//------------------------------------------------------------------------------------------------------------
//! \brief  Indexing function for array of vectors used for 3D sweeps.
//!
//! Indexing function for vectors used by the sweep algorithms. Assumes that the vector to be indexed is
//! stored within an array of vectors stored contiguously in memory, with one vector per angular ordinate
//! arranged in order of increasing ordinate index.
//!
//! \param[in]      context     sweepContext_t object containing the context of the current sweep.
//! \param[in]      q           Index of angular ordinate.
//! \param[in]      a           Degree of test function with respect to \f$ x \f$.
//! \param[in]      b           Degree of test function with respect to \f$ y \f$.
//! \param[in]      c           Degree of test function with respect to \f$ z \f$.
//!
//! \see    SWC_matrixIndexBlock_columnMajor()
//! \see    SWC_matrixIndexBlock_rowMajor()
//------------------------------------------------------------------------------------------------------------
static inline size_t SWC_vectorIndexBlock (

    const sweepContext_t context,
    const int64_t q,
    const int64_t a,
    const int64_t b,
    const int64_t c
) {
    return (
        c + (context->dgOrder + 1)*(
            b + (context->dgOrder + 1)*(
                a + (context->dgOrder + 1)*(
                    q
                )
            )
        )
    );
}


//------------------------------------------------------------------------------------------------------------
//! \brief  Indexing function for array of vectors used in the 3D sweep algorithms.
//!
//! Indexing function for vectors used by the sweep algorithms. Indexes into the vector assuming it is stored
//! as a single vector.
//!
//! \param[in]      context     sweepContext_t object containing the context of the current sweep.
//! \param[in]      a           Degree of test function with respect to \f$ x \f$.
//! \param[in]      b           Degree of test function with respect to \f$ y \f$.
//! \param[in]      c           Degree of test function with respect to \f$ z \f$.
//!
//! \see    SWC_matrixIndexSingle_columnMajor()
//! \see    SWC_matrixIndexSingle_rowMajor()
//------------------------------------------------------------------------------------------------------------
static inline size_t SWC_vectorIndexSingle (

    const sweepContext_t context,
    const int64_t a,
    const int64_t b,
    const int64_t c
) {
    return (
        c + (context->dgOrder + 1)*(
            b + (context->dgOrder + 1)*(
                a
            )
        )
    );
}


//------------------------------------------------------------------------------------------------------------
//! \brief  Initialize sweepContext_t object for use.
//!
//! Sets the parameters for the sweep based on the parameters from the given angularFlux_t object.
//!
//! \param[out]     context     sweepContext_t object to initialize.
//! \param[in]      type        sweepType_t specifying the sweep algorithm to use. If this is SWEEP_NONE
//!                             then g_theInputList is searched for a valid sweepType_t using key
//!                             \e "sweep_type".
//! \param[in]      dgOrder     Maximum degree of DG basis polynomials.
//! \param[in]      spatParams  Contains the parameters of the spatial discretization.
//! \param[in]      ordinateSet ordinateSet_t value specifying type of discrete ordinate quadrature to use.
//! \param[in]      angOrder    Order of discrete ordinate quadrature set of type \pp{ordinateSet} to use.
//!
//! \see    SWC_create()
//! \see    SWC_destroy()
//! \see    SWC_cleanup()
//------------------------------------------------------------------------------------------------------------
void SWC_setup (

    const sweepContext_t context,
    const sweepType_t type,
    const int64_t dgOrder,
    const spatialParams_t * const spatParams,
    const ordinateSet_t ordinateSet,
    const int64_t angOrder
) {

    SWC_cleanup( context );

    // Set values from parameters provided.
    context->dgOrder = dgOrder;

    memcpy( context, spatParams, sizeof(spatialParams_t) );

    context->ordinateSet = ordinateSet;
    context->angOrder = angOrder;

    if ( type == SWEEP_NONE )
        IL_getValue_swt( g_theInputList, "sweep_type", &context->type );
    else
        context->type = type;

    computeOrdinates( ordinateSet, angOrder, &context->nq,
                      &context->xi, &context->eta, &context->mu, &context->w );

    switch ( context->type ) {

    case SWEEP_SEQUENTIAL_LAPACK_LU:
    {
        // Size of data structures containing linear systems.
        context->sizeofA = sizeof(double) * context->nq
                           * (context->dgOrder + 1)*(context->dgOrder + 1)*(context->dgOrder + 1)
                           * (context->dgOrder + 1)*(context->dgOrder + 1)*(context->dgOrder + 1);

        context->sizeofB = sizeof(double) * context->nq
                           * (context->dgOrder + 1)* (context->dgOrder + 1)* (context->dgOrder + 1);

        // Allocate memory for linear systems.
        context->A  = malloc( context->sizeofA );
        context->Aq = malloc( context->sizeofA );
        context->B  = malloc( context->sizeofB );

        // Precompute Aq.
        memset( context->Aq, 0, context->sizeofA );

        # pragma omp parallel for
        for ( int64_t q = 0; q < context->nq; ++q ) {

            // Compute A_q.
            for ( int64_t d = 0; d <= context->dgOrder; ++d ) {
            for ( int64_t e = 0; e <= context->dgOrder; ++e ) {
            for ( int64_t f = 0; f <= context->dgOrder; ++f ) {

                for ( int64_t a = d-1; a >= 0; a -= 2 )
                    IAQ_BLOCK_COL(context,q,d,e,f,a,e,f)
                        -= 2.0 * context->xi[q] * context->dy * context->dz / ( (2*e + 1)*(2*f + 1) );

                for ( int64_t b = e-1; b >= 0; b -= 2 )
                    IAQ_BLOCK_COL(context,q,d,e,f,d,b,f)
                        -= 2.0 * context->eta[q] * context->dx * context->dz / ( (2*d + 1)*(2*f + 1) );

                for ( int64_t c = f-1; c >= 0; c -= 2 )
                    IAQ_BLOCK_COL(context,q,d,e,f,d,e,c)
                        -= 2.0 * context->mu[q] * context->dx * context->dy / ( (2*d + 1)*(2*e + 1) );

                // Components from advection term for ξ angular direction.
                if ( context->xi[q] < 0 ) {

                    for ( int64_t a = 0; a <= context->dgOrder; ++a )
                        IAQ_BLOCK_COL(context,q,d,e,f,a,e,f)
                            += neg1pow(d+a+1) * context->xi[q] * context->dy * context->dz
                               / ( (2*e + 1)*(2*f + 1) );
                } else {

                    for ( int64_t a = 0; a <= context->dgOrder; ++a )
                        IAQ_BLOCK_COL(context,q,d,e,f,a,e,f)
                            += context->xi[q] * context->dy * context->dz / ( (2*e + 1)*(2*f + 1) );
                }

                // Components from advection term for η angular direction.
                if ( context->eta[q] < 0 ) {

                    for ( int64_t b = 0; b <= context->dgOrder; ++b )
                        IAQ_BLOCK_COL(context,q,d,e,f,d,b,f)
                            += neg1pow(e+b+1) * context->eta[q] * context->dx * context->dz
                               / ( (2*d + 1)*(2*f + 1) );
                } else {

                    for ( int64_t b = 0; b <= context->dgOrder; ++b )
                        IAQ_BLOCK_COL(context,q,d,e,f,d,b,f)
                            += context->eta[q] * context->dx * context->dz / ( (2*d + 1)*(2*f + 1) );
                }

                // Components from advection term for μ angular direction.
                if ( context->mu[q] < 0 ) {

                    for ( int64_t c = 0; c <= context->dgOrder; ++c )
                        IAQ_BLOCK_COL(context,q,d,e,f,d,e,c)
                            += neg1pow(f+c+1) * context->mu[q] * context->dx * context->dy
                               / ( (2*d + 1)*(2*e + 1) );
                } else {

                    for ( int64_t c = 0; c <= context->dgOrder; ++c )
                        IAQ_BLOCK_COL(context,q,d,e,f,d,e,c)
                            += context->mu[q] * context->dx * context->dy / ( (2*d + 1)*(2*e + 1) );
                }
            }}}
        }
    } break;

    case SWEEP_SEQUENTIAL_HESSENBERG:
    {
        // Size of data structures containing linear systems.
        context->sizeofA = sizeof(double) * context->nq
                           * (context->dgOrder + 1)*(context->dgOrder + 1)*(context->dgOrder + 1)
                           * (context->dgOrder + 1)*(context->dgOrder + 1)*(context->dgOrder + 1);

        context->sizeofB = sizeof(double) * context->nq
                           * (context->dgOrder + 1)*(context->dgOrder + 1)*(context->dgOrder + 1);

        // Allocate memory for linear systems.
        context->A   = malloc( context->sizeofA );
        context->Aq  = malloc( context->sizeofA );
        context->B   = malloc( context->sizeofB );

        context->P   = malloc( context->sizeofA );
        context->tau = malloc( context->sizeofB );

        // Precompute Aq.
        memset( context->Aq, 0, context->sizeofA );

        #pragma omp parallel for
        for ( int64_t q = 0; q < context->nq; ++q ) {

            // Compute A_q.
            for ( int64_t d = 0; d <= context->dgOrder; ++d ) {
            for ( int64_t e = 0; e <= context->dgOrder; ++e ) {
            for ( int64_t f = 0; f <= context->dgOrder; ++f ) {

                for ( int64_t a = d-1; a >= 0; a -= 2 )
                    IAQ_BLOCK_COL(context,q,d,e,f,a,e,f)
                        -= 2.0 * context->xi[q] * context->dy * context->dz * (2*d + 1);

                for ( int64_t b = e-1; b >= 0; b -= 2 )
                    IAQ_BLOCK_COL(context,q,d,e,f,d,b,f)
                        -= 2.0 * context->eta[q] * context->dx * context->dz * (2*e + 1);

                for ( int64_t c = f-1; c >= 0; c -= 2 )
                    IAQ_BLOCK_COL(context,q,d,e,f,d,e,c)
                        -= 2.0 * context->mu[q] * context->dx * context->dy * (2*f + 1);

                // Components from advection term for ξ angular direction.
                if ( context->xi[q] < 0 ) {

                    for ( int64_t a = 0; a <= context->dgOrder; ++a )
                        IAQ_BLOCK_COL(context,q,d,e,f,a,e,f)
                            += neg1pow(d+a+1) * context->xi[q] * context->dy * context->dz * (2*d + 1);
                } else {

                    for ( int64_t a = 0; a <= context->dgOrder; ++a )
                        IAQ_BLOCK_COL(context,q,d,e,f,a,e,f)
                            += context->xi[q] * context->dy * context->dz * (2*d + 1);
                }

                // Components from advection term for η angular direction.
                if ( context->eta[q] < 0 ) {

                    for ( int64_t b = 0; b <= context->dgOrder; ++b )
                        IAQ_BLOCK_COL(context,q,d,e,f,d,b,f)
                            += neg1pow(e+b+1) * context->eta[q] * context->dx * context->dz * (2*e + 1);
                } else {

                    for ( int64_t b = 0; b <= context->dgOrder; ++b )
                        IAQ_BLOCK_COL(context,q,d,e,f,d,b,f)
                            += context->eta[q] * context->dx * context->dz * (2*e + 1);
                }

                // Components from advection term for μ angular direction.
                if ( context->mu[q] < 0 ) {

                    for ( int64_t c = 0; c <= context->dgOrder; ++c )
                        IAQ_BLOCK_COL(context,q,d,e,f,d,e,c)
                            += neg1pow(f+c+1) * context->mu[q] * context->dx * context->dy * (2*f + 1);
                } else {

                    for ( int64_t c = 0; c <= context->dgOrder; ++c )
                        IAQ_BLOCK_COL(context,q,d,e,f,d,e,c)
                            += context->mu[q] * context->dx * context->dy * (2*f + 1);
                }
            }}}

            //
            // Reduce the matrix A_q to upper Hessenberg form by an orthogonal similarity transformation. This
            // is done through the following steps. We follow the source code of the DGEES function of the
            // LAPACK library, beginning at approximately line 267.
            //
            //   1. Use DGEHRD to compute the reduction to upper Hessenberg form by orthogonal similarity
            //      transformation A = P*H*P'. This routine computes (and stores as a return value) the
            //      coefficients and vectors required to form the Householder reflections which determine Q.
            //
            //   2. Use DLACPY to copy the Householder vectors from the lower triangular portion of the matrix
            //      returned by DGEHRD into the memory for the matrix P.
            //
            //      NOTE: The matrix P is not constructed explicitly. Instead the LAPACK routine DORMHR is
            //            is used to compute the application of P to any vectors or matrices. This routine
            //            applies each Householder transformation which together form P in terms of the
            //            scalars (stored in s_tau) and elementary reflectors (stored in the lower triangular
            //            portion of s_P) which define each Householder transformation.
            //

            char UPLO = 'L';
            int N     = (context->dgOrder + 1)*(context->dgOrder + 1)*(context->dgOrder + 1);
            int M     = N;
            int ILO   = 1;
            int IHI   = N;
            int LDA   = N;
            int LDB   = N;
            int LWORK = N*N;
            int INFO;

            double * WORK = malloc( LWORK * sizeof(double) );

            double * tau = context->tau + q * N;
            double * Aq  = context->Aq  + q * N*N;
            double * P   = context->P   + q * N*N;

            // 1. Compute upper Hessenberg form via orthogonal similarity transformation.
            dgehrd_( &N, &ILO, &IHI, Aq, &LDA, tau, WORK, &LWORK, &INFO );

            // 2. Copy Householder vectors.
            dlacpy_( &UPLO, &M, &N, Aq, &LDA, P, &LDB );

            free( WORK );
        }
    } break;

    case SWEEP_PLANEWAVE_LAPACK_LU:
    {
        // Size of linear systems.
        context->sizeofA = sizeof(double)
                           * (context->dgOrder + 1)*(context->dgOrder + 1)*(context->dgOrder + 1)
                           * (context->dgOrder + 1)*(context->dgOrder + 1)*(context->dgOrder + 1);

        context->sizeofB = sizeof(double)
                           * (context->dgOrder + 1)*(context->dgOrder + 1)*(context->dgOrder + 1);

        // Allocate memory for pre-computed components of matrices.
        context->Aq_wave = malloc( context->nq * sizeof(double *) );

        // Precompute Aq.
        # pragma omp parallel for
        for ( int64_t q = 0; q < context->nq; ++q ) {

            context->Aq_wave[q] = malloc( context->sizeofA );

            memset( context->Aq_wave[q], 0, context->sizeofA );

            // Compute A_q.
            for ( int64_t d = 0; d <= context->dgOrder; ++d ) {
            for ( int64_t e = 0; e <= context->dgOrder; ++e ) {
            for ( int64_t f = 0; f <= context->dgOrder; ++f ) {

                for ( int64_t a = d-1; a >= 0; a -= 2 )
                    IAQ_NONCONTIG_COL(context,q,d,e,f,a,e,f)
                        -= 2.0 * context->xi[q] * context->dy * context->dz / ( (2*e + 1)*(2*f + 1) );

                for ( int64_t b = e-1; b >= 0; b -= 2 )
                    IAQ_NONCONTIG_COL(context,q,d,e,f,d,b,f)
                        -= 2.0 * context->eta[q] * context->dx * context->dz / ( (2*d + 1)*(2*f + 1) );

                for ( int64_t c = f-1; c >= 0; c -= 2 )
                    IAQ_NONCONTIG_COL(context,q,d,e,f,d,e,c)
                        -= 2.0 * context->mu[q] * context->dx * context->dy / ( (2*d + 1)*(2*e + 1) );

                // Components from advection term for ξ angular direction.
                if ( context->xi[q] < 0 ) {

                    for ( int64_t a = 0; a <= context->dgOrder; ++a )
                        IAQ_NONCONTIG_COL(context,q,d,e,f,a,e,f)
                            += neg1pow(d+a+1) * context->xi[q] * context->dy * context->dz
                               / ( (2*e + 1)*(2*f + 1) );
                } else {

                    for ( int64_t a = 0; a <= context->dgOrder; ++a )
                        IAQ_NONCONTIG_COL(context,q,d,e,f,a,e,f) += context->xi[q] * context->dy * context->dz
                                                                    / ( (2*e + 1)*(2*f + 1) );
                }

                // Components from advection term for η angular direction.
                if ( context->eta[q] < 0 ) {

                    for ( int64_t b = 0; b <= context->dgOrder; ++b )
                        IAQ_NONCONTIG_COL(context,q,d,e,f,d,b,f)
                            += neg1pow(e+b+1) * context->eta[q] * context->dx * context->dz
                               / ( (2*d + 1)*(2*f + 1) );
                } else {

                    for ( int64_t b = 0; b <= context->dgOrder; ++b )
                        IAQ_NONCONTIG_COL(context,q,d,e,f,d,b,f)
                            += context->eta[q] * context->dx * context->dz / ( (2*d + 1)*(2*f + 1) );
                }

                // Components from advection term for μ angular direction.
                if ( context->mu[q] < 0 ) {

                    for ( int64_t c = 0; c <= context->dgOrder; ++c )
                        IAQ_NONCONTIG_COL(context,q,d,e,f,d,e,c)
                            += neg1pow(f+c+1) * context->mu[q] * context->dx * context->dy
                               / ( (2*d + 1)*(2*e + 1) );
                } else {

                    for ( int64_t c = 0; c <= context->dgOrder; ++c )
                        IAQ_NONCONTIG_COL(context,q,d,e,f,d,e,c)
                            += context->mu[q] * context->dx * context->dy / ( (2*d + 1)*(2*e + 1) );
                }
            }}}
        }
    } break;

    case SWEEP_PLANEWAVE_HESSENBERG:
    {
        // Size of linear systems.
        context->sizeofA = sizeof(double)
                           * (context->dgOrder + 1)*(context->dgOrder + 1)*(context->dgOrder + 1)
                           * (context->dgOrder + 1)*(context->dgOrder + 1)*(context->dgOrder + 1);

        context->sizeofB = sizeof(double)
                           * (context->dgOrder + 1)*(context->dgOrder + 1)*(context->dgOrder + 1);

        // Allocate memory for pre-computed components of matrices.
        context->Aq_wave  = malloc( context->nq * sizeof(double *) );
        context->P_wave   = malloc( context->nq * sizeof(double *) );
        context->tau_wave = malloc( context->nq * sizeof(double *) );

        // Precompute Aq.
        #pragma omp parallel for
        for ( int64_t q = 0; q < context->nq; ++q ) {

        # if defined (USE_ALIGNED_ALLOC)

            context->Aq_wave[q]  = aligned_alloc( ALIGNED_ALLOC_ALIGNMENT, context->sizeofA );
            context->tau_wave[q] = aligned_alloc( ALIGNED_ALLOC_ALIGNMENT, context->sizeofB );
            context->P_wave[q]   = aligned_alloc( ALIGNED_ALLOC_ALIGNMENT, context->sizeofA );
        # else
            context->Aq_wave[q]  = malloc( context->sizeofA );
            context->tau_wave[q] = malloc( context->sizeofB );
            context->P_wave[q]   = malloc( context->sizeofA );
        # endif

            memset( context->Aq_wave[q], 0, context->sizeofA );

            // Compute A_q.
            for ( int64_t d = 0; d <= context->dgOrder; ++d ) {
            for ( int64_t e = 0; e <= context->dgOrder; ++e ) {
            for ( int64_t f = 0; f <= context->dgOrder; ++f ) {

                for ( int64_t a = d-1; a >= 0; a -= 2 )
                    IAQ_NONCONTIG_COL(context,q,d,e,f,a,e,f)
                        -= 2.0 * context->xi[q] * context->dy * context->dz * (2*d + 1);

                for ( int64_t b = e-1; b >= 0; b -= 2 )
                    IAQ_NONCONTIG_COL(context,q,d,e,f,d,b,f)
                        -= 2.0 * context->eta[q] * context->dx * context->dz * (2*e + 1);

                for ( int64_t c = f-1; c >= 0; c -= 2 )
                    IAQ_NONCONTIG_COL(context,q,d,e,f,d,e,c)
                        -= 2.0 * context->mu[q] * context->dx * context->dy * (2*f + 1);

                // Components from advection term for ξ angular direction.
                if ( context->xi[q] < 0 ) {

                    for ( int64_t a = 0; a <= context->dgOrder; ++a )
                        IAQ_NONCONTIG_COL(context,q,d,e,f,a,e,f)
                            += neg1pow(d+a+1) * context->xi[q] * context->dy * context->dz * (2*d + 1);

                } else {

                    for ( int64_t a = 0; a <= context->dgOrder; ++a )
                        IAQ_NONCONTIG_COL(context,q,d,e,f,a,e,f)
                            += context->xi[q] * context->dy * context->dz * (2*d + 1);
                }

                // Components from advection term for η angular direction.
                if ( context->eta[q] < 0 ) {

                    for ( int64_t b = 0; b <= context->dgOrder; ++b )
                        IAQ_NONCONTIG_COL(context,q,d,e,f,d,b,f)
                            += neg1pow(e+b+1) * context->eta[q] * context->dx * context->dz * (2*e + 1);

                } else {

                    for ( int64_t b = 0; b <= context->dgOrder; ++b )
                        IAQ_NONCONTIG_COL(context,q,d,e,f,d,b,f)
                            += context->eta[q] * context->dx * context->dz * (2*e + 1);
                }

                // Components from advection term for μ angular direction.
                if ( context->mu[q] < 0 ) {

                    for ( int64_t c = 0; c <= context->dgOrder; ++c )
                        IAQ_NONCONTIG_COL(context,q,d,e,f,d,e,c)
                            += neg1pow(f+c+1) * context->mu[q] * context->dx * context->dy * (2*f + 1);

                } else {

                    for ( int64_t c = 0; c <= context->dgOrder; ++c )
                        IAQ_NONCONTIG_COL(context,q,d,e,f,d,e,c)
                            += context->mu[q] * context->dx * context->dy * (2*f + 1);
                }
            }}}

            //
            // Reduce the matrix A_q to upper Hessenberg form by an orthogonal similarity transformation. This
            // is done through the following steps. We follow the source code of the DGEES function of the
            // LAPACK library, beginning at approximately line 267.
            //
            //   1. Use DGEHRD to compute the reduction to upper Hessenberg form by orthogonal similarity
            //      transformation A = P*H*P'. This routine computes (and stores as a return value) the
            //      coefficients and vectors required to form the Householder reflections which determine Q.
            //
            //   2. Use DLACPY to copy the Householder vectors from the lower triangular portion of the matrix
            //      returned by DGEHRD into the memory for the matrix P.
            //
            //      NOTE: The matrix P is not constructed explicitly. Instead the LAPACK routine DORMHR is
            //            is used to compute the application of P to any vectors or matrices. This routine
            //            applies each Householder transformation which together form P in terms of the
            //            scalars (stored in s_tau) and elementary reflectors (stored in the lower triangular
            //            portion of s_P) which define each Householder transformation.
            //

            char UPLO = 'L';
            int N     = (context->dgOrder + 1)*(context->dgOrder + 1)*(context->dgOrder + 1);
            int M     = N;
            int ILO   = 1;
            int IHI   = N;
            int LDA   = N;
            int LDB   = N;
            int LWORK = N*N;
            int INFO;

            double * const restrict WORK = malloc( LWORK * sizeof(double) );

            // 1. Compute upper Hessenberg form via orthogonal similarity transformation.
            dgehrd_( &N, &ILO, &IHI, context->Aq_wave[q], &LDA, context->tau_wave[q], WORK, &LWORK, &INFO );

            // 2. Copy Householder vectors.
            dlacpy_( &UPLO, &M, &N, context->Aq_wave[q], &LDA, context->P_wave[q], &LDB );

            free( WORK );
        }
    } break;

    default:
    {
        PRINT_ERROR( "Invalid sweep type %d in %s.\n", context->type, __func__ )
        SWC_cleanup( context );
    } return;
    }

    PRINT_LOG( "  %-20s  %s\n", "Sweep Type:", g_sweepStrings[context->type] )
}


//============================================================================================================
//=== SEQUENTIAL REDUCED UPPER HESSENBERG ROUTINES ===========================================================
//============================================================================================================


//------------------------------------------------------------------------------------------------------------
//! \brief Computes the elements of the matrices \f$ A \f$ for the Hessenberg sweep.
//!
//! Computes the local matrices solved at the given step of the sweep algorithm. For given spatial indices
//! \pp{i},\pp{j},\pp{k} the local mass matrices are computed for each discrete ordinate angle. For
//! left-moving angles the spatial index is reflected across the spatial domain to reverse the direction of
//! the sweep. In this case \f$ \sigma \f$ is assumed to be piecewise constant in space (and therefore
//! constant across each spatial cell).
//!
//! For solving steady-state problems one should use <tt>\pp{dt} = inf</tt>.
//!
//! \param[in,out]  context     sweepContext_t object containing the context of the current sweep.
//!                             The matrices \f$ A \f$ computed by this function is stored in \pp{context->A}.
//! \param[in]      sigma       crossSection_t object containing the material cross-section \f$ \sigma \f$.
//! \param[in]      dt          Timestep size.
//! \param[in]      i           Index for the current spatial cell in \f$ x \f$ dimension.
//! \param[in]      j           Index for the current spatial cell in \f$ y \f$ dimension.
//! \param[in]      k           Index for the current spatial cell in \f$ z \f$ dimension.
//!
//! \see    sequentialHess_calcB()
//------------------------------------------------------------------------------------------------------------
static void sequentialHess_calcA (

    const sweepContext_t context,
    const crossSection_t sigma,
    const double dt,
    const int64_t i,
    const int64_t j,
    const int64_t k
) {

    // Set value using angular component of matrices.
    memcpy( context->A, context->Aq, context->sizeofA );

    #pragma omp parallel for
    for ( int64_t q = 0; q < context->nq; ++q ) {

        int64_t ii = i;
        int64_t jj = j;
        int64_t kk = k;

        // Reverse direction of sweep for negative angles.
        if ( context->xi[q]  < 0 ) {  ii = context->nx - i + 1;  }
        if ( context->eta[q] < 0 ) {  jj = context->ny - j + 1;  }
        if ( context->mu[q]  < 0 ) {  kk = context->nz - k + 1;  }

        for ( int64_t d = 0; d <= context->dgOrder; ++d ) {
        for ( int64_t e = 0; e <= context->dgOrder; ++e ) {
        for ( int64_t f = 0; f <= context->dgOrder; ++f ) {

            IA_BLOCK_COL(context,q,d,e,f,d,e,f)
                += ( ICSD(sigma,ii,jj,kk,0,0,0) + 1.0 / dt ) * context->dx * context->dy * context->dz;
        }}}
    }
}


//------------------------------------------------------------------------------------------------------------
//! \brief  Computes the elements of the right-hand vectors \f$ b \f$ for the Hessenberg sweep.
//!
//! Computes the known vectors for the local linear systems to be solved and the given step of the sweep
//! algorithm. For given spatial indices \pp{i},\pp{j},\pp{k} the vectors are computed for each discrete
//! ordinate angle. For left-moving angles the spatial index is reflected across the spatial domain to reverse
//! the direction of the sweep.
//!
//! \param[in,out]  context     sweepContext_t object containing the context of the current sweep.
//!                             The vectors \f$ b \f$ computed by this function are stored in \pp{context->B}.
//! \param[in]      Q           angularFlux_t object containing the coefficients of the source term \f$ Q \f$
//!                             for the transport system.
//! \param[in]      i           Index for the current spatial cell in \f$ x \f$ dimension.
//! \param[in]      j           Index for the current spatial cell in \f$ y \f$ dimension.
//! \param[in]      k           Index for the current spatial cell in \f$ z \f$ dimension.
//!
//! \see    sequentialHess_calcA()
//------------------------------------------------------------------------------------------------------------
static void sequentialHess_calcB (

    const sweepContext_t context,
    const angularFlux_t Q,
    const int64_t i,
    const int64_t j,
    const int64_t k
) {
# pragma omp parallel
{

    char SIDE  = 'L';
    char TRANS = 'T';
    int M      = (context->dgOrder + 1)*(context->dgOrder + 1)*(context->dgOrder + 1);
    int N      = 1;
    int INFO;

    double * const restrict WORK = malloc( M * sizeof(double) );

    # pragma omp for
    for ( int64_t q = 0; q < context->nq; ++q ) {

        int64_t ii = i;
        int64_t jj = j;
        int64_t kk = k;

        // Reverse direction of sweep for negative angles.
        if ( context->xi[q]  < 0 ) {  ii = context->nx - i + 1;  }
        if ( context->eta[q] < 0 ) {  jj = context->ny - j + 1;  }
        if ( context->mu[q]  < 0 ) {  kk = context->nz - k + 1;  }

        for ( int64_t d = 0; d <= context->dgOrder; ++d ) {
        for ( int64_t e = 0; e <= context->dgOrder; ++e ) {
        for ( int64_t f = 0; f <= context->dgOrder; ++f ) {

            // Set value from known source term.
            IB_BLOCK(context,q,d,e,f) = context->dx * context->dy * context->dz * IAF(Q,q,ii,jj,kk,d,e,f);

            // Include known components of advection term for ξ angular direction.
            if ( context->xi[q] < 0 ) {

                for ( int64_t a = 0; a <= context->dgOrder; ++a )
                    IB_BLOCK(context,q,d,e,f) -= neg1pow(a) * context->xi[q] * context->dy * context->dz
                                                 * IAF(Q,q,ii+1,jj,kk,a,e,f) * (2*d + 1);
            } else {

                for ( int64_t a = 0; a <= context->dgOrder; ++a )
                    IB_BLOCK(context,q,d,e,f) -= neg1pow(d+1) * context->xi[q] * context->dy * context->dz
                                                 * IAF(Q,q,ii-1,jj,kk,a,e,f) * (2*d + 1);
            }

            // Include known components of advection term for η angular direction.
            if ( context->eta[q] < 0 ) {

                for ( int64_t b = 0; b <= context->dgOrder; ++b )
                    IB_BLOCK(context,q,d,e,f) -= neg1pow(b) * context->eta[q] * context->dx * context->dz
                                                 * IAF(Q,q,ii,jj+1,kk,d,b,f) * (2*e + 1);
            } else {

                for ( int64_t b = 0; b <= context->dgOrder; ++b )
                    IB_BLOCK(context,q,d,e,f) -= neg1pow(e+1) * context->eta[q] * context->dx * context->dz
                                                 * IAF(Q,q,ii,jj-1,kk,d,b,f) * (2*e + 1);
            }

            // Include known components of advection term for μ angular direction.
            if ( context->mu[q] < 0 ) {

                for ( int64_t c = 0; c <= context->dgOrder; ++c )
                    IB_BLOCK(context,q,d,e,f) -= neg1pow(c) * context->mu[q] * context->dx * context->dy
                                                 * IAF(Q,q,ii,jj,kk+1,d,e,c) * (2*f + 1);
            } else {

                for ( int64_t c = 0; c <= context->dgOrder; ++c )
                    IB_BLOCK(context,q,d,e,f) -= neg1pow(f+1) * context->mu[q] * context->dx * context->dy
                                                 * IAF(Q,q,ii,jj,kk-1,d,e,c) * (2*f + 1);
            }
        }}}

        // Pᵀb.
        double * const P   = context->P   + q * M*M;
        double * const tau = context->tau + q * M;
        double * const B   = context->B   + q * M;

        dormhr_( &SIDE, &TRANS, &M, &N, &N, &M, P, &M, tau, B, &M, WORK, &M, &INFO );

    # if defined (STRICT_CHECK)
        if ( INFO ) {

            PRINT_ERROR( "dormhr_ returned INFO = %d.\n", INFO )
            exit( EXIT_FAILURE );
        }
    # endif
    }

    free( WORK );
}}


//============================================================================================================
//=== SEQUENTIAL LU ROUTINES =================================================================================
//============================================================================================================


//------------------------------------------------------------------------------------------------------------
//! \brief Computes the elements of the matrices \f$ A \f$ for the \f$ LU \f$ sweep.
//!
//! Computes the local matrices solved at the given step of the sweep algorithm. For given spatial indices
//! \pp{i},\pp{j},\pp{k} the local mass matrices are computed for each discrete ordinate angle. For
//! left-moving angles the spatial index is reflected across the spatial domain to reverse the direction of
//! the sweep. In this case \f$ \sigma \f$ is assumed to be high order in space, and thus the full tensor of
//! Legendre product integrals is used.
//!
//! For solving steady-state problems one should use <tt>\pp{dt} = inf</tt>.
//!
//! \param[in,out]  context     sweepContext_t object containing the context of the current sweep.
//!                             The matrices \f$ A \f$ computed by this function is stored in \pp{context->A}.
//! \param[in]      sigma       crossSection_t object containing the material cross-section \f$ \sigma \f$.
//! \param[in]      dt          Timestep size.
//! \param[in]      i           Index for the current spatial cell in \f$ x \f$ dimension.
//! \param[in]      j           Index for the current spatial cell in \f$ y \f$ dimension.
//! \param[in]      k           Index for the current spatial cell in \f$ z \f$ dimension.
//!
//! \see    sequentialLU_calcB()
//------------------------------------------------------------------------------------------------------------
static void sequentialLU_calcA (

    const sweepContext_t context,
    const crossSection_t sigma,
    const double dt,
    const int64_t i,
    const int64_t j,
    const int64_t k
) {

    // Set value using angular component of matrices.
    memcpy( context->A, context->Aq, context->sizeofA );

    #pragma omp parallel for
    for ( int64_t q = 0; q < context->nq; ++q ) {

        int64_t ii = i;
        int64_t jj = j;
        int64_t kk = k;

        // Reverse direction of sweep for negative angles.
        if ( context->xi[q]  < 0 ) {  ii = context->nx - i + 1;  }
        if ( context->eta[q] < 0 ) {  jj = context->ny - j + 1;  }
        if ( context->mu[q]  < 0 ) {  kk = context->nz - k + 1;  }

        for ( int64_t d = 0; d <= context->dgOrder; ++d ) {
        for ( int64_t e = 0; e <= context->dgOrder; ++e ) {
        for ( int64_t f = 0; f <= context->dgOrder; ++f ) {

            IA_BLOCK_COL(context,q,d,e,f,d,e,f)
                += context->dx * context->dy * context->dz / ( (2*d + 1)*(2*e + 1)*(2*f + 1) * dt );

        for ( int64_t u = 0; u <= context->dgOrder; ++u ) {
        for ( int64_t v = 0; v <= context->dgOrder; ++v ) {
        for ( int64_t w = 0; w <= context->dgOrder; ++w ) {

            IA_BLOCK_COL(context,q,d,e,f,u,v,w)
                += context->dx * context->dy * context->dz * ICST(sigma,ii,jj,kk,d,e,f,u,v,w) / 8.0;
        }}}}}}
    }
}


//------------------------------------------------------------------------------------------------------------
//! \brief  Computes the elements of the right-hand vectors \f$ b \f$ for the \f$ LU \f$ sweep.
//!
//! Computes the known vectors for the local linear systems to be solved and the given step of the sweep
//! algorithm. For given spatial indices \pp{i},\pp{j},\pp{k} the vectors are computed for each discrete
//! ordinate angle. For left-moving angles the spatial index is reflected across the spatial domain to reverse
//! the direction of the sweep.
//!
//! \param[in,out]  context     sweepContext_t object containing the context of the current sweep.
//!                             The vectors \f$ b \f$ computed by this function are stored in \pp{context->B}.
//! \param[in]      Q           angularFlux_t object containing the coefficients of the source term of the
//!                             transport system.
//! \param[in]      i           Index for the current spatial cell in \f$ x \f$ dimension.
//! \param[in]      j           Index for the current spatial cell in \f$ y \f$ dimension.
//! \param[in]      k           Index for the current spatial cell in \f$ z \f$ dimension.
//!
//! \see    sequentialLU_calcA()
//------------------------------------------------------------------------------------------------------------
static void sequentialLU_calcB (

    const sweepContext_t context,
    const angularFlux_t Q,
    const int64_t i,
    const int64_t j,
    const int64_t k
) {

    #pragma omp parallel for
    for ( int64_t q = 0; q < context->nq; ++q ) {

        int64_t ii = i;
        int64_t jj = j;
        int64_t kk = k;

        // Reverse direction of sweep for negative angles.
        if ( context->xi[q]  < 0 ) {  ii = context->nx - i + 1;  }
        if ( context->eta[q] < 0 ) {  jj = context->ny - j + 1;  }
        if ( context->mu[q]  < 0 ) {  kk = context->nz - k + 1;  }

        for ( int64_t d = 0; d <= context->dgOrder; ++d ) {
        for ( int64_t e = 0; e <= context->dgOrder; ++e ) {
        for ( int64_t f = 0; f <= context->dgOrder; ++f ) {

            // Set value from known source term.
            IB_BLOCK(context,q,d,e,f) = context->dx * context->dy * context->dz * IAF(Q,q,ii,jj,kk,d,e,f)
                                        / ( (2*d + 1)*(2*e + 1)*(2*f + 1) );

            // Include known components of advection term for ξ angular direction.
            if ( context->xi[q] < 0 ) {

                for ( int64_t a = 0; a <= context->dgOrder; ++a )
                    IB_BLOCK(context,q,d,e,f) -= neg1pow(a) * context->xi[q] * context->dy * context->dz
                                                 * IAF(Q,q,ii+1,jj,kk,a,e,f) / ( (2*e + 1)*(2*f + 1) );
            } else {

                for ( int64_t a = 0; a <= context->dgOrder; ++a )
                    IB_BLOCK(context,q,d,e,f) -= neg1pow(d+1) * context->xi[q] * context->dy * context->dz
                                                 * IAF(Q,q,ii-1,jj,kk,a,e,f) / ( (2*e + 1)*(2*f + 1) );
            }

            // Include known components of advection term for η angular direction.
            if ( context->eta[q] < 0 ) {

                for ( int64_t b = 0; b <= context->dgOrder; ++b )
                    IB_BLOCK(context,q,d,e,f) -= neg1pow(b) * context->eta[q] * context->dx * context->dz
                                                 * IAF(Q,q,ii,jj+1,kk,d,b,f) / ( (2*d + 1)*(2*f + 1) );
            } else {

                for ( int64_t b = 0; b <= context->dgOrder; ++b )
                    IB_BLOCK(context,q,d,e,f) -= neg1pow(e+1) * context->eta[q] * context->dx * context->dz
                                                 * IAF(Q,q,ii,jj-1,kk,d,b,f) / ( (2*d + 1)*(2*f + 1) );
            }

            // Include known components of advection term for μ angular direction.
            if ( context->mu[q] < 0 ) {

                for ( int64_t c = 0; c <= context->dgOrder; ++c )
                    IB_BLOCK(context,q,d,e,f) -= neg1pow(c) * context->mu[q] * context->dx * context->dy
                                                 * IAF(Q,q,ii,jj,kk+1,d,e,c) / ( (2*d + 1)*(2*e + 1) );
            } else {

                for ( int64_t c = 0; c <= context->dgOrder; ++c )
                    IB_BLOCK(context,q,d,e,f) -= neg1pow(f+1) * context->mu[q] * context->dx * context->dy
                                                 * IAF(Q,q,ii,jj,kk-1,d,e,c) / ( (2*d + 1)*(2*e + 1) );
            }
        }}}
    }
}


//============================================================================================================
//=== SEQUENTIAL INTERFACE ROUTINES ==========================================================================
//============================================================================================================


//------------------------------------------------------------------------------------------------------------
//! \brief  Calls the appropriate \c sequential*_calcA() routine corresponding to the sequential sweep
//!         algorithm as set by \pp{context->type}.
//!
//! \param[in,out]  context     sweepContext_t object containing the context of the current sweep.
//! \param[in]      sigma       crossSection_t object containing the material cross-section \f$ \sigma \f$.
//! \param[in]      dt          Timestep size.
//! \param[in]      i           Index for the current spatial cell in \f$ x \f$ dimension.
//! \param[in]      j           Index for the current spatial cell in \f$ y \f$ dimension.
//! \param[in]      k           Index for the current spatial cell in \f$ z \f$ dimension.
//!
//! \see    sequential_sweep()
//! \see    sequential_calcB()
//! \see    sequentialHess_calcA()
//! \see    sequentialLU_calcA()
//------------------------------------------------------------------------------------------------------------
static void sequential_calcA (

    const sweepContext_t context,
    const crossSection_t sigma,
    const double dt,
    const int64_t i,
    const int64_t j,
    const int64_t k
) {

    switch ( context->type ) {

    case SWEEP_SEQUENTIAL_LAPACK_LU:

        sequentialLU_calcA( context, sigma, dt, i,j,k );
        return;

    case SWEEP_SEQUENTIAL_HESSENBERG:

        sequentialHess_calcA( context, sigma, dt, i,j,k );
        return;

    default:

        PRINT_ERROR( "Invalid sweep type %d in %s.\n", context->type, __func__ )
        return;
    }
}


//------------------------------------------------------------------------------------------------------------
//! \brief  Calls the appropriate \c sequential*_calcB() routine corresponding to the sequential sweep
//!         algorithm as set by \pp{context->type}.
//!
//! \param[in,out]  context     sweepContext_t object containing the context of the current sweep.
//! \param[in]      Q           angularFlux_t object containing the coefficients of the source term of the
//!                             transport system.
//! \param[in]      i           Index for the current spatial cell in \f$ x \f$ dimension.
//! \param[in]      j           Index for the current spatial cell in \f$ y \f$ dimension.
//! \param[in]      k           Index for the current spatial cell in \f$ z \f$ dimension.
//!
//! \see    sequential_sweep()
//! \see    sequential_calcA()
//! \see    sequentialHess_calcB()
//! \see    sequentialLU_calcB()
//------------------------------------------------------------------------------------------------------------
static void sequential_calcB (

    const sweepContext_t context,
    const angularFlux_t Q,
    const int64_t i,
    const int64_t j,
    const int64_t k
) {

    switch ( context->type ) {

    case SWEEP_SEQUENTIAL_LAPACK_LU:

        sequentialLU_calcB( context, Q, i,j,k );
        return;

    case SWEEP_SEQUENTIAL_HESSENBERG:

        sequentialHess_calcB( context, Q, i,j,k );
        return;

    default:

        PRINT_ERROR( "Invalid sweep type %d in %s.\n", context->type, __func__ )
        return;
    }
}


//------------------------------------------------------------------------------------------------------------
//! \brief  Calls the linear solver corresponding to the sequential sweep algorithm as set by
//!         \pp{context->type}.
//!
//! Calls the proper linear solver to solve the current linear system for each ordinate.
//!
//! \param[in,out]  context     sweepContext_t object containing the context of the current sweep.
//------------------------------------------------------------------------------------------------------------
static void sequential_linearSolve (

    const sweepContext_t context
) {

    switch ( context->type ) {

    case SWEEP_SEQUENTIAL_LAPACK_LU:
    # pragma omp parallel
    {
        int N    = (context->dgOrder + 1)*(context->dgOrder + 1)*(context->dgOrder + 1);
        int NRHS = 1;
        int INFO;

        int * const restrict IPIV = malloc( N * sizeof(int) );

        # pragma omp for
        for ( int64_t q = 0; q < context->nq; ++q ) {

            double * const A = context->A + q * N*N;
            double * const B = context->B + q * N;

            // Call solver from LAPACK.
            dgesv_( &N, &NRHS, A, &N, IPIV, B, &N, &INFO );

        # if defined (STRICT_CHECK)
            if ( INFO ) {

                PRINT_ERROR( "dgesv_ returned INFO = %d.\n", INFO )
                exit( EXIT_FAILURE );
            }
        # endif
        }

        free( IPIV );

    } return;

    case SWEEP_SEQUENTIAL_HESSENBERG:
    # pragma omp parallel
    {
        char SIDE  = 'L';
        char TRANS = 'N';
        int M      = (context->dgOrder + 1)*(context->dgOrder + 1)*(context->dgOrder + 1);
        int N      = 1;
        int INFO;

        void * const restrict WORK = malloc( M * sizeof(double) );

        # pragma omp for
        for ( int64_t q = 0; q < context->nq; ++q ) {

            double * const A   = context->A   + q * M*M;
            double * const P   = context->P   + q * M*M;
            double * const B   = context->B   + q * M;
            double * const tau = context->tau + q * M;

            // Solve for PᵀΨ.
            upperHessenbergSolve_columnMajor( M, A, B, WORK );

            // Solve for Ψ using PᵀΨ.
            dormhr_( &SIDE, &TRANS, &M, &N, &N, &M, P, &M, tau, B, &M, WORK, &M, &INFO );

        # if defined (STRICT_CHECK)
            if ( INFO ) {

                PRINT_ERROR( "dormhr_ returned INFO = %d.\n", INFO )
                exit( EXIT_FAILURE );
            }
        # endif
        }

        free( WORK );

    } return;

    default:
    {
        PRINT_ERROR( "Invalid sweep type %d in %s.\n", context->type, __func__ )
        return;
    }
    }
}


//------------------------------------------------------------------------------------------------------------
//! \brief  Driver function for solving the system \f$ \mathcal{L} \Psi = Q \f$ using sequential sweep
//!         algorithms.
//!
//! Solves local linear systems either using \f$ LU \f$ decompositions or upper Hessenberg reductions.
//! Multi-threaded parallelism is provided over the angular variable only.
//!
//! \param[in]      context     sweepContext_t object containing the context of the current sweep.
//! \param[in,out]  Q           Initially contains the coefficients of the source term \f$ Q \f$ for the
//!                             transport system.
//!                             Upon return, contains the coefficients of \f$ \Psi \f$ which solves the
//!                             transport system.
//! \param[in]      sigma       crossSection_t object containing the material cross-section \f$ \sigma \f$.
//! \param[in]      dt          Timestep size.
//!
//! \see    sequential_calcA()
//! \see    sequential_calcB()
//------------------------------------------------------------------------------------------------------------
int sequential_sweep (

    const sweepContext_t context,
    const angularFlux_t Q,
    const crossSection_t sigma,
    const double dt
) {

    PRINT_STATUS( "Computing transport sweep with %" PRId64 " angles.\n", Q->nq )

    int returnValue = 0;

    // Sequential sweeps are not designed to run with MPI.
    if ( g_MPI_numRanks > 1 ) {

        PRINT_ERROR( "MPI communication not implemented for sequential sweep algorithms.\n" )
        MPI_Abort( MPI_COMM_WORLD, EXIT_FAILURE );
    }

    // Sweep across the spatial mesh.
    for ( int64_t i = 1; i <= context->nx; ++i ) {
    for ( int64_t j = 1; j <= context->ny; ++j ) {
    for ( int64_t k = 1; k <= context->nz; ++k ) {

        // Construct and solve linear systems.
        sequential_calcB( context, Q, i, j, k );
        sequential_calcA( context, sigma, dt, i, j, k );
        sequential_linearSolve( context );

        // Save the solution values.
        #pragma omp parallel for
        for ( int64_t q = 0; q < context->nq; ++q ) {

            int64_t ii = i;
            int64_t jj = j;
            int64_t kk = k;

            // Reverse direction of sweep for negative angles.
            if ( context->xi[q]  < 0 ) {  ii = context->nx - i + 1;  }
            if ( context->eta[q] < 0 ) {  jj = context->ny - j + 1;  }
            if ( context->mu[q]  < 0 ) {  kk = context->nz - k + 1;  }

            for ( int64_t d = 0; d <= context->dgOrder; ++d ) {
            for ( int64_t e = 0; e <= context->dgOrder; ++e ) {
            for ( int64_t f = 0; f <= context->dgOrder; ++f ) {

                IAF(Q,q,ii,jj,kk,d,e,f) = IB_BLOCK(context,q,d,e,f);
            }}}
        }
    }}}

    return returnValue;
}


//============================================================================================================
//=== WAVE-TYPE REDUCED UPPER HESSENBERG ROUTINES ============================================================
//============================================================================================================


//------------------------------------------------------------------------------------------------------------
//! \brief  Computes the elements of the matrix \f$ A \f$ for the linear system to be solved at a given step
//!         of the wavefront Hessenberg sweep.
//!
//! Computes the matrix for the local linear system to be solved by the sweep algorithm for the \pp{q}th
//! ordinate at the spatial cell with indices \pp{i},\pp{j},\pp{k}. In this case \f$ \sigma \f$ is assumed to
//! be piecewise constant in space (and therefore constant across each spatial cell).
//!
//! For solving steady-state problems one should use <tt>\pp{dt} = inf</tt>.
//!
//! \param[out]     A           Pointer to memory location at which to store the computed matrix.
//! \param[in]      context     sweepContext_t object containing the context of the current sweep.
//! \param[in]      sigma       crossSection_t object containing the material cross-section \f$ \sigma \f$.
//! \param[in]      dt          Timestep size.
//! \param[in]      q           Index of angular ordinate to compute matrix for.
//! \param[in]      i           Index for the current spatial cell in \f$ x \f$ dimension.
//! \param[in]      j           Index for the current spatial cell in \f$ y \f$ dimension.
//! \param[in]      k           Index for the current spatial cell in \f$ z \f$ dimension.
//!
//! \see    waveHess_calcB()
//------------------------------------------------------------------------------------------------------------
static void waveHess_calcA (

    double * const restrict A,
    const sweepContext_t context,
    const crossSection_t sigma,
    const double dt,
    const int64_t q,
    const int64_t i,
    const int64_t j,
    const int64_t k
) {

    memcpy( A, context->Aq_wave[q], context->sizeofA );

    for ( int64_t d = 0; d <= context->dgOrder; ++d ) {
    for ( int64_t e = 0; e <= context->dgOrder; ++e ) {
    for ( int64_t f = 0; f <= context->dgOrder; ++f ) {

        A[ SWC_matrixIndexSingle_columnMajor(context,d,e,f,d,e,f) ]
            += ( ICSD(sigma,i,j,k,0,0,0) + 1.0 / dt ) * context->dx * context->dy * context->dz;
    }}}
}


//------------------------------------------------------------------------------------------------------------
//! \brief  Computes the elements of the right-hand vector \f$ b \f$ for the linear system to be solved at a
//!         given step of the wavefront Hessenberg sweep.
//!
//! Computes the known vector for the local linear system for the \pp{q}th ordinate at the spatial cell with
//! indices \pp{i},\pp{j},\pp{k}.
//!
//! \param[out]     B           Pointer to memory location at which to store the computed vector.
//! \param[in,out]  context     sweepContext_t object containing the context of the current sweep.
//! \param[in]      Q           angularFlux_t object containing the coefficients of the source term \f$ Q \f$
//!                             for the transport system.
//! \param[in]      q           Index of angular ordinate to compute vector for.
//! \param[in]      i           Index for the current spatial cell in \f$ x \f$ dimension.
//! \param[in]      j           Index for the current spatial cell in \f$ y \f$ dimension.
//! \param[in]      k           Index for the current spatial cell in \f$ z \f$ dimension.
//! \param[in]      WORK        Pointer to memory location used for working space for LAPACK routines. Assumed
//!                             to be at least \pp{context->sizeofB} Bytes.
//!
//! \see    waveHess_calcA()
//------------------------------------------------------------------------------------------------------------
static void waveHess_calcB (

    double * const restrict B,
    const sweepContext_t context,
    const angularFlux_t Q,
    const int64_t q,
    const int64_t i,
    const int64_t j,
    const int64_t k,
    void * const restrict WORK
) {

    for ( int64_t d = 0; d <= context->dgOrder; ++d ) {
    for ( int64_t e = 0; e <= context->dgOrder; ++e ) {
    for ( int64_t f = 0; f <= context->dgOrder; ++f ) {

        // Set value from known source term.
        B[ SWC_vectorIndexSingle(context,d,e,f) ]
            = context->dx * context->dy * context->dz * IAF(Q,q,i,j,k,d,e,f);

        // Include known components of advection term for ξ angular direction.
        if ( context->xi[q] < 0 ) {

            for ( int64_t a = 0; a <= context->dgOrder; ++a )
                B[ SWC_vectorIndexSingle(context,d,e,f) ]
                    -= neg1pow(a) * context->xi[q] * context->dy * context->dz * (2*d + 1)
                       * IAF(Q,q,i+1,j,k,a,e,f);
        } else {

            for ( int64_t a = 0; a <= context->dgOrder; ++a )
                B[ SWC_vectorIndexSingle(context,d,e,f) ]
                    -= neg1pow(d+1) * context->xi[q] * context->dy * context->dz * (2*d + 1)
                       * IAF(Q,q,i-1,j,k,a,e,f);
        }

        // Include known components of advection term for η angular direction.
        if ( context->eta[q] < 0 ) {

            for ( int64_t b = 0; b <= context->dgOrder; ++b )
                B[ SWC_vectorIndexSingle(context,d,e,f) ]
                    -= neg1pow(b) * context->eta[q] * context->dx * context->dz * (2*e + 1)
                       * IAF(Q,q,i,j+1,k,d,b,f);
        } else {

            for ( int64_t b = 0; b <= context->dgOrder; ++b )
                B[ SWC_vectorIndexSingle(context,d,e,f) ]
                    -= neg1pow(e+1) * context->eta[q] * context->dx * context->dz * (2*e + 1)
                       * IAF(Q,q,i,j-1,k,d,b,f);
        }

        // Include known components of advection term for μ angular direction.
        if ( context->mu[q] < 0 ) {

            for ( int64_t c = 0; c <= context->dgOrder; ++c )
                B[ SWC_vectorIndexSingle(context,d,e,f) ]
                    -= neg1pow(c) * context->mu[q] * context->dx * context->dy * (2*f + 1)
                       * IAF(Q,q,i,j,k+1,d,e,c);
        } else {

            for ( int64_t c = 0; c <= context->dgOrder; ++c )
                B[ SWC_vectorIndexSingle(context,d,e,f) ]
                    -= neg1pow(f+1) * context->mu[q] * context->dx * context->dy * (2*f + 1)
                       * IAF(Q,q,i,j,k-1,d,e,c);
        }
    }}}

    // Pᵀb.
    const int N = (context->dgOrder + 1)*(context->dgOrder + 1)*(context->dgOrder + 1);

    multHouseholderReflections_transpose( B, context->P_wave[q], context->tau_wave[q], N );
}


//============================================================================================================
//=== WAVE-TYPE LU ROUTINES ==================================================================================
//============================================================================================================


//------------------------------------------------------------------------------------------------------------
//! \brief Computes the elements of the matrix \f$ A \f$ to be solved at a given step of the wave-type
//! \f$ LU \f$ sweeps.
//!
//! Computes the local matrix to be solved by the sweep algorithm for the \pp{q}th ordinate at the spatial
//! cell with indices \pp{i},\pp{j},\pp{K}. In this case \f$ \sigma \f$ is assumed to be high order in space,
//! and Thus the full tensor of Legendre product integrals (as stored in \pp{sigma->tensor}) is used.
//!
//! For solving steady-state problems one should use <tt>\pp{dt} = inf</tt>.
//!
//! \param[out]     A           Pointer to memory location to store computed matrix.
//! \param[in]      context     sweepContext_t object containing the context of the current sweep.
//! \param[in]      sigma       crossSection_t object containing the material cross-section \f$ \sigma \f$.
//! \param[in]      dt          Timestep size.
//! \param[in]      q           Index of angular ordinate to compute matrix for.
//! \param[in]      i           Index for the current spatial cell in \f$ x \f$ dimension.
//! \param[in]      j           Index for the current spatial cell in \f$ y \f$ dimension.
//! \param[in]      k           Index for the current spatial cell in \f$ z \f$ dimension.
//!
//! \see    planewave_sweep()
//! \see    wavefront_sweep()
//! \see    waveLU_calcB()
//------------------------------------------------------------------------------------------------------------
static void waveLU_calcA (

    double * const restrict A,
    const sweepContext_t context,
    const crossSection_t sigma,
    const double dt,
    const int64_t q,
    const int64_t i,
    const int64_t j,
    const int64_t k
) {

    memcpy( A, context->Aq_wave[q], context->sizeofA );

    for ( int64_t d = 0; d <= context->dgOrder; ++d ) {
    for ( int64_t e = 0; e <= context->dgOrder; ++e ) {
    for ( int64_t f = 0; f <= context->dgOrder; ++f ) {

        A[ SWC_matrixIndexSingle_columnMajor(context,d,e,f,d,e,f) ]
            += context->dx * context->dy * context->dz / ( (2*d + 1)*(2*e + 1)*(2*f + 1) * dt );

    for ( int64_t u = 0; u <= context->dgOrder; ++u ) {
    for ( int64_t v = 0; v <= context->dgOrder; ++v ) {
    for ( int64_t w = 0; w <= context->dgOrder; ++w ) {

        A[ SWC_matrixIndexSingle_columnMajor(context,d,e,f,u,v,f) ]
            += context->dx * context->dy * context->dz * ICST(sigma,i,j,k,d,e,f,u,v,w) / 8.0;
    }}}}}}
}


//------------------------------------------------------------------------------------------------------------
//! \brief  Computes the elements of the right-hand vector \f$ b \f$ to be used for solving the linear system
//!         at the given step of the wave-type \f$ LU \f$ sweeps.
//!
//! Computes the known vector for the local linear system for the \pp{q}th ordinate at the spatial cell with
//! indices \pp{i},\pp{j},\pp{k}.
//!
//! \param[out]     B           Pointer to memory location to store computed vector.
//! \param[in]      context     sweepContext_t object containing the context of the current sweep.
//! \param[in]      Q           angularFlux_t object containing the coefficients of the source term of the
//!                             transport system.
//! \param[in]      q           Index of angular ordinate to compute vector for.
//! \param[in]      i           Index for the current spatial cell in \f$ x \f$ dimension.
//! \param[in]      j           Index for the current spatial cell in \f$ y \f$ dimension.
//! \param[in]      k           Index for the current spatial cell in \f$ z \f$ dimension.
//!
//! \see    planewave_sweep()
//! \see    wavefront_sweep()
//! \see    waveLU_calcA()
//------------------------------------------------------------------------------------------------------------
static void waveLU_calcB (

    double * const restrict B,
    const sweepContext_t context,
    const angularFlux_t Q,
    const int64_t q,
    const int64_t i,
    const int64_t j,
    const int64_t k
) {

    for ( int64_t d = 0; d <= context->dgOrder; ++d ) {
    for ( int64_t e = 0; e <= context->dgOrder; ++e ) {
    for ( int64_t f = 0; f <= context->dgOrder; ++f ) {

        // Set value from known source term.
        B[ SWC_vectorIndexSingle(context,d,e,f) ]
            = context->dx * context->dy * context->dz * IAF(Q,q,i,j,k,d,e,f)
              / ( (2*d + 1)*(2*e + 1)*(2*f + 1) );

        // Include known components of advection term for ξ angular direction.
        if ( context->xi[q] < 0 ) {

            for ( int64_t a = 0; a <= context->dgOrder; ++a )
                B[ SWC_vectorIndexSingle(context,d,e,f) ]
                    -= neg1pow(a) * context->xi[q] * context->dy * context->dz * IAF(Q,q,i+1,j,k,a,e,f)
                       / ( (2*e + 1)*(2*f + 1) );
        } else {

            for ( int64_t a = 0; a <= context->dgOrder; ++a )
                B[ SWC_vectorIndexSingle(context,d,e,f) ]
                    -= neg1pow(d+1) * context->xi[q] * context->dy * context->dz * IAF(Q,q,i-1,j,k,a,e,f)
                       / ( (2*e + 1)*(2*f + 1) );
        }

        // Include known components of advection term for η angular direction.
        if ( context->eta[q] < 0 ) {

            for ( int64_t b = 0; b <= context->dgOrder; ++b )
                B[ SWC_vectorIndexSingle(context,d,e,f) ]
                    -= neg1pow(b) * context->eta[q] * context->dx * context->dz * IAF(Q,q,i,j+1,k,d,b,f)
                       / ( (2*d + 1)*(2*f + 1) );
        } else {

            for ( int64_t b = 0; b <= context->dgOrder; ++b )
                B[ SWC_vectorIndexSingle(context,d,e,f) ]
                    -= neg1pow(e+1) * context->eta[q] * context->dx * context->dz * IAF(Q,q,i,j-1,k,d,b,f)
                       / ( (2*d + 1)*(2*f + 1) );
        }

        // Include known components of advection term for μ angular direction.
        if ( context->mu[q] < 0 ) {

            for ( int64_t c = 0; c <= context->dgOrder; ++c )
                B[ SWC_vectorIndexSingle(context,d,e,f) ]
                    -= neg1pow(c) * context->mu[q] * context->dx * context->dy * IAF(Q,q,i,j,k+1,d,e,c)
                       / ( (2*d + 1)*(2*e + 1) );
        } else {

            for ( int64_t c = 0; c <= context->dgOrder; ++c )
                B[ SWC_vectorIndexSingle(context,d,e,f) ]
                    -= neg1pow(f+1) * context->mu[q] * context->dx * context->dy * IAF(Q,q,i,j,k-1,d,e,c)
                       / ( (2*d + 1)*(2*e + 1) );
        }
    }}}
}


//============================================================================================================
//=== WAVE-TYPE INTERFACE ROUTINES ===========================================================================
//============================================================================================================


//------------------------------------------------------------------------------------------------------------
//! \brief  Calls the appropriate \c wave*_calcA() routine corresponding to the wave-type sweep algorithm as
//!         set by \pp{context->type}.
//!
//! \param[out]     A           Pointer to memory location at which to store the computed matrix.
//! \param[in]      context     sweepContext_t object containing the context of the current sweep.
//! \param[in]      sigma       crossSection_t object containing the material cross-section \f$ \sigma \f$.
//! \param[in]      dt          Timestep size.
//! \param[in]      q           Index of angular ordinate to compute matrix for.
//! \param[in]      i           Index for the current spatial cell in \f$ x \f$ dimension.
//! \param[in]      j           Index for the current spatial cell in \f$ y \f$ dimension.
//! \param[in]      k           Index for the current spatial cell in \f$ z \f$ dimension.
//!
//! \see    wave_calcB()
//! \see    wave_linearSolve()
//------------------------------------------------------------------------------------------------------------
static void wave_calcA (

    double * const restrict A,
    const sweepContext_t context,
    const crossSection_t sigma,
    const double dt,
    const int64_t q,
    const int64_t i,
    const int64_t j,
    const int64_t k
) {

# if defined (DO_SWEEP_SUBTIMING)
    TMR_start( &g_TMR_calcA );
# endif

    switch ( context->type ) {

    case SWEEP_PLANEWAVE_LAPACK_LU:
    case SWEEP_WAVEFRONT_LAPACK_LU:

        waveLU_calcA( A, context, sigma, dt, q, i,j,k );
        goto cleanup;

    case SWEEP_PLANEWAVE_GAUSSIAN_ELIMINATION:
    case SWEEP_WAVEFRONT_GAUSSIAN_ELIMINATION:

//        waveGE_calcA( A, context, sigma, dt, q, i,j,k );
        goto cleanup;

    case SWEEP_PLANEWAVE_HESSENBERG:
    case SWEEP_WAVEFRONT_HESSENBERG:

        waveHess_calcA( A, context, sigma, dt, q, i,j,k );
        goto cleanup;

    default:

        PRINT_ERROR( "Invalid sweep type %d in %s.\n", context->type, __func__ )
        return;
    }

cleanup:

# if defined (DO_SWEEP_SUBTIMING)
    TMR_stop( &g_TMR_calcA );
# endif

    return;
}


//------------------------------------------------------------------------------------------------------------
//! \brief  Calls the appropriate \c wave*_calcB() routine corresponding to the wave-type sweep algorithm as
//!         set by \pp{context->type}.
//!
//! \param[out]     B           Pointer to memory location at which to store the computed vector.
//! \param[in,out]  context     sweepContext_t object containing the context of the current sweep.
//! \param[in]      Q           angularFlux_t object containing the coefficients of the source term \f$ Q \f$
//!                             for the transport system.
//! \param[in]      q           Index of angular ordinate to compute vector for.
//! \param[in]      i           Index for the current spatial cell in \f$ x \f$ dimension.
//! \param[in]      j           Index for the current spatial cell in \f$ y \f$ dimension.
//! \param[in]      k           Index for the current spatial cell in \f$ z \f$ dimension.
//! \param[in]      WORK        Pointer to memory location used for working space for LAPACK routines. Should
//!                             be at least \pp{context->sizeofB} Bytes.
//!                             If \pp{context->type} is \c SWEEP_PLANEWAVE_LAPACK_LU or
//!                             \c SWEEP_WAVEFRONT_LAPACK_LU then this pointer is not used and may be null.
//!
//! \see    wave_calcA()
//! \see    wave_linearSolve()
//------------------------------------------------------------------------------------------------------------
static void wave_calcB (

    double * const restrict B,
    const sweepContext_t context,
    const angularFlux_t Q,
    const int64_t q,
    const int64_t i,
    const int64_t j,
    const int64_t k,
    void * const restrict WORK
) {

# if defined (DO_SWEEP_SUBTIMING)
    TMR_start( &g_TMR_calcB );
# endif

    switch ( context->type ) {

    case SWEEP_PLANEWAVE_LAPACK_LU:
    case SWEEP_WAVEFRONT_LAPACK_LU:
    case SWEEP_PLANEWAVE_GAUSSIAN_ELIMINATION:
    case SWEEP_WAVEFRONT_GAUSSIAN_ELIMINATION:

        waveLU_calcB( B, context, Q, q, i,j,k );
        goto cleanup;

    case SWEEP_PLANEWAVE_HESSENBERG:
    case SWEEP_WAVEFRONT_HESSENBERG:

        waveHess_calcB( B, context, Q, q, i,j,k, WORK );
        goto cleanup;

    default:

        PRINT_ERROR( "Invalid sweep type %d in %s.\n", context->type, __func__ )
        return;
    }

cleanup:

# if defined (DO_SWEEP_SUBTIMING)
    TMR_stop( &g_TMR_calcB );
# endif

    return;
}


//------------------------------------------------------------------------------------------------------------
//! \brief  Calls the linear solver corresponding to the wave-type sweep algorithm as set by
//!         \pp{context->type}.
//!
//! Calls the proper linear solver to solve the current linear system for each ordinate.
//!
//! \param[in]      context     sweepContext_t object containing the context of the current sweep.
//! \param[in]      A           Pointer to memory location containing matrix for linear system to be solved.
//! \param[in,out]  B           Pointer to memory location containing RHS vector for linear system to be
//!                             solved.
//!                             Upon return, contains the solution of the linear system.
//! \param[in]      q           Index of angular ordinate corresponding to current linear system.
//! \param[in]      WORK        Pointer to memory location used for working space for LAPACK routines. Should
//!                             be at least \pp{context->sizeofB} Bytes.
//!                             If \pp{context->type} is \c SWEEP_PLANEWAVE_LAPACK_LU or
//!                             \c SWEEP_WAVEFRONT_LAPACK_LU then this pointer is not used and may be null.
//!
//! \see    wave_calcA()
//! \see    wave_calcB()
//------------------------------------------------------------------------------------------------------------
static void wave_linearSolve (

    const sweepContext_t context,
    double * const restrict A,
    double * const restrict B,
    const int64_t q,
    void * const restrict WORK
) {

# if defined (DO_SWEEP_SUBTIMING)
    TMR_start( &g_TMR_linearSolve );
# endif

    switch ( context->type ) {

    case SWEEP_PLANEWAVE_LAPACK_LU:
    case SWEEP_WAVEFRONT_LAPACK_LU:
    {
        int N = (context->dgOrder + 1)*(context->dgOrder + 1)*(context->dgOrder + 1);
        int NRHS = 1;
        int INFO;

        dgesv_( &N, &NRHS, A, &N, WORK, B, &N, &INFO );

    # if defined (STRICT_CHECK)
        if ( INFO ) {

            PRINT_ERROR( "dgesv_ returned INFO = %d.\n", INFO )
            exit( EXIT_FAILURE );
        }
    # endif

    } goto cleanup;

    case SWEEP_PLANEWAVE_GAUSSIAN_ELIMINATION:
    case SWEEP_WAVEFRONT_GAUSSIAN_ELIMINATION:
    {
        const int N = (context->dgOrder + 1)*(context->dgOrder + 1)*(context->dgOrder + 1);

        GaussianElimination_rowMajor( N, A, B );

    } goto cleanup;

    case SWEEP_PLANEWAVE_HESSENBERG:
    case SWEEP_WAVEFRONT_HESSENBERG:
    {
        const int N = (context->dgOrder + 1)*(context->dgOrder + 1)*(context->dgOrder + 1);

        // Solve for PᵀΨ.
        upperHessenbergSolve_columnMajor( N, A, B, WORK );

        // Compute Ψ from PᵀΨ.
        multHouseholderReflections( B, context->P_wave[q], context->tau_wave[q], N );

    } goto cleanup;

    default:
    {
        PRINT_ERROR( "Invalid sweep type %d in %s.\n", context->type, __func__ )
    } return;
    }

cleanup:

# if defined (DO_SWEEP_SUBTIMING)
    TMR_stop( &g_TMR_linearSolve );
# endif

    return;
}


//============================================================================================================
//=== PLANEWAVE INTERFACE ROUTINES ===========================================================================
//============================================================================================================


//------------------------------------------------------------------------------------------------------------
//! \brief  Performs a planewave sweep for all ordinates in the first octant.
//!
//! Executes a planewave sweep across the entire spatial mesh for all angles in the first octant, i.e.,
//! all angles such that \f$ \xi, \eta, \mu > 0 \f$.
//!
//! \param[in]      context     sweepContext_t object containing the context of the current sweep.
//! \param[in,out]  Q           Initially contains the coefficients of the source term \f$ Q \f$ for the
//!                             transport system.
//!                             Upon return, contains the coefficients of \f$ \Psi \f$ which solves the
//!                             transport system.
//! \param[in]      sigma       Pointer to array of \f$ \sigma \f$ cross-sections.
//! \param[in]      dt          Timestep size.
//!
//! \see    planewave_sweep_octantTwo()
//! \see    planewave_sweep_octantThree()
//! \see    planewave_sweep_octantFour()
//! \see    planewave_sweep_octantFive()
//! \see    planewave_sweep_octantSix()
//! \see    planewave_sweep_octantSeven()
//! \see    planewave_sweep_octantEight()
//!
//! \see    wave_calcA()
//! \see    wave_calcB()
//! \see    wave_linearSolve()
//! \see    planewave_sweep()
//------------------------------------------------------------------------------------------------------------
static void planewave_sweep_octantOne (

    const sweepContext_t context,
    const angularFlux_t Q,
    const crossSection_t sigma,
    const double dt
) {

    const int64_t numDiags = context->nx + context->ny - 1;

    // --- Parallel work region. ---------------------------------------------------------------------- //

    # pragma omp parallel
    {

    # if defined (USE_ALIGNED_ALLOC)
        double * const restrict A    = aligned_alloc( ALIGNED_ALLOC_ALIGNMENT, context->sizeofA );
        double * const restrict B    = aligned_alloc( ALIGNED_ALLOC_ALIGNMENT, context->sizeofB );
        void   * const restrict WORK = aligned_alloc( ALIGNED_ALLOC_ALIGNMENT, context->sizeofB );
    # else
        double * const restrict A    = malloc( context->sizeofA );
        double * const restrict B    = malloc( context->sizeofB );
        void   * const restrict WORK = malloc( context->sizeofB );
    # endif

        for ( int64_t k = 1; k <= context->nz; ++k ) {
        for ( int64_t diag = 1; diag <= numDiags; ++diag ) {

            // --- Sweep diagonal - First quadrant: ξ,η > 0. ------------------------------------------ //

            int64_t diagLength = MIN( diag, MIN( numDiags - diag + 1, MIN( context->nx, context->ny ) ) );

            # pragma omp for collapse(2)
            for ( int64_t cell = 1; cell <= diagLength; ++cell ) {
            for ( int64_t q = 0; q < context->nq / 8; ++q ) {

                const int64_t i = MIN( diag, context->nx ) - cell + 1;
                const int64_t j = diag - i + 1;

                wave_calcA( A, context, sigma, dt, q, i,j,k );
                wave_calcB( B, context, Q, q, i,j,k, WORK );
                wave_linearSolve( context, A, B, q, WORK );

                for ( int64_t d = 0; d <= context->dgOrder; ++d ) {
                for ( int64_t e = 0; e <= context->dgOrder; ++e ) {
                for ( int64_t f = 0; f <= context->dgOrder; ++f ) {

                    IAF(Q,q,i,j,k,d,e,f) = B[ SWC_vectorIndexSingle(context,d,e,f) ];
                }}}
            }}
        }}

        free( A );
        free( B );
        free( WORK );
    }
}


//------------------------------------------------------------------------------------------------------------
//! \brief  Performs a wavefront sweep for all ordinates in the second octant.
//!
//! Executes a planewave sweep across the entire spatial mesh for all angles in the second octant, i.e.,
//! all angles such that \f$ \xi < 0 \f$ and \f$ \eta, \mu > 0 \f$.
//!
//! \param[in]      context     sweepContext_t object containing the context of the current sweep.
//! \param[in,out]  Q           Initially contains the coefficients of the source term \f$ Q \f$ for the
//!                             transport system.
//!                             Upon return, contains the coefficients of \f$ \Psi \f$ which solves the
//!                             transport system.
//! \param[in]      sigma       Pointer to array of \f$ \sigma \f$ cross-sections.
//! \param[in]      dt          Timestep size.
//!
//! \see    planewave_sweep_octantOne()
//! \see    planewave_sweep_octantThree()
//! \see    planewave_sweep_octantFour()
//! \see    planewave_sweep_octantFive()
//! \see    planewave_sweep_octantSix()
//! \see    planewave_sweep_octantSeven()
//! \see    planewave_sweep_octantEight()
//!
//! \see    wave_calcA()
//! \see    wave_calcB()
//! \see    wave_linearSolve()
//! \see    planewave_sweep()
//------------------------------------------------------------------------------------------------------------
static void planewave_sweep_octantTwo (

    const sweepContext_t context,
    const angularFlux_t Q,
    const crossSection_t sigma,
    const double dt
) {

    const int64_t numDiags = context->nx + context->ny - 1;

    // --- Parallel work region. ---------------------------------------------------------------------- //

    # pragma omp parallel
    {

    # if defined (USE_ALIGNED_ALLOC)
        double * const restrict A    = aligned_alloc( ALIGNED_ALLOC_ALIGNMENT, context->sizeofA );
        double * const restrict B    = aligned_alloc( ALIGNED_ALLOC_ALIGNMENT, context->sizeofB );
        void   * const restrict WORK = aligned_alloc( ALIGNED_ALLOC_ALIGNMENT, context->sizeofB );
    # else
        double * const restrict A    = malloc( context->sizeofA );
        double * const restrict B    = malloc( context->sizeofB );
        void   * const restrict WORK = malloc( context->sizeofB );
    # endif

        for ( int64_t k = 1; k <= context->nz; ++k ) {
        for ( int64_t diag = 1; diag <= numDiags; ++diag ) {

            // --- Sweep diagonal - Second quadrant: ξ < 0, η > 0. ------------------------------------ //

            int64_t diagLength = MIN( diag, MIN( numDiags - diag + 1, MIN( context->nx, context->ny ) ) );

            # pragma omp for collapse(2)
            for ( int64_t cell = 1; cell <= diagLength; ++cell ) {
            for ( int64_t q = context->nq / 2; q < 5 * context->nq / 8; ++q ) {

                const int64_t i = context->nx - MIN( diag, context->nx ) + cell;
                const int64_t j = diag - context->nx + i;

                wave_calcA( A, context, sigma, dt, q, i,j,k );
                wave_calcB( B, context, Q, q, i,j,k, WORK );
                wave_linearSolve( context, A, B, q, WORK );

                for ( int64_t d = 0; d <= context->dgOrder; ++d ) {
                for ( int64_t e = 0; e <= context->dgOrder; ++e ) {
                for ( int64_t f = 0; f <= context->dgOrder; ++f ) {

                    IAF(Q,q,i,j,k,d,e,f) = B[ SWC_vectorIndexSingle(context,d,e,f) ];
                }}}
            }}
        }}

        free( A );
        free( B );
        free( WORK );
    }
}


//------------------------------------------------------------------------------------------------------------
//! \brief  Performs a wavefront sweep for all ordinates in the third octant.
//!
//! Executes a wavefront sweep across the entire spatial mesh for all angles in the third octant, i.e.,
//! all angles such that \f$ \xi, \eta < 0 \f$ and \f$ \mu > 0 \f$.
//!
//! \param[in]      context     sweepContext_t object containing the context of the current sweep.
//! \param[in,out]  Q           Initially contains the coefficients of the source term \f$ Q \f$ for the
//!                             transport system.
//!                             Upon return, contains the coefficients of \f$ \Psi \f$ which solves the
//!                             transport system.
//! \param[in]      sigma       Pointer to array of \f$ \sigma \f$ cross-sections.
//! \param[in]      dt          Timestep size.
//!
//! \see    planewave_sweep_octantOne()
//! \see    planewave_sweep_octantTwo()
//! \see    planewave_sweep_octantFour()
//! \see    planewave_sweep_octantFive()
//! \see    planewave_sweep_octantSix()
//! \see    planewave_sweep_octantSeven()
//! \see    planewave_sweep_octantEight()
//!
//! \see    wave_calcA()
//! \see    wave_calcB()
//! \see    wave_linearSolve()
//! \see    planewave_sweep()
//------------------------------------------------------------------------------------------------------------
static void planewave_sweep_octantThree (

    const sweepContext_t context,
    const angularFlux_t Q,
    const crossSection_t sigma,
    const double dt
) {

    const int64_t numDiags = context->nx + context->ny - 1;

    // --- Parallel work region. ---------------------------------------------------------------------- //

    # pragma omp parallel
    {

    # if defined (USE_ALIGNED_ALLOC)
        double * const restrict A    = aligned_alloc( ALIGNED_ALLOC_ALIGNMENT, context->sizeofA );
        double * const restrict B    = aligned_alloc( ALIGNED_ALLOC_ALIGNMENT, context->sizeofB );
        void   * const restrict WORK = aligned_alloc( ALIGNED_ALLOC_ALIGNMENT, context->sizeofB );
    # else
        double * const restrict A    = malloc( context->sizeofA );
        double * const restrict B    = malloc( context->sizeofB );
        void   * const restrict WORK = malloc( context->sizeofB );
    # endif

        for ( int64_t k = 1; k <= context->nz; ++k ) {
        for ( int64_t diag = 1; diag <= numDiags; ++diag ) {

            // --- Sweep diagonal - Third quadrant: ξ,η < 0. ------------------------------------------ //

            int64_t diagLength = MIN( diag, MIN( numDiags - diag + 1, MIN( context->nx, context->ny ) ) );

            # pragma omp for collapse(2)
            for ( int64_t cell = 1; cell <= diagLength; ++cell ) {
            for ( int64_t q = 3 * context->nq / 4; q < 7 * context->nq / 8; ++q ) {

                const int64_t i = context->nx - MIN( diag, context->nx ) + cell;
                const int64_t j = numDiags - diag - i + 2;

                wave_calcA( A, context, sigma, dt, q, i,j,k );
                wave_calcB( B, context, Q, q, i,j,k, WORK );
                wave_linearSolve( context, A, B, q, WORK );

                for ( int64_t d = 0; d <= context->dgOrder; ++d ) {
                for ( int64_t e = 0; e <= context->dgOrder; ++e ) {
                for ( int64_t f = 0; f <= context->dgOrder; ++f ) {

                    IAF(Q,q,i,j,k,d,e,f) = B[ SWC_vectorIndexSingle(context,d,e,f) ];
                }}}
            }}
        }}

        free( A );
        free( B );
        free( WORK );
    }
}


//------------------------------------------------------------------------------------------------------------
//! \brief  Performs a wavefront sweep for all ordinates in the fourth octant.
//!
//! Executes a wavefront sweep across the entire spatial mesh for all angles in the fourth octant, i.e.,
//! all angles such that \f$ \xi, \mu > 0 \f$ and \f$ \eta < 0 \f$.
//!
//! \param[in]      context     sweepContext_t object containing the context of the current sweep.
//! \param[in,out]  Q           Initially contains the coefficients of the source term \f$ Q \f$ for the
//!                             transport system.
//!                             Upon return, contains the coefficients of \f$ \Psi \f$ which solves the
//!                             transport system.
//! \param[in]      sigma       Pointer to array of \f$ \sigma \f$ cross-sections.
//! \param[in]      dt          Timestep size.
//!
//! \see    planewave_sweep_octantOne()
//! \see    planewave_sweep_octantTwo()
//! \see    planewave_sweep_octantThree()
//! \see    planewave_sweep_octantFive()
//! \see    planewave_sweep_octantSix()
//! \see    planewave_sweep_octantSeven()
//! \see    planewave_sweep_octantEight()
//!
//! \see    wave_calcA()
//! \see    wave_calcB()
//! \see    wave_linearSolve()
//! \see    planewave_sweep()
//------------------------------------------------------------------------------------------------------------
static void planewave_sweep_octantFour (

    const sweepContext_t context,
    const angularFlux_t Q,
    const crossSection_t sigma,
    const double dt
) {

    const int64_t numDiags = context->nx + context->ny - 1;

    // --- Parallel work region. ---------------------------------------------------------------------- //

    # pragma omp parallel
    {

    # if defined (USE_ALIGNED_ALLOC)
        double * const restrict A    = aligned_alloc( ALIGNED_ALLOC_ALIGNMENT, context->sizeofA );
        double * const restrict B    = aligned_alloc( ALIGNED_ALLOC_ALIGNMENT, context->sizeofB );
        void   * const restrict WORK = aligned_alloc( ALIGNED_ALLOC_ALIGNMENT, context->sizeofB );
    # else
        double * const restrict A    = malloc( context->sizeofA );
        double * const restrict B    = malloc( context->sizeofB );
        void   * const restrict WORK = malloc( context->sizeofB );
    # endif

        for ( int64_t k = 1; k <= context->nz; ++k ) {
        for ( int64_t diag = 1; diag <= numDiags; ++diag ) {

            // --- Sweep diagonal - Fourth quadrant: ξ > 0, η < 0. ------------------------------------ //

            int64_t diagLength = MIN( diag, MIN( numDiags - diag + 1, MIN( context->nx, context->ny ) ) );

            # pragma omp for collapse(2)
            for ( int64_t cell = 1; cell <= diagLength; ++cell ) {
            for ( int64_t q = context->nq / 4; q < 3 * context->nq / 8; ++q ) {

                const int64_t i = MIN( diag, context->nx ) - cell + 1;
                const int64_t j = context->ny - diag + i;

                wave_calcA( A, context, sigma, dt, q, i,j,k );
                wave_calcB( B, context, Q, q, i,j,k, WORK );
                wave_linearSolve( context, A, B, q, WORK );

                // Save result.
                for ( int64_t d = 0; d <= context->dgOrder; ++d ) {
                for ( int64_t e = 0; e <= context->dgOrder; ++e ) {
                for ( int64_t f = 0; f <= context->dgOrder; ++f ) {

                    IAF(Q,q,i,j,k,d,e,f) = B[ SWC_vectorIndexSingle(context,d,e,f) ];
                }}}
            }}
        }}

        free( A );
        free( B );
        free( WORK );
    }
}


//------------------------------------------------------------------------------------------------------------
//! \brief  Performs a planewave sweep for all ordinates in the fifth octant.
//!
//! Executes a wavefront sweep across the entire spatial mesh for all angles in the fifth octant, i.e.,
//! all angles such that \f$ \xi, \eta > 0 \f$ and \f$ \mu < 0 \f$.
//!
//! \param[in]      context     sweepContext_t object containing the context of the current sweep.
//! \param[in,out]  Q           Initially contains the coefficients of the source term \f$ Q \f$ for the
//!                             transport system.
//!                             Upon return, contains the coefficients of \f$ \Psi \f$ which solves the
//!                             transport system.
//! \param[in]      sigma       Pointer to array of \f$ \sigma \f$ cross-sections.
//! \param[in]      dt          Timestep size.
//!
//! \see    planewave_sweep_octantOne()
//! \see    planewave_sweep_octantTwo()
//! \see    planewave_sweep_octantThree()
//! \see    planewave_sweep_octantFour()
//! \see    planewave_sweep_octantSix()
//! \see    planewave_sweep_octantSeven()
//! \see    planewave_sweep_octantEight()
//!
//! \see    wave_calcA()
//! \see    wave_calcB()
//! \see    wave_linearSolve()
//! \see    planewave_sweep()
//------------------------------------------------------------------------------------------------------------
static void planewave_sweep_octantFive (

    const sweepContext_t context,
    const angularFlux_t Q,
    const crossSection_t sigma,
    const double dt
) {

    const int64_t numDiags = context->nx + context->ny - 1;


    // --- Parallel work region. ---------------------------------------------------------------------- //

    # pragma omp parallel
    {

    # if defined (USE_ALIGNED_ALLOC)
        double * const restrict A    = aligned_alloc( ALIGNED_ALLOC_ALIGNMENT, context->sizeofA );
        double * const restrict B    = aligned_alloc( ALIGNED_ALLOC_ALIGNMENT, context->sizeofB );
        void   * const restrict WORK = aligned_alloc( ALIGNED_ALLOC_ALIGNMENT, context->sizeofB );
    # else
        double * const restrict A    = malloc( context->sizeofA );
        double * const restrict B    = malloc( context->sizeofB );
        void   * const restrict WORK = malloc( context->sizeofB );
    # endif

        for ( int64_t k = context->nz; k >= 1; --k ) {
        for ( int64_t diag = 1; diag <= numDiags; ++diag ) {

            // --- Sweep diagonal - First quadrant: ξ,η > 0. ------------------------------------------ //

            int64_t diagLength = MIN( diag, MIN( numDiags - diag + 1, MIN( context->nx, context->ny ) ) );

            # pragma omp for collapse(2)
            for ( int64_t cell = 1; cell <= diagLength; ++cell ) {
            for ( int64_t q = context->nq / 8; q <  context->nq / 4; ++q ) {

                const int64_t i = MIN( diag, context->nx ) - cell + 1;
                const int64_t j = diag - i + 1;

                wave_calcA( A, context, sigma, dt, q, i,j,k );
                wave_calcB( B, context, Q, q, i,j,k, WORK );
                wave_linearSolve( context, A, B, q, WORK );

                for ( int64_t d = 0; d <= context->dgOrder; ++d ) {
                for ( int64_t e = 0; e <= context->dgOrder; ++e ) {
                for ( int64_t f = 0; f <= context->dgOrder; ++f ) {

                    IAF(Q,q,i,j,k,d,e,f) = B[ SWC_vectorIndexSingle(context,d,e,f) ];
                }}}
            }}
        }}

        free( A );
        free( B );
        free( WORK );
    }
}


//------------------------------------------------------------------------------------------------------------
//! \brief  Performs a wavefront sweep for all ordinates in the sixth octant.
//!
//! Executes a planewave sweep across the entire spatial mesh for all angles in the sixth octant, i.e.,
//! all angles such that \f$ \xi, \mu < 0 \f$ and \f$ \eta > 0 \f$.
//!
//! \param[in]      context     sweepContext_t object containing the context of the current sweep.
//! \param[in,out]  Q           Initially contains the coefficients of the source term \f$ Q \f$ for the
//!                             transport system.
//!                             Upon return, contains the coefficients of \f$ \Psi \f$ which solves the
//!                             transport system.
//! \param[in]      sigma       Pointer to array of \f$ \sigma \f$ cross-sections.
//! \param[in]      dt          Timestep size.
//!
//! \see    planewave_sweep_octantOne()
//! \see    planewave_sweep_octantTwo()
//! \see    planewave_sweep_octantThree()
//! \see    planewave_sweep_octantFour()
//! \see    planewave_sweep_octantFive()
//! \see    planewave_sweep_octantSeven()
//! \see    planewave_sweep_octantEight()
//!
//! \see    wave_calcA()
//! \see    wave_calcB()
//! \see    wave_linearSolve()
//! \see    planewave_sweep()
//------------------------------------------------------------------------------------------------------------
static void planewave_sweep_octantSix (

    const sweepContext_t context,
    const angularFlux_t Q,
    const crossSection_t sigma,
    const double dt
) {

    const int64_t numDiags = context->nx + context->ny - 1;

    // --- Parallel work region. ---------------------------------------------------------------------- //

    # pragma omp parallel
    {

    # if defined (USE_ALIGNED_ALLOC)
        double * const restrict A    = aligned_alloc( ALIGNED_ALLOC_ALIGNMENT, context->sizeofA );
        double * const restrict B    = aligned_alloc( ALIGNED_ALLOC_ALIGNMENT, context->sizeofB );
        void   * const restrict WORK = aligned_alloc( ALIGNED_ALLOC_ALIGNMENT, context->sizeofB );
    # else
        double * const restrict A    = malloc( context->sizeofA );
        double * const restrict B    = malloc( context->sizeofB );
        void   * const restrict WORK = malloc( context->sizeofB );
    # endif

        for ( int64_t k = context->nz; k >= 1; --k ) {
        for ( int64_t diag = 1; diag <= numDiags; ++diag ) {

            // --- Sweep diagonal - Second quadrant: ξ < 0, η > 0. ------------------------------------ //

            int64_t diagLength = MIN( diag, MIN( numDiags - diag + 1, MIN( context->nx, context->ny ) ) );

            # pragma omp for collapse(2)
            for ( int64_t cell = 1; cell <= diagLength; ++cell ) {
            for ( int64_t q = 5 * context->nq / 8; q < 3 * context->nq / 4; ++q ) {

                const int64_t i = context->nx - MIN( diag, context->nx ) + cell;
                const int64_t j = diag - context->nx + i;

                wave_calcA( A, context, sigma, dt, q, i,j,k );
                wave_calcB( B, context, Q, q, i,j,k, WORK );
                wave_linearSolve( context, A, B, q, WORK );

                for ( int64_t d = 0; d <= context->dgOrder; ++d ) {
                for ( int64_t e = 0; e <= context->dgOrder; ++e ) {
                for ( int64_t f = 0; f <= context->dgOrder; ++f ) {

                    IAF(Q,q,i,j,k,d,e,f) = B[ SWC_vectorIndexSingle(context,d,e,f) ];
                }}}
            }}
        }}

        free( A );
        free( B );
        free( WORK );
    }
}


//------------------------------------------------------------------------------------------------------------
//! \brief  Performs a wavefront sweep for all ordinates in the seventh octant.
//!
//! Executes a wavefront sweep across the entire spatial mesh for all angles in the seventh octant, i.e.,
//! all angles such that \f$ \xi, \eta, \mu < 0 \f$.
//!
//! \param[in]      context     sweepContext_t object containing the context of the current sweep.
//! \param[in,out]  Q           Initially contains the coefficients of the source term \f$ Q \f$ for the
//!                             transport system.
//!                             Upon return, contains the coefficients of \f$ \Psi \f$ which solves the
//!                             transport system.
//! \param[in]      sigma       Pointer to array of \f$ \sigma \f$ cross-sections.
//! \param[in]      dt          Timestep size.
//!
//! \see    planewave_sweep_octantOne()
//! \see    planewave_sweep_octantTwo()
//! \see    planewave_sweep_octantThree()
//! \see    planewave_sweep_octantFour()
//! \see    planewave_sweep_octantFive()
//! \see    planewave_sweep_octantSix()
//! \see    planewave_sweep_octantEight()
//!
//! \see    wave_calcA()
//! \see    wave_calcB()
//! \see    wave_linearSolve()
//! \see    planewave_sweep()
//------------------------------------------------------------------------------------------------------------
static void planewave_sweep_octantSeven (

    const sweepContext_t context,
    const angularFlux_t Q,
    const crossSection_t sigma,
    const double dt
) {

    const int64_t numDiags = context->nx + context->ny - 1;

    // --- Parallel work region. ---------------------------------------------------------------------- //

    # pragma omp parallel
    {

    # if defined (USE_ALIGNED_ALLOC)
        double * const restrict A    = aligned_alloc( ALIGNED_ALLOC_ALIGNMENT, context->sizeofA );
        double * const restrict B    = aligned_alloc( ALIGNED_ALLOC_ALIGNMENT, context->sizeofB );
        void   * const restrict WORK = aligned_alloc( ALIGNED_ALLOC_ALIGNMENT, context->sizeofB );
    # else
        double * const restrict A    = malloc( context->sizeofA );
        double * const restrict B    = malloc( context->sizeofB );
        void   * const restrict WORK = malloc( context->sizeofB );
    # endif

        for ( int64_t k = context->nz; k >= 1; --k ) {
        for ( int64_t diag = 1; diag <= numDiags; ++diag ) {

            // --- Sweep diagonal - Third quadrant: ξ,η < 0. ------------------------------------------ //

            int64_t diagLength = MIN( diag, MIN( numDiags - diag + 1, MIN( context->nx, context->ny ) ) );

            # pragma omp for collapse(2)
            for ( int64_t cell = 1; cell <= diagLength; ++cell ) {
            for ( int64_t q = 7 * context->nq / 8; q < context->nq; ++q ) {

                const int64_t i = context->nx - MIN( diag, context->nx ) + cell;
                const int64_t j = numDiags - diag - i + 2;

                wave_calcA( A, context, sigma, dt, q, i,j,k );
                wave_calcB( B, context, Q, q, i,j,k, WORK );
                wave_linearSolve( context, A, B, q, WORK );

                for ( int64_t d = 0; d <= context->dgOrder; ++d ) {
                for ( int64_t e = 0; e <= context->dgOrder; ++e ) {
                for ( int64_t f = 0; f <= context->dgOrder; ++f ) {

                    IAF(Q,q,i,j,k,d,e,f) = B[ SWC_vectorIndexSingle(context,d,e,f) ];
                }}}
            }}
        }}

        free( A );
        free( B );
        free( WORK );
    }
}


//------------------------------------------------------------------------------------------------------------
//! \brief  Performs a wavefront sweep for all ordinates in the eighth octant.
//!
//! Executes a wavefront sweep across the entire spatial mesh for all angles in the eighth octant, i.e.,
//! all angles such that \f$ \xi > 0 \f$ and \f$ \eta, \mu < 0 \f$.
//!
//! \param[in]      context     sweepContext_t object containing the context of the current sweep.
//! \param[in,out]  Q           Initially contains the coefficients of the source term \f$ Q \f$ for the
//!                             transport system.
//!                             Upon return, contains the coefficients of \f$ \Psi \f$ which solves the
//!                             transport system.
//! \param[in]      sigma       Pointer to array of \f$ \sigma \f$ cross-sections.
//! \param[in]      dt          Timestep size.
//!
//! \see    planewave_sweep_octantOne()
//! \see    planewave_sweep_octantTwo()
//! \see    planewave_sweep_octantThree()
//! \see    planewave_sweep_octantFour()
//! \see    planewave_sweep_octantFive()
//! \see    planewave_sweep_octantSix()
//! \see    planewave_sweep_octantSeven()
//!
//! \see    wave_calcA()
//! \see    wave_calcB()
//! \see    wave_linearSolve()
//! \see    planewave_sweep()
//------------------------------------------------------------------------------------------------------------
static void planewave_sweep_octantEight (

    const sweepContext_t context,
    const angularFlux_t Q,
    const crossSection_t sigma,
    const double dt
) {

    const int64_t numDiags = context->nx + context->ny - 1;

    // --- Parallel work region. ---------------------------------------------------------------------- //

    # pragma omp parallel
    {

    # if defined (USE_ALIGNED_ALLOC)
        double * const restrict A    = aligned_alloc( ALIGNED_ALLOC_ALIGNMENT, context->sizeofA );
        double * const restrict B    = aligned_alloc( ALIGNED_ALLOC_ALIGNMENT, context->sizeofB );
        void   * const restrict WORK = aligned_alloc( ALIGNED_ALLOC_ALIGNMENT, context->sizeofB );
    # else
        double * const restrict A    = malloc( context->sizeofA );
        double * const restrict B    = malloc( context->sizeofB );
        void   * const restrict WORK = malloc( context->sizeofB );
    # endif

        for ( int64_t k = context->nz; k >= 1; --k ) {
        for ( int64_t diag = 1; diag <= numDiags; ++diag ) {

            // --- Sweep diagonal - Fourth quadrant: ξ > 0, η < 0. ------------------------------------ //

            int64_t diagLength = MIN( diag, MIN( numDiags - diag + 1, MIN( context->nx, context->ny ) ) );

            # pragma omp for collapse(2)
            for ( int64_t cell = 1; cell <= diagLength; ++cell ) {
            for ( int64_t q = 3 * context->nq / 8; q < context->nq / 2; ++q ) {

                const int64_t i = MIN( diag, context->nx ) - cell + 1;
                const int64_t j = context->ny - diag + i;

                wave_calcA( A, context, sigma, dt, q, i,j,k );
                wave_calcB( B, context, Q, q, i,j,k, WORK );
                wave_linearSolve( context, A, B, q, WORK );

                // Save result.
                for ( int64_t d = 0; d <= context->dgOrder; ++d ) {
                for ( int64_t e = 0; e <= context->dgOrder; ++e ) {
                for ( int64_t f = 0; f <= context->dgOrder; ++f ) {

                    IAF(Q,q,i,j,k,d,e,f) = B[ SWC_vectorIndexSingle(context,d,e,f) ];
                }}}
            }}
        }}

        free( A );
        free( B );
        free( WORK );
    }
}


//------------------------------------------------------------------------------------------------------------
//! \brief  Interface function for performing planewave sweeps.
//!
//! \param[in]      context     sweepContext_t object containing the context of the current sweep.
//! \param[in,out]  Q           Initially contains the coefficients of the source term \f$ Q \f$ for the
//!                             transport system.
//!                             Upon return, contains the coefficients of \f$ \Psi \f$ which solves the
//!                             transport system.
//! \param[in]      sigma       Pointer to array of \f$ \sigma \f$ cross-sections.
//! \param[in]      dt          Timestep size.
//!
//! \see    planewave_sweep_octantOne()
//! \see    planewave_sweep_octantTwo()
//! \see    planewave_sweep_octantThree()
//! \see    planewave_sweep_octantFour()
//! \see    planewave_sweep_octantFive()
//! \see    planewave_sweep_octantSix()
//! \see    planewave_sweep_octantSeven()
//! \see    planewave_sweep_octantEight()
//------------------------------------------------------------------------------------------------------------
int planewave_sweep (

    const sweepContext_t context,
    const angularFlux_t Q,
    const crossSection_t sigma,
    const double dt
) {

    planewave_sweep_octantOne( context, Q, sigma, dt );
    planewave_sweep_octantTwo( context, Q, sigma, dt );
    planewave_sweep_octantThree( context, Q, sigma, dt );
    planewave_sweep_octantFour( context, Q, sigma, dt );

    planewave_sweep_octantFive( context, Q, sigma, dt );
    planewave_sweep_octantSix( context, Q, sigma, dt );
    planewave_sweep_octantSeven( context, Q, sigma, dt );
    planewave_sweep_octantEight( context, Q, sigma, dt );
}


//============================================================================================================
//=== WAVEFRONT INTERFACE ROUTINES ===========================================================================
//============================================================================================================


//------------------------------------------------------------------------------------------------------------
//! \brief  Interface function for performing wavefront sweeps.
//!
//! \param[in]      context     sweepContext_t object containing the context of the current sweep.
//! \param[in,out]  Q           Initially contains the coefficients of the source term \f$ Q \f$ for the
//!                             transport system.
//!                             Upon return, contains the coefficients of \f$ \Psi \f$ which solves the
//!                             transport system.
//! \param[in]      sigma       Pointer to array of \f$ \sigma \f$ cross-sections.
//! \param[in]      dt          Timestep size.
//!
//------------------------------------------------------------------------------------------------------------
int wavefront_sweep (

    const sweepContext_t context,
    const angularFlux_t Q,
    const crossSection_t sigma,
    const double dt
) {

}


# endif // if SPACE_DIMS == 3
