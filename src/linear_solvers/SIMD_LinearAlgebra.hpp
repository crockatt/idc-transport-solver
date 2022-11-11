//------------------------------------------------------------------------------------------------------------
//! \file   linear_solvers/SIMD_LinearAlgebra.hpp
//! \brief  Contains templates of SIMD blocked linear algebra routines used by transport sweeps.
//!
//! \author Michael M. Crockatt
//! \date   March 2018
//------------------------------------------------------------------------------------------------------------

//
// Clang cannot find cstdalign, but will find stdalign.h. Prefer C++ style cstdalign for other compilers.
//
# if defined (__clang__)
    # include <stdalign.h>
# else
    # include <cstdalign>
# endif

# include "utils/SIMD.hpp"


# ifndef MAX_GE_BLOCK_DIM
    # define MAX_GE_BLOCK_DIM 16
# endif

# ifndef MAX_HESS_BLOCK_DIM
    # define MAX_HESS_BLOCK_DIM 16
# endif

# ifdef MIN
    # undef MIN
# endif
# define MIN(a,b) ( (a) < (b) ? (a) : (b) )


//------------------------------------------------------------------------------------------------------------
//! \brief  Solves a set of stacked linear systems \f$ Ax = b \f$ using Gaussian elimination without pivoting.
//!
//! \attention  The matrices \f$ A \f$ are assumed to be stored in _column-major_ ordering.
//!
//! \param[in]      N           Dimension of the systems to be solved.
//! \param[in]      A           Pointer to the matrices..
//! \param[in,out]  b           Initially contains the known vectors \f$ b \f$ for the linear systems.
//!                             Upon return contains the solution vectors \f$ x \f$ of the linear systems.
//------------------------------------------------------------------------------------------------------------
template< int SIMD_length >
void SIMD_GaussianElimination_CM (

    const int N,
    double * const restrict A,
    double * const restrict b
) {

//
// \brief  Macro for indexing into SIMD matrix arrays. Assumes column-major storage.
//
// \param[in]   i   Row index of element in matrix.
// \param[in]   j   Column index of element in matrix.
// \param[in]   l   Index for system in SIMD block. In [ 0, SIMD_length ).
//
# define IMAT(i,j,l) ( (l) + SIMD_length*( (i) + N*(j) ) )

//
// \brief  Macro for indexing into SIMD vector arrays.
//
// \param[in]  i   Index of element in vector.
// \param[in]  l   Index for system in SIMD block. IN [ 0, SIMD_length ).
//
# define IVEC(i,l) ( (l) + SIMD_length*(i) )

# if defined (GE_BLOCK)

    if ( N <= MAX_GE_BLOCK_DIM )

# endif // if defined (GE_BLOCK)

    {
        // Forward elimination.
        for ( int j = 0; j < N; ++j ) {

            // Compute L_{i,j} and apply reductions to RHS vector.
            for ( int i = j+1; i < N; ++i ) {

                # pragma omp simd aligned( A : SIMD_length * SIZEOF_DOUBLE )
                for ( int l = 0; l < SIMD_length; ++l )
                    A[ IMAT(i,j,l) ] /= - A[ IMAT(j,j,l) ];

                # pragma omp simd aligned( A, b : SIMD_length * SIZEOF_DOUBLE )
                for ( int l = 0; l < SIMD_length; ++l )
                    b[ IVEC(i,l) ] += A[ IMAT(i,j,l) ] * b[ IVEC(j,l) ];
            }

            // Reduce rows.
            for ( int k = j+1; k < N; ++k ) {
            for ( int i = j+1; i < N; ++i ) {
            # pragma omp simd aligned( A : SIMD_length * SIZEOF_DOUBLE )
            for ( int l = 0; l < SIMD_length; ++l ) {

                A[ IMAT(i,k,l) ] += A[ IMAT(i,j,l) ] * A[ IMAT(j,k,l) ];
            }}}
        }

        // Backward substitution.
        for ( int i = N-1; i >= 0; --i ) {

            # pragma omp simd aligned( A, b : SIMD_length * SIZEOF_DOUBLE )
            for ( int l = 0; l < SIMD_length; ++l )
                b[ IVEC(i,l) ] /= A[ IMAT(i,i,l) ];

            for ( int j = 0; j < i; ++j ) {
            # pragma omp simd aligned( A, b : SIMD_length * SIZEOF_DOUBLE )
            for ( int l = 0; l < SIMD_length; ++l ) {

                b[ IVEC(j,l) ] -= A[ IMAT(j,i,l) ] * b[ IVEC(i,l) ];
            }}
        }

        return;
    }

# if defined (GE_BLOCK)

    else {

        const int num_blocks = 1 + ((N - 1) / MAX_GE_BLOCK_DIM);   // if N != 0.
        const int block_dim = N / num_blocks + ( N % num_blocks ? 1 : 0 );

        // --- Forward elimination. ----------------------------------------------------------------------- //

        for ( int blk_i = 0; blk_i < num_blocks; ++blk_i ) {

            const int blk_i_min = blk_i * block_dim;
            const int blk_i_max = MIN( N, (blk_i + 1) * block_dim );

            // Compute LU factorization of diagonal block.
            for ( int i = blk_i_min; i < blk_i_max; ++i ) {

                // Compute L_{j,i} and apply reductions to RHS vector.
                for ( int j = i+1; j < blk_i_max; ++j ) {

                    # pragma omp simd aligned( A : SIMD_length * SIZEOF_DOUBLE )
                    for ( int l = 0; l < SIMD_length; ++l )
                        A[ IMAT(j,i,l) ] /= - A[ IMAT(i,i,l) ];

                    // Apply reductions to RHS vector.
                    # pragma omp simd aligned( A, b : SIMD_length * SIZEOF_DOUBLE )
                    for ( int l = 0; l < SIMD_length; ++l )
                        b[ IVEC(j,l) ] += A[ IMAT(j,i,l) ] * b[ IVEC(i,l) ];
                }

                // Reduce rows.
                for ( int k = i+1; k < blk_i_max; ++k ) {
                for ( int j = i+1; j < blk_i_max; ++j ) {
                # pragma omp simd aligned( A : SIMD_length * SIZEOF_DOUBLE )
                for ( int l = 0; l < SIMD_length; ++l ) {

                    A[ IMAT(j,k,l) ] += A[ IMAT(j,i,l) ] * A[ IMAT(i,k,l) ];
                }}}
            }

            // Reduce block row.
            for ( int blk_k = blk_i +1; blk_k < num_blocks; ++blk_k ) {

                const int blk_k_min = blk_k * block_dim;
                const int blk_k_max = MIN( N, (blk_k + 1) * block_dim );

                for ( int k = blk_k_min; k < blk_k_max; ++k ) {
                for ( int i = blk_i_min; i < blk_i_max; ++i ) {
                for ( int j = i+1;       j < blk_i_max; ++j ) {
                # pragma omp simd aligned( A : SIMD_length * SIZEOF_DOUBLE )
                for ( int l = 0; l < SIMD_length; ++l ) {

                    A[ IMAT(j,k,l) ] += A[ IMAT(j,i,l) ] * A[ IMAT(i,k,l) ];
                }}}}
            }

            // Reduce block column.
            for ( int blk_j = blk_i +1; blk_j < num_blocks; ++blk_j ) {

                const int blk_j_min = blk_j * block_dim;
                const int blk_j_max = MIN( N, (blk_j + 1) * block_dim );

                for ( int i = blk_i_min; i < blk_i_max; ++i ) {

                    // Compute L_{j,i} and apply reductions to RHS vector.
                    for ( int j = blk_j_min; j < blk_j_max; ++j ) {

                        # pragma omp simd aligned( A : SIMD_length * SIZEOF_DOUBLE )
                        for ( int l = 0; l < SIMD_length; ++l )
                            A[ IMAT(j,i,l) ] /= - A[ IMAT(i,i,l) ];

                        # pragma omp simd aligned( A, b : SIMD_length * SIZEOF_DOUBLE )
                        for ( int l = 0; l < SIMD_length; ++l )
                            b[ IVEC(j,l) ] += A[ IMAT(j,i,l) ] * b[ IVEC(i,l) ];
                    }

                    // Reduce row.
                    for ( int k = i+1; k < blk_i_max; ++k ) {
                    for ( int j = blk_j_min; j < blk_j_max; ++j ) {
                    # pragma omp simd aligned( A : SIMD_length * SIZEOF_DOUBLE )
                    for ( int l = 0; l < SIMD_length; ++l ) {

                        A[ IMAT(j,k,l) ] += A[ IMAT(j,i,l) ] * A[ IMAT(i,k,l) ];
                    }}}
                }
            }

            // Compute Schur complement.
            for ( int blk_k = blk_i +1; blk_k < num_blocks; ++blk_k ) {

                const int blk_k_min = blk_k * block_dim;
                const int blk_k_max = MIN( N, (blk_k + 1) * block_dim );

                for ( int blk_j = blk_i +1; blk_j < num_blocks; ++blk_j ) {

                    const int blk_j_min = blk_j * block_dim;
                    const int blk_j_max = MIN( N, (blk_j + 1) * block_dim );

                    for ( int k = blk_k_min; k < blk_k_max; ++k ) {
                    for ( int i = blk_i_min; i < blk_i_max; ++i ) {
                    for ( int j = blk_j_min; j < blk_j_max; ++j ) {
                    # pragma omp simd aligned( A : SIMD_length * SIZEOF_DOUBLE )
                    for ( int l = 0; l < SIMD_length; ++l ) {

                        A[ IMAT(j,k,l) ] += A[ IMAT(j,i,l) ] * A[ IMAT(i,k,l) ];
                    }}}}
                }
            }
        }

        // --- Backward substitution. --------------------------------------------------------------------- //

        // For each block.
        for ( int blk_i = num_blocks -1; blk_i >= 0; --blk_i ) {

            const int blk_i_min = blk_i * block_dim;
            const int blk_i_max = MIN( N, (blk_i + 1) * block_dim );

            // Perform backward substitution on block.
            for ( int i = blk_i_max -1; i >= blk_i_min; --i ) {

                # pragma omp simd aligned( A, b : SIMD_length * SIZEOF_DOUBLE )
                for ( int l = 0; l < SIMD_length; ++l )
                    b[ IVEC(i,l) ] /= A[ IMAT(i,i,l) ];

                for ( int j = blk_i_min; j < i; ++j ) {
                # pragma omp simd aligned( A, b : SIMD_length * SIZEOF_DOUBLE )
                for ( int l = 0; l < SIMD_length; ++l ) {

                    b[ IVEC(j,l) ] -= A[ IMAT(j,i,l) ] * b[ IVEC(i,l) ];
                }}

            }

            // Reduce column block of A.
            for ( int blk_k = 0; blk_k < blk_i; ++blk_k ) {
            for ( int i = blk_i_min; i < blk_i_max; ++i ) {
            for ( int k = blk_k * block_dim; k < (blk_k + 1) * block_dim; ++k ) {

                # pragma omp simd aligned( A, b : SIMD_length * SIZEOF_DOUBLE )
                for ( int l = 0; l < SIMD_length; ++l )
                    b[ IVEC(k,l) ] -= A[ IMAT(k,i,l) ] * b[ IVEC(i,l) ];
            }}}
        }
    }

# endif // if defined (GE_BLOCK)

# undef IMAT
# undef IVEC
}


//------------------------------------------------------------------------------------------------------------
//! \brief  Solves a set of stacked linear systems of the form \f$ QHQ^Tx = b \f$ where \f$ H \f$ is upper
//!         Hessenberg and \f$ Q \f$ is an orthogonal matrix composed of Householder reflections.
//!
//! Uses a reduced Gaussian elimination algorithm. Partial row pivoting is enabled if the macro #HESS_PIVOT is
//! defined during compilation.
//!
//! \param[in]      N           The dimension of the linear system.
//! \param[in]      H           The upper Hessenberg matrix.
//! \param[in,out]  B           The right hand size \f$ b \f$ of the linear system.
//!                             Upon return contains the vector \f$ x \f$ which solves the linear system.
//! \param[in]      tau         Pointer to vector containing the scalars which define each Householder
//!                             reflection.
//! \param[in,out]  piv         Pointer to array used to store pivot indices. Expected to be at least
//!                             <tt>\pp{N} * sizeof(int)</tt> bytes. If #HESS_PIVOT is not defined during
//!                             compilation then this pointer is not dereferenced.
//------------------------------------------------------------------------------------------------------------
template< int SIMD_length >
void SIMD_UpperHessenbergSolve_CM (

    const int N,
    double * const restrict H,
    double * const restrict B,
    const double * const restrict tau,
    int * const restrict piv
) {

# if ! defined (HESS_PIVOT)
    (void) piv;
# endif

//
// \brief  Macro for indexing into SIMD matrix arrays. Assumes column-major storage.
//
// \param[in]   i   Row index of element in matrix.
// \param[in]   j   Column index of element in matrix.
// \param[in]   l   Index for system in SIMD block. In [ 0, SIMD_length ).
//
# define IMAT(i,j,l) ( (l) + SIMD_length*( (i) + N*(j) ) )

//
// \brief  Macro for indexing into SIMD vector arrays.
//
// \param[in]  i   Index of element in vector.
// \param[in]  l   Index for system in SIMD block. IN [ 0, SIMD_length ).
//
# define IVEC(i,l) ( (l) + SIMD_length*(i) )

    // --- PᵀB ---------------------------------------------------------------------------------------- //

    for ( int i = 0; i < N-2; ++i ) {

        alignas( SIMD_length * sizeof(double) ) double scal [SIMD_length];

        # pragma omp simd aligned( B, scal : SIMD_length * SIZEOF_DOUBLE )
        for ( int l = 0; l < SIMD_length; ++l )
            scal[l] = B[ IVEC(i+1,l) ];

        for ( int j = i+2; j < N; ++j ) {
        # pragma omp simd aligned( B, H, scal : SIMD_length * SIZEOF_DOUBLE )
        for ( int l = 0; l < SIMD_length; ++l ) {

            scal[l] += B[ IVEC(j,l) ] * H[ IMAT(j,i,l) ];
        }}

        # pragma omp simd aligned( tau, scal : SIMD_length * SIZEOF_DOUBLE )
        for ( int l = 0; l < SIMD_length; ++l )
            scal[l] *= -tau[ IVEC(i,l) ];

        # pragma omp simd aligned( B, scal : SIMD_length * SIZEOF_DOUBLE )
        for ( int l = 0; l < SIMD_length; ++l )
            B[ IVEC(i+1,l) ] += scal[l];

        for ( int j = i+2; j < N; ++j ) {
        # pragma omp simd aligned( B, H, scal : SIMD_length * SIZEOF_DOUBLE )
        for ( int l = 0; l < SIMD_length; ++l ) {

            B[ IVEC(j,l) ] += scal[l] * H[ IMAT(j,i,l) ];
        }}
    }

    // --- H⁻¹PᵀB ------------------------------------------------------------------------------------- //

# if defined (HESS_BLOCK)

    if ( N <= MAX_HESS_BLOCK_DIM )

# endif // if defined (HESS_BLOCK)

    {
        // Zero sub-diagonal.
        for ( int i = 0; i < N-1; ++i ) {

            # pragma omp simd aligned( H : SIMD_length * SIZEOF_DOUBLE )
            for ( int l = 0; l < SIMD_length; ++l )
                H[ IMAT(i+1,i,l) ] /= - H[ IMAT(i,i,l) ];

            # pragma omp simd aligned( B, H : SIMD_length * SIZEOF_DOUBLE )
            for ( int l = 0; l < SIMD_length; ++l )
                B[ IVEC(i+1,l) ] += H[ IMAT(i+1,i,l) ] * B[ IVEC(i,l) ];

            for ( int j = i+1; j < N; ++j ) {
            # pragma omp simd aligned( H : SIMD_length * SIZEOF_DOUBLE )
            for ( int l = 0; l < SIMD_length; ++l ) {

                H[ IMAT(i+1,j,l) ] += H[ IMAT(i+1,i,l) ] * H[ IMAT(i,j,l) ];
            }}
        }

        // Backward substitution.
        for ( int i = N-1; i >= 0; --i ) {

            # pragma omp simd aligned( B, H : SIMD_length * SIZEOF_DOUBLE )
            for ( int l = 0; l < SIMD_length; ++l )
                B[ IVEC(i,l) ] /= H[ IMAT(i,i,l) ];

            for ( int j = 0; j < i; ++j ) {
            # pragma omp simd aligned( B, H : SIMD_length * SIZEOF_DOUBLE )
            for ( int l = 0; l < SIMD_length; ++l ) {

                B[ IVEC(j,l) ] -= H[ IMAT(j,i,l) ] * B[ IVEC(i,l) ];
            }}
        }
    }

# if defined (HESS_BLOCK)

    else {

        const int num_blocks = 1 + ((N - 1) / MAX_HESS_BLOCK_DIM);  // if N != 0.
        const int block_dim = N / num_blocks + ( N % num_blocks ? 1 : 0 );

        // --- Zero sub-diagonal. --------------------------------------------------------------------- /

        for ( int blk_i = 0; blk_i < num_blocks; ++blk_i ) {

            const int blk_i_min = blk_i * block_dim;
            const int blk_i_max = MIN( N-1, (blk_i + 1) * block_dim );

            // Zero subdiagonal and reduce diagonal block.
            for ( int i = blk_i_min; i < blk_i_max; ++i ) {

                # pragma omp simd aligned( H : SIMD_length * SIZEOF_DOUBLE )
                for ( int l = 0; l < SIMD_length; ++l )
                    H[ IMAT(i+1,i,l) ] /= - H[ IMAT(i,i,l) ];

                # pragma omp simd aligned( B, H : SIMD_length * SIZEOF_DOUBLE )
                for ( int l = 0; l < SIMD_length; ++l )
                    B[ IVEC(i+1,l) ] += H[ IMAT(i+1,i,l) ] * B[ IVEC(i,l) ];

                for ( int j = i+1; j < blk_i_max; ++j ) {
                # pragma omp simd aligned( H : SIMD_length * SIZEOF_DOUBLE )
                for ( int l = 0; l < SIMD_length; ++l ) {

                    H[ IMAT(i+1,j,l) ] += H[ IMAT(i+1,i,l) ] * H[ IMAT(i,j,l) ];
                }}
            }

            // Reduce remainder of block row.
            for ( int j = blk_i_max; j < N; ++j ) {
            for ( int i = blk_i_min; i < blk_i_max; ++i ) {
            # pragma omp simd aligned( H : SIMD_length * SIZEOF_DOUBLE )
            for ( int l = 0; l < SIMD_length; ++l ) {

                H[ IMAT(i+1,j,l) ] += H[ IMAT(i+1,i,l) ] * H[ IMAT(i,j,l) ];
            }}}
        }

        // --- Backward substitution. ----------------------------------------------------------------- //

        // For each block.
        for ( int blk_i = num_blocks -1; blk_i >= 0; --blk_i ) {

            const int blk_i_min = blk_i * block_dim;
            const int blk_i_max = MIN( N, (blk_i + 1) * block_dim );

            // Perform backward substitution on block.
            for ( int i = blk_i_max -1; i >= blk_i_min; --i ) {

                # pragma omp simd aligned( H, B : SIMD_length * SIZEOF_DOUBLE )
                for ( int l = 0; l < SIMD_length; ++l )
                    B[ IVEC(i,l) ] /= H[ IMAT(i,i,l) ];

                for ( int j = blk_i_min; j < i; ++j ) {
                # pragma omp simd aligned( H, B : SIMD_length * SIZEOF_DOUBLE )
                for ( int l = 0; l < SIMD_length; ++l ) {

                    B[ IVEC(j,l) ] -= H[ IMAT(j,i,l) ] * B[ IVEC(i,l) ];
                }}

            }

            // Reduce column block of A.
            for ( int blk_k = 0; blk_k < blk_i; ++blk_k ) {
            for ( int i = blk_i_min; i < blk_i_max; ++i ) {
            for ( int k = blk_k * block_dim; k < (blk_k + 1) * block_dim; ++k ) {

                # pragma omp simd aligned( H, B : SIMD_length * SIZEOF_DOUBLE )
                for ( int l = 0; l < SIMD_length; ++l )
                    B[ IVEC(k,l) ] -= H[ IMAT(k,i,l) ] * B[ IVEC(i,l) ];
            }}}
        }
    }

# endif // if defined (HESS_BLOCK)

    // --- PH⁻¹PᵀB ------------------------------------------------------------------------------------ //

    for ( int i = N-3; i >= 0; --i ) {

        alignas( SIMD_length * sizeof(double) ) double scal[SIMD_length];

        # pragma omp simd aligned( B, scal : SIMD_length * SIZEOF_DOUBLE )
        for ( int l = 0; l < SIMD_length; ++l )
            scal[l] = B[ IVEC(i+1,l) ];

        for ( int j = i+2; j < N; ++j ) {
        # pragma omp simd aligned( B, H, scal : SIMD_length * SIZEOF_DOUBLE )
        for ( int l = 0; l < SIMD_length; ++l ) {

            scal[l] += B[ IVEC(j,l) ] * H[ IMAT(j,i,l) ];
        }}

        # pragma omp simd aligned( tau, scal : SIMD_length * SIZEOF_DOUBLE )
        for ( int l = 0; l < SIMD_length; ++l )
            scal[l] *= -tau[ IVEC(i,l) ];

        # pragma omp simd aligned( B, scal : SIMD_length * SIZEOF_DOUBLE )
        for ( int l = 0; l < SIMD_length; ++l )
            B[ IVEC(i+1,l) ] += scal[l];

        for ( int j = i+2; j < N; ++j ) {
        # pragma omp simd aligned( B, H, scal : SIMD_length * SIZEOF_DOUBLE )
        for ( int l = 0; l < SIMD_length; ++l ) {

            B[ IVEC(j,l) ] += scal[l] * H[ IMAT(j,i,l) ];
        }}
    }
}
