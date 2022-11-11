//------------------------------------------------------------------------------------------------------------
//! \file   linear_solvers/LinearAlgebra.c
//! \brief  Contains implementations of linear algebra routines used by transport sweeps.
//!
//! \author Michael M. Crockatt
//! \date   December 2017
//------------------------------------------------------------------------------------------------------------

# include <math.h>

# include "LinearAlgebra.h"


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


//============================================================================================================
//=== SINGLE SYSTEM ROUTINES =================================================================================
//============================================================================================================


//------------------------------------------------------------------------------------------------------------
//! \brief  Solves the given linear system \f$ Ax = b \f$ using Gaussian elimination without pivoting.
//!
//! \attention  The matrix \f$ A \f$ is assumed to be stored in _column-major_ ordering.
//!
//! \param[in]      N           Dimension of the system to be solved.
//! \param[in]      A           Pointer to the matrix.
//! \param[in,out]  b           Initially contains the known vector \f$ b \f$ for the linear system.
//!                             Upon return contains the solution vector \f$ x \f$ of the linear system.
//!
//! Total flop count: \f$ \frac{4}{6} N^3 + \frac{9}{6} N^2 - \frac{7}{6} N \f$.
//!
//! Flop composition without \c FMA instructions:
//!
//! | Operation | Count |
//! | :-------- | ----: |
//! | \c DIV    | \f$ \frac{1}{2} N^2 + \frac{1}{2} N \f$ |
//! | \c MUL    | \f$ \frac{1}{3} N^3 + \frac{1}{2} N^2 - \frac{5}{6} N \f$ |
//! | \c ADD    | \f$ \frac{1}{3} N^3 + \frac{1}{2} N^2 - \frac{5}{6} N \f$ |
//!
//! Flop composition with \c FMA instructions:
//!
//! | Operation | Count |
//! | :-------- | ----: |
//! | \c DIV    | \f$ \frac{1}{2} N^2 + \frac{1}{2} N \f$ |
//! | \c FMA    | \f$ \frac{1}{3} N^3 + \frac{1}{2} N^2 - \frac{5}{6} N \f$ |
//!
//------------------------------------------------------------------------------------------------------------
void GaussianElimination_CM (

    const int N,
    double * const restrict A,
    double * const restrict b
) {

//
// \brief  Macro for indexing into matrix arrays. Assumes column-major storage.
//
// \param[in]   i   Row index of element in matrix.
// \param[in]   j   Column index of element in matrix.
//
# define IMAT(i,j) ( (i) + N*(j) )

    // Forward elimination.
    for ( int j = 0; j < N; ++j ) {

        // Compute L_{i,j} and apply reductions to RHS vector.
        for ( int i = j+1; i < N; ++i ) {
        /*
         *  Total executions of loop body: (1/2) N^2 - (1/2) N
         *  with:
         *          1 DIV   |   1 DIV
         *          1 MUL   |   1 FMA
         *          1 ADD   |
         */
            A[ IMAT(i,j) ] /= - A[ IMAT(j,j) ]; // 1 DIV.

            b[i] += A[ IMAT(i,j) ] * b[j];      // 1 MUL, 1 ADD = 1 FMA.
        }

        // Reduce rows.
        for ( int k = j+1; k < N; ++k ) {
        for ( int i = j+1; i < N; ++i ) {
        /*
         *  Total executions of loop body: (1/3) N^3 - (1/2) N^2 + (1/6) N
         *  with:
         *          1 MUL   |   1 FMA
         *          1 ADD   |
         */
            A[ IMAT(i,k) ] += A[ IMAT(i,j) ] * A[ IMAT(j,k) ];  // 1 MUL, 1 ADD = 1 FMA.
        }}
    }

    // Backward substitution.
    for ( int i = N-1; i >= 0; --i ) {
    /*
     *  Total executions of loop body: N
     *  with:
     *          1 DIV
     */
        b[i] /= A[ IMAT(i,i) ]; // 1 DIV.

        for ( int j = 0; j < i; ++j ) {
        /*
         *  Total executions of loop body: (1/2) N^2 - (1/2) N
         *  with:
         *          1 MUL   |   1 FMA
         *          1 ADD   |
         */
            b[j] -= A[ IMAT(j,i) ] * b[i];  // 1 MUL, 1 ADD = 1 FMA.
        }
    }

# undef IMAT
}


//------------------------------------------------------------------------------------------------------------
//! \brief  Solves a linear system of the form \f$ Hx = b \f$ where \f$ H \f$ is upper Hessenberg.
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
//!
//! Total flop count: \f$ 6 N^2 - 6 N - 3 \f$.
//!
//! Flop composition without \c FMA instructions:
//!
//! | Operation | Count |
//! | :-------- | ----: |
//! | \c DIV    | \f$ 2N - 1 \f$ |
//! | \c MUL    | \f$ 3 N^2 - 4 N - 1 \f$ |
//! | \c ADD    | \f$ 3 N^2 - 4 N - 1 \f$ |
//!
//! Flop composition with \c FMA instructions:
//!
//! | Operation | Count |
//! | :-------- | ----: |
//! | \c DIV    | \f$ 2N - 1 \f$ |
//! | \c FMA    | \f$ 3 N^2 - 6 N + 3 \f$ |
//! | \c MUL    | \f$ 2 N - 4 \f$ |
//! | \c ADD    | \f$ 2 N - 4 \f$ |
//!
//------------------------------------------------------------------------------------------------------------
void UpperHessenbergSolve_CM (

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
// \brief  Macro for indexing into matrix arrays. Assumes column-major storage.
//
// \param[in]   i   Row index of element in matrix.
// \param[in]   j   Column index of element in matrix.
//
# define IMAT(i,j) ( (i) + N*(j) )

    // --- PᵀB. --------------------------------------------------------------------------------------- //
    /*
     *  Total flop count: 2 N^2 - 4 N
     *
     *  Flop composition without FMA instructions:
     *
     *      MUL     N^2 - 2 N
     *      ADD     N^2 - 2 N
     *
     *  Flop composition with FMA instructions:
     *
     *      FMA     N^2 - 3 N + 2
     *      MUL     N - 2
     *      ADD     N - 2
     */
    for ( int i = 0; i < N-2; ++i ) {
    /*
     *  Total executions of loop body: N - 2
     *  with:
     *          1 MUL
     *          1 ADD
     */
        double scal = B[i+1];

        for ( int j = i+2; j < N; ++j ) {
        /*
         *  Total executions of loop body: (1/2) N^2 - (3/2) N + 1
         *  with:
         *          1 MUL   |   1 FMA
         *          1 ADD   |
         */
            scal += B[j] * H[ IMAT(j,i) ];  // 1 MUL, 1 ADD = 1 FMA.
        }

        scal *= -tau[i];    // 1 MUL.
        B[i+1] += scal;     // 1 ADD.

        for ( int j = i+2; j < N; ++j ) {
        /*
         *  Total executions of loop body: (1/2) N^2 - (3/2) N + 1
         *  with:
         *          1 MUL   |   1 FMA
         *          1 ADD   |
         */
            B[j] += scal * H[ IMAT(j,i) ];  // 1 MUL, 1 ADD = 1 FMA.
        }
    }

    // --- H⁻¹PᵀB. ------------------------------------------------------------------------------------ //
    /*
     *  Total flop count: 2 N^2 + 2 N - 3
     *
     *  Flop composition without FMA instructions:
     *
     *      DIV     2 N - 1
     *      MUL     N^2 - 1
     *      ADD     N^2 - 1
     *
     *  Flop composition with FMA instructions:
     *
     *      DIV     2 N - 1
     *      FMA     N^2 - 1
     */
    // Zero sub-diagonal.
    for ( int i = 0; i < N-1; ++i ) {
    /*
     *  Total executions of loop body: N - 1
     *  with:
     *          1 DIV   |   1 DIV
     *          1 MUL   |   1 FMA
     *          1 ADD   |
     */
        H[ IMAT(i+1,i) ] /= - H[ IMAT(i,i) ];   // 1 DIV.

        B[i+1] += H[ IMAT(i+1,i) ] * B[i];      // 1 MUL, 1 ADD = 1 FMA.

        for ( int j = i+1; j < N; ++j ) {
        /*
         *  Total executions of loop body: (1/2) N^2 - (1/2) N
         *  with:
         *          1 MUL   |   1 FMA
         *          1 ADD   |
         */
            H[ IMAT(i+1,j) ] += H[ IMAT(i+1,i) ] * H[ IMAT(i,j) ];  // 1 MUL, 1 ADD = 1 FMA.
        }
    }

    // Backward substitution.
    for ( int i = N-1; i >= 0; --i ) {
    /*
     *  Total executions of loop body: N
     *  with:
     *          1 DIV
     */
        B[i] /= H[ IMAT(i,i) ]; // 1 DIV.

        for ( int j = 0; j < i; ++j ) {
        /*
         *  Total executions of loop body: (1/2) N^2 - (1/2) N
         *  with:
         *          1 MUL   |   1 FMA
         *          1 ADD
         */
            B[j] -= H[ IMAT(j,i) ] * B[i];  // 1 MUL, 1 ADD = 1 FMA.
        }
    }

    // --- PH⁻¹PᵀB. ----------------------------------------------------------------------------------- //
    /*
     *  Total flop count: 2 N^2 - 4 N
     *
     *  Flop composition without FMA instructions:
     *
     *      MUL     N^2 - 2 N
     *      ADD     N^2 - 2 N
     *
     *  Flop composition with FMA instructions:
     *
     *      FMA     N^2 - 3 N + 2
     *      MUL     N - 2
     *      ADD     N - 2
     */
    for ( int i = N-3; i >= 0; --i ) {
    /*
     *  Total executions of loop body: N - 2
     *  with:
     *          1 MUL
     *          1 ADD
     */
        double scal = B[i+1];

        for ( int j = i+2; j < N; ++j ) {
        /*
         *  Total executions of loop body: (1/2) N^2 - (3/2) N + 1
         *  with:
         *          1 MUL   |   1 FMA
         *          1 ADD   |
         */
            scal += B[j] * H[ IMAT(j,i) ];  // 1 MUL, 1 ADD = 1 FMA.
        }

        scal *= -tau[i];    // 1 MUL.
        B[i+1] += scal;     // 1 ADD.

        for ( int j = i+2; j < N; ++j ) {
        /*
         *  Total executions of loop body: (1/2) N^2 - (3/2) N + 1
         *  with:
         *          1 MUL   |   1 FMA
         *          1 ADD   |
         */
            B[j] += scal * H[ IMAT(j,i) ];  // 1 MUL, 1 ADD = 1 FMA.
        }
    }

# undef IMAT
}
