//------------------------------------------------------------------------------------------------------------
//! \file   operators/RelabelOperator.cpp
//! \brief  Implementation of RelabelOperator class.
//!
//! \author Michael M. Crockatt
//! \date   December 2017
//------------------------------------------------------------------------------------------------------------


# include <cstring>
# include <cstdlib>
# include <limits>
# include <stdexcept>

# include "linear_solvers/LinearAlgebra.h"
# include "operators/RelabelOperator.hpp"
# include "utils/CLog.hpp"
# include "utils/global.hpp"


using namespace Quadrule;


//============================================================================================================
//=== STATIC DEFINITIONS =====================================================================================
//============================================================================================================


//------------------------------------------------------------------------------------------------------------
//! \brief  Used to convert RelabelOperator::RelabelType values to their corresponding string descriptors.
//------------------------------------------------------------------------------------------------------------
std::map< RelabelOperator::RelabelType, std::string > RelabelOperator::RelabelType_to_String = {

    { RelabelOperator::RelabelType::None,               "none"              }

# if SPACE_DIMS == 1

,   { RelabelOperator::RelabelType::PolyInterp,         "poly-interp"       }
,   { RelabelOperator::RelabelType::DoublePolyInterp,   "dbl-poly-interp"   }

# elif SPACE_DIMS >= 2

,   { RelabelOperator::RelabelType::SphHyperinterp,     "sph-hyperinterp"   }
,   { RelabelOperator::RelabelType::SphLeastSquares,    "sph-ls"            }
,   { RelabelOperator::RelabelType::PWC_Hierarchical,   "pwc-hierarchical"  }

# endif
};


//------------------------------------------------------------------------------------------------------------
//! \brief  Used to convert string descriptors to RelabelOperator::RelabelType values.
//------------------------------------------------------------------------------------------------------------
std::map< std::string, RelabelOperator::RelabelType > RelabelOperator::String_to_RelabelType = {

    { "none",               RelabelOperator::RelabelType::None,             }

# if SPACE_DIMS == 1

,   { "poly-interp",        RelabelOperator::RelabelType::PolyInterp        }
,   { "dbl-poly-interp",    RelabelOperator::RelabelType::DoublePolyInterp  }

# elif SPACE_DIMS >= 2

,   { "sph-hyperinterp",    RelabelOperator::RelabelType::SphHyperinterp    }
,   { "sph-ls",             RelabelOperator::RelabelType::SphLeastSquares   }
,   { "pwc-hierarchical",   RelabelOperator::RelabelType::PWC_Hierarchical  }

# endif
};


//------------------------------------------------------------------------------------------------------------
//! \brief  Used to convert RelabelOperator::RelabelBLASOp values to their corresponding string descriptors.
//------------------------------------------------------------------------------------------------------------
std::map< RelabelOperator::RelabelBLASOp, std::string > RelabelOperator::RelabelBLASOp_to_String = {

    { RelabelOperator::RelabelBLASOp::None,             "none"          }
,   { RelabelOperator::RelabelBLASOp::OneMat_GEMV,      "one-mat-gemv"  }
,   { RelabelOperator::RelabelBLASOp::OneMat_GEMM,      "one-mat-gemm"  }

# if SPACE_DIMS >= 2

,   { RelabelOperator::RelabelBLASOp::TwoMat_GEMV,      "two-mat-gemv"  }
,   { RelabelOperator::RelabelBLASOp::TwoMat_GEMM,      "two-mat-gemm"  }

# endif
};


//------------------------------------------------------------------------------------------------------------
//! \brief  Used to convert string descriptors to RelabelOperator::RelabelBLASOp values.
//------------------------------------------------------------------------------------------------------------
std::map< std::string, RelabelOperator::RelabelBLASOp > RelabelOperator::String_to_RelabelBLASOp = {

    { "none",           RelabelOperator::RelabelBLASOp::None,       }
,   { "one-mat-gemv",   RelabelOperator::RelabelBLASOp::OneMat_GEMV }
,   { "one-mat-gemm",   RelabelOperator::RelabelBLASOp::OneMat_GEMM }

# if SPACE_DIMS >= 2

,   { "two-mat-gemv",   RelabelOperator::RelabelBLASOp::TwoMat_GEMV }
,   { "two-mat-gemm",   RelabelOperator::RelabelBLASOp::TwoMat_GEMM }

# endif
};


//============================================================================================================
//=== CONSTRUCTORS AND DESTRUCTOR ============================================================================
//============================================================================================================


# if SPACE_DIMS == 1

//------------------------------------------------------------------------------------------------------------
//! \brief  Construct a RelabelOperator object using the given parameters.
//!
//! \param[in]  src_ordinate_type       Type of discrete ordinate quadrature to map from.
//! \param[in]  src_ang_order           Order of discrete ordinate quadrature to map from.
//! \param[in]  src_symmetric_reduce    Specifies whether the number of ordinates should be reduced by
//!                                     symmetry in the discrete ordinate quadrature to map from.
//! \param[in]  dst_ordinate_type       Type of discrete ordinate quadrature to map onto.
//! \param[in]  dst_ang_order           Order of discrete ordinate quadrature to map onto.
//! \param[in]  dst_symmetric_reduce    Specifies whether the number of ordinates should be reduced by
//!                                     symmetry in the discrete ordinate quadrature to map to.
//------------------------------------------------------------------------------------------------------------
RelabelOperator::RelabelOperator (

    const OrdinateSet::OrdinateType src_ordinate_type,
    const int64_t src_ang_order,
    const bool src_symmetric_reduce,
    const OrdinateSet::OrdinateType dst_ordinate_type,
    const int64_t dst_ang_order,
    const bool dst_symmetric_reduce
) :
    relabel_type {},
    src( src_ang_order, src_symmetric_reduce, src_ordinate_type ), src_basis_ptr{ nullptr },
    dst( dst_ang_order, dst_symmetric_reduce, dst_ordinate_type ), dst_basis_ptr{ nullptr },
    stride(0),
    R_mat_ptr{ nullptr }
{

    Global::input_list.GetValue( "relabel_type", this->relabel_type, String_to_RelabelType );

    // Determine RelabelBLASOp to use when applying relabel operator. Use default value if none found.
    try         {  Global::input_list.GetValue( "relabel_blas_op", blas_op, String_to_RelabelBLASOp );  }
    catch (...) {  blas_op = RelabelBLASOp::OneMat_GEMM;                                                }

    PRINT_LOG( "  %-20s  %s\n", "Relabel type:", RelabelType_to_String.at( relabel_type ).c_str() )
    PRINT_LOG( "  %-20s  %s\n", "Relabel BLAS Op:", RelabelBLASOp_to_String.at( blas_op ).c_str() )

    // --- Allocate memory and compute matrix. -------------------------------------------------------- //

    R_mat_ptr = (double *)
        # if defined (USE_ALIGNED_ALLOC)
            aligned_alloc( ALIGNED_ALLOC_ALIGNMENT, src.nq() * dst.nq() * sizeof(double) );
        # else
            std::malloc( src.nq() * dst.nq() * sizeof(double) );
        # endif

    std::memset( R_mat_ptr, 0, src.nq() * dst.nq() * sizeof(double) );

    switch ( this->relabel_type ) {

        //
        // Construct relabeling operator using Lagrange interpolating polynomial over [-1,1].
        //
        case RelabelType::PolyInterp:
        {
            for ( int64_t src_q = 0; src_q < src.nq(); ++src_q ) {
            for ( int64_t dst_q = 0; dst_q < dst.nq(); ++dst_q ) {

                RMat(dst_q,src_q) = 1.0;

                for ( int64_t k = 0; k < src.nq(); ++k ) {

                    if ( k == src_q ) {  continue;  }

                    RMat(dst_q,src_q) *= ( dst.xi(dst_q) - src.xi(k) ) / ( src.xi(src_q) - src.xi(k) );
                }
            }}
        } break;

        //
        // Construct relabeling operator using Lagrange interpolating polynomials over [-1,0] and [0,1]
        // separately. Assumes that both the source and destination quadratures contain the same number of
        // points in the negative and positive subintervals.
        //
        case RelabelType::DoublePolyInterp:
        {
            // Negative angles.
            for ( int64_t src_q = 0; src_q < src.nq() /2; ++src_q ) {
            for ( int64_t dst_q = 0; dst_q < dst.nq() /2; ++dst_q ) {

                RMat(dst_q,src_q) = 1.0;

                for ( int64_t k = 0; k < src.nq() /2; ++k ) {

                    if ( k == src_q ) {  continue;  }

                    RMat(dst_q,src_q) *= ( dst.xi(dst_q) - src.xi(k) ) / ( src.xi(src_q) - src.xi(k) );
                }
            }}

            // Positive angles.
            for ( int64_t src_q = src.nq() /2; src_q < src.nq(); ++src_q ) {
            for ( int64_t dst_q = dst.nq() /2; dst_q < dst.nq(); ++dst_q ) {

                RMat(dst_q,src_q) = 1.0;

                for ( int64_t k = src.nq() /2; k < src.nq(); ++k ) {

                    if ( k == src_q ) {  continue;  }

                    RMat(dst_q,src_q) *= ( dst.xi(dst_q) - src.xi(k) ) / ( src.xi(src_q) - src.xi(k) );
                }
            }}
        } break;

        default:
        {   std::string error_message = "Invalid OrdinateType '"
                                        + OrdinateSet::OrdinateType_to_String.at( src.GetOrdinateType() )
                                        + "' for source OrdinateSet in '" + __func__ + "'.\n";

            PRINT_ERROR( error_message.c_str() )
            throw std::invalid_argument( error_message );
        }
    }
}

# endif // if SPACE_DIMS == 1


# if SPACE_DIMS >= 2 || defined (DOXYCOMPILE)

//------------------------------------------------------------------------------------------------------------
//! \brief  Construct a RelabelOperator object using the given parameters.
//!
//! \param[in]  src_ordinate_type       Type of discrete ordinate quadrature to map from.
//! \param[in]  src_ang_order           Order of discrete ordinate quadrature to map from.
//! \param[in]  src_symmetric_reduce    Specifies whether the number of ordinates should be reduced by
//!                                     symmetry in the discrete ordinate quadrature to map from.
//! \param[in]  dst_ordinate_type       Type of discrete ordinate quadrature to map onto.
//! \param[in]  dst_ang_order           Order of discrete ordinate quadrature to map onto.
//! \param[in]  dst_symmetric_reduce    Specifies whether the number of ordinates should be reduced by
//!                                     symmetry in the discrete ordinate quadrature to map to.
//------------------------------------------------------------------------------------------------------------
RelabelOperator::RelabelOperator (

    const OrdinateSet::OrdinateType src_ordinate_type,
    const int64_t src_ang_order,
    const bool src_symmetric_reduce,
    const OrdinateSet::OrdinateType dst_ordinate_type,
    const int64_t dst_ang_order,
    const bool dst_symmetric_reduce
) :
    relabel_type {},
    src( src_ang_order, src_symmetric_reduce, src_ordinate_type ), src_basis_ptr{ nullptr },
    dst( dst_ang_order, dst_symmetric_reduce, dst_ordinate_type ), dst_basis_ptr{ nullptr },
    interpolant_degree {},
    stride {},
    R_mat_ptr{ nullptr }
{

    Global::input_list.GetValue( "relabel_type", relabel_type, String_to_RelabelType );

    // Determine RelabelBLASOp to use when applying relabel operator. Use default value if none found.
    try         {  Global::input_list.GetValue( "relabel_blas_op", blas_op, String_to_RelabelBLASOp );  }
    catch (...) {  blas_op = RelabelBLASOp::TwoMat_GEMM;                                                }

    switch ( relabel_type ) {

        //
        // Relabel from collided quadrature via hyperinterpolation by spherical harmonics.
        //
        case RelabelType::SphHyperinterp:
        {
            // Determine degree of spherical harmonics basis.
            switch ( src.GetOrdinateType() ) {

                case OrdinateSet::OrdinateType::ChebyshevLegendre:

                    interpolant_degree = src.GetAngOrder() - 1;
                    break;

            # if defined (ENABLE_LEBEDEV)

                case OrdinateSet::OrdinateType::Lebedev:

                    interpolant_degree = src.GetAngOrder() / 2;
                    break;

            # endif // if defined (ENABLE_LEBEDEV)

                default:
                {   std::string error_message =
                          "Invalid OrdinateType '"
                        + OrdinateSet::OrdinateType_to_String.at( src.GetOrdinateType() )
                        + "' for constructing relabeling operator with RelabelType '"
                        + RelabelType_to_String.at( relabel_type )
                        + "'.\n";

                    PRINT_ERROR( error_message.c_str() )
                    throw std::invalid_argument( error_message );
                }
            }

            stride = (interpolant_degree + 1)*(interpolant_degree + 1);

            // Allocate arrays for storing matrices.
            src_basis_ptr = (double *)
                # if defined (USE_ALIGNED_ALLOC)
                    aligned_alloc( ALIGNED_ALLOC_ALIGNMENT, src.nq() * stride * sizeof(double) );
                # else
                    std::malloc( src.nq() * stride * sizeof(double) );
                # endif

            dst_basis_ptr = (double *)
                # if defined (USE_ALIGNED_ALLOC)
                    aligned_alloc( ALIGNED_ALLOC_ALIGNMENT, dst.nq() * stride * sizeof(double) );
                # else
                    std::malloc( dst.nq() * stride * sizeof(double) );
                # endif

            //
            // Compute maps used for relabeling.
            //
            // SOURCE MAP:
            //
            //      The first map is the discrete-to-moment operator that maps the values at the source
            //      quadrature nodes to a set of moments with respect to a basis of unit-normalized spherical
            //      harmonics.
            //

            for ( int64_t q = 0; q < src.nq(); ++q ) {

                for ( int64_t n = 0; n <= interpolant_degree; ++n ) {
                for ( int64_t l = 1; l <  2*(n+1);            ++l ) {

                    const int64_t k = n*n + l - 1;

                # if SPACE_DIMS == 2

                    if ( src.GetOrdinateSymmetry() ) {

                        SrcMat(k,q) = src.w(q) /2
                                      * ( SphericalHarmonic( n, l, src.mu(q), src.OrdinatePhi(q) )
                                          + SphericalHarmonic( n, l, -src.mu(q), src.OrdinatePhi(q) ) );
                    } else
                # endif // if SPACE_DIMS == 2
                    {
                        SrcMat(k,q) = src.w(q) * SphericalHarmonic( n, l, src.mu(q), src.OrdinatePhi(q) );
                    }
                }}
            }

            //
            // DESTINATION MAP:
            //
            //      The second map is a moment-to-discrete operator that maps the moments with respect to the
            //      basis of unit-normalized spherical harmonics to values at the destination quadrature nodes
            //      by sampling the expansion of spherical harmonics at each quadrature node.
            //

            for ( int64_t q = 0; q < dst.nq(); ++q ) {

                for ( int64_t n = 0; n <= interpolant_degree; ++n ) {
                for ( int64_t l = 1; l <  2*(n+1);            ++l ) {

                    const int64_t k = n*n + l - 1;

                # if SPACE_DIMS == 2

                    if ( dst.GetOrdinateSymmetry() ) {

                        DstMat(k,q) = 0.5 * ( SphericalHarmonic( n, l, dst.mu(q), dst.OrdinatePhi(q) )
                                              + SphericalHarmonic( n, l, -dst.mu(q), dst.OrdinatePhi(q) ) );
                    } else
                # endif // if SPACE_DIMS == 2
                    {
                        DstMat(k,q) = SphericalHarmonic( n, l, dst.mu(q), dst.OrdinatePhi(q) );
                    }
                }}
            }

            // Compute multiplied-out matrix R.
            R_mat_ptr = (double *)
                # if defined (USE_ALIGNED_ALLOC)
                    aligned_alloc( ALIGNED_ALLOC_ALIGNMENT, src.nq() * dst.nq() * sizeof(double) );
                # else
                    std::malloc( src.nq() * dst.nq() * sizeof(double) );
                # endif

            std::memset( R_mat_ptr, 0, src.nq() * dst.nq() * sizeof(double) );

            for ( int64_t dst_q = 0; dst_q < dst.nq(); ++dst_q ) {
            for ( int64_t src_q = 0; src_q < src.nq(); ++src_q ) {
            for ( int64_t k     = 0; k     < stride;   ++k     ) {

                RMat(dst_q,src_q) += SrcMat(k,src_q) * DstMat(k,dst_q);
            }}}

            break;
        }

        //
        // Relabel from collided quadrature via least-squares reconstruction onto a basis of spherical
        // harmonics.
        //
        // The degree of the spherical harmonic basis used is chosen such that the total number of spherical
        // harmonics used is less than or equal to the number of quadrature nodes. This yields an
        // over-determined least-squares problem.
        //
        case RelabelType::SphLeastSquares:
        {
            // Determine degree of spherical harmonic basis.
            interpolant_degree = ((int64_t)(std::sqrt( src.nq() ))) - 1;

            stride = (interpolant_degree + 1)*(interpolant_degree + 1);

            // Allocate arrays for storing matrices.
            src_basis_ptr = (double *)
                # if defined (USE_ALIGNED_ALLOC)
                    aligned_alloc( ALIGNED_ALLOC_ALIGNMENT, src.nq() * stride * sizeof(double) );
                # else
                    std::malloc( src.nq() * stride * sizeof(double) );
                # endif

            dst_basis_ptr = (double *)
                # if defined (USE_ALIGNED_ALLOC)
                    aligned_alloc( ALIGNED_ALLOC_ALIGNMENT, dst.nq() * stride * sizeof(double) );
                # else
                    std::malloc( dst.nq() * stride * sizeof(double) );
                # endif

            //
            // Compute maps used for relabeling.
            //
            // SOURCE MAP:
            //
            //      The first map is a discrete-to-moment operator that maps the values at the source
            //      quadrature nodes to a set of moments with respect to a basis of unit-normalized spherical
            //      harmonics. This operator is constructed by taking the Moore-Penrose pseudoinverse of the
            //      matrix that defines the associated linear least-squares "interpolation" problem.
            //

            // Temporary variables used for LAPACK calls.
            int M = src.nq();
            int N = stride;

            double * A = new double[ M*N ];
            double * U = new double[ M*M ];
            double * S = new double[ std::min(M,N) ];
            double * V = new double[ N*N ];

            int LWORK = -1;
            double * WORK = new double[1];

            int * IWORK = new int[ 8 * std::min(M,N) ];

            // Construct least-squares matrix.
            for ( int64_t q = 0; q < src.nq(); ++q ) {

                for ( int64_t n = 0; n <= interpolant_degree; ++n ) {
                for ( int64_t l = 1; l <  2*(n+1);            ++l ) {

                    const int64_t k = n*n + l - 1;

                    A[ q + M*k ] = SphericalHarmonic( n, l, src.mu(q), src.OrdinatePhi(q) );
                }}
            }

            //
            // Compute norm of matrix and determine tolerance for computing numerical rank.
            //
            char NORM = 'F';    // Compute Frobenius norm of matrix.

            double F_norm = dlange_( &NORM,     // NORM:    Return Frobenius norm.
                                     &M,        // M:       Number of rows of the input matrix A.
                                     &N,        // N:       Number of columns of the input matrix A.
                                     A,         // A:       Pointer to the matrix A.
                                     &M,        // LDA:     Leading dimension of the matrix A.
                                     WORK );    // WORK:    Not referenced for Frobenius norm computation.

//             const double TOL = std::numeric_limits<double>::epsilon() * F_norm;
            const double TOL = 1.0e-4 * F_norm;

            PRINT_STATUS( "TOL = %.4e\n\n", TOL )

            //
            // Compute SVD of least-squares matrix. First call determines space required for WORK, second
            // call performs SVD computation.
            //
            char JOBZ = 'A';    // Return all M columns of U and all N rows of V^T.
            int INFO = 0;

            dgesdd_( &JOBZ,     // JOBZ:    Return all M columns of U and all N rows of V^T.
                     &M,        // M:       Number of rows of the input matrix A.
                     &N,        // N:       Number of columns of the input matrix A.
                     A,         // A:       Pointer to the matrix A.
                     &M,        // LDA:     Leading dimension of the matrix A.
                     S,         // S:       Pointer to the matrix Σ.
                     U,         // U:       Pointer to the matrix U.
                     &M,        // LDU:     Leading dimension of the matrix U.
                     V,         // VT:      Pointer to the matrix V^T.
                     &N,        // LDVT:    Leading dimension of the matrix V^T.
                     WORK,      //
                     &LWORK,    // LWORK:   Dimension of the array WORK.
                     IWORK,     //
                     &INFO );   //

            LWORK = WORK[0];
            delete [] WORK;
            WORK = new double[ LWORK ];

            dgesdd_( &JOBZ,     // JOBZ:    Return all M columns of U and all N rows of V^T.
                     &M,        // M:       Number of rows of the input matrix A.
                     &N,        // N:       Number of columns of the input matrix A.
                     A,         // A:       Pointer to the matrix A.
                     &M,        // LDA:     Leading dimension of the matrix A.
                     S,         // S:       Pointer to the matrix Σ.
                     U,         // U:       Pointer to the matrix U.
                     &M,        // LDU:     Leading dimension of the matrix U.
                     V,         // VT:      Pointer to the matrix V^T.
                     &N,        // LDVT:    Leading dimension of the matrix V^T.
                     WORK,      //
                     &LWORK,    // LWORK:   Dimension of the array WORK.
                     IWORK,     //
                     &INFO );   //

            //
            // Compute pseudoinverse using SVD.
            //

            // First invert matrix of singular values, ignoring singular values below tolerance.
            {
                int64_t i = 0;

                for ( ; i < std::min(M,N) && S[i] > TOL; ++i ) {

                    PRINT_STATUS( "S[%4" PRId64 "] = %.4e\n", i, S[i] )
                    S[i] = 1.0 / S[i];
                }

                PRINT_STATUS( "\n" )

                for ( ; i < std::min(M,N); ++i ) {

                    PRINT_STATUS( "S[%4" PRId64 "] = %.4e\n", i, S[i] )
                    S[i] = 0.0;
                }
            }

            // Then apply inverted singular values to U.
            for ( int64_t i = 0; i < std::min(M,N); ++i ) {
            for ( int64_t j = 0; j < M;             ++j ) {

                U[ j + M*i ] *= S[i];
            }}

            // Compute pinv(A) = (V^T)^T * (UΣ^{-T})^T.
            char TRANST = 'T';  // Specifies op(*) = A**T.
            double ZERO = 0.0;
            double ONE  = 1.0;

            dgemm_( &TRANST,        // TRANSA:  op(A) = A**T.
                    &TRANST,        // TRANSB:  op(B) = B**T.
                    &N,             // M:       Number of rows of the matrix op(A).
                    &M,             // N:       Number of columns of the matrix op(B).
                    &N,             // K:       Number of columns in op(A) and number of rows in op(B).
                    &ONE,           // ALPHA:   The scalar alpha.
                    V,              // A:       Pointer to the matrix A.
                    &N,             // LDA:     Leading dimension of the matrix A.
                    U,              // B:       Pointer to the matrix B.
                    &M,             // LDB:     Leading dimension of the matrix B.
                    &ZERO,          // BETA:    The scalar beta.
                    src_basis_ptr,  // C:       Pointer to the matrix C.
                    &N );           // LDC:     Leading dimension of the matrix C.

            // Cleanup.
            delete [] A;
            delete [] U;
            delete [] S;
            delete [] V;
            delete [] WORK;
            delete [] IWORK;

            //
            // DESTINATION MAP:
            //
            //      The second map is a moment-to-discrete operator that maps the moments with respect to the
            //      basis of unit-normalized spherical harmonics to values at the destination quadrature nodes
            //      by sampling the expansion of spherical harmonics at each quadrature node.
            //

            for ( int64_t q = 0; q < dst.nq(); ++q ) {

                for ( int64_t n = 0; n <= interpolant_degree; ++n ) {
                for ( int64_t l = 1; l <  2*(n+1);            ++l ) {

                    const int64_t k = n*n + l - 1;

                    DstMat(k,q) = SphericalHarmonic( n, l, dst.mu(q), dst.OrdinatePhi(q) );
                }}
            }

            // Compute multiplied-out matrix R.
            R_mat_ptr = (double *)
                # if defined (USE_ALIGNED_ALLOC)
                    aligned_alloc( ALIGNED_ALLOC_ALIGNMENT, src.nq() * dst.nq() * sizeof(double) );
                # else
                    std::malloc( src.nq() * dst.nq() * sizeof(double) );
                # endif

            std::memset( R_mat_ptr, 0, src.nq() * dst.nq() * sizeof(double) );

            for ( int64_t dst_q = 0; dst_q < dst.nq(); ++dst_q ) {
            for ( int64_t src_q = 0; src_q < src.nq(); ++src_q ) {
            for ( int64_t k     = 0; k     < stride;   ++k     ) {

                RMat(dst_q,src_q) += SrcMat(k,src_q) * DstMat(k,dst_q);
            }}}

            break;
        }

        //
        // Relabel between spherical-triangle quadrature sets using hierarchical decomposition.
        //
        case RelabelType::PWC_Hierarchical:
        {
            // Check that angular quadratures are compatible with given RelabelType.
            if (    src.GetOrdinateType() != OrdinateSet::OrdinateType::SphericalTriangle
                 || dst.GetOrdinateType() != OrdinateSet::OrdinateType::SphericalTriangle
                 || !IsPowerOfTwo( src.GetAngOrder() )
                 || !IsPowerOfTwo( dst.GetAngOrder() )
                 || dst.GetAngOrder() < src.GetAngOrder()
            ) {
                std::string error_message =
                      "Invalid RelabelType '"
                    + RelabelType_to_String.at( relabel_type )
                    + "' for constructing relabeling operator to map from OrdinateSet '"
                    + OrdinateSet::OrdinateType_to_String.at( src.GetOrdinateType() )
                    + ":" + std::to_string( src.GetAngOrder() )
                    + "' to OrdinateSet '"
                    + OrdinateSet::OrdinateType_to_String.at( dst.GetOrdinateType() )
                    + ":" + std::to_string( dst.GetAngOrder() )
                    + "'.\n";

                PRINT_ERROR( error_message.c_str() )
                throw std::invalid_argument( error_message );
            }

            // Set BLAS op to None to prevent allocation of temporary memory in Relabel routines.
            this->blas_op = RelabelBLASOp::None;

            /* No setup necessary for this RelabelType. */

            break;
        }

        default:
        {   std::string error_message = "Invalid RelabelType '"
                                        + RelabelType_to_String.at( relabel_type )
                                        + "' for constructing relabeling operator.\n";

            PRINT_ERROR( error_message.c_str() )
            throw std::invalid_argument( error_message );
        }
    }

# if 0
    // Compute SVD of relabel operator and print singular values.
    char JOBZN  = 'N';  // Do not compute singular vectors.

    double * S  = new double[ src_nq ];
    int * IWORK = new int[ 8 * dst_nq ];
    int INFO    = -1;
    int LWORK   = -1;

    double * WORK = new double;

    dgesdd_( &JOBZN,            // = 'N' --> Do not compute singular vectors.
             &dst.nq(),         // M     = number of rows in the matrix A.
             &src.nq(),         // N     = number of columns in the matrix A.
             R_mat_ptr,         // A     = pointer to the matrix A.
             &dst.nq(),         // LDA   = leading dimension of the matrix A.
             S,                 // S     = pointer to array for storing singular values.
             nullptr,           // U     --> if JOBZ = 'N' then U is not referenced.
             &dst.nq(),         // LDU   = leading dimension of the matrix U.
             nullptr,           // VT    --> if JOBZ = 'N' then VT is not referenced.
             &src.nq(),         // LDVT  = leading dimension of the matrix VT.
             WORK,              // WORK  = double precision workspace array.
             &LWORK,            // LWORK = dimension of the array WORK.
             IWORK,             // IWORK = integer array workspace.
             &INFO );           // INFO  = 0 if successful.

    LWORK = *WORK;
    delete WORK;
    WORK = new double[ LWORK ];

    dgesdd_( &JOBZN,            // = 'N' --> Do not compute singular vectors.
             &dst.nq(),         // M     = number of rows in the matrix A.
             &src.nq(),         // N     = number of columns in the matrix A.
             R_mat_ptr,         // A     = pointer to the matrix A.
             &dst.nq(),         // LDA   = leading dimension of the matrix A.
             S,                 // S     = pointer to array for storing singular values.
             nullptr,           // U     --> if JOBZ = 'N' then U is not referenced.
             &dst.nq(),         // LDU   = leading dimension of the matrix U.
             nullptr,           // VT    --> if JOBZ = 'N' then VT is not referenced.
             &src.nq(),         // LDVT  = leading dimension of the matrix VT.
             WORK,              // WORK  = double precision workspace array.
             &LWORK,            // LWORK = dimension of the array WORK.
             IWORK,             // IWORK = integer array workspace.
             &INFO );           // INFO  = 0 if successful.

    PRINT_LOG( "%d\n", INFO )

    PRINT_LOG( "\nSingular values of relabeling operator:\n\n" )
    for ( int64_t i = 0; i < src_nq; ++i )
        PRINT_LOG( "%4d  %.8e\n", i, S[i] )

    delete [] S;
    delete [] IWORK;
    delete [] WORK;
# endif // if 0
}

# endif // if SPACE_DIMS >= 2


//------------------------------------------------------------------------------------------------------------
//! \brief  Construct a RelabelOperator object using the given parameters.
//!
//! Delegates to RelabelOperator( const Quadrule::OrdinateSet::OrdinateType, const int64_t,
//!                               const Quadrule::OrdinateSet::OrdinateType, const int64_t ).
//!
//! \param[in]  src_in  OrdinateSet to relabel from.
//! \param[in]  dst_in  OrdinateSet to relabel onto.
//------------------------------------------------------------------------------------------------------------
RelabelOperator::RelabelOperator (

    const OrdinateSet & src_in,
    const OrdinateSet & dst_in
) : RelabelOperator( src_in.GetOrdinateType(), src_in.GetAngOrder(), src_in.GetOrdinateSymmetry(),
                     dst_in.GetOrdinateType(), dst_in.GetAngOrder(), dst_in.GetOrdinateSymmetry() )
{}



//------------------------------------------------------------------------------------------------------------
//! \brief  Destructor for class RelabelOperator.
//------------------------------------------------------------------------------------------------------------
RelabelOperator::~RelabelOperator( void ) {

# if SPACE_DIMS >= 2

    std::free( src_basis_ptr );
    std::free( dst_basis_ptr );

# endif // if SPACE_DIMS >= 2

    std::free( R_mat_ptr );
}


//============================================================================================================
//=== ADDITIONAL MEMBER FUNCTIONS ============================================================================
//============================================================================================================


//------------------------------------------------------------------------------------------------------------
//! \param[in]  prefix  String prepended to each line of output, intended to be used for controlling
//!                     indentation of output.
//------------------------------------------------------------------------------------------------------------
void RelabelOperator::Print (

    const std::string & prefix // = "  "

) const {

    PRINT_LOG( "%s%-*s  %s\n", prefix.c_str(),
               Global::col_width, "Relabel type:",
               RelabelType_to_String.at( relabel_type ).c_str() )

    PRINT_LOG( "%s%-*s  %s\n", prefix.c_str(),
               Global::col_width, "Relabel BLAS Op:",
               RelabelBLASOp_to_String.at( blas_op ).c_str() )

# if LOGLEVEL >= 2

    PRINT_LOG( "%sFrom:\n", prefix.c_str() )
    this->src.Print( prefix + "    " );

    PRINT_LOG( "%sTo:\n", prefix.c_str() )
    this->dst.Print( prefix + "    " );

# endif // if LOGLEVEL >= 2
}


//============================================================================================================
//=== RELABEL ROUTINES =======================================================================================
//============================================================================================================


//------------------------------------------------------------------------------------------------------------
//! \brief  Computes \f$ \vec{y} \gets \vec{y} + \alpha \mathcal{R} \vec{x} \f$.
//!
//! Computes \f$ \vec{y} \gets \vec{y} + \alpha \mathcal{R} \vec{x} \f$ where \f$ \mathcal{R} \f$ is the
//! angular flux relabeling operator which maps \f$ \vec{x} \f$ to the discrete ordinate quadrature used by
//! \f$ \vec{y} \f$.
//!
//! \todo   Use OpDomain specification in RelabelOperator application routines.
//!
//! \param[in]      alpha       Scalar to augment \f$ \mathcal{R} \f$.
//! \param[in]      x           Angular flux to relabel and scale.
//! \param[in,out]  y           Angular flux to which the scaled and relabeled flux is added to.
//------------------------------------------------------------------------------------------------------------
void RelabelOperator::Relabel (

    double alpha,
    const RKDG::OrdinateFlux & x,
    RKDG::OrdinateFlux & y
) {

    PRINT_STATUS( "Relabeling angular flux from %" PRId64 " angles to %" PRId64 " angles.\n",
                  src.nq(), dst.nq() )

    // Perform consistency checks between input objects and this.
    if ( !OrdinateSet::AreMatching( this->src, x ) ) {

        std::string error_message
            = "Function '" + std::string(__func__) + "' given OrdinateFlux using OrdinateSet '"
            + OrdinateSet::OrdinateType_to_String.at( x.GetOrdinateType() ) + ":"
            + std::to_string( x.GetAngOrder() )
            + "' but RelabelOperator is configured to map from OrdinateSet '"
            + OrdinateSet::OrdinateType_to_String.at( this->src.GetOrdinateType() ) + ":"
            + std::to_string( this->src.GetAngOrder() ) + "'.\n";

        PRINT_ERROR( error_message.c_str() )
        throw std::invalid_argument( error_message );
    }

    if ( !OrdinateSet::AreMatching( this->dst, y ) ) {

        std::string error_message
            = "Function '" + std::string(__func__) + "' given OrdinateFlux using OrdinateSet '"
            + OrdinateSet::OrdinateType_to_String.at( y.GetOrdinateType() ) + ":"
            + std::to_string( y.GetAngOrder() )
            + "' but RelabelOperator is configured to map to OrdinateSet '"
            + OrdinateSet::OrdinateType_to_String.at( this->dst.GetOrdinateType() ) + ":"
            + std::to_string( this->dst.GetAngOrder() ) + "'.\n";

        PRINT_ERROR( error_message.c_str() )
        throw std::invalid_argument( error_message );
    }

    Global::TMR_AF_Rmv.Start();

    // Number of spatial degrees of freedom per spatial cell.
    int DGO
        # if SPACE_DIMS == 1
            = y.DG_degree + 1;
        # elif SPACE_DIMS == 2
            = (y.DG_degree + 1)*(y.DG_degree + 1);
        # elif SPACE_DIMS == 3
            = (y.DG_degree + 1)*(y.DG_degree + 1)*(y.DG_degree + 1);
        # endif

    size_t sizeof_interp_coeffs = 0;

    switch ( this->blas_op ) {

    # if SPACE_DIMS >= 2

        case RelabelBLASOp::TwoMat_GEMV:

            sizeof_interp_coeffs = stride * sizeof(double);
            break;

        case RelabelBLASOp::TwoMat_GEMM:

            sizeof_interp_coeffs = stride * (y.DG_degree + 1)*(y.DG_degree + 1) * sizeof(double);
            break;

    # endif // if SPACE_DIMS >= 2

        default: break;
    }

    // --- Apply relabel operator. -------------------------------------------------------------------- //

    # pragma omp parallel
    {
        double * interp_coeffs = nullptr;

        // Allocate temporary memory if necessary.
        if ( sizeof_interp_coeffs != 0 ) {

        # if defined (ENABLE_HWLOC) && defined (USE_HWLOC_ALLOC)

            int tid
                # if defined (_OPENMP)
                    = omp_get_thread_num();
                # else
                    = 0;
                # endif

            interp_coeffs = (double *)
                hwloc_alloc_membind( Global::machine_topology, sizeof_interp_coeffs, Global::thread_masks[tid],
                                     HWLOC_MEMBIND_BIND, HWLOC_MEMBIND_THREAD );
        # else

            interp_coeffs = (double *)
                # if defined (USE_ALIGNED_ALLOC)
                    aligned_alloc( ALIGNED_ALLOC_ALIGNMENT, sizeof_interp_coeffs );
                # else
                    std::malloc( sizeof_interp_coeffs );
                # endif

        # endif // if defined (ENABLE_HWLOC) && defined (USE_HWLOC_ALLOC)
        }

        // Compute offsets for BLAS routines.
    # if SPACE_DIMS == 1

        int src_INCX = x.Index(1,0,0) - x.Index(0,0,0);
        int dst_INCY = y.Index(1,0,0) - y.Index(0,0,0);

    # elif SPACE_DIMS == 2

        int src_INCX = x.Index(1,0,0,0,0) - x.Index(0,0,0,0,0);
        int dst_INCY = y.Index(1,0,0,0,0) - y.Index(0,0,0,0,0);

    # elif SPACE_DIMS == 3

        int src_INCX = x.Index(1,0,0,0,0,0,0) - x.Index(0,0,0,0,0,0,0);
        int dst_INCY = y.Index(1,0,0,0,0,0,0) - y.Index(0,0,0,0,0,0,0);

    # endif // if SPACE_DIMS == ?

        // Additional auxiliary variables for BLAS routines.
        double DBL_ONE  = 1.0;

    # if SPACE_DIMS >= 2
        double DBL_ZERO = 0.0;
        int INT_ONE     = 1;
    # endif

        char TRANSN = 'N';
        char TRANST = 'T';

        int src_nq = src.nq();
        int dst_nq = dst.nq();


        switch ( this->relabel_type ) {

        # if SPACE_DIMS == 1
            case RelabelType::PolyInterp:
            case RelabelType::DoublePolyInterp:
        # elif SPACE_DIMS >= 2
            case RelabelType::SphHyperinterp:
            case RelabelType::SphLeastSquares:
        # endif
            {
                switch ( this->blas_op ) {

                # if SPACE_DIMS == 1

                    case RelabelBLASOp::OneMat_GEMV:
                    {
                        # pragma omp for
                        for ( int64_t i = 0; i <= y.nx(0) + 1; ++i ) {
                        for ( int64_t d = 0; d <= y.DG_degree; ++d ) {

                            dgemv_( &TRANSN,                // = 'N' --> A*x
                                    &dst_nq,                // M     = number of rows in the matrix A.
                                    &src_nq,                // N     = number of columns in the matrix A.
                                    &alpha,                 // ALPHA = scalar for A*x.
                                    R_mat_ptr,              // A     = pointer to the matrix A.
                                    &dst_nq,                // LDA   = leading dimension of the matrix A.
                                    x.PointerAt(0,i,d),     // X     = vector x to apply matrix A to.
                                    &src_INCX,              // INVX  = increment for the elements of x.
                                    &DBL_ONE,               // BETA  = scalar for Y.
                                    y.PointerAt(0,i,d),     // Y     = vector to store result.
                                    &dst_INCY );            // INCY  = increment for the elements of Y.
                        }}
                    } break;

                    case RelabelBLASOp::OneMat_GEMM:
                    {
                        # pragma omp for
                        for ( int64_t i = 0; i <= y.nx(0) + 1; ++i ) {

                            dgemm_( &TRANSN,                // = 'N' --> op(A) = A
                                    &TRANST,                // = 'T' --> op(B) = B**T
                                    &DGO,                   // M     = number of rows in the matrix A and C.
                                    &dst_nq,                // N     = number of columns in the matrix C and B**T.
                                    &src_nq,                // K     = number of columns of A, rows of B**T.
                                    &alpha,                 // ALPHA = scalar for ABᵀ
                                    x.PointerAt(0,i,0),     // A     = pointer to matrix A.
                                    &DGO,                   // LDA   = leading dimension of the matrix A.
                                    R_mat_ptr,              // B     = pointer to matrix B.
                                    &dst_nq,                // LDB   = leading dimension of the matrix B.
                                    &DBL_ONE,               // BETA  = scalar for C.
                                    y.PointerAt(0,i,0),     // C     = pointer to matrix C.
                                    &DGO );                 // LDC   = leading dimension of the matrix C.
                        }
                    } break;

                # elif SPACE_DIMS == 2

                    case RelabelBLASOp::OneMat_GEMV:
                    {
                        const int64_t (& nx) [SPACE_DIMS] = y.nx();

                        # pragma omp for collapse(2)
                        for ( int64_t i = 0; i <= nx[0] + 1;   ++i ) {
                        for ( int64_t j = 0; j <= nx[1] + 1;   ++j ) {
                        for ( int64_t d = 0; d <= y.DG_degree; ++d ) {
                        for ( int64_t e = 0; e <= y.DG_degree; ++e ) {

                            dgemv_( &TRANSN,                // = 'N' --> A*x
                                    &dst_nq,                // M     = number of rows in the matrix A.
                                    &src_nq,                // N     = number of columns in the matrix A.
                                    &alpha,                 // ALPHA = scalar for A*x.
                                    R_mat_ptr,              // A     = pointer to the matrix A.
                                    &dst_nq,                // LDA   = leading dimension of the matrix A.
                                    x.PointerAt(0,i,j,d,e), // X     = vector x to apply matrix A to.
                                    &src_INCX,              // INVX  = increment for the elements of x.
                                    &DBL_ONE,               // BETA  = scalar for Y.
                                    y.PointerAt(0,i,j,d,e), // Y     = vector to store result.
                                    &dst_INCY );            // INCY  = increment for the elements of Y.
                        }}}}
                    } break;

                    case RelabelBLASOp::OneMat_GEMM:
                    {
                        const int64_t (& nx) [SPACE_DIMS] = y.nx();

                        # pragma omp for collapse(2)
                        for ( int64_t i = 0; i <= nx[0] + 1; ++i ) {
                        for ( int64_t j = 0; j <= nx[1] + 1; ++j ) {

                            dgemm_( &TRANSN,                // = 'N' --> op(A) = A
                                    &TRANST,                // = 'T' --> op(B) = B**T
                                    &DGO,                   // M     = number of rows in the matrix A and C.
                                    &dst_nq,                // N     = number of columns in the matrix C and B**T.
                                    &src_nq,                // K     = number of columns of A, rows of B**T.
                                    &alpha,                 // ALPHA = scalar for ABᵀ
                                    x.PointerAt(0,i,j,0,0), // A     = pointer to matrix A.
                                    &DGO,                   // LDA   = leading dimension of the matrix A.
                                    R_mat_ptr,              // B     = pointer to matrix B.
                                    &dst_nq,                // LDB   = leading dimension of the matrix B.
                                    &DBL_ONE,               // BETA  = scalar for C.
                                    y.PointerAt(0,i,j,0,0), // C     = pointer to matrix C.
                                    &DGO );                 // LDC   = leading dimension of the matrix C.
                        }}
                    } break;

                    case RelabelBLASOp::TwoMat_GEMV:
                    {
                        const int64_t (& nx) [SPACE_DIMS] = y.nx();

                        # pragma omp for collapse(2)
                        for ( int64_t i = 0; i <= nx[0] + 1;   ++i ) {
                        for ( int64_t j = 0; j <= nx[1] + 1;   ++j ) {
                        for ( int64_t d = 0; d <= y.DG_degree; ++d ) {
                        for ( int64_t e = 0; e <= y.DG_degree; ++e ) {

                            dgemv_( &TRANSN,                // = 'N' --> A*x
                                    &stride,                // M     = number of rows in the matrix A.
                                    &src_nq,                // N     = number of columns in the matrix A.
                                    &DBL_ONE,               // ALPHA = scalar for A*x.
                                    src_basis_ptr,          // A     = pointer to the matrix A.
                                    &stride,                // LDA   = leading dimension of the matrix A.
                                    x.PointerAt(0,i,j,d,e), // X     = vector x to apply matrix A to.
                                    &src_INCX,              // INVX  = increment for the elements of x.
                                    &DBL_ZERO,              // BETA  = scalar for Y.
                                    interp_coeffs,          // Y     = vector to store result.
                                    &INT_ONE );             // INCY  = increment for the elements of Y.

                            dgemv_( &TRANST,                // = 'T' --> Aᵀ*x
                                    &stride,                // M     = number of rows in the matrix A.
                                    &dst_nq,                // N     = number of columns in the matrix A.
                                    &alpha,                 // ALPHA = scalar for Aᵀ*x
                                    dst_basis_ptr,          // A     = pointer to the matrix A.
                                    &stride,                // LDA   = leading dimension of the matrix A.
                                    interp_coeffs,          // X     = vector x to apply matrix Aᵀ to.
                                    &INT_ONE,               // INCX  = increment for the elements of x.
                                    &DBL_ONE,               // BETA  = scalar for Y.
                                    y.PointerAt(0,i,j,d,e), // Y     = vector to store result.
                                    &dst_INCY );            // INCY  = increment for the elements of Y.
                        }}}}
                    } break;

                    case RelabelBLASOp::TwoMat_GEMM:
                    {
                        const int64_t (& nx) [SPACE_DIMS] = y.nx();

                        # pragma omp for collapse(2)
                        for ( int64_t i = 0; i <= nx[0] + 1; ++i ) {
                        for ( int64_t j = 0; j <= nx[1] + 1; ++j ) {

                            dgemm_( &TRANSN,                // = 'N' --> op(A) = A
                                    &TRANST,                // = 'T' --> op(B) = B**T
                                    &DGO,                   // M     = number of rows in the matrix A and C.
                                    &stride,                // N     = number of columns in the matrix C and B**T.
                                    &src_nq,                // K     = number of columns of A, rows of B**T.
                                    &DBL_ONE,               // ALPHA = scalar for ABᵀ
                                    x.PointerAt(0,i,j,0,0), // A     = pointer to matrix A.
                                    &DGO,                   // LDA   = leading dimension of the matrix A.
                                    src_basis_ptr,          // B     = pointer to matrix B.
                                    &stride,                // LDB   = leading dimension of the matrix B.
                                    &DBL_ZERO,              // BETA  = scalar for C.
                                    interp_coeffs,          // C     = pointer to matrix C.
                                    &DGO );                 // LDC   = leading dimension of the matrix C.

                            dgemm_( &TRANSN,                // = 'N' --> op(A) = A
                                    &TRANSN,                // = 'N' --> op(B) = B
                                    &DGO,                   // M     = number of rows in the matrix A and C.
                                    &dst_nq,                // N     = number of columns in the matrix C and B.
                                    &stride,                // K     = number of columns of A, rows of B.
                                    &alpha,                 // ALPHA = scalar for ABᵀ
                                    interp_coeffs,          // A     = pointer to matrix A.
                                    &DGO,                   // LDA   = leading dimension of the matrix A.
                                    dst_basis_ptr,          // B     = pointer to matrix B.
                                    &stride,                // LDB   = leading dimension of the matrix B.
                                    &DBL_ONE,               // BETA  = scalar for C.
                                    y.PointerAt(0,i,j,0,0), // C     = pointer to matrix C.
                                    &DGO );                 // LDC   = leading dimension of the matrix C.
                        }}
                    } break;

                # elif SPACE_DIMS == 3

                    case RelabelBLASOp::OneMat_GEMV:
                    {
                        const int64_t (& nx) [SPACE_DIMS] = y.nx();

                        # pragma omp for collapse(3)
                        for ( int64_t i = 0; i <= nx[0] + 1;   ++i ) {
                        for ( int64_t j = 0; j <= nx[1] + 1;   ++j ) {
                        for ( int64_t k = 0; k <= nx[2] + 1;   ++k ) {
                        for ( int64_t d = 0; d <= y.DG_degree; ++d ) {
                        for ( int64_t e = 0; e <= y.DG_degree; ++e ) {
                        for ( int64_t f = 0; f <= y.DG_degree; ++f ) {

                            dgemv_( &TRANSN,                    // = 'N' --> A*x
                                    &dst_nq,                    // M     = number of rows in the matrix A.
                                    &src_nq,                    // N     = number of columns in the matrix A.
                                    &alpha,                     // ALPHA = scalar for A*x.
                                    R_mat_ptr,                  // A     = pointer to the matrix A.
                                    &dst_nq,                    // LDA   = leading dimension of the matrix A.
                                    x.PointerAt(0,i,j,k,d,e,f), // X     = vector x to apply matrix A to.
                                    &src_INCX,                  // INVX  = increment for the elements of x.
                                    &DBL_ONE,                   // BETA  = scalar for Y.
                                    y.PointerAt(0,i,j,k,d,e,f), // Y     = vector to store result.
                                    &dst_INCY );                // INCY  = increment for the elements of Y.
                        }}}}}}
                    } break;

                    case RelabelBLASOp::OneMat_GEMM:
                    {
                        const int64_t (& nx) [SPACE_DIMS] = y.nx();

                        # pragma omp for collapse(3)
                        for ( int64_t i = 0; i <= nx[0] + 1; ++i ) {
                        for ( int64_t j = 0; j <= nx[1] + 1; ++j ) {
                        for ( int64_t k = 0; k <= nx[2] + 1; ++k ) {

                            dgemm_( &TRANSN,                    // = 'N' --> op(A) = A
                                    &TRANST,                    // = 'T' --> op(B) = B**T
                                    &DGO,                       // M     = number of rows in the matrix A and C.
                                    &dst_nq,                    // N     = number of columns in the matrix C and B**T.
                                    &src_nq,                    // K     = number of columns of A, rows of B**T.
                                    &alpha,                     // ALPHA = scalar for ABᵀ
                                    x.PointerAt(0,i,j,k,0,0,0), // A     = pointer to matrix A.
                                    &DGO,                       // LDA   = leading dimension of the matrix A.
                                    R_mat_ptr,                  // B     = pointer to matrix B.
                                    &dst_nq,                    // LDB   = leading dimension of the matrix B.
                                    &DBL_ONE,                   // BETA  = scalar for C.
                                    y.PointerAt(0,i,j,k,0,0,0), // C     = pointer to matrix C.
                                    &DGO );                     // LDC   = leading dimension of the matrix C.
                        }}}
                    } break;

                    case RelabelBLASOp::TwoMat_GEMV:
                    {
                        const int64_t (& nx) [SPACE_DIMS] = y.nx();

                        # pragma omp for collapse(3)
                        for ( int64_t i = 0; i <= nx[0] + 1;   ++i ) {
                        for ( int64_t j = 0; j <= nx[1] + 1;   ++j ) {
                        for ( int64_t k = 0; k <= nx[2] + 1;   ++k ) {
                        for ( int64_t d = 0; d <= y.DG_degree; ++d ) {
                        for ( int64_t e = 0; e <= y.DG_degree; ++e ) {
                        for ( int64_t f = 0; f <= y.DG_degree; ++f ) {

                            dgemv_( &TRANSN,                    // = 'N' --> A*x
                                    &stride,                    // M     = number of rows in the matrix A.
                                    &src_nq,                    // N     = number of columns in the matrix A.
                                    &DBL_ONE,                   // ALPHA = scalar for A*x.
                                    src_basis_ptr,              // A     = pointer to the matrix A.
                                    &stride,                    // LDA   = leading dimension of the matrix A.
                                    x.PointerAt(0,i,j,k,d,e,f), // X     = vector x to apply matrix A to.
                                    &src_INCX,                  // INVX  = increment for the elements of x.
                                    &DBL_ZERO,                  // BETA  = scalar for Y.
                                    interp_coeffs,              // Y     = vector to store result.
                                    &INT_ONE );                 // INCY  = increment for the elements of Y.

                            dgemv_( &TRANST,                    // = 'T' --> Aᵀ*x
                                    &stride,                    // M     = number of rows in the matrix A.
                                    &dst_nq,                    // N     = number of columns in the matrix A.
                                    &alpha,                     // ALPHA = scalar for Aᵀ*x
                                    dst_basis_ptr,              // A     = pointer to the matrix A.
                                    &stride,                    // LDA   = leading dimension of the matrix A.
                                    interp_coeffs,              // X     = vector x to apply matrix Aᵀ to.
                                    &INT_ONE,                   // INCX  = increment for the elements of x.
                                    &DBL_ONE,                   // BETA  = scalar for Y.
                                    y.PointerAt(0,i,j,k,d,e,f), // Y     = vector to store result.
                                    &dst_INCY );                // INCY  = increment for the elements of Y.
                        }}}}}}
                    } break;

                    case RelabelBLASOp::TwoMat_GEMM:
                    {
                        const int64_t (& nx) [SPACE_DIMS] = y.nx();

                        # pragma omp for collapse(3)
                        for ( int64_t i = 0; i <= nx[0] + 1; ++i ) {
                        for ( int64_t j = 0; j <= nx[1] + 1; ++j ) {
                        for ( int64_t k = 0; k <= nx[2] + 1; ++k ) {

                            dgemm_( &TRANSN,                    // = 'N' --> op(A) = A
                                    &TRANST,                    // = 'T' --> op(B) = B**T
                                    &DGO,                       // M     = number of rows in the matrix A and C.
                                    &stride,                    // N     = number of columns in the matrix C and B**T.
                                    &src_nq,                    // K     = number of columns of A, rows of B**T.
                                    &DBL_ONE,                   // ALPHA = scalar for ABᵀ
                                    x.PointerAt(0,i,j,k,0,0,0), // A     = pointer to matrix A.
                                    &DGO,                       // LDA   = leading dimension of the matrix A.
                                    src_basis_ptr,              // B     = pointer to matrix B.
                                    &stride,                    // LDB   = leading dimension of the matrix B.
                                    &DBL_ZERO,                  // BETA  = scalar for C.
                                    interp_coeffs,              // C     = pointer to matrix C.
                                    &DGO );                     // LDC   = leading dimension of the matrix C.

                            dgemm_( &TRANSN,                    // = 'N' --> op(A) = A
                                    &TRANSN,                    // = 'N' --> op(B) = B
                                    &DGO,                       // M     = number of rows in the matrix A and C.
                                    &dst_nq,                    // N     = number of columns in the matrix C and B.
                                    &stride,                    // K     = number of columns of A, rows of B.
                                    &alpha,                     // ALPHA = scalar for ABᵀ
                                    interp_coeffs,              // A     = pointer to matrix A.
                                    &DGO,                       // LDA   = leading dimension of the matrix A.
                                    dst_basis_ptr,              // B     = pointer to matrix B.
                                    &stride,                    // LDB   = leading dimension of the matrix B.
                                    &DBL_ONE,                   // BETA  = scalar for C.
                                    y.PointerAt(0,i,j,k,0,0,0), // C     = pointer to matrix C.
                                    &DGO );                     // LDC   = leading dimension of the matrix C.
                        }}}
                    } break;

                # endif // if SPACE_DIMS == ?

                    default:
                        # pragma omp master
                    {   std::string error_message = "Invalid RelabelBLASOp '"
                                                    + RelabelBLASOp_to_String.at( this->blas_op )
                                                    + "' in '" + std::string(__func__) + "'.\n";

                        PRINT_ERROR( error_message.c_str() )
                        throw std::invalid_argument( error_message );
                    }
                }
            } break;

        # if SPACE_DIMS >= 2
            case RelabelType::PWC_Hierarchical:
        # endif
            {

            # if SPACE_DIMS == 2

                const int64_t (& nx) [SPACE_DIMS] = y.nx();

                # pragma omp for collapse(2)
                for ( int64_t i = 0; i <= nx[0] + 1; ++i ) {
                for ( int64_t j = 0; j <= nx[1] + 1; ++j ) {

                    for ( int64_t qc = 0; qc < x.nq();          ++qc ) {
                    for ( int64_t l  = 0; l  < y.nq() / x.nq(); ++l  ) {

                        const int64_t qu = l + qc * ( y.nq() / x.nq() );

                        for ( int64_t d = 0; d <= y.DG_degree; ++d ) {
                        for ( int64_t e = 0; e <= y.DG_degree; ++e ) {

                            y(qu,i,j,d,e) += alpha * x(qc,i,j,d,e);
                        }}
                    }}
                }}

            # elif SPACE_DIMS == 3

                # pragma omp parallel for collapse(3)
                for ( int64_t i = 0; i <= nx[0] + 1; ++i ) {
                for ( int64_t j = 0; j <= nx[1] + 1; ++j ) {
                for ( int64_t k = 0; k <= nx[2] + 1; ++k ) {

                    for ( int64_t qc = 0; qc < x.nq();          ++qc ) {
                    for ( int64_t l  = 0; l  < y.nq() / x.nq(); ++l  ) {

                        const int64_t qu = l + qc * ( y.nq() / x.nq() );

                        for ( int64_t d = 0; d <= y.DG_degree; ++d ) {
                        for ( int64_t e = 0; e <= y.DG_degree; ++e ) {
                        for ( int64_t f = 0; f <= y.DG_degree; ++f ) {

                            y(qu,i,j,k,d,e,f) += alpha * x(qc,i,j,k,d,e,f);
                        }}}
                    }}
                }}}

            # endif // if SPACE_DIMS == ?
            } break;

            default:
                # pragma omp master
            {   std::string error_message = "Invalid RelabelType '"
                                            + RelabelType_to_String.at( this->relabel_type )
                                            + "' in '" + std::string(__func__) + "'.\n";

                PRINT_ERROR( error_message.c_str() )
                throw std::invalid_argument( error_message );
            }
        }

        if ( sizeof_interp_coeffs != 0 ) {

        # if defined (ENABLE_HWLOC) && defined (USE_HWLOC_ALLOC)
            hwloc_free( Global::machine_topology, interp_coeffs, sizeof_interp_coeffs );
        # else
            std::free( interp_coeffs );
        # endif
        }
    }

    y.halo_cells_dirty |= x.halo_cells_dirty;
    y.upwind_halo_cells_dirty |= x.halo_cells_dirty;

    Global::TMR_AF_Rmv.Stop();
}


//------------------------------------------------------------------------------------------------------------
//! \brief  Computes \f$ \vec{y} \gets \vec{y} + \alpha \mathcal{R} \vec{x} \f$.
//!
//! Computes \f$ \vec{y} \gets \vec{y} + \alpha \mathcal{R} \vec{x} \f$ where \f$ \mathcal{R} \f$ is the
//! angular flux relabeling operator which maps \f$ \vec{x} \f$ to the discrete ordinate quadrature used by
//! \f$ \vec{y} \f$.
//!
//! The STDG::OrdinateFlux \pp{x} is evaluated at the endpoint \f$ t_{n+1} \f$ of the timestep interval
//! before relabeling into the RKDG::OrdinateFlux \pp{y}.
//!
//! \param[in]      alpha       Scalar to augment \f$ \mathcal{R} \f$.
//! \param[in]      x           Angular flux to relabel and scale.
//! \param[in,out]  y           Angular flux to which the scaled and relabeled flux is added to.
//------------------------------------------------------------------------------------------------------------
void RelabelOperator::Relabel (

    double alpha,
    const STDG::OrdinateFlux & x,
    RKDG::OrdinateFlux & y
) {

    PRINT_STATUS( "Relabeling angular flux from %" PRId64 " angles to %" PRId64 " angles.\n",
                  src.nq(), dst.nq() )

    // Perform consistency checks between input objects and this.
    if ( !OrdinateSet::AreMatching( this->src, x ) ) {

        std::string error_message
            = "Function '" + std::string(__func__) + "' given OrdinateFlux using OrdinateSet '"
            + OrdinateSet::OrdinateType_to_String.at( x.GetOrdinateType() ) + ":"
            + std::to_string( x.GetAngOrder() )
            + "' but RelabelOperator is configured to map from OrdinateSet '"
            + OrdinateSet::OrdinateType_to_String.at( this->src.GetOrdinateType() ) + ":"
            + std::to_string( this->src.GetAngOrder() ) + "'.\n";

        PRINT_ERROR( error_message.c_str() )
        throw std::invalid_argument( error_message );
    }

    if ( !OrdinateSet::AreMatching( this->dst, y ) ) {

        std::string error_message
            = "Function '" + std::string(__func__) + "' given OrdinateFlux using OrdinateSet '"
            + OrdinateSet::OrdinateType_to_String.at( y.GetOrdinateType() ) + ":"
            + std::to_string( y.GetAngOrder() )
            + "' but RelabelOperator is configured to map to OrdinateSet '"
            + OrdinateSet::OrdinateType_to_String.at( this->dst.GetOrdinateType() ) + ":"
            + std::to_string( this->dst.GetAngOrder() ) + "'.\n";

        PRINT_ERROR( error_message.c_str() )
        throw std::invalid_argument( error_message );
    }

    Global::TMR_AF_Rmv.Start();

    // Number of spatial degrees of freedom per spatial cell.
    int DGO
        # if SPACE_DIMS == 1
            = y.DG_degree + 1;
        # elif SPACE_DIMS == 2
            = (y.DG_degree + 1)*(y.DG_degree + 1);
        # elif SPACE_DIMS == 3
            = (y.DG_degree + 1)*(y.DG_degree + 1)*(y.DG_degree + 1);
        # endif

    size_t sizeof_interp_coeffs = 0;
    size_t sizeof_ST_coeffs = 0;

    switch ( blas_op ) {

        case RelabelBLASOp::OneMat_GEMV:

            sizeof_ST_coeffs = src.nq() * sizeof(double);
            break;

        case RelabelBLASOp::OneMat_GEMM:

            sizeof_ST_coeffs = src.nq() * DGO * sizeof(double);
            break;

    # if SPACE_DIMS >= 2

        case RelabelBLASOp::TwoMat_GEMV:

            sizeof_ST_coeffs = src.nq() * sizeof(double);
            sizeof_interp_coeffs = stride * sizeof(double);
            break;

        case RelabelBLASOp::TwoMat_GEMM:

            sizeof_ST_coeffs = src.nq() * DGO * sizeof(double);
            sizeof_interp_coeffs = stride * DGO * sizeof(double);
            break;

    # endif // if SPACE_DIMS >= 2

        default: break;
    }

    // --- Apply relabeling operator. ----------------------------------------------------------------- //

    # pragma omp parallel
    {
        // Pointers for thread-local buffers.
        double * ST_coeffs = nullptr;
        double * interp_coeffs = nullptr;

        // Allocate buffer ST_coeffs.
    # if defined (ENABLE_HWLOC) && defined (USE_HWLOC_ALLOC)

        int tid
            # if defined (_OPENMP)
                = omp_get_thread_num();
            # else
                = 0;
            # endif

        ST_coeffs = (double *)
            hwloc_alloc_membind( Global::machine_topology, sizeof_ST_coeffs, Global::thread_masks[tid],
                                 HWLOC_MEMBIND_BIND, HWLOC_MEMBIND_THREAD );
    # else

        ST_coeffs = (double *)
            # if defined (USE_ALIGNED_ALLOC)
                aligned_alloc( ALIGNED_ALLOC_ALIGNMENT, sizeof_ST_coeffs );
            # else
                std::malloc( sizeof_ST_coeffs );
            # endif

    # endif // if defined (ENABLE_HWLOC) && defined (USE_HWLOC_ALLOC)

        // Allocate buffer interp_coeffs for two-matrix relabel.
        if ( sizeof_interp_coeffs != 0 ) {

        # if defined (ENABLE_HWLOC) && defined (USE_HWLOC_ALLOC)

            interp_coeffs = (double *)
                hwloc_alloc_membind( Global::machine_topology, sizeof_interp_coeffs, Global::thread_masks[tid],
                                     HWLOC_MEMBIND_BIND, HWLOC_MEMBIND_THREAD );
        # else

            interp_coeffs = (double *)
                # if defined (USE_ALIGNED_ALLOC)
                    aligned_alloc( ALIGNED_ALLOC_ALIGNMENT, sizeof_interp_coeffs );
                # else
                    std::malloc( sizeof_interp_coeffs );
                # endif

        # endif // if defined (ENABLE_HWLOC) && defined (USE_HWLOC_ALLOC)
        }

        // Compute offsets for BLAS routines, and define macros for indexing into temporary arrays for STDG
        // coefficients.
    # if SPACE_DIMS == 1

        int dst_INCY = y.Index(1,0,0) - y.Index(0,0,0);

        # define IST_coeffs(x,q,d) ( (d) + ((x).DG_degree_x + 1)*(q) )

    # elif SPACE_DIMS == 2

        int dst_INCY = y.Index(1,0,0,0,0) - y.Index(0,0,0,0,0);

        # define IST_coeffs(x,q,d,e) ( (e) + ((x).DG_degree_x + 1)*( (d) + ((x).DG_degree_x + 1)*(q) ) )

    # elif SPACE_DIMS == 3

        int dst_INCY = y.Index(1,0,0,0,0,0,0) - y.Index(0,0,0,0,0,0,0);

        # define IST_coeffs(x,q,d,e,f) ( (f) + ((x).DG_degree_x + 1)*( (e) + ((x).DG_degree_x + 1)*( (d) + ((x).DG_degree_x + 1)*(q) ) ) )

    # endif // if SPACE_DIMS == ?

        // Additional auxiliary variables for BLAS routines.
        double DBL_ONE  = 1.0;
    # if SPACE_DIMS >= 2
        double DBL_ZERO = 0.0;
    # endif
        int INT_ONE     = 1;

        char TRANSN = 'N';
        char TRANST = 'T';

        int src_nq = src.nq();
        int dst_nq = dst.nq();


        switch ( this->relabel_type ) {

        # if SPACE_DIMS == 1
            case RelabelType::PolyInterp:
            case RelabelType::DoublePolyInterp:
        # elif SPACE_DIMS >= 2
            case RelabelType::SphHyperinterp:
            case RelabelType::SphLeastSquares:
        # endif
            {
                switch ( this->blas_op ) {

                # if SPACE_DIMS == 1

                    case RelabelBLASOp::OneMat_GEMV:
                    {
                        # pragma omp for
                        for ( int64_t i = 0; i <= y.nx(0) + 1; ++i ) {
                        for ( int64_t d = 0; d <= y.DG_degree; ++d ) {

                            std::memset( ST_coeffs, 0, sizeof_ST_coeffs );

                            for ( int64_t q = 0; q <  x.nq();        ++q ) {
                            for ( int64_t s = 0; s <= x.DG_degree_t; ++s ) {

                                ST_coeffs[q] += x(q,i,d,s);
                            }}

                            dgemv_( &TRANSN,                // = 'N' --> A*x
                                    &dst_nq,                // M     = number of rows in the matrix A.
                                    &src_nq,                // N     = number of columns in the matrix A.
                                    &alpha,                 // ALPHA = scalar for A*x.
                                    R_mat_ptr,              // A     = pointer to the matrix A.
                                    &dst_nq,                // LDA   = leading dimension of the matrix A.
                                    ST_coeffs,              // X     = vector x to apply matrix A to.
                                    &INT_ONE,               // INVX  = increment for the elements of x.
                                    &DBL_ONE,               // BETA  = scalar for Y.
                                    y.PointerAt(0,i,d),     // Y     = vector to store result.
                                    &dst_INCY );            // INCY  = increment for the elements of Y.
                        }}
                    } break;

                    case RelabelBLASOp::OneMat_GEMM:
                    {
                        # pragma omp for
                        for ( int64_t i = 0; i <= y.nx(0) + 1;  ++i ) {

                            std::memset( ST_coeffs, 0, sizeof_ST_coeffs );

                            for ( int64_t q = 0; q <  x.nq();        ++q ) {
                            for ( int64_t d = 0; d <= x.DG_degree_x; ++d ) {
                            for ( int64_t s = 0; s <= x.DG_degree_t; ++s ) {

                                ST_coeffs[ IST_coeffs(x,q,d) ] += x(q,i,d,s);
                            }}}

                            dgemm_( &TRANSN,                // = 'N' --> op(A) = A
                                    &TRANST,                // = 'T' --> op(B) = B**T
                                    &DGO,                   // M     = number of rows in the matrix A and C.
                                    &dst_nq,                // N     = number of columns in the matrix C and B**T.
                                    &src_nq,                // K     = number of columns of A, rows of B**T.
                                    &alpha,                 // ALPHA = scalar for ABᵀ
                                    ST_coeffs,              // A     = pointer to matrix A.
                                    &DGO,                   // LDA   = leading dimension of the matrix A.
                                    R_mat_ptr,              // B     = pointer to matrix B.
                                    &dst_nq,                // LDB   = leading dimension of the matrix B.
                                    &DBL_ONE,               // BETA  = scalar for C.
                                    y.PointerAt(0,i,0),     // C     = pointer to matrix C.
                                    &DGO );                 // LDC   = leading dimension of the matrix C.
                        }
                    } break;

                # elif SPACE_DIMS == 2

                    case RelabelBLASOp::OneMat_GEMV:
                    {
                        const int64_t (& nx) [SPACE_DIMS] = y.nx();

                        # pragma omp for collapse(2)
                        for ( int64_t i = 0; i <= nx[0] + 1;   ++i ) {
                        for ( int64_t j = 0; j <= nx[1] + 1;   ++j ) {
                        for ( int64_t d = 0; d <= y.DG_degree; ++d ) {
                        for ( int64_t e = 0; e <= y.DG_degree; ++e ) {

                            std::memset( ST_coeffs, 0, sizeof_ST_coeffs );

                            for ( int64_t q = 0; q <  x.nq();        ++q ) {
                            for ( int64_t s = 0; s <= x.DG_degree_t; ++s ) {

                                ST_coeffs[q] += x(q,i,j,d,e,s);
                            }}

                            dgemv_( &TRANSN,                // = 'N' --> A*x
                                    &dst_nq,                // M     = number of rows in the matrix A.
                                    &src_nq,                // N     = number of columns in the matrix A.
                                    &alpha,                 // ALPHA = scalar for A*x.
                                    R_mat_ptr,              // A     = pointer to the matrix A.
                                    &dst_nq,                // LDA   = leading dimension of the matrix A.
                                    ST_coeffs,              // X     = vector x to apply matrix A to.
                                    &INT_ONE,               // INVX  = increment for the elements of x.
                                    &DBL_ONE,               // BETA  = scalar for Y.
                                    y.PointerAt(0,i,j,d,e), // Y     = vector to store result.
                                    &dst_INCY );            // INCY  = increment for the elements of Y.
                        }}}}
                    } break;

                    case RelabelBLASOp::OneMat_GEMM:
                    {
                        const int64_t (& nx) [SPACE_DIMS] = y.nx();

                        # pragma omp for collapse(2)
                        for ( int64_t i = 0; i <= nx[0] + 1; ++i ) {
                        for ( int64_t j = 0; j <= nx[1] + 1; ++j ) {

                            std::memset( ST_coeffs, 0, sizeof_ST_coeffs );

                            for ( int64_t q = 0; q <  x.nq();        ++q ) {
                            for ( int64_t d = 0; d <= x.DG_degree_x; ++d ) {
                            for ( int64_t e = 0; e <= x.DG_degree_x; ++e ) {
                            for ( int64_t s = 0; s <= x.DG_degree_t; ++s ) {

                                ST_coeffs[ IST_coeffs(x,q,d,e) ] += x(q,i,j,d,e,s);
                            }}}}

                            dgemm_( &TRANSN,                // = 'N' --> op(A) = A
                                    &TRANST,                // = 'T' --> op(B) = B**T
                                    &DGO,                   // M     = number of rows in the matrix A and C.
                                    &dst_nq,                // N     = number of columns in the matrix C and B**T.
                                    &src_nq,                // K     = number of columns of A, rows of B**T.
                                    &alpha,                 // ALPHA = scalar for ABᵀ
                                    ST_coeffs,              // A     = pointer to matrix A.
                                    &DGO,                   // LDA   = leading dimension of the matrix A.
                                    R_mat_ptr,              // B     = pointer to matrix B.
                                    &dst_nq,                // LDB   = leading dimension of the matrix B.
                                    &DBL_ONE,               // BETA  = scalar for C.
                                    y.PointerAt(0,i,j,0,0), // C     = pointer to matrix C.
                                    &DGO );                 // LDC   = leading dimension of the matrix C.
                        }}
                    } break;

                    case RelabelBLASOp::TwoMat_GEMV:
                    {
                        const int64_t (& nx) [SPACE_DIMS] = y.nx();

                        # pragma omp for collapse(2)
                        for ( int64_t i = 0; i <= nx[0] + 1;   ++i ) {
                        for ( int64_t j = 0; j <= nx[1] + 1;   ++j ) {
                        for ( int64_t d = 0; d <= y.DG_degree; ++d ) {
                        for ( int64_t e = 0; e <= y.DG_degree; ++e ) {

                            std::memset( ST_coeffs, 0, sizeof_ST_coeffs );

                            for ( int64_t q = 0; q <  x.nq();        ++q ) {
                            for ( int64_t s = 0; s <= x.DG_degree_t; ++s ) {

                                ST_coeffs[q] += x(q,i,j,d,e,s);
                            }}

                            dgemv_( &TRANSN,                // = 'N' --> A*x
                                    &stride,                // M     = number of rows in the matrix A.
                                    &src_nq,                // N     = number of columns in the matrix A.
                                    &DBL_ONE,               // ALPHA = scalar for A*x.
                                    src_basis_ptr,          // A     = pointer to the matrix A.
                                    &stride,                // LDA   = leading dimension of the matrix A.
                                    ST_coeffs,              // X     = vector x to apply matrix A to.
                                    &INT_ONE,               // INVX  = increment for the elements of x.
                                    &DBL_ZERO,              // BETA  = scalar for Y.
                                    interp_coeffs,          // Y     = vector to store result.
                                    &INT_ONE );             // INCY  = increment for the elements of Y.

                            dgemv_( &TRANST,                // = 'T' --> Aᵀ*x
                                    &stride,                // M     = number of rows in the matrix A.
                                    &dst_nq,                // N     = number of columns in the matrix A.
                                    &alpha,                 // ALPHA = scalar for Aᵀ*x
                                    dst_basis_ptr,          // A     = pointer to the matrix A.
                                    &stride,                // LDA   = leading dimension of the matrix A.
                                    interp_coeffs,          // X     = vector x to apply matrix Aᵀ to.
                                    &INT_ONE,               // INCX  = increment for the elements of x.
                                    &DBL_ONE,               // BETA  = scalar for Y.
                                    y.PointerAt(0,i,j,d,e), // Y     = vector to store result.
                                    &dst_INCY );            // INCY  = increment for the elements of Y.
                        }}}}
                    } break;

                    case RelabelBLASOp::TwoMat_GEMM:
                    {
                        const int64_t (& nx) [SPACE_DIMS] = y.nx();

                        # pragma omp for collapse(2)
                        for ( int64_t i = 0; i <= nx[0] + 1; ++i ) {
                        for ( int64_t j = 0; j <= nx[1] + 1; ++j ) {

                            std::memset( ST_coeffs, 0, sizeof_ST_coeffs );

                            for ( int64_t q = 0; q <  x.nq();        ++q ) {
                            for ( int64_t d = 0; d <= x.DG_degree_x; ++d ) {
                            for ( int64_t e = 0; e <= x.DG_degree_x; ++e ) {
                            for ( int64_t s = 0; s <= x.DG_degree_t; ++s ) {

                                ST_coeffs[ IST_coeffs(x,q,d,e) ] += x(q,i,j,d,e,s);
                            }}}}

                            dgemm_( &TRANSN,                // = 'N' --> op(A) = A
                                    &TRANST,                // = 'T' --> op(B) = B**T
                                    &DGO,                   // M     = number of rows in the matrix A and C.
                                    &stride,                // N     = number of columns in the matrix C and B**T.
                                    &src_nq,                // K     = number of columns of A, rows of B**T.
                                    &DBL_ONE,               // ALPHA = scalar for ABᵀ
                                    ST_coeffs,              // A     = pointer to matrix A.
                                    &DGO,                   // LDA   = leading dimension of the matrix A.
                                    src_basis_ptr,          // B     = pointer to matrix B.
                                    &stride,                // LDB   = leading dimension of the matrix B.
                                    &DBL_ZERO,              // BETA  = scalar for C.
                                    interp_coeffs,          // C     = pointer to matrix C.
                                    &DGO );                 // LDC   = leading dimension of the matrix C.

                            dgemm_( &TRANSN,                // = 'N' --> op(A) = A
                                    &TRANSN,                // = 'N' --> op(B) = B
                                    &DGO,                   // M     = number of rows in the matrix A and C.
                                    &dst_nq,                // N     = number of columns in the matrix C and B.
                                    &stride,                // K     = number of columns of A, rows of B.
                                    &alpha,                 // ALPHA = scalar for ABᵀ
                                    interp_coeffs,          // A     = pointer to matrix A.
                                    &DGO,                   // LDA   = leading dimension of the matrix A.
                                    dst_basis_ptr,          // B     = pointer to matrix B.
                                    &stride,                // LDB   = leading dimension of the matrix B.
                                    &DBL_ONE,               // BETA  = scalar for C.
                                    y.PointerAt(0,i,j,0,0), // C     = pointer to matrix C.
                                    &DGO );                 // LDC   = leading dimension of the matrix C.
                        }}
                    } break;

                # elif SPACE_DIMS == 3

                    case RelabelBLASOp::OneMat_GEMV:
                    {
                        const int64_t (& nx) [SPACE_DIMS] = y.nx();

                        # pragma omp for collapse(3)
                        for ( int64_t i = 0; i <= nx[0] + 1;   ++i ) {
                        for ( int64_t j = 0; j <= nx[1] + 1;   ++j ) {
                        for ( int64_t k = 0; k <= nx[2] + 1;   ++k ) {
                        for ( int64_t d = 0; d <= y.DG_degree; ++d ) {
                        for ( int64_t e = 0; e <= y.DG_degree; ++e ) {
                        for ( int64_t f = 0; f <= y.DG_degree; ++f ) {

                            std::memset( ST_coeffs, 0, sizeof_ST_coeffs );

                            for ( int64_t q = 0; q <  x.nq();        ++q ) {
                            for ( int64_t s = 0; s <= x.DG_degree_t; ++s ) {

                                ST_coeffs[q] += x(q,i,j,k,d,e,f,s);
                            }}

                            dgemv_( &TRANSN,                    // = 'N' --> A*x
                                    &dst_nq,                    // M     = number of rows in the matrix A.
                                    &src_nq,                    // N     = number of columns in the matrix A.
                                    &alpha,                     // ALPHA = scalar for A*x.
                                    R_mat_ptr,                  // A     = pointer to the matrix A.
                                    &dst_nq,                    // LDA   = leading dimension of the matrix A.
                                    ST_coeffs,                  // X     = vector x to apply matrix A to.
                                    &INT_ONE,                   // INVX  = increment for the elements of x.
                                    &DBL_ONE,                   // BETA  = scalar for Y.
                                    y.PointerAt(0,i,j,k,d,e,f), // Y = vector to store result.
                                    &dst_INCY );                // INCY  = increment for the elements of Y.
                        }}}}}}
                    } break;

                    case RelabelBLASOp::OneMat_GEMM:
                    {
                        const int64_t (& nx) [SPACE_DIMS] = y.nx();

                        # pragma omp for collapse(3)
                        for ( int64_t i = 0; i <= nx[0] + 1; ++i ) {
                        for ( int64_t j = 0; j <= nx[1] + 1; ++j ) {
                        for ( int64_t k = 0; k <= nx[2] + 1; ++k ) {

                            std::memset( ST_coeffs, 0, sizeof_ST_coeffs );

                            for ( int64_t q = 0; q <  x.nq();        ++q ) {
                            for ( int64_t d = 0; d <= x.DG_degree_x; ++d ) {
                            for ( int64_t e = 0; e <= x.DG_degree_x; ++e ) {
                            for ( int64_t f = 0; f <= x.DG_degree_x; ++f ) {
                            for ( int64_t s = 0; s <= x.DG_degree_t; ++s ) {

                                ST_coeffs[ IST_coeffs(x,q,d,e,f) ] += x(q,i,j,k,d,e,f,s);
                            }}}}}

                            dgemm_( &TRANSN,                    // = 'N' --> op(A) = A
                                    &TRANST,                    // = 'T' --> op(B) = B**T
                                    &DGO,                       // M     = number of rows in the matrix A and C.
                                    &dst_nq,                    // N     = number of columns in the matrix C and B**T.
                                    &src_nq,                    // K     = number of columns of A, rows of B**T.
                                    &alpha,                     // ALPHA = scalar for ABᵀ
                                    ST_coeffs,                  // A     = pointer to matrix A.
                                    &DGO,                       // LDA   = leading dimension of the matrix A.
                                    R_mat_ptr,                  // B     = pointer to matrix B.
                                    &dst_nq,                    // LDB   = leading dimension of the matrix B.
                                    &DBL_ONE,                   // BETA  = scalar for C.
                                    y.PointerAt(0,i,j,k,0,0,0), // C = pointer to matrix C.
                                    &DGO );                     // LDC   = leading dimension of the matrix C.
                        }}}
                    } break;

                    case RelabelBLASOp::TwoMat_GEMV:
                    {
                        const int64_t (& nx) [SPACE_DIMS] = y.nx();

                        # pragma omp for collapse(3)
                        for ( int64_t i = 0; i <= nx[0] + 1;   ++i ) {
                        for ( int64_t j = 0; j <= nx[1] + 1;   ++j ) {
                        for ( int64_t k = 0; k <= nx[2] + 1;   ++k ) {
                        for ( int64_t d = 0; d <= y.DG_degree; ++d ) {
                        for ( int64_t e = 0; e <= y.DG_degree; ++e ) {
                        for ( int64_t f = 0; f <= y.DG_degree; ++f ) {

                            std::memset( ST_coeffs, 0, sizeof_ST_coeffs );

                            for ( int64_t q = 0; q <  x.nq();        ++q ) {
                            for ( int64_t s = 0; s <= x.DG_degree_x; ++s ) {

                                ST_coeffs[q] += x(q,i,j,k,d,e,f,s);
                            }}

                            dgemv_( &TRANSN,                    // = 'N' --> A*x
                                    &stride,                    // M     = number of rows in the matrix A.
                                    &src_nq,                    // N     = number of columns in the matrix A.
                                    &DBL_ONE,                   // ALPHA = scalar for A*x.
                                    src_basis_ptr,              // A     = pointer to the matrix A.
                                    &stride,                    // LDA   = leading dimension of the matrix A.
                                    ST_coeffs,                  // X     = vector x to apply matrix A to.
                                    &INT_ONE,                   // INVX  = increment for the elements of x.
                                    &DBL_ZERO,                  // BETA  = scalar for Y.
                                    interp_coeffs,              // Y     = vector to store result.
                                    &INT_ONE );                 // INCY  = increment for the elements of Y.

                            dgemv_( &TRANST,                    // = 'T' --> Aᵀ*x
                                    &stride,                    // M     = number of rows in the matrix A.
                                    &dst_nq,                    // N     = number of columns in the matrix A.
                                    &alpha,                     // ALPHA = scalar for Aᵀ*x
                                    dst_basis_ptr,              // A     = pointer to the matrix A.
                                    &stride,                    // LDA   = leading dimension of the matrix A.
                                    interp_coeffs,              // X     = vector x to apply matrix Aᵀ to.
                                    &INT_ONE,                   // INCX  = increment for the elements of x.
                                    &DBL_ONE,                   // BETA  = scalar for Y.
                                    y.PointerAt(0,i,j,k,d,e,f), // Y     = vector to store result.
                                    &dst_INCY );                // INCY  = increment for the elements of Y.
                        }}}}}}
                    } break;

                    case RelabelBLASOp::TwoMat_GEMM:
                    {
                        const int64_t (& nx) [SPACE_DIMS] = y.nx();

                        # pragma omp for collapse(3)
                        for ( int64_t i = 0; i <= nx[0] + 1; ++i ) {
                        for ( int64_t j = 0; j <= nx[1] + 1; ++j ) {
                        for ( int64_t k = 0; k <= nx[2] + 1; ++k ) {

                            std::memset( ST_coeffs, 0, sizeof_ST_coeffs );

                            for ( int64_t q = 0; q <  x.nq();        ++q ) {
                            for ( int64_t d = 0; d <= x.DG_degree_x; ++d ) {
                            for ( int64_t e = 0; e <= x.DG_degree_x; ++e ) {
                            for ( int64_t f = 0; f <= x.DG_degree_x; ++f ) {
                            for ( int64_t s = 0; s <= x.DG_degree_t; ++s ) {

                                ST_coeffs[ IST_coeffs(x,q,d,e,f) ] += x(q,i,j,k,d,e,f,s);
                            }}}}}

                            dgemm_( &TRANSN,                    // = 'N' --> op(A) = A
                                    &TRANST,                    // = 'T' --> op(B) = B**T
                                    &DGO,                       // M     = number of rows in the matrix A and C.
                                    &stride,                    // N     = number of columns in the matrix C and B**T.
                                    &src_nq,                    // K     = number of columns of A, rows of B**T.
                                    &DBL_ONE,                   // ALPHA = scalar for ABᵀ
                                    ST_coeffs,                  // A     = pointer to matrix A.
                                    &DGO,                       // LDA   = leading dimension of the matrix A.
                                    src_basis_ptr,              // B     = pointer to matrix B.
                                    &stride,                    // LDB   = leading dimension of the matrix B.
                                    &DBL_ZERO,                  // BETA  = scalar for C.
                                    interp_coeffs,              // C     = pointer to matrix C.
                                    &DGO );                     // LDC   = leading dimension of the matrix C.

                            dgemm_( &TRANSN,                    // = 'N' --> op(A) = A
                                    &TRANSN,                    // = 'N' --> op(B) = B
                                    &DGO,                       // M     = number of rows in the matrix A and C.
                                    &dst_nq,                    // N     = number of columns in the matrix C and B.
                                    &stride,                    // K     = number of columns of A, rows of B.
                                    &alpha,                     // ALPHA = scalar for ABᵀ
                                    interp_coeffs,              // A     = pointer to matrix A.
                                    &DGO,                       // LDA   = leading dimension of the matrix A.
                                    dst_basis_ptr,              // B     = pointer to matrix B.
                                    &stride,                    // LDB   = leading dimension of the matrix B.
                                    &DBL_ONE,                   // BETA  = scalar for C.
                                    y.PointerAt(0,i,j,k,0,0,0), // C = pointer to matrix C.
                                    &DGO );                     // LDC   = leading dimension of the matrix C.
                        }}}
                    } break;

                # endif // if SPACE_DIMS == ?

                    default:
                        # pragma omp master
                    {   std::string error_message = "Invalid RelabelBLASOp '"
                                                    + RelabelBLASOp_to_String.at( this->blas_op )
                                                    + "' in '" + std::string(__func__) + "'.\n";

                        PRINT_ERROR( error_message.c_str() )
                        throw std::invalid_argument( error_message );
                    }
                }
            } break;

        # if SPACE_DIMS >= 2
            case RelabelType::PWC_Hierarchical:
        # endif
            {
            # if SPACE_DIMS == 2

                const int64_t (& nx) [SPACE_DIMS] = y.nx();

                # pragma omp for collapse(2)
                for ( int64_t i = 0; i <= nx[0] + 1; ++i ) {
                for ( int64_t j = 0; j <= nx[1] + 1; ++j ) {

                    for ( int64_t qc = 0; qc < x.nq();          ++qc ) {
                    for ( int64_t l  = 0; l  < y.nq() / x.nq(); ++l  ) {

                        const int64_t qu = l + qc * ( y.nq() / x.nq() );

                        for ( int64_t d = 0; d <= x.DG_degree_x; ++d ) {
                        for ( int64_t e = 0; e <= x.DG_degree_x; ++e ) {
                        for ( int64_t s = 0; s <= x.DG_degree_t; ++s ) {

                            y(qu,i,j,d,e) += alpha * x(qc,i,j,d,e,s);
                        }}}
                    }}
                }}

            # elif SPACE_DIMS == 3

                # pragma omp for collapse(3)
                for ( int64_t i = 0; i <= nx[0] + 1; ++i ) {
                for ( int64_t j = 0; j <= nx[1] + 1; ++j ) {
                for ( int64_t k = 0; k <= nx[2] + 1; ++k ) {

                    for ( int64_t qc = 0; qc < x.nq();          ++qc ) {
                    for ( int64_t l  = 0; l  < y.nq() / x.nq(); ++l  ) {

                        const int64_t qu = l + qc * ( y.nq() / x.nq() );

                        for ( int64_t d = 0; d <= x.DG_degree_x; ++d ) {
                        for ( int64_t e = 0; e <= x.DG_degree_x; ++e ) {
                        for ( int64_t f = 0; f <= x.DG_degree_x; ++f ) {
                        for ( int64_t s = 0; s <= x.DG_degree_t; ++s ) {

                            y(qu,i,j,k,d,e,f) += alpha * x(qc,i,j,k,d,e,f,s);
                        }}}}
                    }}
                }}}

            # endif // if SPACE_DIMS == ?
            } break;

            default:
                # pragma omp master
            {   std::string error_message = "Invalid RelabelType '"
                                            + RelabelType_to_String.at( this->relabel_type )
                                            + "' in '" + std::string(__func__) + "'.\n";

                PRINT_ERROR( error_message.c_str() )
                throw std::invalid_argument( error_message );
            }
        }


    # undef IST_coeffs

    # if defined (ENABLE_HWLOC) && defined (USE_HWLOC_ALLOC)
        hwloc_free( Global::machine_topology, ST_coeffs, sizeof_ST_coeffs );
    # else
        std::free( ST_coeffs );
    # endif

        if ( sizeof_interp_coeffs != 0 ) {

        # if defined (ENABLE_HWLOC) && defined (USE_HWLOC_ALLOC)
            hwloc_free( Global::machine_topology, interp_coeffs, sizeof_interp_coeffs );
        # else
            std::free( interp_coeffs );
        # endif
        }
    }

    y.halo_cells_dirty |= x.halo_cells_dirty;
    y.upwind_halo_cells_dirty |= x.halo_cells_dirty;

    Global::TMR_AF_Rmv.Stop();
}
