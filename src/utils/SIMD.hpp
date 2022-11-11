//------------------------------------------------------------------------------------------------------------
//! \file   utils/SIMD.hpp
//! \brief  Contains some useful definitions for implementing SIMD blocked algorithms.
//!
//! \author Michael M. Crockatt
//! \date   January 2018
//------------------------------------------------------------------------------------------------------------

# ifndef __SIMD_HPP__
# define __SIMD_HPP__


# include <cstdint>


# if defined (ENABLE_SIMD_BLOCKING)

    // Set default vector lengths if not defined.
    # ifndef SIMD_LEN
        # define SIMD_LEN 4
    # endif

    # ifndef SWEEP_SIMD_LEN
        # define SWEEP_SIMD_LEN SIMD_LEN
    # endif

    # ifndef LMV_SIMD_LEN
        # define LMV_SIMD_LEN SIMD_LEN
    # endif

# else // if defined (ENABLE_SIMD_BLOCKING)

    /*  If blocking is disabled, then all vector lengths are forced to be defined to 1. */

    # if defined (SIMD_LEN)
        # warning "Undefining previously defined macro SIMD_LEN."
        # undef SIMD_LEN
    # endif
    # define SIMD_LEN 1

    # if defined (SWEEP_SIMD_LEN)
        # warning "Undefining previously defined macro SWEEP_SIMD_LEN."
        # undef SWEEP_SIMD_LEN
    # endif
    # define SWEEP_SIMD_LEN 1

    # if defined (LMV_SIMD_LEN)
        # warning "Undefining previously defined macro LMV_SIMD_LEN."
        # undef LMV_SIMD_LEN
    # endif
    # define LMV_SIMD_LEN 1

# endif // if defined (ENABLE_SIMD_BLOCKING)


//
// It seems that Intel's compilers refuse to allow the sizeof operator to be used within an aligned clause of
// an OpenMP simd construct; e.g.,
//
//      # pragma omp simd aligned( 4 * sizeof(double) )
//
// This hack seeks to circumvent this limitation. The value SIZEOF_DOUBLE is by default replaced by
// sizeof(double). If Intel's compilers are used, then SIZEOF_DOUBLE is replaced by 8.
//
# if defined (__INTEL_COMPILER)
    # define SIZEOF_DOUBLE 8
# else
    # define SIZEOF_DOUBLE sizeof(double)
# endif


//------------------------------------------------------------------------------------------------------------
//! \brief  Specifies several groupings of indices in SIMD blocks for which optimizations of certain routines
//!         can be performed.
//------------------------------------------------------------------------------------------------------------
enum class BlockType {

    //!
    //! \brief  Most general case. Specifies that spatial cell and ordinate indices may all differ between
    //!         components of the SIMD block. It is assumed that all ordinate indices come from the same
    //!         quadrant (octant) for 2D (3D) implementations.
    //!
    None,

    //!
    //! \brief  Specifies that spatial cell indices are the same for all elements of the SIMD block, and
    //!         that the ordinate indices are sequential (specifically, monotonically increasing).
    //!
    SameCell_ContiguousAngles,

    //!
    //! \brief  Specifies that all elements of the SIMD block have the same ordinate index, but may have
    //!         different spatial cell indices.
    //!
    SameAngle
};


//------------------------------------------------------------------------------------------------------------
//! \brief  Used to store index data for elements in blocks for SIMD blocked methods.
//!
//! This struct stores independent indices for each element of the SIMD block. This allows blocking of
//! elements in arbitrary arrangements. Special arrangements may be specified using the BlockType enumeration.
//------------------------------------------------------------------------------------------------------------
template< int64_t simd_len >
struct SIMD_BlkIdx {

    alignas( simd_len * sizeof(int64_t) ) int64_t q [simd_len];
    alignas( simd_len * sizeof(int64_t) ) int64_t i [simd_len];

# if SPACE_DIMS >= 2

    alignas( simd_len * sizeof(int64_t) ) int64_t j [simd_len];

# if SPACE_DIMS == 3

    alignas( simd_len * sizeof(int64_t) ) int64_t k [simd_len];

# endif // if SPACE_DIMS == 3
# endif // if SPACE_DIMS >= 2

    //!
    //! \brief  Specifies the number of valid elements stored in the SIMD block. Used to determine what
    //!         values should be stored after computation.
    //!
    int64_t len;

# if SPACE_DIMS == 2 || defined (DOXYCOMPILE)

    //!
    //! \brief  Specifies the angular quadrant that the indices stored in SIMD_BlkIdx:q belong to.
    //!
    int64_t quadrant;

# elif SPACE_DIMS == 3

    //!
    //! \brief  Specifies the angular octant that the indices stored in SIMD_BlkIdx:q belong to.
    //!
    int64_t octant;

# endif // if SPACE_DIMS == ?

    //!
    //! \brief  Specifies whether some common structure exists between elements of the SIMD block.
    //!
    BlockType type = BlockType::None;
};


//------------------------------------------------------------------------------------------------------------
//! \brief  Used to store index data for elements in blocks for SIMD blocked methods.
//!
//! This specialization is used for simplification and reducing the struct size for cases where the indices
//! for all elements of the block can be determined directly from the indices of the first element of the
//! block. The documentation for the members of this specialization reflect the fact that this specialization
//! is typically used for blocking sequential ordinate indices with the same spatial cell index.
//------------------------------------------------------------------------------------------------------------
template<>
struct SIMD_BlkIdx<0> {

    int64_t q;      //!< Index of first ordinate in SIMD block.
    int64_t i;      //!< Index of spatial cell in the \f$ x_1 \f$ dimension for all elements in SIMD block.

# if SPACE_DIMS >= 2 || defined (DOXYCOMPILE)

    int64_t j;      //!< Index of spatial cell in the \f$ x_2 \f$ dimension for all elements in SIMD block.

# if SPACE_DIMS == 3

    int64_t k;      //!< Index of spatial cell in the \f$ x_3 \f$ dimension for all elements in SIMD block.

# endif // if SPACE_DIMS == 3
# endif // if SPACE_DIMS >= 2

    int64_t len;    //!< Number of sequential ordinates in SIMD block.

# if SPACE_DIMS == 2 || defined (DOXYCOMPILE)

    //!
    //! \brief  Specifies the angular quadrant that the indices stored in SIMD_BlkIdx:q belong to.
    //!
    int64_t quadrant;

# elif SPACE_DIMS == 3

    //!
    //! \brief  Specifies the angular octant that the indices stored in SIMD_BlkIdx:q belong to.
    //!
    int64_t octant;

# endif // if SPACE_DIMS == ?
};


# endif // ifndef __SIMD_HPP__
