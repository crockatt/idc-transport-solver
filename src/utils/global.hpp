//------------------------------------------------------------------------------------------------------------
//! \file   utils/global.hpp
//! \brief  Header file for declaring global functions and variables.
//!
//! \author Michael M. Crockatt, Kris Garrett
//! \date   January 2018
//------------------------------------------------------------------------------------------------------------

# ifndef __GLOBAL_HPP__
# define __GLOBAL_HPP__


# include <cinttypes>
# include <cmath>
# include <cstdint>
# include <string>
# include <type_traits>

# if defined (ENABLE_MPI) || defined (DOXYCOMPILE)
    # include <mpi.h>
# endif

# if defined (ENABLE_HWLOC) || defined (DOXYCOMPILE)
    # include <hwloc.h>
# endif

# if defined (_OPENMP)
    # include <omp.h>
# endif


//
// Define alignment macro used in calls to aligned_alloc.
//
# if defined (USE_ALIGNED_ALLOC) || defined (DOXYCOMPILE)

    //!
    //! \brief  Memory alignment (in bytes) used by calls to \c aligned_alloc.
    //!
    # if !defined (ALIGNED_ALLOC_ALIGNMENT)
        # define ALIGNED_ALLOC_ALIGNMENT 4096
    # endif

# endif // if defined (USE_ALIGNED_ALLOC)


# include "utils/bitmask_operators.hpp"
# include "utils/ParameterList.hpp"
# include "utils/Walltimer.hpp"


//------------------------------------------------------------------------------------------------------------
//! \brief  Template used to check that the input variable has integral type.
//!
//! \param[in]  input   Variable to check type of.
//------------------------------------------------------------------------------------------------------------
template<class T>
void assert_integral_type ( T input ) {

    static_assert( std::is_integral<T>::value, "Integral value required." );
    (void) input;
}


//------------------------------------------------------------------------------------------------------------
//! \brief  Macro for checking if a dimension index is valid for the current compilation environment.
//!
//! \param[in]  dim     Dimension index to check.
//------------------------------------------------------------------------------------------------------------
# define CHECK_DIM_VALUE(dim) \
\
    assert_integral_type(dim); \
\
    if ( dim < 0 || dim >= SPACE_DIMS ) { \
\
        std::string error_message =   "Spatial dimension " \
                                    + std::to_string(dim) \
                                    + " not in valid range [0," \
                                    + std::to_string( SPACE_DIMS ) \
                                    + ") for " \
                                    + std::string(__func__) \
                                    + ".\n"; \
\
        PRINT_ERROR( error_message.c_str() ) \
        throw std::out_of_range( error_message ); \
    }


//------------------------------------------------------------------------------------------------------------
//! \brief  Flags used to specify which edges of the domain boundary a condition may be satisfied or to which
//!         an operation should be applied.
//------------------------------------------------------------------------------------------------------------
enum class BoundaryEdge {

    //!
    //! Null value. Specifies no boundaries.
    //!
    None = 0x00

    //!
    //! Denotes the boundary at which the spatial coordinate in the first dimension is fixed at its minimum
    //! value.
    //!
,   X_Min = 0x01

    //!
    //! Denotes the boundary at which the spatial coordinate in the first dimension is fixed at its maximum
    //! value.
    //!
,   X_Max = 0x02

# if SPACE_DIMS >= 2 || defined (DOXYCOMPILE)

    //!
    //! Denotes the boundary at which the spatial coordinate in the second dimension is fixed at its minimum
    //! value.
    //!
,   Y_Min = 0x04

    //!
    //! Denotes the boundary at which the spatial coordinate in the second dimension is fixed at its maximum
    //! value.
    //!
,   Y_Max = 0x08

# if SPACE_DIMS == 3 || defined (DOXYCOMPILE)

    //!
    //! Denotes the boundary at which the spatial coordinate in the third dimension is fixed at its minimum
    //! value.
    //!
,   Z_Min = 0x10

    //!
    //! Denotes the boundary at which the spatial coordinate in the third dimension is fixed at its maximum
    //! value.
    //!
,   Z_Max = 0x20

# endif // if SPACE_DIMS == 3
# endif // if SPACE_DIMS >= 2

    //!
    //! Combined bitmask for all boundary edges in all spatial dimensions.
    //!
,   All =
        # if SPACE_DIMS == 1
            0x03
        # elif SPACE_DIMS == 2 || defined (DOXYCOMPILE)
            0x0F
        # elif SPACE_DIMS == 3
            0x3F
        # endif
};

ENABLE_BITMASK_OPERATORS( BoundaryEdge )


namespace ExitCode {

//------------------------------------------------------------------------------------------------------------
//! \brief  Specifies exit codes for program termination.
//------------------------------------------------------------------------------------------------------------
enum ExitCode_enum {

    //!
    //! \brief  Successful completion of execution.
    //!
    Success = EXIT_SUCCESS,

    //!
    //! \brief  Generic error code.
    //!
    Failure = EXIT_FAILURE,

    //!
    //! \brief  Indicates that program terminated by checkpoint.
    //!
    Checkpoint = 150

};

}


//============================================================================================================
//=== GLOBAL VARIABLES =======================================================================================
//============================================================================================================


//!
//! \brief  A hack so global variables only need to be declared in this file.
//!
# ifdef NO_EXTERN_GLOBAL
    # define EXTERN_GLOBAL
# else
    # define EXTERN_GLOBAL extern
# endif


//!
//! \brief  Contains global variables.
//!
namespace Global {


    //!
    //! \brief  Used to map from a (dimension, edge) index pair to the associated member of the BoundaryEdge
    //!         enum class.
    //!
    //! Used to determine whether certain edges should be included in operations using a dimension-independent
    //! implementation.
    //!
# if ! defined (NO_EXTERN_GLOBAL)
    EXTERN_GLOBAL const BoundaryEdge DimEdgeIdx_to_BoundaryEdge [SPACE_DIMS][2];
# endif


# if defined (ENABLE_MPI) || defined (DOXYCOMPILE)

    //!
    //! \brief  Holds the total number of MPI ranks in the given simulation.
    //!
    EXTERN_GLOBAL int MPI_num_ranks;

    //!
    //! \brief  Holds the MPI rank of the current process in the given simulation.
    //!
    EXTERN_GLOBAL int MPI_rank;

    //!
    //! \brief  MPI communicator containing topology information of spatial domain decomposition.
    //!
    EXTERN_GLOBAL MPI_Comm MPI_cart_comm;

# if defined (ENABLE_KBA_WAIT_TIMING) || defined (DOXYCOMPILE)

    //!
    //! \brief  Records the time the communication team spends waiting for computation to complete in
    //!         SweepPattern::KBA.
    //!
    EXTERN_GLOBAL Walltimer TMR_KBA_comm_wait;

    //!
    //! \brief  Records the time the computation team spends waiting for communication to complete in
    //!         SweepPattern::KBA.
    //!
    EXTERN_GLOBAL Walltimer TMR_KBA_comp_wait;

# endif // if defined (ENABLE_KBA_WAIT_TIMING)
# endif // if defined (ENABLE_MPI)


//!
//! \brief  Contains the prefix used for all output files from the current run.
//!
EXTERN_GLOBAL std::string file_prefix;


//!
//! \brief  Specifies the width of columns used in text outputs.
//!
# if ! defined (NO_EXTERN_GLOBAL) || defined (DOXYCOMPILE)
    extern const int col_width;
# endif


//!
//! \brief  Maximum partial degree of DG basis functions.
//!
EXTERN_GLOBAL int64_t DG_degree;


//!
//! \brief  Pointer to array containing coefficients of the Legendre polynomial triple product tensor.
//!
//! \see    \c #ITPI
//! \see    ComputeTPI()
//!
EXTERN_GLOBAL double * TPI;


//!
//! \brief  Determines if the spatial domain is periodic or not.
//!
EXTERN_GLOBAL bool periodic;


//!
//! \brief  Specifies edges of the domain boundary with reflecting boundary conditions.
//!
EXTERN_GLOBAL BoundaryEdge reflecting_boundaries;


//!
//! \brief  Contains the input parameters.
//!
EXTERN_GLOBAL ParameterList input_list;


//!
//! \brief  Timer for AF_Lsv().
//!
EXTERN_GLOBAL Walltimer TMR_AF_Lsv;


//!
//! \brief  Timer for AF_Lmv().
//!
EXTERN_GLOBAL Walltimer TMR_AF_Lmv;


//!
//! \brief  Timer for AF_Rmv().
//!
EXTERN_GLOBAL Walltimer TMR_AF_Rmv;


//!
//! \brief  Timer for Pmv().
//!
EXTERN_GLOBAL Walltimer TMR_Pmv;


//!
//! \brief  Timer for AF_Smv().
//!
EXTERN_GLOBAL Walltimer TMR_AF_Smv;


//!
//! \brief  Timer for AF_axpy().
//!
EXTERN_GLOBAL Walltimer TMR_AF_axpy;


//!
//! \brief  Timer for AF_scal().
//!
EXTERN_GLOBAL Walltimer TMR_AF_scal;


//!
//! \brief  Timer for AF_copy().
//!
EXTERN_GLOBAL Walltimer TMR_AF_copy;


//!
//! \brief  Timer for AF_zeroCoefficients().
//!
EXTERN_GLOBAL Walltimer TMR_AF_zero;


//!
//! \brief  Timer for OrdinateFlux::ReflectBoundaries().
//!
EXTERN_GLOBAL Walltimer TMR_AF_boundary;


//!
//! \brief  Timer for SD_zeroCoefficients().
//!
EXTERN_GLOBAL Walltimer TMR_SD_zero;


//!
//! \brief  Timer for DensityFunction::AXPY().
//!
EXTERN_GLOBAL Walltimer TMR_SD_axpy;


//!
//! \brief  Timer for SD_packPETScVec() and SD_unpackPETScVec().
//!
EXTERN_GLOBAL Walltimer TMR_SD_pack;


//!
//! \brief  Timer for iterative solver internals.
//!
//! For ImplicitSolverOLD::SolveType::Gmres, times the execution of PETSc internal routines.
//!
//! For ImplicitSolverOLD::SolveType::SourceIteration, times the execution of vector norms.
//!
EXTERN_GLOBAL Walltimer TMR_PETSc;


//!
//! \brief  Timer for halo exchange routines.
//!
//! \see    OrdinateFlux::SynchronizeHalos().
//!
EXTERN_GLOBAL Walltimer TMR_SynchronizeHalos;


//!
//! \brief  Timer for assembly of diffusion matrix in DSA preconditioner.
//!
EXTERN_GLOBAL Walltimer TMR_DSA_Assemble;


//!
//! \brief  Timer for setting up and solving diffusion equations in DSA preconditioner.
//!
EXTERN_GLOBAL Walltimer TMR_DSA_Solve;


//!
//! \brief  Timer for collided components of hybrid methods.
//!
EXTERN_GLOBAL Walltimer TMR_collided;


//!
//! \brief  Timer for uncollided components of hybrid methods.
//!
EXTERN_GLOBAL Walltimer TMR_uncollided;


//!
//! \brief  Timer for relabel or reconstruction step of hybrid methods.
//!
EXTERN_GLOBAL Walltimer TMR_relabel;


//!
//! \brief  Timer for disk IO operations.
//!
EXTERN_GLOBAL Walltimer TMR_io;


# if defined (DO_SWEEP_SUBTIMING)

    //!
    //! \brief  Timer for CalcA routines of sweeps.
    //!
    EXTERN_GLOBAL Walltimer * TMR_calcA;


    //!
    //! \brief  Timer for CalcB routines of Sweeps.
    //!
    EXTERN_GLOBAL Walltimer * TMR_calcB;


    //!
    //! \brief  Timer for linear solves of sweeps.
    //!
    EXTERN_GLOBAL Walltimer * TMR_linearSolve;

# endif // if defined (DO_SWEEP_SUBTIMING)


# if defined (ENABLE_HWLOC) || defined (DOXYCOMPILE)

    //!
    //! \brief  Holds the topology information for the machine on which the program is running.
    //!
    EXTERN_GLOBAL hwloc_topology_t machine_topology;

    //!
    //! \brief  Pointer to array of CPU set bitmaps specifying the core binding of each OpenMP thread.
    //!
    EXTERN_GLOBAL hwloc_cpuset_t * thread_masks;

    //!
    //! \brief  Contains the logical "or" of the core binding bit masks of all OpenMP threads.
    //!
    EXTERN_GLOBAL hwloc_bitmap_t active_core_mask;

# endif // if defined (ENABLE_HWLOC)


} // namespace Global


//============================================================================================================
//=== ADDITIONAL SIMPLE ROUTINES =============================================================================
//============================================================================================================


//------------------------------------------------------------------------------------------------------------
//! \brief  Returns the value of \f$ (-1)^{\texttt{pow}} \f$.
//!
//! \param[in]      pow         Integral exponent for -1.
//!
//! \return     Returns -1 to the power \pp{pow}.
//------------------------------------------------------------------------------------------------------------
inline int64_t neg1pow( const int64_t pow ) {  return ( pow % 2 ? -1 : 1 );  }
// inline int64_t neg1pow( const int64_t pow ) {  return 1 - ((pow & 1) << 1 );  }


//------------------------------------------------------------------------------------------------------------
//! \brief  Returns \f$ x^n \f$ for \f$ n \f$ nonnegative.
//!
//! If #STRICT_CHECK is defined, an exception is thrown in \f$ n \f$ is negative.
//!
//! \param[in]      x           Value to raise to the specified power.
//! \param[in]      n           Integral value denoting the exponent.
//!
//! \return     Returns \f$ x^n \f$.
//------------------------------------------------------------------------------------------------------------
template< typename x_type, typename n_type = int64_t >
x_type IntPow (

    x_type x,
    n_type n
) {
    static_assert( std::is_integral<n_type>::value, "Integral type required for exponent." );

# if defined (STRICT_CHECK)

    if ( n < 0 ) {

        std::string error_message =   "Negative value '" + std::to_string(n) + "' is invalid in '"
                                    + std::string(__func__) + "'.\n";

        throw std::invalid_argument( error_message );
    }

# endif // if defined (STRICT_CHECK)

    x_type result = 1;

    while ( n != 0 ) {

        if ( n & 1 ) {  result *= x;  }

        n >>= 1, x *= x;
    }

    return result;
}


//------------------------------------------------------------------------------------------------------------
//! \brief  Computes the integer component of the base-2 logarithm of an integer.
//!
//! \attention  The input \pp{x} must be of integral type.
//!
//! \param[in]  x   Value to compute logarithm for.
//!
//! \return     Returns the integer component of the base-2 logarithm of \pp{x}.
//------------------------------------------------------------------------------------------------------------
template< typename x_type >
constexpr x_type IntLog2 (

    x_type x
) {
    // Make sure input is of integral type.
    static_assert( std::is_integral<x_type>::value,
                   "Integral type required for input argument to IntLog2" );

# if 0
    // Make sure input is non-negative.
    if ( x < 0 ) {

        std::string error_message =   "Negative value '"
                                    + std::to_string(x)
                                    + "' is invalid in '"
                                    + std::string(__func__)
                                    + "'.\n";

        throw std::invalid_argument( error_message );
    }
# endif // if 0

    // Compute integer component of logarithm.
    x_type log_x = 0;

    while ( x >>= 1 ) {  ++log_x;  }

    return log_x;
}


//------------------------------------------------------------------------------------------------------------
//! \brief  Reads a null-terminated C-style string from a binary file at the given file pointer into an
//!         std::string object.
//------------------------------------------------------------------------------------------------------------
void ReadStringFromFile( std::string &,
                    # if defined (ENABLE_MPI)
                         MPI_File &
                    # else
                         std::FILE * &
                    # endif
                         );


//------------------------------------------------------------------------------------------------------------
//! \brief  Writes the contents of an std::string object to a binary file at the given file pointer as a
//!         null-terminated C-style string.
//------------------------------------------------------------------------------------------------------------
void WriteStringToFile( const std::string &,
                    # if defined (ENABLE_MPI)
                        MPI_File &
                    # else
                        std::FILE * &
                    # endif
                        );


//------------------------------------------------------------------------------------------------------------
//! \brief  Returns the degree of the DG approximation for the spatial variables based on the contents of the
//!         given ParameterList.
//------------------------------------------------------------------------------------------------------------
int64_t GetDGDegreeX( const ParameterList & );


//------------------------------------------------------------------------------------------------------------
//! \brief  Returns the degree of the DG approximation in time based on the contents of the given
//!         ParameterList.
//------------------------------------------------------------------------------------------------------------
int64_t GetDGDegreeT( const ParameterList & );


//!
//! \brief  Forward declaration for OutputSolution function definition.
//!
namespace RKDG {

    class OrdinateFlux;
}


//============================================================================================================
//=== MACROS AND ROUTINES FOR TENSOR OF LEGENDRE TRIPLE PRODUCT INTEGRALS ====================================
//============================================================================================================


//------------------------------------------------------------------------------------------------------------
//! \brief  Allocates memory for and computes the full tensor of triple product integrals of Legendre
//!         polynomials.
//------------------------------------------------------------------------------------------------------------
void ComputeTPI();


//------------------------------------------------------------------------------------------------------------
//! \brief  Macro for indexing the array of Legendre triple product integrals #Global::TPI.
//!
//! \see    #Global::TPI
//------------------------------------------------------------------------------------------------------------
# define ITPI(a,b,c) ( (a) + (Global::DG_degree + 1)*( (b) + (Global::DG_degree + 1)*(c) ) )


# endif // ifndef __GLOBAL_HPP__
