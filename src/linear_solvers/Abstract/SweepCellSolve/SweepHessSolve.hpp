//------------------------------------------------------------------------------------------------------------
//! \file   linear_solvers/Abstract/SweepCellSolve/SweepHessSolve.hpp
//! \brief  Header file containing declaration of SweepHessSolve class.
//!
//! \author Michael M. Crockatt
//! \date   March 2018
//------------------------------------------------------------------------------------------------------------

# ifndef __ABSTRACT__SWEEP_HESS_SOLVE_HPP__
# define __ABSTRACT__SWEEP_HESS_SOLVE_HPP__


# include <cstring>
# include <limits>

# include "linear_solvers/Abstract/SweepOperator.hpp"
# include "linear_solvers/Abstract/SweepCellSolve/SweepLinearSolveSpec.hpp"


//------------------------------------------------------------------------------------------------------------
//! \brief  Declaration of SweepHessSolve class.
//!
//! Solves local linear systems in each element using a Hessenberg reduction.
//!
//! \attention  This solve method assumes that cross section values are zone-wise constant with respect to
//!             space.
//!
//! The idea behind the construction of this solver is as follows. If the total cross section is independent
//! of angle and constant within a spatial cell, then the matrix \f$ A \f$ for the local systems within that
//! cell can be written as \f$ A = A_q + \sigma_{\mathrm{t}} I \f$, where \f$ A_q \f$ depends only on the
//! angular ordinate. Given a Hessenberg decomposition of \f$ A_q \f$ of the form
//! \f$ A_q = Q_q H_q Q_q^T \f$, the local system \f$ Ax = b \f$ can be written as
//! \f[ Q_q \left( H_q + \sigma_{\mathrm{t}} I \right) Q_q^T x = b . \f]
//! The solution is then given by
//! \f[ x = Q_q \left( H_q + \sigma_{\mathrm{t}} I \right)^{-1} Q_q^T b . \f]
//! Because \f$ H_q \f$ is a Hessenberg matrix, so is \f$ H_q + \sigma_{\mathrm{t}} I \f$: Therefore \f$ x \f$
//! can be computed in \f$ O(n^2) \f$ operations, where \f$ n \f$ is the dimension of the linear system,
//! if the Hessenberg decomposition of \f$ A_q \f$ is known.
//!
//! The Hessenberg decomposition of the angular components of the matrices is computed by
//! SweepHessSolve::ComputeMatrices.
//!
//! This class template implements:
//!     - SweepCellSolveBase::SolveSystem
//!
//------------------------------------------------------------------------------------------------------------
template< class OrdinateFlux, int64_t SIMD_length >
class SweepHessSolve :
    public SweepLinearSolveSpec< OrdinateFlux, SIMD_length >
{

public:

    //========================================================================================================
    //=== CONSTRUCTORS AND DESTRUCTOR ========================================================================
    //========================================================================================================

    //--------------------------------------------------------------------------------------------------------
    //! \brief  Delete default constructor.
    //--------------------------------------------------------------------------------------------------------
    SweepHessSolve( void ) = delete;

    //--------------------------------------------------------------------------------------------------------
    //! \brief  Delete copy constructor.
    //--------------------------------------------------------------------------------------------------------
    SweepHessSolve( const SweepHessSolve< OrdinateFlux, SIMD_length > & ) = delete;

    //--------------------------------------------------------------------------------------------------------
    //! \brief  Construct solver object from given parameters.
    //--------------------------------------------------------------------------------------------------------
    SweepHessSolve (
        const SweepOperator<OrdinateFlux> & enclosing,
        const ParameterList & input_list,
        const double dt_in = std::numeric_limits<double>::infinity()
    );

    //--------------------------------------------------------------------------------------------------------
    //! \brief  Virtual destructor.
    //--------------------------------------------------------------------------------------------------------
    virtual ~SweepHessSolve( void );


    //========================================================================================================
    //=== OPERATOR OVERLOADS =================================================================================
    //========================================================================================================

    //--------------------------------------------------------------------------------------------------------
    //! \brief  Delete copy assignment operator.
    //--------------------------------------------------------------------------------------------------------
    SweepHessSolve< OrdinateFlux, SIMD_length > &
        operator=( const SweepHessSolve< OrdinateFlux, SIMD_length > & ) = delete;


    //========================================================================================================
    //=== INTERFACE OVERRIDES ================================================================================
    //========================================================================================================

    //--------------------------------------------------------------------------------------------------------
    //! \brief  Returns the string descriptor of the class type.
    //--------------------------------------------------------------------------------------------------------
    static const std::string Descriptor( void );

    //--------------------------------------------------------------------------------------------------------
    //! \brief  Returns the string descriptor for an objects type.
    //--------------------------------------------------------------------------------------------------------
    const std::string GetDescriptor( void ) const override;

    //--------------------------------------------------------------------------------------------------------
    //! \brief  Prints a summary of the object's configuration to the logging interface.
    //--------------------------------------------------------------------------------------------------------
    void Print( const std::string & = "  " ) const override;

protected:

    //--------------------------------------------------------------------------------------------------------
    //! \brief  Computes the elements of the matrix/matrices for the linear system(s) to be solved at the
    //!         given step of a sweep.
    //--------------------------------------------------------------------------------------------------------
    void CalcA (
        double * const A,
        const SIMD_BlkIdx<0> & idx,
        const RKDG::CrossSection & sigma,
        void * const work_ptr
    ) const override;

    //--------------------------------------------------------------------------------------------------------
    //! \brief  Solves the system for the updated coefficients.
    //--------------------------------------------------------------------------------------------------------
    void SolveSystem (
        double * const A,
        double * const B,
        void * const work_ptr
    ) const override;


    //========================================================================================================
    //=== PROTECTED HELPER FUNCTIONS =========================================================================
    //========================================================================================================

    //--------------------------------------------------------------------------------------------------------
    //! \brief  Zeros the arrays for storing precomputed values (e.g., angular component of matrices).
    //--------------------------------------------------------------------------------------------------------
    void ZeroMatrices( void ) const override;

    //--------------------------------------------------------------------------------------------------------
    //! \brief  Computes Hessenberg reductions of the angular components of the matrices that define the
    //!         linear systems.
    //--------------------------------------------------------------------------------------------------------
    void ComputeMatrices( void ) const override;


    //========================================================================================================
    //=== DIMENSION-INDEPENDENT INDEXING FUNCTIONS FOR LOCAL SYSTEMS =========================================
    //========================================================================================================

# if defined (STRICT_CHECK) || defined (DOXYCOMPILE)

    //--------------------------------------------------------------------------------------------------------
    //! \brief  Performs bounds checking on indices.
    //--------------------------------------------------------------------------------------------------------
    template< int64_t len >
    inline void CheckBoundsDI (
        const int64_t i,
        const int64_t j,
        const int64_t l,
        const char * const func
    ) const;

# endif // if defined (STRICT_CHECK)

    //--------------------------------------------------------------------------------------------------------
    //! \brief  Column-major dimension-independent indexing function for matrices used in the sweep
    //!         algorithms.
    //--------------------------------------------------------------------------------------------------------
    template< int64_t len >
    inline size_t MatIdxDI_CM (
        const int64_t i,
        const int64_t j,
        const int64_t l = 0
    ) const;

    //--------------------------------------------------------------------------------------------------------
    //! \brief  Dimension-independent indexing function for vectors used in the sweep algorithms.
    //--------------------------------------------------------------------------------------------------------
    template< int64_t len >
    inline size_t VecIdxDI (
        const int64_t i,
        const int64_t l = 0
    ) const;


    //========================================================================================================
    //=== INTERNAL VARIABLES =================================================================================
    //========================================================================================================

    //!
    //! \brief  Pointer for storing scalar components of Householder reflectors for Hessenberg reductions of
    //!         angular components of matrices.
    //!
    double ** tau = nullptr;

    int64_t dimof_tau;      //!< Number of tau matrices that are stored.

};


//============================================================================================================
//=== DEFINITIONS OF INLINE FUNCTIONS ========================================================================
//============================================================================================================

# if defined (STRICT_CHECK) || defined (DOXYCOMPILE)

//------------------------------------------------------------------------------------------------------------
//! \brief  Performs bounds checking on indices.
//!
//! An exception of type \c std::out_of_range is thrown if an out-of-bounds index is detected.
//!
//! This function is intended to be called as:
//! \code
//!     CheckBoundsDI(i,j,l,__func__)
//! \endcode
//!
//! \tparam         len     Number of stacked systems for indexing function. Note that this may be different
//!                         from SIMD_length.
//! \param[in]      i       Row of matrix.
//! \param[in]      j       Column of matrix.
//! \param[in]      l       Index for system in stack. In [ 0, len ).
//! \param[in]      func    Calling function.
//!
//! \see    MatIdxDI_CM()
//! \see    VecIdxDI()
//------------------------------------------------------------------------------------------------------------
template< class OrdinateFlux, int64_t SIMD_length >
template< int64_t len >
inline void SweepHessSolve< OrdinateFlux, SIMD_length >
::CheckBoundsDI (

    const int64_t i,
    const int64_t j,
    const int64_t l,
    const char * const func

) const {

    if ( i < 0 || i >= this->system_dim ) {

        std::string error_message
            = "Value " + std::to_string(i) + " for row index i in " + std::string(func)
                + " outsize permissible range [ 0, " + std::to_string( this->system_dim ) + " ).\n";

        PRINT_ERROR( error_message.c_str() )
        throw std::out_of_range( error_message );
    }

    if ( j < 0 || j >= this->system_dim ) {

        std::string error_message
            = "Value " + std::to_string(j) + " for column index j in " + std::string(func)
                + " outsize permissible range [ 0, " + std::to_string( this->system_dim ) + " ).\n";

        PRINT_ERROR( error_message.c_str() )
        throw std::out_of_range( error_message );
    }

    if ( l < 0 || l >= len ) {

        std::string error_message
            = "Value " + std::to_string(l) + " for SIMD index l in " + std::string(func)
                + " outsize permissible range [ 0, " + std::to_string( len ) + " ).\n";

        PRINT_ERROR( error_message.c_str() )
        throw std::out_of_range( error_message );
    }
}

# endif // if defined (STRICT_CHECK)


//------------------------------------------------------------------------------------------------------------
//! \brief  Column-major dimension-independent indexing function for matrices used in the sweep algorithms.
//!
//! \tparam         len     Number of stacked systems for indexing function. Note that this may be different
//!                         from SIMD_length.
//! \param[in]      i       Row of matrix.
//! \param[in]      j       Column of matrix.
//! \param[in]      l       Index for system in stack. In [ 0, len ).
//!
//! \see    MatIdxDI_RM()
//! \see    VecIdxDI()
//------------------------------------------------------------------------------------------------------------
template< class OrdinateFlux, int64_t SIMD_length >
template< int64_t len >
inline size_t SweepHessSolve< OrdinateFlux, SIMD_length >
::MatIdxDI_CM (

    const int64_t i,
    const int64_t j,
    const int64_t l // = 0

) const {

# if defined (STRICT_CHECK)

    CheckBoundsDI<len>(i,j,l, __func__ );

# endif // if defined (STRICT_CHECK)

    return (
        l + (len)*
        (
            i + (this->system_dim)*
            (
                j
            )
        )
    );
}


//------------------------------------------------------------------------------------------------------------
//! \brief  Dimension-independent indexing function for vectors used in the sweep algorithms.
//!
//! \tparam         len     Number of stacked systems for indexing function. Note that this may be different
//!                         from SIMD_length.
//! \param[in]      i       Index of element of vector.
//! \param[in]      l       Index for system in stack. In [ 0, len ).
//!
//! \see    VecIdxDI()
//------------------------------------------------------------------------------------------------------------
template< class OrdinateFlux, int64_t SIMD_length >
template< int64_t len >
inline size_t SweepHessSolve< OrdinateFlux, SIMD_length >
::VecIdxDI (

    const int64_t i,
    const int64_t l // = 0

) const {

# if defined (STRICT_CHECK)

    CheckBoundsDI<len>(i,i,l, __func__ );

# endif // if defined (STRICT_CHECK)

    return (
        l + (len)*
        (
            i
        )
    );
}


# endif // ifndef __ABSTRACT__SWEEP_HESS_SOLVE_HPP__
