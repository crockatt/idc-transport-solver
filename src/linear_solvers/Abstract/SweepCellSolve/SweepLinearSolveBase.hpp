//------------------------------------------------------------------------------------------------------------
//! \file   linear_solvers/Abstract/SweepCellSolve/SweepLinearSolveBase.hpp
//! \brief  Header file containing declaration of SweepLinearSolveBase abstract class.
//!
//! \author Michael M. Crockatt
//! \date   March 2018
//------------------------------------------------------------------------------------------------------------

# ifndef __ABSTRACT__SWEEP_LINEAR_SOLVE_BASE_HPP__
# define __ABSTRACT__SWEEP_LINEAR_SOLVE_BASE_HPP__


# include "linear_solvers/Abstract/SweepCellSolve/SweepCellSolve.hpp"


//------------------------------------------------------------------------------------------------------------
//! \brief  Declaration of SweepLinearSolveBase abstract class.
//!
//!
//------------------------------------------------------------------------------------------------------------
template< class OrdinateFlux, int64_t SIMD_length >
class SweepLinearSolveBase :
    public SweepCellSolve< OrdinateFlux, SIMD_length >
{

public:

    //========================================================================================================
    //=== CONSTRUCTORS AND DESTRUCTOR ========================================================================
    //========================================================================================================

    //--------------------------------------------------------------------------------------------------------
    //! \brief  Delete default constructor.
    //--------------------------------------------------------------------------------------------------------
    SweepLinearSolveBase( void ) = delete;

    //--------------------------------------------------------------------------------------------------------
    //! \brief  Delete copy constructor.
    //--------------------------------------------------------------------------------------------------------
    SweepLinearSolveBase( const SweepLinearSolveBase< OrdinateFlux, SIMD_length > & ) = delete;

    //--------------------------------------------------------------------------------------------------------
    //! \brief  Construct solver object from given parameters.
    //--------------------------------------------------------------------------------------------------------
    SweepLinearSolveBase (
        const SweepOperator<OrdinateFlux> & enclosing,
        const ParameterList & input_list,
        const double dt_in = std::numeric_limits<double>::infinity()
    );

    //--------------------------------------------------------------------------------------------------------
    //! \brief  Virtual destructor.
    //--------------------------------------------------------------------------------------------------------
    virtual ~SweepLinearSolveBase( void );


    //========================================================================================================
    //=== OPERATOR OVERLOADS =================================================================================
    //========================================================================================================

    //--------------------------------------------------------------------------------------------------------
    //! \brief  Delete copy assignment operator.
    //--------------------------------------------------------------------------------------------------------
    SweepLinearSolveBase< OrdinateFlux, SIMD_length > &
        operator=( const SweepLinearSolveBase< OrdinateFlux, SIMD_length > & ) = delete;


    //========================================================================================================
    //=== OVERRIDES OF INTERFACE FUNCTIONS ===================================================================
    //========================================================================================================

    //--------------------------------------------------------------------------------------------------------
    //! \brief  Returns the dimension of the array to be allocated for A.
    //--------------------------------------------------------------------------------------------------------
    size_t GetDimofA ( void ) const override;

    //--------------------------------------------------------------------------------------------------------
    //! \brief  Returns the dimension of the array to be allocated for B.
    //--------------------------------------------------------------------------------------------------------
    size_t GetDimofB ( void ) const override;

    //--------------------------------------------------------------------------------------------------------
    //! \brief  Returns the dimension of the array to be allocated for extra work space.
    //--------------------------------------------------------------------------------------------------------
    size_t GetDimofWork ( void ) const override;


    //========================================================================================================
    //=== PROTECTED HELPER FUNCTIONS =========================================================================
    //========================================================================================================

    //--------------------------------------------------------------------------------------------------------
    //! \brief  Zeros the arrays for storing precomputed values (e.g., angular component of matrices).
    //--------------------------------------------------------------------------------------------------------
    virtual void ZeroMatrices( void ) const;

    //--------------------------------------------------------------------------------------------------------
    //! \brief  Computes the angular component of the matrices that define the linear systems.
    //--------------------------------------------------------------------------------------------------------
    virtual void ComputeMatrices( void ) const = 0;


    //========================================================================================================
    //=== INTERNAL VARIABLES =================================================================================
    //========================================================================================================

    double ** Aq;               //!< Pointer for storing angular component of matrices.

    int64_t dimof_A;            //!< Number of elements in each matrix.
    int64_t dimof_B;            //!< Number of elements in each vector.

    size_t sizeof_A;            //!< Number of bytes allocated for each matrix.
    size_t sizeof_B;            //!< Number of bytes allocated for each vector.

    bool force_zwc;             //!< Force SweepLinearSolveBase::CalcA to assume zone-wise constant cross sections.

};


# endif // ifndef __ABSTRACT__SWEEP_LINEAR_SOLVE_BASE_HPP__
