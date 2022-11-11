//------------------------------------------------------------------------------------------------------------
//! \file   linear_solvers/STDG/SweepCellSolve/SweepCellSolve.hpp
//! \brief  Header file containing declaration of specialization of SweepCellSolve for STDG::OrdinateFlux
//!         objects.
//!
//! \author Michael M. Crockatt
//! \date   February 2018
//------------------------------------------------------------------------------------------------------------

# ifndef __STDG__SWEEP_CELL_SOLVE_HPP__
# define __STDG__SWEEP_CELL_SOLVE_HPP__


# include "linear_solvers/Abstract/SweepCellSolve/SweepCellSolve.hpp"
# include "objects/STDG/OrdinateFlux.hpp"


//------------------------------------------------------------------------------------------------------------
//! \brief  Declaration of SweepCellSolve<STDG::OrdinateFlux> specialization.
//!
//! This specialization implements:
//!     - Indexing functions for local systems.
//!     - SweepCellSolveBase::StoreResult.
//!
//------------------------------------------------------------------------------------------------------------
template< int64_t SIMD_length >
class SweepCellSolve< STDG::OrdinateFlux, SIMD_length > :
    public SweepCellSolveBase<STDG::OrdinateFlux>
{

public:

    //========================================================================================================
    //=== CONSTRUCTORS AND DESTRUCTOR ========================================================================
    //========================================================================================================

    //--------------------------------------------------------------------------------------------------------
    //! \brief  Delete default constructor.
    //--------------------------------------------------------------------------------------------------------
    SweepCellSolve( void ) = delete;

    //--------------------------------------------------------------------------------------------------------
    //! \brief  Delete copy constructor.
    //--------------------------------------------------------------------------------------------------------
    SweepCellSolve( const SweepCellSolve< STDG::OrdinateFlux, SIMD_length > & ) = delete;

    //--------------------------------------------------------------------------------------------------------
    //! \brief  Construct solver object from given parameters.
    //--------------------------------------------------------------------------------------------------------
    SweepCellSolve (
        const SweepOperator<STDG::OrdinateFlux> & enclosing,
        const ParameterList & input_list,
        const double dt_in = std::numeric_limits<double>::infinity()
    );

    //--------------------------------------------------------------------------------------------------------
    //! \brief  Virtual destructor.
    //--------------------------------------------------------------------------------------------------------
    virtual ~SweepCellSolve( void );


    //========================================================================================================
    //=== OPERATOR OVERLOADS =================================================================================
    //========================================================================================================

    //--------------------------------------------------------------------------------------------------------
    //! \brief  Delete copy assignment operator.
    //--------------------------------------------------------------------------------------------------------
    SweepCellSolve< STDG::OrdinateFlux, SIMD_length > &
        operator=( const SweepCellSolve< STDG::OrdinateFlux, SIMD_length > & ) = delete;


    //========================================================================================================
    //=== OVERRIDES OF INTERFACE FUNCTIONS ===================================================================
    //========================================================================================================

protected:

    //--------------------------------------------------------------------------------------------------------
    //! \brief  Stores the coefficients obtained from the solve into the output object.
    //--------------------------------------------------------------------------------------------------------
    void StoreResult (
        const double * const B,
        const SIMD_BlkIdx<0> & idx,
        const STDG::OrdinateFlux & result
    ) const override;


    //========================================================================================================
    //=== INDEXING FUNCTIONS FOR LOCAL SYSTEMS ===============================================================
    //========================================================================================================

# if defined (STRICT_CHECK) || defined (DOXYCOMPILE)

    //--------------------------------------------------------------------------------------------------------
    //! \brief  Performs bounds checking on indices.
    //--------------------------------------------------------------------------------------------------------
    template< int64_t len >
    inline void CheckBounds (
            const int64_t a
    # if SPACE_DIMS >= 2 || defined (DOXYCOMPILE)
        ,   const int64_t b
    # endif
    # if SPACE_DIMS == 3
        ,   const int64_t c
    # endif
        ,   const int64_t r
        ,   const int64_t d
    # if SPACE_DIMS >= 2 || defined (DOXYCOMPILE)
        ,   const int64_t e
    # endif
    # if SPACE_DIMS == 3
        ,   const int64_t f
    # endif
        ,   const int64_t s
        ,   const int64_t l
        ,   const char * const func
    ) const;

# endif // if defined (STRICT_CHECK)

    //--------------------------------------------------------------------------------------------------------
    //! \brief  Column-major indexing function for matrices used in the STDG sweep algorithms.
    //--------------------------------------------------------------------------------------------------------
    template< int64_t len >
    inline size_t MatIdx_CM (
            const int64_t a
    # if SPACE_DIMS >= 2 || defined (DOXYCOMPILE)
        ,   const int64_t b
    # endif
    # if SPACE_DIMS == 3
        ,   const int64_t c
    # endif
        ,   const int64_t r
        ,   const int64_t d
    # if SPACE_DIMS >= 2 || defined (DOXYCOMPILE)
        ,   const int64_t e
    # endif
    # if SPACE_DIMS == 3
        ,   const int64_t f
    # endif
        ,   const int64_t s
        ,   const int64_t l = 0
    ) const;

    //--------------------------------------------------------------------------------------------------------
    //! \brief  Indexing function for vectors used in the STDG sweep algorithms.
    //--------------------------------------------------------------------------------------------------------
    template< int64_t len >
    inline size_t VecIdx (
            const int64_t a
    # if SPACE_DIMS >= 2 || defined (DOXYCOMPILE)
        ,   const int64_t b
    # endif
    # if SPACE_DIMS == 3
        ,   const int64_t c
    # endif
        ,   const int64_t r
        ,   const int64_t l = 0
    ) const;


    //========================================================================================================
    //=== INTERNAL VARIABLES =================================================================================
    //========================================================================================================

    const int64_t DG_degree_x;      //! Degree of DG spatial discretization.
    const int64_t DG_degree_t;      //! Degree of DG temporal discretization.

};


//============================================================================================================
//=== DEFINITIONS OF INLINED INDEXING FUNCTIONS FOR LOCAL SYSTEMS ============================================
//============================================================================================================

# if defined (STRICT_CHECK) || defined (DOXYCOMPILE)

//------------------------------------------------------------------------------------------------------------
//! \brief  Performs bounds checking on indices.
//!
//! An exception of type \c std::out_of_range is thrown if an out-of-bounds index is detected.
//!
//! This function is intended to be called as:
//! \code
//!     CheckBounds(a,b,d,e,__func__)
//! \endcode
//!
//! \param[in]      a           Degree of test function with respect to \f$ x_1 \f$.
//! \param[in]      b           Degree of test function with respect to \f$ x_2 \f$.
//! \param[in]      r           Degree of test function with respect to time.
//! \param[in]      d           Degree of basis function with respect to \f$ x_1 \f$.
//! \param[in]      e           Degree of basis function with respect to \f$ x_2 \f$.
//! \param[in]      s           Degree of basis function with respect to time.
//! \param[in]      l           Index for system in SIMD block. In [ 0, len ).
//! \param[in]      func        Calling function.
//!
//! \see    MatIdx_RM()
//! \see    VecIdx()
//------------------------------------------------------------------------------------------------------------
template< int64_t SIMD_length >
template< int64_t len >
inline void SweepCellSolve< STDG::OrdinateFlux, SIMD_length >::CheckBounds (

        const int64_t a
# if SPACE_DIMS >= 2 || defined (DOXYCOMPILE)
    ,   const int64_t b
# endif
# if SPACE_DIMS == 3
    ,   const int64_t c
# endif

    ,   const int64_t r

    ,   const int64_t d
# if SPACE_DIMS >= 2 || defined (DOXYCOMPILE)
    ,   const int64_t e
# endif
# if SPACE_DIMS == 3
    ,   const int64_t f
# endif

    ,   const int64_t s
    ,   const int64_t l
    ,   const char * const func

) const {

    if ( a < 0 || a > this->DG_degree_x ) {

        std::string error_message
            = "Value " + std::to_string(a) + " for degree a in " + std::string(func)
                + " outsize permissible range [ 0, " + std::to_string( this->DG_degree_x ) + " ].\n";

        PRINT_ERROR( error_message.c_str() )
        throw std::out_of_range( error_message );
    }

    if ( r < 0 || r > this->DG_degree_t ) {

        std::string error_message
            = "Value " + std::to_string(r) + " for degree r in " + std::string(func)
                + " outsize permissible range [ 0, " + std::to_string( this->DG_degree_t ) + " ].\n";

        PRINT_ERROR( error_message.c_str() )
        throw std::out_of_range( error_message );
    }

    if ( d < 0 || d > this->DG_degree_x ) {

        std::string error_message
            = "Value " + std::to_string(d) + " for degree d in " + std::string(func)
                + " outsize permissible range [ 0, " + std::to_string( this->DG_degree_x ) + " ].\n";

        PRINT_ERROR( error_message.c_str() )
        throw std::out_of_range( error_message );
    }

    if ( s < 0 || s > this->DG_degree_t ) {

        std::string error_message
            = "Value " + std::to_string(s) + " for degree s in " + std::string(func)
                + " outsize permissible range [ 0, " + std::to_string( this->DG_degree_t ) + " ].\n";

        PRINT_ERROR( error_message.c_str() )
        throw std::out_of_range( error_message );
    }

    if ( l < 0 || l >= len ) {

        std::string error_message
            = "Value " + std::to_string(l) + " for SIMD index l in " + std::string(func)
                + " outsize permissible range [ 0, " + std::to_string( len ) + " ].\n";

        PRINT_ERROR( error_message.c_str() )
        throw std::out_of_range( error_message );
    }

# if SPACE_DIMS >= 2

    if ( b < 0 || b > this->DG_degree_x ) {

        std::string error_message
            = "Value " + std::to_string(b) + " for degree b in " + std::string(func)
                + " outsize permissible range [ 0, " + std::to_string( this->DG_degree_x ) + " ].\n";

        PRINT_ERROR( error_message.c_str() )
        throw std::out_of_range( error_message );
    }

    if ( e < 0 || e > this->DG_degree_x ) {

        std::string error_message
            = "Value " + std::to_string(e) + " for degree e in " + std::string(func)
                + " outsize permissible range [ 0, " + std::to_string( this->DG_degree_x ) + " ].\n";

        PRINT_ERROR( error_message.c_str() )
        throw std::out_of_range( error_message );
    }

# endif // if SPACE_DIMS >= 2
# if SPACE_DIMS == 3

    if ( c < 0 || c > this->DG_degree_x ) {

        std::string error_message
            = "Value " + std::to_string(c) + " for degree c in " + std::string(func)
                + " outsize permissible range [ 0, " + std::to_string( this->DG_degree_x ) + " ].\n";

        PRINT_ERROR( error_message.c_str() )
        throw std::out_of_range( error_message );
    }

    if ( f < 0 || f > this->DG_degree_x ) {

        std::string error_message
            = "Value " + std::to_string(f) + " for degree f in " + std::string(func)
                + " outsize permissible range [ 0, " + std::to_string( this->DG_degree_x ) + " ].\n";

        PRINT_ERROR( error_message.c_str() )
        throw std::out_of_range( error_message );
    }

# endif // if SPACE_DIMS == 3
}

# endif // if defined (STRICT_CHECK)


//------------------------------------------------------------------------------------------------------------
//! \brief  Column-major indexing function for matrices used in the STDG sweep algorithms.
//!
//! \param[in]      a           Degree of test function with respect to \f$ x_1 \f$.
//! \param[in]      b           Degree of test function with respect to \f$ x_2 \f$.
//! \param[in]      r           Degree of test function with respect to time.
//! \param[in]      d           Degree of basis function with respect to \f$ x_1 \f$.
//! \param[in]      e           Degree of basis function with respect to \f$ x_2 \f$.
//! \param[in]      s           Degree of basis function with respect to time.
//! \param[in]      l           Index for system in SIMD block. In [ 0, len ).
//!
//! \see    VecIdx()
//------------------------------------------------------------------------------------------------------------
template< int64_t SIMD_length >
template< int64_t len >
inline size_t SweepCellSolve< STDG::OrdinateFlux, SIMD_length >::MatIdx_CM (

        const int64_t a
# if SPACE_DIMS >= 2 || defined (DOXYCOMPILE)
    ,   const int64_t b
# endif
# if SPACE_DIMS == 3
    ,   const int64_t c
# endif

    ,   const int64_t r

    ,   const int64_t d
# if SPACE_DIMS >= 2 || defined (DOXYCOMPILE)
    ,   const int64_t e
# endif
# if SPACE_DIMS == 3
    ,   const int64_t f
# endif

    ,   const int64_t s
    ,   const int64_t l // = 0

) const {

# if defined (STRICT_CHECK)

    CheckBounds<len>(
            # if SPACE_DIMS == 1
                a,r,d,s
            # elif SPACE_DIMS == 2
                a,b,r,d,e,s
            # elif SPACE_DIMS == 3
                a,b,c,r,d,e,f,s
            # endif
                , l
                , __func__ );

# endif // if defined (STRICT_CHECK)

    return (
        l + (len)*
        (
            r + (this->DG_degree_t + 1)*
            (
            # if SPACE_DIMS == 3
                c + (this->DG_degree_x + 1)*
            # endif
                (
                # if SPACE_DIMS >= 2
                    b + (this->DG_degree_x + 1)*
                # endif
                    (
                        a + (this->DG_degree_x + 1)*
                        (
                            s + (this->DG_degree_t + 1)*
                            (
                            # if SPACE_DIMS == 3
                                f + (this->DG_degree_x + 1)*
                            # endif
                                (
                                # if SPACE_DIMS >= 2
                                    e + (this->DG_degree_x + 1)*
                                # endif
                                    (
                                        d
                                    )
                                )
                            )
                        )
                    )
                )
            )
        )
    );
}


//------------------------------------------------------------------------------------------------------------
//! \brief  Indexing function for vectors used in the STDG sweep algorithms.
//!
//! \param[in]      a           Degree of test function with respect to \f$ x_1 \f$.
//! \param[in]      b           Degree of test function with respect to \f$ x_2 \f$.
//! \param[in]      r           Degree of test function with respect to time.
//! \param[in]      l           Index for system in SIMD block. In [ 0, len ).
//!
//! \see    MatIdx_CM()
//------------------------------------------------------------------------------------------------------------
template< int64_t SIMD_length >
template< int64_t len >
inline size_t SweepCellSolve< STDG::OrdinateFlux, SIMD_length >::VecIdx (

        const int64_t a
# if SPACE_DIMS >= 2 || defined (DOXYCOMPILE)
    ,   const int64_t b
# endif
# if SPACE_DIMS == 3
    ,   const int64_t c
# endif

    ,   const int64_t r
    ,   const int64_t l // = 0

) const {

# if defined (STRICT_CHECK)

    CheckBounds<len>(
            # if SPACE_DIMS == 1
                a,r,a,r
            # elif SPACE_DIMS == 2
                a,b,r,a,b,r
            # elif SPACE_DIMS == 3
                a,b,c,r,a,b,c,r
            # endif
                , l
                , __func__ );

# endif // if defined (STRICT_CHECK)

    return (
        l + (len)*
        (
            r + (this->DG_degree_t + 1)*
            (
            # if SPACE_DIMS == 3
                c + (this->DG_degree_x + 1)*
            # endif
                (
                # if SPACE_DIMS >= 2
                    b + (this->DG_degree_x + 1)*
                # endif
                    (
                        a
                    )
                )
            )
        )
    );
}


# endif // ifndef __STDG__SWEEP_CELL_SOLVE_HPP__
