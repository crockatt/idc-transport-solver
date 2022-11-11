//------------------------------------------------------------------------------------------------------------
//! \file   objects/STDG/DensityFunction.hpp
//! \brief  Header for STDG::DensityFunction class.
//!
//! \author Michael M. Crockatt
//! \date   February 2018
//------------------------------------------------------------------------------------------------------------

# ifndef __STDG__DENSITY_FUNCTION_HPP__
# define __STDG__DENSITY_FUNCTION_HPP__


# if defined (ENABLE_AZTEC_EPETRA)
    # include <Epetra_MultiVector.h>
# endif

# include "utils/CLog.hpp"
# include "utils/global.hpp"
# include "objects/Abstract/DensityFunction.hpp"

# if defined (ENABLE_BELOS_TPETRA)
    # include "utils/BelosTpetraInterface.hpp"
# endif



namespace STDG {


//------------------------------------------------------------------------------------------------------------
//! \brief  Class used for storing spatially and time-varying scalar-valued density functions with a
//!         space-time discontinuous Galerkin approximation.
//!
//! Density functions use meshes extruded in time and are represented as a space-time slab with a width of one
//! cell in the temporal dimension.
//!
//!
//! \attention  Virtual destructor inherited from Abstract::DensityFunction class.
//!
//!
//------------------------------------------------------------------------------------------------------------
class DensityFunction : public virtual Abstract::DensityFunction {

public:

    //!
    //! Maximum degree of DG basis polynomials in spatial dimensions.
    //!
    int64_t DG_degree_x;

    //!
    //! Maximum degree of DG basis polynomials in temporal dimension.
    //!
    int64_t DG_degree_t;


    //========================================================================================================
    //=== CONSTRUCTORS, DESTRUCTOR, AND RECONFIGURATION ROUTINES =============================================
    //========================================================================================================

    //--------------------------------------------------------------------------------------------------------
    //! \brief  Constructs an empty STDG::DensityFunction object with zero-initialized parameters.
    //--------------------------------------------------------------------------------------------------------
    DensityFunction( void );

    //--------------------------------------------------------------------------------------------------------
    //! \brief  Delete copy constructor.
    //--------------------------------------------------------------------------------------------------------
    DensityFunction( const DensityFunction & ) = delete;

    //--------------------------------------------------------------------------------------------------------
    //! \brief  Constructs an empty RKDG::DensityFunction object with zero-initialized parameters.
    //--------------------------------------------------------------------------------------------------------
    DensityFunction(
        const DomainDecomposition &,
        const int64_t
    );

    //--------------------------------------------------------------------------------------------------------
    //! \brief  Constructs an empty RKDG::DensityFunction object with zero-initialized parameters.
    //--------------------------------------------------------------------------------------------------------
    DensityFunction(
        const DomainDecomposition &,
        const int64_t,
        const int64_t
    );

    //
    // Acknowledge use of Reconfigure routines from base class.
    //
    using Abstract::DensityFunction::Reconfigure;

    //--------------------------------------------------------------------------------------------------------
    //! \brief  Reconfigures an STDG::DensityFunction object with a different set of parameters.
    //--------------------------------------------------------------------------------------------------------
    virtual DensityFunction & Reconfigure( const int64_t );

    //--------------------------------------------------------------------------------------------------------
    //! \brief  Reconfigures an STDG::DensityFunction object with a different set of parameters.
    //--------------------------------------------------------------------------------------------------------
    virtual DensityFunction & Reconfigure( const int64_t, const int64_t );

    //--------------------------------------------------------------------------------------------------------
    //! \brief  Reconfigures an STDG::DensityFunction object with a different set of parameters.
    //--------------------------------------------------------------------------------------------------------
    virtual DensityFunction & Reconfigure(
        const DomainDecomposition &,
        const int64_t
    );

    //--------------------------------------------------------------------------------------------------------
    //! \brief  Reconfigures an STDG::DensityFunction object with a different set of parameters.
    //--------------------------------------------------------------------------------------------------------
    virtual DensityFunction & Reconfigure(
        const DomainDecomposition &,
        const int64_t,
        const int64_t
    );


    //========================================================================================================
    //=== INDEXING MEMBER FUNCTIONS ==========================================================================
    //========================================================================================================

# if defined (STRICT_CHECK) || defined (DOXYCOMPILE)

    //--------------------------------------------------------------------------------------------------------
    //! \brief  Performs bounds checking on indices.
    //--------------------------------------------------------------------------------------------------------
    inline void CheckBounds (
            const int64_t i
    # if SPACE_DIMS >= 2 || defined (DOXYCOMPILE)
        ,   const int64_t j
    # endif
    # if SPACE_DIMS == 3
        ,   const int64_t k
    # endif
        ,   const int64_t d
    # if SPACE_DIMS >= 2 || defined (DOXYCOMPILE)
        ,   const int64_t e
    # endif
    # if SPACE_DIMS == 3
        ,   const int64_t f
    # endif
        ,   const int64_t s
        ,   const char * const func
    ) const;

# endif // if defined (STRICT_CHECK)

    //--------------------------------------------------------------------------------------------------------
    //! \brief  Returns the one-dimensional index for the specified DG coefficient in the array
    //!         DensityFunction::density.
    //--------------------------------------------------------------------------------------------------------
    inline size_t Index (
            const int64_t i
    # if SPACE_DIMS >= 2 || defined (DOXYCOMPILE)
        ,   const int64_t j
    # endif
    # if SPACE_DIMS == 3
        ,   const int64_t k
    # endif
        ,   const int64_t d
    # if SPACE_DIMS >= 2 || defined (DOXYCOMPILE)
        ,   const int64_t e
    # endif
    # if SPACE_DIMS == 3
        ,   const int64_t f
    # endif
        ,   const int64_t s
    ) const;

    //--------------------------------------------------------------------------------------------------------
    //! \brief  Returns a reference to the specified DG coefficient inside the array DensityFunction::density.
    //--------------------------------------------------------------------------------------------------------
    inline double & operator() (
            const int64_t i
    # if SPACE_DIMS >= 2 || defined (DOXYCOMPILE)
        ,   const int64_t j
    # endif
    # if SPACE_DIMS == 3
        ,   const int64_t k
    # endif
        ,   const int64_t d
    # if SPACE_DIMS >= 2 || defined (DOXYCOMPILE)
        ,   const int64_t e
    # endif
    # if SPACE_DIMS == 3
        ,   const int64_t f
    # endif
        ,   const int64_t s
    ) const;


    //========================================================================================================
    //=== OPERATOR OVERLOADS =================================================================================
    //========================================================================================================

    //--------------------------------------------------------------------------------------------------------
    //! \brief  Delete copy assignment operator.
    //--------------------------------------------------------------------------------------------------------
    DensityFunction & operator=( const DensityFunction & ) = delete;


    //========================================================================================================
    //=== STATIC MEMBER FUNCTIONS ============================================================================
    //========================================================================================================

    //--------------------------------------------------------------------------------------------------------
    //! \brief  Creates and returns a DensityFunction object.
    //--------------------------------------------------------------------------------------------------------
    static DensityFunction * Create( const DomainDecomposition &, const ParameterList & );


    //========================================================================================================
    //=== ADDITIONAL PUBLIC MEMBER FUNCTIONS =================================================================
    //========================================================================================================

    //--------------------------------------------------------------------------------------------------------
    //! \brief  Prints a summary of the object configuration to the logging interface.
    //--------------------------------------------------------------------------------------------------------
    void Print( const std::string & = "  " ) const override;

    //--------------------------------------------------------------------------------------------------------
    //! \brief  Deallocates memory at all internal pointers and zeros all parameters of an
    //!         STDG::DensityFunction object.
    //--------------------------------------------------------------------------------------------------------
    DensityFunction & Zero( void ) override;

    //--------------------------------------------------------------------------------------------------------
    //! \brief  Returns the stride in the global array between subsequent spatial cells, moving in the
    //!         specified dimension.
    //--------------------------------------------------------------------------------------------------------
    int64_t CellStride( const int64_t ) const override;

    //--------------------------------------------------------------------------------------------------------
    //! \brief  Returns a pointer to the location within DensityFunction::density where the coefficients for
    //!         the given spatial cell begin.
    //--------------------------------------------------------------------------------------------------------
    double * PointerAtCell(
            const int64_t i
    # if SPACE_DIMS >= 2
        ,   const int64_t j
    # endif
    # if SPACE_DIMS == 3
        ,   const int64_t k
    # endif
    ) const override;


protected:

    //========================================================================================================
    //=== PROTECTED MEMBER FUNCTIONS =========================================================================
    //========================================================================================================

    //--------------------------------------------------------------------------------------------------------
    //! \brief  Sets appropriate values for object parameters.
    //--------------------------------------------------------------------------------------------------------
    void SetDensityDimensions( void ) override;

};


//============================================================================================================
//=== DEFINITIONS OF INLINED INDEXING FUNCTIONS ==============================================================
//============================================================================================================

# if defined (STRICT_CHECK) || defined (DOXYCOMPILE)

//------------------------------------------------------------------------------------------------------------
//! \brief  Performs bounds checking on indices.
//!
//! An exception of type \c std::out_of_range is thrown if an out-of-bounds index is detected.
//!
//! This function is intended to be called as:
//! \code
//!     CheckBounds(i,j,d,e,s,__func__)
//! \endcode
//!
//! \param[in]      i           Spatial cell index (including ghost cells) in the \f$ x_1 \f$ dimension.
//! \param[in]      j           Spatial cell index (including ghost cells) in the \f$ x_2 \f$ dimension.
//! \param[in]      d           Degree of spatial basis polynomial wrt. \f$ x_1 \f$.
//! \param[in]      e           Degree of spatial basis polynomial wrt. \f$ x_2 \f$.
//! \param[in]      s           Degree of temporal basis polynomial.
//! \param[in]      func        Calling function.
//!
//------------------------------------------------------------------------------------------------------------
inline void DensityFunction::CheckBounds (

        const int64_t i
# if SPACE_DIMS >= 2 || defined (DOXYCOMPILE)
    ,   const int64_t j
# endif
# if SPACE_DIMS == 3
    ,   const int64_t k
# endif

    ,   const int64_t d
# if SPACE_DIMS >= 2 || defined (DOXYCOMPILE)
    ,   const int64_t e
# endif
# if SPACE_DIMS == 3
    ,   const int64_t f
# endif

    ,   const int64_t s
    ,   const char * const func

) const {

    if ( i < 0 || i >= this->nx(0) +2 ) {

        std::string error_message
            = "Value " + std::to_string(i) + " for index i in " + std::string(func)
                + " outsize permissible range [ 0, " + std::to_string( this->nx(0) +2 ) + " ).\n";

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

# if SPACE_DIMS >= 2

    if ( j < 0 || j >= this->nx(1) +2 ) {

        std::string error_message
            = "Value " + std::to_string(j) + " for index j in " + std::string(func)
                + " outsize permissible range [ 0, " + std::to_string( this->nx(1) +2 ) + " ).\n";

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

    if ( k < 0 || k >= this->nx(2) +2 ) {

        std::string error_message
            = "Value " + std::to_string(k) + " for index k in " + std::string(func)
                + " outsize permissible range [ 0, " + std::to_string( this->nx(2) +2 ) + " ).\n";

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
//! \brief  Returns the one-dimensional index for the specified DG coefficient in the array
//!         DensityFunction::density.
//!
//! \param[in]      i           Spatial cell index (including ghost cells) in the \f$ x_1 \f$ dimension.
//! \param[in]      j           Spatial cell index (including ghost cells) in the \f$ x_2 \f$ dimension.
//! \param[in]      d           Degree of spatial basis polynomial wrt. \f$ x_1 \f$.
//! \param[in]      e           Degree of spatial basis polynomial wrt. \f$ x_2 \f$.
//! \param[in]      s           Degree of temporal basis polynomial.
//!
//! \return     Appropriate index in DensityFunction::density array for specified coefficient.
//!
//! \see    STDG::DensityFunction::operator()()
//------------------------------------------------------------------------------------------------------------
inline size_t DensityFunction::Index (

        const int64_t i
# if SPACE_DIMS >= 2 || defined (DOXYCOMPILE)
    ,   const int64_t j
# endif
# if SPACE_DIMS == 3
    ,   const int64_t k
# endif

    ,   const int64_t d
# if SPACE_DIMS >= 2 || defined (DOXYCOMPILE)
    ,   const int64_t e
# endif
# if SPACE_DIMS == 3
    ,   const int64_t f
# endif

    ,   const int64_t s

) const {

# if defined (STRICT_CHECK)

    CheckBounds(
            # if SPACE_DIMS == 1
                i,d,s
            # elif SPACE_DIMS == 2
                i,j,d,e,s
            # elif SPACE_DIMS == 3
                i,j,k,d,e,f,s
            # endif
                , __func__ );

# endif // if defined (STRICT_CHECK)

    return (
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
                    d + (this->DG_degree_x + 1)*
                    (
                    # if SPACE_DIMS == 3
                        k + (this->nx(2) + 2)*
                    # endif
                        (
                        # if SPACE_DIMS >= 2
                            j + (this->nx(1) + 2)*
                        # endif
                            (
                                i
                            )
                        )
                    )
                )
            )
        )
    );
}


//------------------------------------------------------------------------------------------------------------
//! \brief  Returns a reference to the specified DG coefficient inside the array DensityFunction::density.
//!
//! \param[in]      i           Spatial cell index (including ghost cells) in the \f$ x_1 \f$ dimension.
//! \param[in]      j           Spatial cell index (including ghost cells) in the \f$ x_2 \f$ dimension.
//! \param[in]      d           Degree of spatial basis polynomial wrt. \f$ x_1 \f$.
//! \param[in]      e           Degree of spatial basis polynomial wrt. \f$ x_2 \f$.
//! \param[in]      s           Degree of temporal basis polynomial.
//!
//! \return     Reference to the specified coefficient.
//!
//! \see    STDG::DensityFunction::Index()
//------------------------------------------------------------------------------------------------------------
inline double & DensityFunction::operator() (

        const int64_t i
# if SPACE_DIMS >= 2 || defined (DOXYCOMPILE)
    ,   const int64_t j
# endif
# if SPACE_DIMS == 3
    ,   const int64_t k
# endif

    ,   const int64_t d
# if SPACE_DIMS >= 2 || defined (DOXYCOMPILE)
    ,   const int64_t e
# endif
# if SPACE_DIMS == 3
    ,   const int64_t f
# endif

    ,   const int64_t s

) const {

# if defined (STRICT_CHECK)

    if ( this->density == nullptr ) {

        std::string error_message = "Pointer this->density is NULL in '" + std::string(__func__) + "'.\n";

        PRINT_ERROR( error_message.c_str() )
        throw std::runtime_error( error_message );
    }

# endif // if defined (STRICT_CHECK)

    return this->density[ Index
                            # if SPACE_DIMS == 1
                                (i,d,s)
                            # elif SPACE_DIMS == 2
                                (i,j,d,e,s)
                            # elif SPACE_DIMS == 3
                                (i,j,k,d,e,f,s)
                            # endif
                        ];
}


} // namespace STDG


# endif // ifndef __STDG__DENSITY_FUNCTION_HPP__
