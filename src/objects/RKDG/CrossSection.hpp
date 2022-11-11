//------------------------------------------------------------------------------------------------------------
//! \file   objects/RKDG/CrossSection.hpp
//! \brief  Header for RKDG::CrossSection class.
//!
//! \author Michael M. Crockatt
//! \date   February 2018
//------------------------------------------------------------------------------------------------------------

# ifndef __RKDG__CROSS_SECTION_HPP__
# define __RKDG__CROSS_SECTION_HPP__


# include <cstddef>
# include <cstdint>

# include "utils/global.hpp"
# include "objects/RKDG/DensityFunction.hpp"


namespace RKDG {


//------------------------------------------------------------------------------------------------------------
//! \brief  Class used for storing spatially-varying density functions representing material cross sections
//!         using a discontinuous Galerkin approximation of the spatial variables.
//!
//!
//! \attention  Virtual destructor inherited from Abstract::DensityFunction class.
//!
//! \todo   Look at updating tensor indexing function to align with memory access pattern of updated sweeps.
//!
//------------------------------------------------------------------------------------------------------------
class CrossSection : public RKDG::DensityFunction {

public:

    //========================================================================================================
    //=== CONSTRUCTORS, DESTRUCTOR, AND RECONFIGURATION ROUTINES =============================================
    //========================================================================================================

    //--------------------------------------------------------------------------------------------------------
    //! \brief  Constructs an empty RKDG::CrossSection object with zero-initialized parameters.
    //--------------------------------------------------------------------------------------------------------
    CrossSection( void );

    //--------------------------------------------------------------------------------------------------------
    //! \brief  Delete copy constructor.
    //--------------------------------------------------------------------------------------------------------
    CrossSection( const CrossSection & ) = delete;

    //--------------------------------------------------------------------------------------------------------
    //! \brief  Inherit constructor RKDG::DensityFunction( const DomainDecomposition &, const int64_t ).
    //--------------------------------------------------------------------------------------------------------
    using DensityFunction::DensityFunction;

    //--------------------------------------------------------------------------------------------------------
    //! \brief  Destructor for class RKDG::CrossSection.
    //--------------------------------------------------------------------------------------------------------
    ~CrossSection( void ) override;


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
        ,   const int64_t a
    # if SPACE_DIMS >= 2 || defined (DOXYCOMPILE)
        ,   const int64_t b
    # endif
    # if SPACE_DIMS == 3
        ,   const int64_t c
    # endif
        ,   const int64_t d
    # if SPACE_DIMS >= 2 || defined (DOXYCOMPILE)
        ,   const int64_t e
    # endif
    # if SPACE_DIMS == 3
        ,   const int64_t f
    # endif
        ,   const char * const func
    ) const;

# endif // if defined (STRICT_CHECK)

    //--------------------------------------------------------------------------------------------------------
    //! \brief  Computes the one-dimensional index for the specified DG coefficient in the array
    //!         RKDG::CrossSection::tensor.
    //--------------------------------------------------------------------------------------------------------
    inline size_t TensorIndex (
            const int64_t i
    # if SPACE_DIMS >= 2 || defined (DOXYCOMPILE)
        ,   const int64_t j
    # endif
    # if SPACE_DIMS == 3
        ,   const int64_t k
    # endif
        ,   const int64_t a
    # if SPACE_DIMS >= 2 || defined (DOXYCOMPILE)
        ,   const int64_t b
    # endif
    # if SPACE_DIMS == 3
        ,   const int64_t c
    # endif
        ,   const int64_t d
    # if SPACE_DIMS >= 2 || defined (DOXYCOMPILE)
        ,   const int64_t e
    # endif
    # if SPACE_DIMS == 3
        ,   const int64_t f
    # endif
    ) const;

    //--------------------------------------------------------------------------------------------------------
    //! \brief  Returns a reference to the specified coefficient inside the array RKDG::CrossSection::tensor.
    //--------------------------------------------------------------------------------------------------------
    inline double & operator() (
            const int64_t i
    # if SPACE_DIMS >= 2 || defined (DOXYCOMPILE)
        ,   const int64_t j
    # endif
    # if SPACE_DIMS == 3
        ,   const int64_t k
    # endif
        ,   const int64_t a
    # if SPACE_DIMS >= 2 || defined (DOXYCOMPILE)
        ,   const int64_t b
    # endif
    # if SPACE_DIMS == 3
        ,   const int64_t c
    # endif
        ,   const int64_t d
    # if SPACE_DIMS >= 2 || defined (DOXYCOMPILE)
        ,   const int64_t e
    # endif
    # if SPACE_DIMS == 3
        ,   const int64_t f
    # endif
    ) const;


    //========================================================================================================
    //=== OPERATOR OVERLOADS =================================================================================
    //========================================================================================================

    //--------------------------------------------------------------------------------------------------------
    //! \brief  Delete copy assignment operator.
    //--------------------------------------------------------------------------------------------------------
    CrossSection & operator=( const CrossSection & ) = delete;

    //--------------------------------------------------------------------------------------------------------
    //! \brief  Inherit RKDG::DensityFunction::operator()().
    //--------------------------------------------------------------------------------------------------------
    using DensityFunction::operator();


    //========================================================================================================
    //=== ADDITIONAL PUBLIC MEMBER FUNCTIONS =================================================================
    //========================================================================================================

    //--------------------------------------------------------------------------------------------------------
    //! \brief  Deallocates memory at all internal pointers and zeros all parameters of an
    //!         RKDG::CrossSection object.
    //--------------------------------------------------------------------------------------------------------
    CrossSection & Zero( void ) override;

    //--------------------------------------------------------------------------------------------------------
    //! \brief  Computes the coefficients of the tensor contraction between the coefficients at
    //!         DensityFunction::density and the 3-tensor of Legendre triple product integrals, storing the
    //!         resulting 2-tensor at RKDG::CrossSection::tensor.
    //--------------------------------------------------------------------------------------------------------
    void ComputeTensor( void );


protected:

    //========================================================================================================
    //=== INTERNAL VARIABLES =================================================================================
    //========================================================================================================

    //!
    //! \brief  Pointer for storing application of TPI tensor to coefficient vector.
    //!
    double * tensor;

    //!
    //! \brief  Number of bytes allocated at the pointer RKDG::CrossSection::tensor.
    //!
    size_t sizeof_tensor;

    //!
    //! \brief  Number of doubles allocated at the pointer RKDG::CrossSection::tensor.
    //!
    int64_t dimof_tensor;


    //========================================================================================================
    //=== PROTECTED MEMBER FUNCTIONS =========================================================================
    //========================================================================================================

    //--------------------------------------------------------------------------------------------------------
    //! \brief  Sets parameters and allocates memory at internal pointers.
    //--------------------------------------------------------------------------------------------------------
    void Allocate( void ) override;

    //--------------------------------------------------------------------------------------------------------
    //! \brief  Deallocates memory at all internal pointers.
    //--------------------------------------------------------------------------------------------------------
    void Deallocate( void ) override;

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
//!     CheckBounds(i,j,a,b,d,e,__func__)
//! \endcode
//!
//! \param[in]  i   Spatial cell index (including ghost cells) in the \f$ x_1 \f$ dimension.
//! \param[in]  j   Spatial cell index (including ghost cells) in the \f$ x_2 \f$ dimension.
//! \param[in]  a   Degree of first basis polynomial in contracted tensor in the \f$ x_1 \f$ dimension.
//! \param[in]  b   Degree of first basis polynomial in contracted tensor in the \f$ x_2 \f$ dimension.
//! \param[in]  d   Degree of second basis polynomial in contracted tensor in the \f$ x_1 \f$ dimension.
//! \param[in]  e   Degree of second basis polynomial in contracted tensor in the \f$ x_2 \f$ dimension.
//! \param[in]      func        Calling function.
//!
//------------------------------------------------------------------------------------------------------------
inline void CrossSection::CheckBounds (

        const int64_t i
# if SPACE_DIMS >= 2 || defined (DOXYCOMPILE)
    ,   const int64_t j
# endif
# if SPACE_DIMS == 3
    ,   const int64_t k
# endif

    ,   const int64_t a
# if SPACE_DIMS >= 2 || defined (DOXYCOMPILE)
    ,   const int64_t b
# endif
# if SPACE_DIMS == 3
    ,   const int64_t c
# endif

    ,   const int64_t d
# if SPACE_DIMS >= 2 || defined (DOXYCOMPILE)
    ,   const int64_t e
# endif
# if SPACE_DIMS == 3
    ,   const int64_t f
# endif

    ,   const char * const func

) const {

    if ( i < 0 || i >= this->nx(0) +2 ) {

        std::string error_message
            = "Value " + std::to_string(i) + " for index i in " + std::string(func)
                + " outsize permissible range [ 0, " + std::to_string( this->nx(0) +2 ) + " ).\n";

        PRINT_ERROR( error_message.c_str() )
        throw std::out_of_range( error_message );
    }

    if ( a < 0 || a > this->DG_degree ) {

        std::string error_message
            = "Value " + std::to_string(a) + " for degree a in " + std::string(func)
                + " outsize permissible range [ 0, " + std::to_string( this->DG_degree ) + " ].\n";

        PRINT_ERROR( error_message.c_str() )
        throw std::out_of_range( error_message );
    }

    if ( d < 0 || d > this->DG_degree ) {

        std::string error_message
            = "Value " + std::to_string(d) + " for degree d in " + std::string(func)
                + " outsize permissible range [ 0, " + std::to_string( this->DG_degree ) + " ].\n";

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

    if ( b < 0 || b > this->DG_degree ) {

        std::string error_message
            = "Value " + std::to_string(b) + " for degree b in " + std::string(func)
                + " outsize permissible range [ 0, " + std::to_string( this->DG_degree ) + " ].\n";

        PRINT_ERROR( error_message.c_str() )
        throw std::out_of_range( error_message );
    }

    if ( e < 0 || e > this->DG_degree ) {

        std::string error_message
            = "Value " + std::to_string(e) + " for degree e in " + std::string(func)
                + " outsize permissible range [ 0, " + std::to_string( this->DG_degree ) + " ].\n";

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

    if ( c < 0 || c > this->DG_degree ) {

        std::string error_message
            = "Value " + std::to_string(c) + " for degree c in " + std::string(func)
                + " outsize permissible range [ 0, " + std::to_string( this->DG_degree ) + " ].\n";

        PRINT_ERROR( error_message.c_str() )
        throw std::out_of_range( error_message );
    }

    if ( f < 0 || f > this->DG_degree ) {

        std::string error_message
            = "Value " + std::to_string(f) + " for degree f in " + std::string(func)
                + " outsize permissible range [ 0, " + std::to_string( this->DG_degree ) + " ].\n";

        PRINT_ERROR( error_message.c_str() )
        throw std::out_of_range( error_message );
    }

# endif // if SPACE_DIMS == 3
}

# endif // if defined (STRICT_CHECK)


//------------------------------------------------------------------------------------------------------------
//! \brief  Computes the one-dimensional index for the specified DG coefficient in the array
//!         RKDG::CrossSection::tensor.
//!
//! \param[in]  i   Spatial cell index (including ghost cells) in the \f$ x_1 \f$ dimension.
//! \param[in]  j   Spatial cell index (including ghost cells) in the \f$ x_2 \f$ dimension.
//! \param[in]  a   Degree of first basis polynomial in contracted tensor in the \f$ x_1 \f$ dimension.
//! \param[in]  b   Degree of first basis polynomial in contracted tensor in the \f$ x_2 \f$ dimension.
//! \param[in]  d   Degree of second basis polynomial in contracted tensor in the \f$ x_1 \f$ dimension.
//! \param[in]  e   Degree of second basis polynomial in contracted tensor in the \f$ x_2 \f$ dimension.
//!
//! \return     Returns the appropriate index in the RKDG::CrossSection::tensor array for the specified
//!             coefficient.
//!
//! \see    RKDG::CrossSection::operator()
//------------------------------------------------------------------------------------------------------------
inline size_t CrossSection::TensorIndex (

        const int64_t i
# if SPACE_DIMS >= 2 || defined (DOXYCOMPILE)
    ,   const int64_t j
# endif
# if SPACE_DIMS == 3
    ,   const int64_t k
# endif

    ,   const int64_t a
# if SPACE_DIMS >= 2 || defined (DOXYCOMPILE)
    ,   const int64_t b
# endif
# if SPACE_DIMS == 3
    ,   const int64_t c
# endif

    ,   const int64_t d
# if SPACE_DIMS >= 2 || defined (DOXYCOMPILE)
    ,   const int64_t e
# endif
# if SPACE_DIMS == 3
    ,   const int64_t f
# endif

) const {

# if defined (STRICT_CHECK)

    CheckBounds(
            # if SPACE_DIMS == 1
                i,a,d
            # elif SPACE_DIMS == 2
                i,j,a,b,d,e
            # elif SPACE_DIMS == 3
                i,j,k,a,b,c,d,e,f
            # endif
                , __func__ );

# endif // if defined (STRICT_CHECK)

    return (
    # if SPACE_DIMS == 3
        c + (this->DG_degree + 1)*
    # endif
        (
        # if SPACE_DIMS >= 2
            b + (this->DG_degree + 1)*
        # endif
            (
                a + (this->DG_degree + 1)*
                (
                # if SPACE_DIMS == 3
                    f + (this->DG_degree + 1)*
                # endif
                    (
                    # if SPACE_DIMS >= 2
                        e + (this->DG_degree + 1)*
                    # endif
                        (
                            d + (this->DG_degree + 1)*
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
            )
        )
    );
}

//------------------------------------------------------------------------------------------------------------
//! \brief  Returns a reference to the specified coefficient inside the array RKDG::CrossSection::tensor.
//!
//! \param[in]  i   Spatial cell index (including ghost cells) in the \f$ x_1 \f$ dimension.
//! \param[in]  j   Spatial cell index (including ghost cells) in the \f$ x_2 \f$ dimension.
//! \param[in]  a   Degree of first basis polynomial in contracted tensor in the \f$ x_1 \f$ dimension.
//! \param[in]  b   Degree of first basis polynomial in contracted tensor in the \f$ x_2 \f$ dimension.
//! \param[in]  d   Degree of second basis polynomial in contracted tensor in the \f$ x_1 \f$ dimension.
//! \param[in]  e   Degree of second basis polynomial in contracted tensor in the \f$ x_2 \f$ dimension.
//!
//! \return     Returns a reference to the specified coefficient.
//!
//! \see    RKDG::CrossSection::TensorIndex()
//------------------------------------------------------------------------------------------------------------
inline double & CrossSection::operator() (

        const int64_t i
# if SPACE_DIMS >= 2 || defined (DOXYCOMPILE)
    ,   const int64_t j
# endif
# if SPACE_DIMS == 3
    ,   const int64_t k
# endif

    ,   const int64_t a
# if SPACE_DIMS >= 2 || defined (DOXYCOMPILE)
    ,   const int64_t b
# endif
# if SPACE_DIMS == 3
    ,   const int64_t c
# endif

    ,   const int64_t d
# if SPACE_DIMS >= 2 || defined (DOXYCOMPILE)
    ,   const int64_t e
# endif
# if SPACE_DIMS == 3
    ,   const int64_t f
# endif

) const {

# if defined (STRICT_CHECK)

    if ( this->tensor == nullptr ) {

        std::string error_message = "Pointer this->tensor is NULL in '" + std::string(__func__) + "'.\n";

        PRINT_ERROR( error_message.c_str() )
        throw std::runtime_error( error_message );
    }

# endif // if defined (STRICT_CHECK)

    return this->tensor[ TensorIndex
                                # if SPACE_DIMS == 1
                                    (i,a,d)
                                # elif SPACE_DIMS == 2
                                    (i,j,a,b,d,e)
                                # elif SPACE_DIMS == 3
                                    (i,j,k,a,b,c,d,e,f)
                                # endif
                        ];
}



} // namespace RKDG


# endif // ifndef __RKDG__CROSS_SECTION_HPP__
