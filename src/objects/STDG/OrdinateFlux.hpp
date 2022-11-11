//------------------------------------------------------------------------------------------------------------
//! \file   objects/STDG/OrdinateFlux.hpp
//! \brief  Header for STDG::OrdinateFlux objects.
//!
//! \author Michael M. Crockatt
//! \date   February 2018
//------------------------------------------------------------------------------------------------------------

# ifndef __STDG__ORDINATE_FLUX_HPP__
# define __STDG__ORDINATE_FLUX_HPP__


# include "utils/Quadrule/Quadrule.hpp"
# include "objects/Abstract/OrdinateFlux.hpp"
# include "objects/RKDG/DensityFunction.hpp"


namespace STDG {


//------------------------------------------------------------------------------------------------------------
//! \brief  Class used for storing spatially-varying scalar-valued angular flux density functions with a
//!         discontinuous Galerkin approximation of the spatial variables and a discrete ordinates
//!         approximation of the angular variable.
//!
//!
//! \attention  Virtual destructor inherited from Abstract::DensityFunction class.
//!
//!
//------------------------------------------------------------------------------------------------------------
class OrdinateFlux : public Abstract::OrdinateFlux, public virtual Abstract::DensityFunction {

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
    //! \brief  Constructs an empty STDG::OrdinateFlux object with zero-initialized parameters.
    //--------------------------------------------------------------------------------------------------------
    OrdinateFlux( void );

    //--------------------------------------------------------------------------------------------------------
    //! \brief  Delete copy constructor.
    //--------------------------------------------------------------------------------------------------------
    OrdinateFlux( const OrdinateFlux & ) = delete;

    //--------------------------------------------------------------------------------------------------------
    //! \brief  Constructs an STDG::OrdinateFlux object with the given set of parameters.
    //--------------------------------------------------------------------------------------------------------
    OrdinateFlux(
        const DomainDecomposition & spatial_params,
        const int64_t DG_degree,
        const int64_t ang_order,
        const bool symmetric_reduce = false,
        const OrdinateType ordinate_type
            # if SPACE_DIMS == 1
                = OrdinateType::GaussLegendre
            # elif SPACE_DIMS >= 2 || defined (DOXYCOMPILE)
                = OrdinateType::ChebyshevLegendre
            # endif
    );

    //--------------------------------------------------------------------------------------------------------
    //! \brief  Constructs an STDG::OrdinateFlux object with the given set of parameters.
    //--------------------------------------------------------------------------------------------------------
    OrdinateFlux(
        const DomainDecomposition & spatial_params,
        const int64_t DG_degree_x_in,
        const int64_t DG_degree_t_in,
        const int64_t ang_order,
        const bool symmetric_reduce = false,
        const OrdinateType ordinate_type
            # if SPACE_DIMS == 1
                = OrdinateType::GaussLegendre
            # elif SPACE_DIMS >= 2 || defined (DOXYCOMPILE)
                = OrdinateType::ChebyshevLegendre
            # endif
    );

    //--------------------------------------------------------------------------------------------------------
    //! \brief  Constructs an STDG::OrdinateFlux object with the given set of parameters.
    //--------------------------------------------------------------------------------------------------------
    OrdinateFlux(
        const DomainDecomposition &,
        const OrdinateSet &,
        const int64_t
    );

    //--------------------------------------------------------------------------------------------------------
    //! \brief  Constructs an STDG::OrdinateFlux object with the given set of parameters.
    //--------------------------------------------------------------------------------------------------------
    OrdinateFlux(
        const DomainDecomposition &,
        const OrdinateSet &,
        const int64_t DG_degree_x_in,
        const int64_t DG_degree_t_in
    );

    //
    // Acknowledge use of Reconfigure routines from base classes.
    //
    using Abstract::OrdinateFlux::Reconfigure;
    using Abstract::DensityFunction::Reconfigure;

    //--------------------------------------------------------------------------------------------------------
    //! \brief  Reconfigures an STDG::OrdinateFlux object with the given set of parameters.
    //--------------------------------------------------------------------------------------------------------
    OrdinateFlux & Reconfigure(
        const int64_t ang_order,
        const bool symmetric_reduce = false,
        const OrdinateType ordinate_type
            # if SPACE_DIMS == 1
                = OrdinateType::GaussLegendre
            # elif SPACE_DIMS >= 2 || defined (DOXYCOMPILE)
                = OrdinateType::ChebyshevLegendre
            # endif
    ) override;

    //--------------------------------------------------------------------------------------------------------
    //! \brief  Reconfigures an STDG::OrdinateFlux object with the parameters of a given OrdinateSet object.
    //--------------------------------------------------------------------------------------------------------
    OrdinateFlux & Reconfigure( const OrdinateSet & ) override;

    //--------------------------------------------------------------------------------------------------------
    //! \brief  Reconfigures an STDG::OrdinateFlux object with the given set of parameters.
    //--------------------------------------------------------------------------------------------------------
    virtual OrdinateFlux & Reconfigure(
        const DomainDecomposition & spatial_params,
        const int64_t DG_degree,
        const int64_t ang_order,
        const bool symmetric_reduce = false,
        const OrdinateType ordinate_type
            # if SPACE_DIMS == 1
                = OrdinateType::GaussLegendre
            # elif SPACE_DIMS >= 2 || defined (DOXYCOMPILE)
                = OrdinateType::ChebyshevLegendre
            # endif
    );

    //--------------------------------------------------------------------------------------------------------
    //! \brief  Reconfigures an STDG::OrdinateFlux object with the given set of parameters.
    //--------------------------------------------------------------------------------------------------------
    virtual OrdinateFlux & Reconfigure(
        const DomainDecomposition & spatial_params,
        const int64_t DG_degree_x_in,
        const int64_t DG_degree_t_in,
        const int64_t ang_order,
        const bool symmetric_reduce = false,
        const OrdinateType ordinate_type
            # if SPACE_DIMS == 1
                = OrdinateType::GaussLegendre
            # elif SPACE_DIMS >= 2 || defined (DOXYCOMPILE)
                = OrdinateType::ChebyshevLegendre
            # endif
    );

    //--------------------------------------------------------------------------------------------------------
    //! \brief  Reconfigures an STDG::OrdinateFlux object with the given set of parameters.
    //--------------------------------------------------------------------------------------------------------
    virtual OrdinateFlux & Reconfigure(
        const DomainDecomposition &,
        const OrdinateSet &,
        const int64_t
    );

    //--------------------------------------------------------------------------------------------------------
    //! \brief  Reconfigures an STDG::OrdinateFlux object with the given set of parameters.
    //--------------------------------------------------------------------------------------------------------
    virtual OrdinateFlux & Reconfigure(
        const DomainDecomposition & spatial_params,
        const OrdinateSet & ordinate_set,
        const int64_t DG_degree_x_in,
        const int64_t DG_degree_t_in
    );


    //========================================================================================================
    //=== INDEXING MEMBER FUNCTIONS ==========================================================================
    //========================================================================================================

# if defined (STRICT_CHECK) || defined (DOXYCOMPILE)

    //--------------------------------------------------------------------------------------------------------
    //! \brief  Performs bounds checking on indices.
    //--------------------------------------------------------------------------------------------------------
    inline void CheckBounds (
            const int64_t q
        ,   const int64_t i
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
    //! \brief  Computes the one-dimensional index for the specified DG coefficient in the array
    //!         OrdinateFlux::density.
    //--------------------------------------------------------------------------------------------------------
    inline size_t Index (
            const int64_t q
        ,   const int64_t i
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
            const int64_t q
        ,   const int64_t i
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
    //! \brief  Returns a pointer to the specified DG coefficient inside the array DensityFunction::density.
    //--------------------------------------------------------------------------------------------------------
    inline double * PointerAt (
            const int64_t q
        ,   const int64_t i
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
    OrdinateFlux & operator=( const OrdinateFlux & ) = delete;

    //--------------------------------------------------------------------------------------------------------
    //! \brief  Inherit overloaded equality operator from Abstract::OrdinateFlux.
    //--------------------------------------------------------------------------------------------------------
    using Abstract::OrdinateFlux::operator==;

    //--------------------------------------------------------------------------------------------------------
    //! \brief  Inherit overloaded inequality operator from Abstract::OrdinateFlux.
    //--------------------------------------------------------------------------------------------------------
    using Abstract::OrdinateFlux::operator!=;


    //========================================================================================================
    //=== STATIC MEMBER FUNCTIONS ============================================================================
    //========================================================================================================

    //--------------------------------------------------------------------------------------------------------
    //! \brief  Creates and returns an OrdinateFlux object.
    //--------------------------------------------------------------------------------------------------------
    static OrdinateFlux * Create(
        const DomainDecomposition &,
        const OrdinateSet &,
        const ParameterList &
    );

    //--------------------------------------------------------------------------------------------------------
    //! \brief  Two STDG::OrdinateFlux objects "are matching" if all of their parameters are equal.
    //--------------------------------------------------------------------------------------------------------
    static bool AreMatching( const OrdinateFlux &, const OrdinateFlux & );

    //--------------------------------------------------------------------------------------------------------
    //! \brief  Swaps the contents of two STDG::OrdinateFlux objects.
    //--------------------------------------------------------------------------------------------------------
    static void swap( OrdinateFlux &, OrdinateFlux & );

    //--------------------------------------------------------------------------------------------------------
    //! \brief  \f$ \vec{y} \gets \beta \vec{y} + \alpha \vec{x} \f$.
    //--------------------------------------------------------------------------------------------------------
    static void AXPY(
        const double,
        const double,
        const OrdinateFlux &,
        RKDG::OrdinateFlux &
    );

    //--------------------------------------------------------------------------------------------------------
    //! \brief  \f$ \vec{y} \gets \beta \vec{y} + \alpha \vec{x} \f$.
    //--------------------------------------------------------------------------------------------------------
    static void AXPY(
        const double,
        const double,
        const RKDG::OrdinateFlux &,
        OrdinateFlux &
    );


    //========================================================================================================
    //=== ADDITIONAL PUBLIC MEMBER FUNCTIONS =================================================================
    //========================================================================================================

    //--------------------------------------------------------------------------------------------------------
    //! \brief  Prints a summary of the object configuration to the logging interface.
    //--------------------------------------------------------------------------------------------------------
    void Print( const std::string & = "  " ) const override;

    //--------------------------------------------------------------------------------------------------------
    //! \brief  Updates ghost cells on reflecting boundaries.
    //--------------------------------------------------------------------------------------------------------
    OrdinateFlux & ReflectBoundaries( void ) override;

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

    //--------------------------------------------------------------------------------------------------------
    //! \brief  Returns a pointer to the location within DensityFunction::density where the coefficients for
    //!         the given spatial cell begin.
    //--------------------------------------------------------------------------------------------------------
    double * PointerAtOrdinateCell(
            const int64_t q
        ,   const int64_t i
    # if SPACE_DIMS >= 2
        ,   const int64_t j
    # endif
    # if SPACE_DIMS == 3
        ,   const int64_t k
    # endif
    ) const override;

# if SPACE_DIMS == 2

    //--------------------------------------------------------------------------------------------------------
    //! \brief  Returns a pointer to the location within DensityFunction::density where the coefficients for
    //!         the given spatial cell begin.
    //--------------------------------------------------------------------------------------------------------
    double * PointerAtQuadrantCell(
        const int64_t quad,
        const int64_t i,
        const int64_t j
    ) const override;

# elif SPACE_DIMS == 3

    //--------------------------------------------------------------------------------------------------------
    //! \brief  Returns a pointer to the location within DensityFunction::density where the coefficients for
    //!         the given spatial cell begin.
    //--------------------------------------------------------------------------------------------------------
    double * PointerAtOctantCell(
        const int64_t oct,
        const int64_t i,
        const int64_t j,
        const int64_t k
    ) const override;

# endif // if SPACE_DIMS == ?


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
//!     CheckBounds(q,i,j,d,e,s,__func__)
//! \endcode
//!
//! \param[in]      q           Index of discrete ordinate.
//! \param[in]      i           Spatial cell index (including ghost cells) in the \f$ x_1 \f$ dimension.
//! \param[in]      j           Spatial cell index (including ghost cells) in the \f$ x_2 \f$ dimension.
//! \param[in]      d           Degree of spatial basis polynomial wrt. \f$ x_1 \f$.
//! \param[in]      e           Degree of spatial basis polynomial wrt. \f$ x_2 \f$.
//! \param[in]      s           Degree of temporal basis polynomial.
//! \param[in]      func        Calling function.
//!
//------------------------------------------------------------------------------------------------------------
inline void OrdinateFlux::CheckBounds (

        const int64_t q

    ,   const int64_t i
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

    if ( q < 0 || q >= this->nq() ) {

        std::string error_message
            = "Value " + std::to_string(q) + " for index q in " + std::string(func)
                + " outsize permissible range [ 0, " + std::to_string( this->nq() ) + " ).\n";

        PRINT_ERROR( error_message.c_str() )
        throw std::out_of_range( error_message );
    }

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
//! \brief  Computes the one-dimensional index for the specified DG coefficient in the array
//!         OrdinateFlux::density.
//!
//! \param[in]      q           Index of discrete ordinate.
//! \param[in]      i           Spatial cell index (including ghost cells) in the \f$ x_1 \f$ dimension.
//! \param[in]      j           Spatial cell index (including ghost cells) in the \f$ x_2 \f$ dimension.
//! \param[in]      d           Degree of spatial basis polynomial wrt. \f$ x_1 \f$.
//! \param[in]      e           Degree of spatial basis polynomial wrt. \f$ x_2 \f$.
//! \param[in]      s           Degree of temporal basis polynomial.
//!
//! \return     Returns the appropriate index in the OrdinateFlux::density array for the specified
//!             coefficient.
//!
//! \see    STDG::OrdinateFlux::operator()()
//------------------------------------------------------------------------------------------------------------
inline size_t OrdinateFlux::Index (

        const int64_t q

    ,   const int64_t i
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
                q,i,d,s
            # elif SPACE_DIMS == 2
                q,i,j,d,e,s
            # elif SPACE_DIMS == 3
                q,i,j,k,d,e,f,s
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
                        q + this->nq()*
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
    );
}


//------------------------------------------------------------------------------------------------------------
//! \brief  Returns a reference to the specified DG coefficient inside the array DensityFunction::density.
//!
//! \param[in]      q           Index of discrete ordinate.
//! \param[in]      i           Spatial cell index (including ghost cells) in the \f$ x_1 \f$ dimension.
//! \param[in]      j           Spatial cell index (including ghost cells) in the \f$ x_2 \f$ dimension.
//! \param[in]      d           Degree of spatial basis polynomial wrt. \f$ x_1 \f$.
//! \param[in]      e           Degree of spatial basis polynomial wrt. \f$ x_2 \f$.
//! \param[in]      s           Degree of temporal basis polynomial.
//!
//! \return     Reference to the specified coefficient.
//!
//! \see    STDG::OrdinateFlux::Index()
//------------------------------------------------------------------------------------------------------------
inline double & OrdinateFlux::operator() (

        const int64_t q

    ,   const int64_t i
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
                                (q,i,d,s)
                            # elif SPACE_DIMS == 2
                                (q,i,j,d,e,s)
                            # elif SPACE_DIMS == 3
                                (q,i,j,k,d,e,f,s)
                            # endif
                        ];
}

//------------------------------------------------------------------------------------------------------------
//! \brief  Returns a pointer to the specified DG coefficient inside the array DensityFunction::density.
//!
//! \param[in]      q           Index of discrete ordinate.
//! \param[in]      i           Spatial cell index (including ghost cells) in the \f$ x_1 \f$ dimension.
//! \param[in]      j           Spatial cell index (including ghost cells) in the \f$ x_2 \f$ dimension.
//! \param[in]      d           Degree of spatial basis polynomial wrt. \f$ x_1 \f$.
//! \param[in]      e           Degree of spatial basis polynomial wrt. \f$ x_2 \f$.
//! \param[in]      s           Degree of temporal basis polynomial.
//!
//! \see    STDG::OrdinateFlux::Index()
//! \see    STDG::OrdinateFlux::operator()()
//------------------------------------------------------------------------------------------------------------
inline double * OrdinateFlux::PointerAt (

        const int64_t q

    ,   const int64_t i
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

    return this->density + Index
                            # if SPACE_DIMS == 1
                                (q,i,d,s)
                            # elif SPACE_DIMS == 2
                                (q,i,j,d,e,s)
                            # elif SPACE_DIMS == 3
                                (q,i,j,k,d,e,f,s)
                            # endif
                            ;
}


} // namespace STDG


# endif // ifndef __STDG__ORDINATE_FLUX_HPP__
