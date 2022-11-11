//------------------------------------------------------------------------------------------------------------
//! \file   objects/RKDG/OrdinateFlux.hpp
//! \brief  Header for RKDG::OrdinateFlux class.
//!
//! \author Michael M. Crockatt
//! \date   January 2018
//------------------------------------------------------------------------------------------------------------

# ifndef __RKDG__ORDINATE_FLUX_HPP__
# define __RKDG__ORDINATE_FLUX_HPP__


# include <vector>

# include "objects/Abstract/OrdinateFlux.hpp"
# include "objects/RKDG/DensityFunction.hpp"


namespace RKDG {


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
    //! Maximum degree of DG basis polynomials.
    //!
    int64_t DG_degree;


    //========================================================================================================
    //=== CONSTRUCTORS, DESTRUCTOR, AND RECONFIGURATION ROUTINES =============================================
    //========================================================================================================

    //--------------------------------------------------------------------------------------------------------
    //! \brief  Constructs an empty RKDG::OrdinateFlux object with zero-initialized parameters.
    //--------------------------------------------------------------------------------------------------------
    OrdinateFlux( void );

    //--------------------------------------------------------------------------------------------------------
    //! \brief  Delete copy constructor.
    //--------------------------------------------------------------------------------------------------------
    OrdinateFlux( const OrdinateFlux & ) = delete;

    //--------------------------------------------------------------------------------------------------------
    //! \brief  Constructs an RKDG::OrdinateFlux object with the given set of parameters.
    //--------------------------------------------------------------------------------------------------------
    OrdinateFlux(
        const DomainDecomposition & spatial_params,
        const int64_t DG_degree_in,
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
    //! \brief  Constructs an RKDG::OrdinateFlux object with the given set of parameters.
    //--------------------------------------------------------------------------------------------------------
    OrdinateFlux(
        const DomainDecomposition &,
        const OrdinateSet &,
        const int64_t
    );

    //
    // Acknowledge use of Reconfigure routines from base classes.
    //
    using Abstract::OrdinateFlux::Reconfigure;
    using Abstract::DensityFunction::Reconfigure;

    //--------------------------------------------------------------------------------------------------------
    //! \brief  Reconfigures an RKDG::OrdinateFlux object with the given set of parameters.
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
    //! \brief  Reconfigures an RKDG::OrdinateFlux object with the parameters of a given OrdinateSet object.
    //--------------------------------------------------------------------------------------------------------
    OrdinateFlux & Reconfigure( const OrdinateSet & ) override;

    //--------------------------------------------------------------------------------------------------------
    //! \brief  Reconfigures an RKDG::OrdinateFlux object with the given set of parameters.
    //--------------------------------------------------------------------------------------------------------
    virtual OrdinateFlux & Reconfigure(
        const DomainDecomposition & spatial_params,
        const int64_t DG_degree_in,
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
    //! \brief  Reconfigures an RKDG::OrdinateFlux object with the given set of parameters.
    //--------------------------------------------------------------------------------------------------------
    virtual OrdinateFlux & Reconfigure(
        const DomainDecomposition &,
        const OrdinateSet &,
        const int64_t
    );


    //========================================================================================================
    //=== INDEXING AND EVALUATION MEMBER FUNCTIONS ===========================================================
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
    ) const;

    //--------------------------------------------------------------------------------------------------------
    //! \brief  Evaluates the density function stored in an RKDG::OrdinateFlux object at the specified
    //!         ordinate and spatial point.
    //--------------------------------------------------------------------------------------------------------
    double EvaluateAtPoint(
            const double x
    # if SPACE_DIMS >= 2 || defined (DOXYCOMPILE)
        ,   const double y
    # endif
    # if SPACE_DIMS == 3
        ,   const double z
    # endif
        ,   const int64_t q
    ) const;

    //--------------------------------------------------------------------------------------------------------
    //! \brief  Evaluates the density function stored in an RKDG::OrdinateFlux object at the specified
    //!         ordinate and spatial point.
    //--------------------------------------------------------------------------------------------------------
    double EvaluateAtLocalPoint(
            const double x
    # if SPACE_DIMS >= 2 || defined (DOXYCOMPILE)
        ,   const double y
    # endif
    # if SPACE_DIMS == 3
        ,   const double z
    # endif
        ,   const int64_t q
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
    //! \brief  Two RKDG::OrdinateFlux objects "are matching" if all of their parameters are equal.
    //--------------------------------------------------------------------------------------------------------
    static bool AreMatching( const OrdinateFlux &, const OrdinateFlux & );

    //--------------------------------------------------------------------------------------------------------
    //! \brief  Swaps the contents of two RKDG::OrdinateFlux objects.
    //--------------------------------------------------------------------------------------------------------
    static void swap( OrdinateFlux &, OrdinateFlux & );

    //--------------------------------------------------------------------------------------------------------
    //! \brief  Computes the difference between two RKDG::OrdinateFlux objects.
    //--------------------------------------------------------------------------------------------------------
    static void ComputeDifference( const OrdinateFlux &, const OrdinateFlux &, OrdinateFlux & );


# if SPACE_DIMS == 1 || defined (DOXYCOMPILE)

    //--------------------------------------------------------------------------------------------------------
    //! \brief  Computes the distances (with respect to the spatial measure) between the angular moments of
    //!         two OrdinateFlux objects.
    //--------------------------------------------------------------------------------------------------------
    static std::vector<double> ComputeMomentDistances(
        const OrdinateFlux &,
        const OrdinateFlux &,
        double (RKDG::DensityFunction::*norm)() const = &RKDG::DensityFunction::L2Norm
    );

# endif // if SPACE_DIMS == 1


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

    //--------------------------------------------------------------------------------------------------------
    //! \brief  Computes the total mass of the density function stored in an RKDG::OrdinateFlux object.
    //--------------------------------------------------------------------------------------------------------
    virtual double TotalMass( void ) const;

    //--------------------------------------------------------------------------------------------------------
    //! \brief  Computes the \f$ L^1( d\vec{x} ) \f$ norm of the density function stored in an
    //!         RKDG::OrdinateFlux object.
    //--------------------------------------------------------------------------------------------------------
    virtual double L1Norm( void ) const;

    //--------------------------------------------------------------------------------------------------------
    //! \brief  Computes the \f$ L^2( d\vec{x} ) \f$ norm of the density function stored in an
    //!         RKDG::OrdinateFlux object.
    //--------------------------------------------------------------------------------------------------------
    virtual double L2Norm( void ) const;

    //--------------------------------------------------------------------------------------------------------
    //! \brief  Computes the \f$ L^{\infty}( d\vec{x} ) \f$ norm of the density function stored in an
    //!         RKDG::OrdinateFlux object.
    //--------------------------------------------------------------------------------------------------------
    virtual double LinfNorm( const double tol = 1e-3 ) const;

    //--------------------------------------------------------------------------------------------------------
    //! \brief  Outputs a grid of sample points from the density function stored in an RKDG::OrdinateFlux
    //!         object for plotting.
    //--------------------------------------------------------------------------------------------------------
    void OutputPlot( const std::string filename ) const;

    //--------------------------------------------------------------------------------------------------------
    //! \brief  Routine to output a grid of points sampling the difference between two OrdinateFlux objects
    //!         for plotting.
    //--------------------------------------------------------------------------------------------------------
    void OutputDiffPlot(
        const OrdinateFlux &,
        const OrdinateFlux &,
        const std::string filename
    ) const;

    //--------------------------------------------------------------------------------------------------------
    //! \brief  Reconstructs an RKDG::OrdinateFlux object using data from a binary file on disk.
    //--------------------------------------------------------------------------------------------------------
    void ReadFromDisk(
    # if defined (ENABLE_MPI)
        MPI_File & fp
    # else
        std::FILE * & fp
    # endif
    );

    //--------------------------------------------------------------------------------------------------------
    //! \brief  Reconstructs an RKDG::OrdinateFlux object using data from a binary file on disk.
    //--------------------------------------------------------------------------------------------------------
    void ReadFromDisk(
        const std::string filename
    );

    //--------------------------------------------------------------------------------------------------------
    //! \brief  Writes the contents of an RKDG::OrdinateFlux object to a binary file on disk.
    //--------------------------------------------------------------------------------------------------------
    void WriteToDisk(
    # if defined (ENABLE_MPI)
        MPI_File & fp
    # else
        std::FILE * & fp
    # endif
    ) const;

    //--------------------------------------------------------------------------------------------------------
    //! \brief  Writes the contents of an RKDG::OrdinateFlux object to a binary file on disk.
    //--------------------------------------------------------------------------------------------------------
    void WriteToDisk( const std::string filename ) const;

# if SPACE_DIMS == 1 || defined (DOXYCOMPILE)

    //--------------------------------------------------------------------------------------------------------
    //! \brief  Computes the angular moments of the ordinate flux. The moments are returned as
    //!         RKDG::DensityFunction elements in an std::vector container, in order of moment degree.
    //--------------------------------------------------------------------------------------------------------
    std::vector<RKDG::DensityFunction> ComputeAngularMoments( void ) const;

# endif // if SPACE_DIMS == 1


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
//!     CheckBounds(q,i,j,d,e,__func__)
//! \endcode
//!
//! \param[in]      q           Index of discrete ordinate.
//! \param[in]      i           Spatial cell index (including ghost cells) in the \f$ x_1 \f$ dimension.
//! \param[in]      j           Spatial cell index (including ghost cells) in the \f$ x_2 \f$ dimension.
//! \param[in]      d           Degree of spatial basis polynomial wrt. \f$ x_1 \f$.
//! \param[in]      e           Degree of spatial basis polynomial wrt. \f$ x_2 \f$.
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
//!         OrdinateFlux::density.
//!
//! \param[in]      q           Index of discrete ordinate.
//! \param[in]      i           Spatial cell index (including ghost cells) in the \f$ x_1 \f$ dimension.
//! \param[in]      j           Spatial cell index (including ghost cells) in the \f$ x_2 \f$ dimension.
//! \param[in]      d           Degree of spatial basis polynomial wrt. \f$ x_1 \f$.
//! \param[in]      e           Degree of spatial basis polynomial wrt. \f$ x_2 \f$.
//!
//! \return     Returns the appropriate index in the OrdinateFlux::density array for the specified
//!             coefficient.
//!
//! \see    RKDG::OrdinateFlux::operator()()
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

) const {

# if defined (STRICT_CHECK)

    CheckBounds(
            # if SPACE_DIMS == 1
                q,i,d
            # elif SPACE_DIMS == 2
                q,i,j,d,e
            # elif SPACE_DIMS == 3
                q,i,j,k,d,e,f
            # endif
                , __func__ );

# endif // if defined (STRICT_CHECK)

    return (
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
//!
//! \return     Reference to the specified coefficient.
//!
//! \see    RKDG::OrdinateFlux::Index()
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
                                (q,i,d)
                            # elif SPACE_DIMS == 2
                                (q,i,j,d,e)
                            # elif SPACE_DIMS == 3
                                (q,i,j,k,d,e,f)
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
//!
//! \see    RKDG::OrdinateFlux::Index()
//! \see    RKDG::OrdinateFlux::operator()()
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
                            (q,i,d)
                        # elif SPACE_DIMS == 2
                            (q,i,j,d,e)
                        # elif SPACE_DIMS == 3
                            (q,i,j,k,d,e,f)
                        # endif
                        ;
}


} // namespace RKDG


# endif // ifndef __RKDG__ORDINATE_FLUX_HPP__
