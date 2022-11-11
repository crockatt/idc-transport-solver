//------------------------------------------------------------------------------------------------------------
//! \file   objects/RKDG/DensityFunction.hpp
//! \brief  Header for RKDG::DensityFunction class.
//!
//! \author Michael M. Crockatt
//! \date   February 2018
//------------------------------------------------------------------------------------------------------------

# ifndef __RKDG__DENSITY_FUNCTION_HPP__
# define __RKDG__DENSITY_FUNCTION_HPP__


# include "utils/CLog.hpp"
# include "utils/global.hpp"
# include "objects/Abstract/DensityFunction.hpp"


namespace RKDG {


//------------------------------------------------------------------------------------------------------------
//! \brief  Class used for storing spatially-varying scalar-valued density functions with a discontinuous
//!         Galerkin approximation of the spatial variables.
//!
//!
//! \attention  Virtual destructor inherited from Abstract::DensityFunction class.
//!
//!
//------------------------------------------------------------------------------------------------------------
class DensityFunction : public virtual Abstract::DensityFunction {

public:

    //!
    //! Maximum degree of DG basis polynomials.
    //!
    int64_t DG_degree;


    //========================================================================================================
    //=== CONSTRUCTORS, DESTRUCTOR, AND RECONFIGURATION ROUTINES =============================================
    //========================================================================================================

    //--------------------------------------------------------------------------------------------------------
    //! \brief  Constructs an empty RKDG::DensityFunction object with zero-initialized parameters.
    //--------------------------------------------------------------------------------------------------------
    DensityFunction( void );

    //--------------------------------------------------------------------------------------------------------
    //! \brief  Delete copy constructor.
    //--------------------------------------------------------------------------------------------------------
    DensityFunction( const DensityFunction & ) = delete;

    //--------------------------------------------------------------------------------------------------------
    //! \brief  Constructs an RKDG::DensityFunction object with the given set of parameters.
    //--------------------------------------------------------------------------------------------------------
    DensityFunction(
        const DomainDecomposition & spatial_params,
        const int64_t DG_degree_in
    );

    //
    // Acknowledge use of Reconfigure routines from base class.
    //
    using Abstract::DensityFunction::Reconfigure;

    //--------------------------------------------------------------------------------------------------------
    //! \brief  Reconfigures an RKDG::DensityFunction object with a different set of parameters.
    //--------------------------------------------------------------------------------------------------------
    virtual DensityFunction & Reconfigure( const int64_t DG_degree_in );

    //--------------------------------------------------------------------------------------------------------
    //! \brief  Reconfigures an RKDG::DensityFunction object with a different set of parameters.
    //--------------------------------------------------------------------------------------------------------
    virtual DensityFunction & Reconfigure(
        const DomainDecomposition & spatial_params,
        const int64_t DG_degree
    );


    //========================================================================================================
    //=== INDEXING AND EVALUATION MEMBER FUNCTIONS ===========================================================
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
    ) const;

    //--------------------------------------------------------------------------------------------------------
    //! \brief  Returns a pointer to the specified DG coefficient inside the array DensityFunction::density.
    //--------------------------------------------------------------------------------------------------------
    inline double * PointerAt (
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
    ) const;

    //--------------------------------------------------------------------------------------------------------
    //! \brief  Evaluates the density function stored in an RKDG::DensityFunction object at the specified
    //!         spatial point.
    //--------------------------------------------------------------------------------------------------------
    double EvaluateAtPoint(
            const double x
    # if SPACE_DIMS >= 2 || defined (DOXYCOMPILE)
        ,   const double y
    # endif
    # if SPACE_DIMS == 3
        ,   const double z
    # endif
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

    //--------------------------------------------------------------------------------------------------------
    //! \brief  Two RKDG::DensityFunction objects "are matching" if all of their parameters are equal.
    //--------------------------------------------------------------------------------------------------------
    static bool AreMatching( const DensityFunction &, const DensityFunction & );

    //--------------------------------------------------------------------------------------------------------
    //! \brief  Outputs a grid of points sampling the difference between two density functions for plotting.
    //--------------------------------------------------------------------------------------------------------
    static void OutputDiffPlot( const DensityFunction &, const DensityFunction &, const std::string );

    //--------------------------------------------------------------------------------------------------------
    //! \brief  Computes the difference between two RKDG::DensityFunction objects.
    //--------------------------------------------------------------------------------------------------------
    static void ComputeDifference( const DensityFunction &, const DensityFunction &, DensityFunction & );


    //========================================================================================================
    //=== ADDITIONAL PUBLIC MEMBER FUNCTIONS =================================================================
    //========================================================================================================

    //--------------------------------------------------------------------------------------------------------
    //! \brief  Prints a summary of the object configuration to the logging interface.
    //--------------------------------------------------------------------------------------------------------
    void Print( const std::string & = "  " ) const override;

    //--------------------------------------------------------------------------------------------------------
    //! \brief  Deallocates memory at all internal pointers and zeros all parameters of an
    //!         RKDG::DensityFunction object.
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

    //--------------------------------------------------------------------------------------------------------
    //! \brief  Computes the total mass of the density function stored in an RKDG::DensityFunction object.
    //--------------------------------------------------------------------------------------------------------
    virtual double TotalMass( void ) const;

    //--------------------------------------------------------------------------------------------------------
    //! \brief  Computes the \f$ L^1( d\vec{x} ) \f$ norm of the density function stored in an
    //!         RKDG::DensityFunction object.
    //--------------------------------------------------------------------------------------------------------
    virtual double L1Norm( void ) const;

    //--------------------------------------------------------------------------------------------------------
    //! \brief  Computes the \f$ L^2( d\vec{x} ) \f$ norm of the density function stored in an
    //!         RKDG::DensityFunction object.
    //--------------------------------------------------------------------------------------------------------
    virtual double L2Norm( void ) const;

    //--------------------------------------------------------------------------------------------------------
    //! \brief  Computes the \f$ L^{\infty}( d\vec{x} ) \f$ norm of the density function stored in an
    //!         RKDG::DensityFunction object.
    //--------------------------------------------------------------------------------------------------------
    virtual double LinfNorm( const double tol = 1e-3 ) const;

    //--------------------------------------------------------------------------------------------------------
    //! \brief  Outputs a grid of sample points from the density function stored in an RKDG::DensityFunction
    //!         object for plotting.
    //--------------------------------------------------------------------------------------------------------
    virtual void OutputPlot( const std::string filename ) const;

# ifndef ENABLE_MPI

    //--------------------------------------------------------------------------------------------------------
    //! \brief  Outputs cell averages of the density function to a binary file.
    //--------------------------------------------------------------------------------------------------------
    virtual void OutputCellAveragesToBin( const std::string filename ) const;

# endif // ifndef ENABLE_MPI

    //--------------------------------------------------------------------------------------------------------
    //! \brief  Reconstructs an RKDG::DensityFunction object using data from a binary file on disk.
    //--------------------------------------------------------------------------------------------------------
    void ReadFromDisk(
    # if defined (ENABLE_MPI)
        MPI_File & fp,
    # else
        std::FILE * & fp,
    # endif
        const int64_t version = 0
    );

    //--------------------------------------------------------------------------------------------------------
    //! \brief  Reconstructs an RKDG::DensityFunction object using data from a binary file on disk.
    //--------------------------------------------------------------------------------------------------------
    void ReadFromDisk(
        const std::string filename,
        const int64_t version = 0
    );

    //--------------------------------------------------------------------------------------------------------
    //! \brief  Writes the contents of an RKDG::DensityFunction object to a binary file on disk.
    //--------------------------------------------------------------------------------------------------------
    void WriteToDisk(
    # if defined (ENABLE_MPI)
        MPI_File & fp
    # else
        std::FILE * & fp
    # endif
    ) const;

    //--------------------------------------------------------------------------------------------------------
    //! \brief  Writes the contents of an RKDG::DensityFunction object to a binary file on disk.
    //--------------------------------------------------------------------------------------------------------
    void WriteToDisk( const std::string filename ) const;


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

//--------------------------------------------------------------------------------------------------------
//! \brief  Performs bounds checking on indices.
//!
//! An exception of type \c std::out_of_range is thrown if an out-of-bounds index is detected.
//!
//! This function is intended to be called as:
//! \code
//!     CheckBounds(i,j,d,e,__func__)
//! \endcode
//!
//! \param[in]      i           Spatial cell index (including ghost cells) in the \f$ x_1 \f$ dimension.
//! \param[in]      j           Spatial cell index (including ghost cells) in the \f$ x_2 \f$ dimension.
//! \param[in]      d           Degree of spatial basis polynomial wrt. \f$ x_1 \f$.
//! \param[in]      e           Degree of spatial basis polynomial wrt. \f$ x_2 \f$.
//! \param[in]      func        Calling function.
//!
//--------------------------------------------------------------------------------------------------------
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

    ,   const char * const func

) const {

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
//! \brief  Returns the one-dimensional index for the specified DG coefficient in the array
//!         DensityFunction::density.
//!
//! \param[in]      i           Spatial cell index (including ghost cells) in the \f$ x_1 \f$ dimension.
//! \param[in]      j           Spatial cell index (including ghost cells) in the \f$ x_2 \f$ dimension.
//! \param[in]      d           Degree of spatial basis polynomial wrt. \f$ x_1 \f$.
//! \param[in]      e           Degree of spatial basis polynomial wrt. \f$ x_2 \f$.
//!
//! \return     Appropriate index in DensityFunction::density array for specified coefficient.
//!
//! \see    RKDG::DensityFunction::IndexNoGC()
//! \see    RKDG::DensityFunction::operator()()
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

) const {

# if defined (STRICT_CHECK)

    CheckBounds(
            # if SPACE_DIMS == 1
                i,d
            # elif SPACE_DIMS == 2
                i,j,d,e
            # elif SPACE_DIMS == 3
                i,j,k,d,e,f
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
    );
}


//------------------------------------------------------------------------------------------------------------
//! \brief  Returns a reference to the specified DG coefficient inside the array DensityFunction::density.
//!
//! \param[in]      i           Spatial cell index (including ghost cells) in the \f$ x_1 \f$ dimension.
//! \param[in]      j           Spatial cell index (including ghost cells) in the \f$ x_2 \f$ dimension.
//! \param[in]      d           Degree of spatial basis polynomial wrt. \f$ x_1 \f$.
//! \param[in]      e           Degree of spatial basis polynomial wrt. \f$ x_2 \f$.
//!
//! \return     Reference to the specified coefficient.
//!
//! \see    RKDG::DensityFunction::Index()
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
                                (i,d)
                            # elif SPACE_DIMS == 2
                                (i,j,d,e)
                            # elif SPACE_DIMS == 3
                                (i,j,k,d,e,f)
                            # endif
                        ];
}


//------------------------------------------------------------------------------------------------------------
//! \brief  Returns a pointer to the specified DG coefficient inside the array DensityFunction::density.
//!
//! \param[in]      i           Spatial cell index (including ghost cells) in the \f$ x_1 \f$ dimension.
//! \param[in]      j           Spatial cell index (including ghost cells) in the \f$ x_2 \f$ dimension.
//! \param[in]      d           Degree of spatial basis polynomial wrt. \f$ x_1 \f$.
//! \param[in]      e           Degree of spatial basis polynomial wrt. \f$ x_2 \f$.
//!
//! \return     Reference to the specified coefficient.
//!
//! \see    RKDG::DensityFunction::Index()
//------------------------------------------------------------------------------------------------------------
inline double * DensityFunction::PointerAt (

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
                                (i,d)
                            # elif SPACE_DIMS == 2
                                (i,j,d,e)
                            # elif SPACE_DIMS == 3
                                (i,j,k,d,e,f)
                            # endif
                            ;
}


} // namespace RKDG


# endif // ifndef __RKDG__DENSITY_FUNCTION_HPP__
