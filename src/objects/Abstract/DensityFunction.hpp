//------------------------------------------------------------------------------------------------------------
//! \file   objects/Abstract/DensityFunction.hpp
//! \brief  Header for Abstract::DensityFunction class.
//!
//! \author Michael M. Crockatt
//! \date   February 2018
//------------------------------------------------------------------------------------------------------------

# ifndef __ABSTRACT__DENSITY_FUNCTION_HPP__
# define __ABSTRACT__DENSITY_FUNCTION_HPP__


# include <cstddef>
# include <cstdint>
# include <initializer_list>

# if defined (ENABLE_PETSC) || defined (DOXYCOMPILE)
    # include "utils/PETScInterface.hpp"
# endif

# if defined (ENABLE_AZTEC_EPETRA) || defined (DOXYCOMPILE)
    # include <Epetra_MultiVector.h>
# endif

# if defined (ENABLE_BELOS_TPETRA) || defined (DOXYCOMPILE)
    # include "utils/BelosTpetraInterface.hpp"
# endif

# include "objects/DomainDecomposition.hpp"
# include "utils/MVWrapper.hpp"


namespace Abstract {


//------------------------------------------------------------------------------------------------------------
//! \brief  Abstract class used as a base for objects used to store numerical representations of
//!         spatially-varying density functions.
//!
//! Objects are assumed to be stored using a spatial mesh with parameters described by a DomainDecomposition
//! object and some number of degrees of freedom per mesh cell. These degrees of freedom may describe
//! different characteristics (e.g., basis coefficients with respect to space, time, or angle) depending on
//! the concrete type of the object. This class provides an abstract interface for interacting with these
//! types of objects when the specifics of what the degrees of freedom in each mesh element represent is not
//! directly relevant. Using this interface, certain operations (e.g., scaling a vector) can be generalized
//! and implemented for any object(s) derived from this class.
//!
//! The majority of simple operations for objects derived from this class are implemented with OpenMP
//! parallelism over the mesh cell indices only; that is, loops over the local degrees of freedom in each mesh
//! element are not included in the worksharing constructs used. This provides a consistent memory access
//! pattern for all such objects that is leveraged in combination with memory pinning through the hwloc
//! library. This has two effects:
//!     1.  This ensures that memory allocated for each object is more evenly distributed across all
//!         available NUMA nodes on the system. By doing so, the aggregate memory bandwidth of systems with
//!         multiple NUMA nodes is more effectively utilized.
//!     2.  When used in combination with thread pinning and a first-touch memory allocation policy, this
//!         increases the locality of memory accesses significantly.
//!
//! \attention  This implementation assumes that the indexing function for any derived object is arranged
//!             such that the local degrees of freedom for each cell are placed contiguously in memory.
//!
//------------------------------------------------------------------------------------------------------------
class DensityFunction : public virtual DomainDecomposition {

public:

    //!
    //! \brief  Used to keep track of the state of overlapping halo cells between MPI ranks.
    //!
    //! This value should be marked \e true whenever the values stored in spatial cells adjacent to the MPI
    //! sub-domain boundary are changed without communicating this change to neighboring MPI ranks.
    //!
    //! Additionally, this value should be checked before reading values from halo cells. If this value is
    //! \e true, then halo cells should be synchronized between MPI ranks by, e.g., calling
    //! DensityFunction::SynchronizeHalos().
    //!
    bool halo_cells_dirty;


    //========================================================================================================
    //=== CONSTRUCTORS, DESTRUCTOR, AND RECONFIGURATION ROUTINES =============================================
    //========================================================================================================

    //--------------------------------------------------------------------------------------------------------
    //! \brief  Constructs an empty Abstract::DensityFunction object with zero-initialized parameters.
    //--------------------------------------------------------------------------------------------------------
    DensityFunction( void );

    //--------------------------------------------------------------------------------------------------------
    //! \brief  Delete copy constructor.
    //--------------------------------------------------------------------------------------------------------
    DensityFunction( const DensityFunction & ) = delete;

    //--------------------------------------------------------------------------------------------------------
    //! \brief  Constructs an Abstract::DensityFunction object with the given set of parameters.
    //--------------------------------------------------------------------------------------------------------
    DensityFunction( const DomainDecomposition & spatial_params );

    //--------------------------------------------------------------------------------------------------------
    //! \brief  Destructor for class Abstract::DensityFunction.
    //--------------------------------------------------------------------------------------------------------
    virtual ~DensityFunction( void );

    //--------------------------------------------------------------------------------------------------------
    //! \brief  Reconfigures an Abstract::DensityFunction object with a different set of parameters.
    //--------------------------------------------------------------------------------------------------------
    virtual DensityFunction & Reconfigure( const DomainDecomposition & spatial_params );


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
    //! \brief  Returns true if the objects have matching DomainDecomposition parameters and the same number
    //!         of degrees of freedom per cell.
    //--------------------------------------------------------------------------------------------------------
    static bool AreCompatible( const DensityFunction &, const DensityFunction & );

    //--------------------------------------------------------------------------------------------------------
    //! \brief  Computes the \f$ \ell_1 \f$ distance between the vectors stored at DensityFunction::density,
    //!         excluding ghost cells.
    //--------------------------------------------------------------------------------------------------------
    static double l1Dist( const DensityFunction &, const DensityFunction & );

    //--------------------------------------------------------------------------------------------------------
    //! \brief  Computes the \f$ \ell_2 \f$ distance between the vectors stored at DensityFunction::density,
    //!         excluding ghost cells.
    //--------------------------------------------------------------------------------------------------------
    static double l2Dist( const DensityFunction &, const DensityFunction & );

    //--------------------------------------------------------------------------------------------------------
    //! \brief  Computes the \f$ \ell_{\infty} \f$ distance between the vectors stored at
    //!         DensityFunction::density, excluding ghost cells.
    //--------------------------------------------------------------------------------------------------------
    static double linfDist( const DensityFunction &, const DensityFunction & );

    //--------------------------------------------------------------------------------------------------------
    //! \brief  Performs a deep copy of one object to another.
    //--------------------------------------------------------------------------------------------------------
    static void Copy(
        DensityFunction &,
        const DensityFunction &,
        const OpDomain = OpDomain::All
    );

    //--------------------------------------------------------------------------------------------------------
    //! \brief  \f$ \vec{y} \gets \beta \vec{y} + \alpha \vec{x} \f$.
    //--------------------------------------------------------------------------------------------------------
    static void AXPY(
        const double,
        const double,
        const DensityFunction &,
        DensityFunction &
    );


    //========================================================================================================
    //=== ADDITIONAL MEMBER FUNCTIONS ========================================================================
    //========================================================================================================

    //--------------------------------------------------------------------------------------------------------
    //! \brief  Prints a summary of the object configuration to the logging interface.
    //--------------------------------------------------------------------------------------------------------
    void Print( const std::string & = "  " ) const override;

    //--------------------------------------------------------------------------------------------------------
    //! \brief  Deallocates memory at all internal pointers and zeros all parameters of an
    //!         Abstract::DensityFunction object.
    //--------------------------------------------------------------------------------------------------------
    DensityFunction & Zero( void ) override;

    //--------------------------------------------------------------------------------------------------------
    //! \brief  Zeros the array of coefficients at the pointer Abstract::DensityFunction::density.
    //--------------------------------------------------------------------------------------------------------
    virtual DensityFunction & ZeroDensity( const OpDomain = OpDomain::All );

    //--------------------------------------------------------------------------------------------------------
    //! \brief  Synchronizes the values of all halo cells between MPI ranks and across periodic boundaries.
    //--------------------------------------------------------------------------------------------------------
    virtual void SynchronizeHalos( void );

    //--------------------------------------------------------------------------------------------------------
    //! \brief  Returns the number of degrees of freedom per spatial mesh element. This includes spatial,
    //!         temporal, and angular degrees of freedom (when applicable).
    //--------------------------------------------------------------------------------------------------------
    virtual int64_t DOFPerCell( void ) const;

    //--------------------------------------------------------------------------------------------------------
    //! \brief  Returns the total number of degrees of freedom used in the underlying local representation of
    //!         the object, either with or without including halo elements.
    //--------------------------------------------------------------------------------------------------------
    virtual int64_t GetLocalArrayDim( const bool with_halo = false ) const;

    //--------------------------------------------------------------------------------------------------------
    //! \brief  Returns the total number of degrees of freedom used in the underlying global representation of
    //!         the object, either with or without including halo elements.
    //--------------------------------------------------------------------------------------------------------
    virtual int64_t GetGlobalArrayDim( const bool with_halo = false ) const;

    //--------------------------------------------------------------------------------------------------------
    //! \brief  Returns the stride in the global array between subsequent spatial cells, moving in the
    //!         specified dimension.
    //--------------------------------------------------------------------------------------------------------
    virtual int64_t CellStride( const int64_t ) const = 0;

    //--------------------------------------------------------------------------------------------------------
    //! \brief  Returns a pointer to the location within Abstract::DensityFunction::density where the
    //!         coefficients for the given spatial cell begin.
    //--------------------------------------------------------------------------------------------------------
    virtual double * PointerAtCell(
            const int64_t i
    # if SPACE_DIMS >= 2
        ,   const int64_t j
    # endif
    # if SPACE_DIMS == 3
        ,   const int64_t k
    # endif
    ) const = 0;

    //--------------------------------------------------------------------------------------------------------
    //! \brief  Returns the state of \c Abstract::DensityFunction::halo_cells_dirty.
    //--------------------------------------------------------------------------------------------------------
    virtual bool HalosAreDirty( void ) const;

    //--------------------------------------------------------------------------------------------------------
    //! \brief  Marks all halo state variables as dirty.
    //--------------------------------------------------------------------------------------------------------
    virtual void MarkAllHalosDirty( void );

    //--------------------------------------------------------------------------------------------------------
    //! \brief  Marks all halo state variables as clean.
    //--------------------------------------------------------------------------------------------------------
    virtual void MarkAllHalosClean( void );

    //--------------------------------------------------------------------------------------------------------
    //! \brief  Scales the coefficients of the object by the given scalar.
    //--------------------------------------------------------------------------------------------------------
    virtual void Scale( const double ) const;

    //--------------------------------------------------------------------------------------------------------
    //! \brief  Computes the \f$ \ell_1 \f$ norm of the vector stored at DensityFunction::Density, excluding
    //!         any ghost and/or halo values.
    //--------------------------------------------------------------------------------------------------------
    virtual double l1Norm( void ) const;

    //--------------------------------------------------------------------------------------------------------
    //! \brief  Computes the \f$ \ell_2 \f$ norm of the vector stored at DensityFunction::Density, excluding
    //!         any ghost and/or halo values.
    //--------------------------------------------------------------------------------------------------------
    virtual double l2Norm( void ) const;

    //--------------------------------------------------------------------------------------------------------
    //! \brief  Computes the \f$ \ell_{\infty} \f$ norm of the vector stored at DensityFunction::Density,
    //!         excluding any ghost and/or halo values.
    //--------------------------------------------------------------------------------------------------------
    virtual double linfNorm( void ) const;

    //--------------------------------------------------------------------------------------------------------
    //! \brief  Performs a deep copy of one object to another.
    //--------------------------------------------------------------------------------------------------------
    void Copy(
        const DensityFunction &,
        const OpDomain = OpDomain::All
    );


    //========================================================================================================
    //=== ROUTINES FOR INTERFACING WITH EXTERNAL LIBRARIES ===================================================
    //========================================================================================================

protected:

    //--------------------------------------------------------------------------------------------------------
    //! \brief  Updates the values in an external object with an operation of the form
    //!         \f$ \vec{y} \gets \beta \vec{y} + \alpha \vec{x} \f$ where \f$ \vec{x} \f$ is \c this.
    //--------------------------------------------------------------------------------------------------------
    template< typename ViewType >
    void PackPointer(
        ViewType &,
        const double,
        const double
    ) const;

    //--------------------------------------------------------------------------------------------------------
    //! \brief  Unpacks the contents of the provided external data structure into the present object.
    //--------------------------------------------------------------------------------------------------------
    template< typename ViewType >
    void UnpackPointer( const ViewType & src );

public:

    //--------------------------------------------------------------------------------------------------------
    //! \brief  Updates the values in an external vector object with an operation of the form
    //!         \f$ \vec{y} \gets \beta \vec{y} + \alpha \vec{x} \f$ where \f$ \vec{x} \f$ is \c this.
    //--------------------------------------------------------------------------------------------------------
    template<
        typename VectorType,
        typename ScalarType = double
    >
    void PackExternalVector(
        VectorType &,
        const ScalarType,
        const ScalarType,
        const size_t = 0
    ) const = delete;

    //--------------------------------------------------------------------------------------------------------
    //! \brief  Unpacks the coefficients stored in an external vector object into the given object.
    //--------------------------------------------------------------------------------------------------------
    template< typename VectorType >
    void UnpackExternalVector(
        const VectorType &,
        const size_t = 0
    ) = delete;


protected:

    //========================================================================================================
    //=== INTERNAL VARIABLES =================================================================================
    //========================================================================================================

    //!
    //! Pointer for storing elements of representation of density function.
    //!
    double * density;

    //!
    //! Number of bytes allocated at the pointer \pp{density}.
    //!
    int64_t sizeof_density;

    //!
    //! Number of doubles allocated at the pointer \pp{density}.
    //!
    int64_t dimof_density;

    //!
    //! \brief  Number of degrees of freedom per mesh cell in object.
    //!
    //! \attention  This variable is owned by the Abstract::DensityFunction class, but it is expected to be
    //!             managed (i.e., properly initialized and updated as necessary) by derived classes.
    //!
    int64_t DOF_per_cell;


    //========================================================================================================
    //=== PROTECTED MEMBER FUNCTIONS =========================================================================
    //========================================================================================================

    //--------------------------------------------------------------------------------------------------------
    //! \brief  Sets appropriate values for object parameters.
    //--------------------------------------------------------------------------------------------------------
    virtual void SetDensityDimensions( void );

    //--------------------------------------------------------------------------------------------------------
    //! \brief  Sets parameters and allocates memory at internal pointers.
    //--------------------------------------------------------------------------------------------------------
    virtual void Allocate( void );

    //--------------------------------------------------------------------------------------------------------
    //! \brief  Deallocates memory at all internal pointers.
    //--------------------------------------------------------------------------------------------------------
    virtual void Deallocate( void );


private:

    //========================================================================================================
    //=== PRIVATE INDEXING FUNCTIONS =========================================================================
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
        ,   const int64_t l
        ,   const char * const func
    ) const;

# endif // if defined (STRICT_CHECK)

    //--------------------------------------------------------------------------------------------------------
    //! \brief  Returns the one-dimensional index for the specified coefficient in the array
    //!         Abstract::DensityFunction::density.
    //--------------------------------------------------------------------------------------------------------
    inline size_t Index (
            const int64_t i
    # if SPACE_DIMS >= 2 || defined (DOXYCOMPILE)
        ,   const int64_t j
    # endif
    # if SPACE_DIMS == 3
        ,   const int64_t k
    # endif

        ,   const int64_t l
    ) const;

// Make IndexNoGC public for computation of diffusion matrix. See PETScPeierlsSolver<>::DSA_ComputeMatrix.
# if SPACE_DIMS == 1
public:
# endif

    //--------------------------------------------------------------------------------------------------------
    //! \brief  Returns the one-dimensional index for the specified coefficient in an array corresponding to
    //!         Abstract::DensityFunction::density, but without ghost cells.
    //--------------------------------------------------------------------------------------------------------
    inline size_t IndexNoGC (
            const int64_t i
    # if SPACE_DIMS >= 2 || defined (DOXYCOMPILE)
        ,   const int64_t j
    # endif
    # if SPACE_DIMS == 3
        ,   const int64_t k
    # endif
        ,   const int64_t l
    ) const;

# if SPACE_DIMS == 1
private:
# endif

    //--------------------------------------------------------------------------------------------------------
    //! \brief  Returns a reference to the coefficient inside the array Abstract::DensityFunction::density.
    //--------------------------------------------------------------------------------------------------------
    inline double & operator() (
            const int64_t i
    # if SPACE_DIMS >= 2 || defined (DOXYCOMPILE)
        ,   const int64_t j
    # endif
    # if SPACE_DIMS == 3
        ,   const int64_t k
    # endif
        ,   const int64_t l
    ) const;

};


//============================================================================================================
//=== DECLARATIONS OF TEMPLATE SPECIALIZATIONS ===============================================================
//============================================================================================================

# if defined (ENABLE_PETSC) || defined (DOXYCOMPILE)


//------------------------------------------------------------------------------------------------------------
//! \brief  Updates the values in PETScVec object with an operation of the form
//!         \f$ \vec{y} \gets \beta \vec{y} + \alpha \vec{x} \f$ where \f$ \vec{x} \f$ is \c this.
//------------------------------------------------------------------------------------------------------------
template<>
void DensityFunction::PackExternalVector< PETScVec, double >(
    PETScVec &,
    const double,
    const double,
    const size_t
) const;

//------------------------------------------------------------------------------------------------------------
//! \brief  Unpacks the coefficients stored in a PETScVec object into the given object.
//------------------------------------------------------------------------------------------------------------
template<>
void DensityFunction::UnpackExternalVector<PETScVec>(
    const PETScVec &,
    const size_t
);


# endif // if defined (ENABLE_PETSC)

# if defined (ENABLE_AZTEC_EPETRA) || defined (DOXYCOMPILE)


//------------------------------------------------------------------------------------------------------------
//! \brief  Updates the values in an Epetra_MultiVector object with an operation of the form
//!         \f$ \vec{y} \gets \beta \vec{y} + \alpha \vec{x} \f$ where \f$ \vec{x} \f$ is \c this.
//------------------------------------------------------------------------------------------------------------
template<>
void DensityFunction::PackExternalVector< Epetra_MultiVector, double >(
    Epetra_MultiVector &,
    const double,
    const double,
    const size_t
) const;

//------------------------------------------------------------------------------------------------------------
//! \brief  Unpacks the coefficients stored in an Epetra_MultiVector into the given object.
//------------------------------------------------------------------------------------------------------------
template<>
void DensityFunction::UnpackExternalVector<Epetra_MultiVector>(
    const Epetra_MultiVector &,
    const size_t
);


# endif // if defined (ENABLE_AZTEC_EPETRA)

# if defined (ENABLE_BELOS_TPETRA) || defined (DOXYCOMPILE)


//------------------------------------------------------------------------------------------------------------
//! \brief  Packs a TpetraVec with the coefficients stored in the given object.
//------------------------------------------------------------------------------------------------------------
template<>
void DensityFunction::PackExternalVector< TpetraVec, TpetraScalar >(
    TpetraVec &,
    const TpetraScalar,
    const TpetraScalar,
    const size_t
) const;

//------------------------------------------------------------------------------------------------------------
//! \brief  Unpacks the coefficients stored in a TpetraVec into the given object.
//------------------------------------------------------------------------------------------------------------
template<>
void DensityFunction::UnpackExternalVector<TpetraVec>(
    const TpetraVec &,
    const size_t
);


# endif // if defined (ENABLE_BELOS_TPETRA)


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
//!     CheckBounds(i,j,l,__func__)
//! \endcode
//!
//! \param[in]      i           Spatial cell index (including ghost cells) in the \f$ x_1 \f$ dimension.
//! \param[in]      j           Spatial cell index (including ghost cells) in the \f$ x_2 \f$ dimension.
//! \param[in]      l           Index of local degree of freedom in given mesh cell.
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

    ,   const int64_t l
    ,   const char * const func

) const {

    if ( i < 0 || i >= this->nx(0) +2 ) {

        std::string error_message
            = "Value " + std::to_string(i) + " for index i in " + std::string(func)
                + " outsize permissible range [ 0, " + std::to_string( this->nx(0) +2 ) + " ).\n";

        PRINT_ERROR( error_message.c_str() )
        throw std::out_of_range( error_message );
    }

    if ( l < 0 || l >= this->DOF_per_cell ) {

        std::string error_message
            = "Value " + std::to_string(l) + " for index l in " + std::string(func)
                + " outsize permissible range [ 0, " + std::to_string( this->DOF_per_cell ) + " ).\n";

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

# endif // if SPACE_DIMS >= 2
# if SPACE_DIMS == 3

    if ( k < 0 || k >= this->nx(2) +2 ) {

        std::string error_message
            = "Value " + std::to_string(k) + " for index k in " + std::string(func)
                + " outsize permissible range [ 0, " + std::to_string( this->nx(2) +2 ) + " ).\n";

        PRINT_ERROR( error_message.c_str() )
        throw std::out_of_range( error_message );
    }

# endif // if SPACE_DIMS == 3
}

# endif // if defined (STRICT_CHECK)


//------------------------------------------------------------------------------------------------------------
//! \brief  Returns the one-dimensional index for the specified coefficient in the array
//!         Abstract::DensityFunction::density.
//!
//! \param[in]      i           Spatial cell index (including ghost cells) in the \f$ x_1 \f$ dimension.
//! \param[in]      j           Spatial cell index (including ghost cells) in the \f$ x_2 \f$ dimension.
//! \param[in]      l           Index of local degree of freedom in given mesh cell.
//!
//! \return     Appropriate index in Abstract::DensityFunction::density array for specified coefficient.
//!
//! \see    Abstract::DensityFunction::IndexNoGC()
//! \see    Abstract::DensityFunction::operator()()
//------------------------------------------------------------------------------------------------------------
inline size_t DensityFunction::Index (

        const int64_t i
# if SPACE_DIMS >= 2 || defined (DOXYCOMPILE)
    ,   const int64_t j
# endif
# if SPACE_DIMS == 3
    ,   const int64_t k
# endif

    ,   const int64_t l

) const {

# if defined (STRICT_CHECK)

    CheckBounds(
            # if SPACE_DIMS == 1
                i,l
            # elif SPACE_DIMS == 2
                i,j,l
            # elif SPACE_DIMS == 3
                i,j,k,l
            # endif
                , __func__ );

# endif // if defined (STRICT_CHECK)

    return (
        l + (this->DOF_per_cell)*
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
    );
}


//------------------------------------------------------------------------------------------------------------
//! \brief  Returns the one-dimensional index for the specified coefficient in an array corresponding to
//!         Abstract::DensityFunction::density, but without ghost cells.
//!
//! Used for packing coefficients into external vectors for solvers.
//!
//! \param[in]      i           Spatial cell index (including ghost cells) in the \f$ x_1 \f$ dimension.
//! \param[in]      j           Spatial cell index (including ghost cells) in the \f$ x_2 \f$ dimension.
//! \param[in]      l           Index of local degree of freedom in given mesh cell.
//!
//! \return     Returns the appropriate index for the specified coefficient.
//!
//! \see    Abstract::DensityFunction::Index()
//------------------------------------------------------------------------------------------------------------
inline size_t DensityFunction::IndexNoGC (

        const int64_t i
# if SPACE_DIMS >= 2 || defined (DOXYCOMPILE)
    ,   const int64_t j
# endif
# if SPACE_DIMS == 3
    ,   const int64_t k
# endif

    ,   const int64_t l

) const {

    return (
        l + (this->DOF_per_cell)*
        (
        # if SPACE_DIMS == 3
            (k - 1) + (this->nx(2))*
        # endif
            (
            # if SPACE_DIMS >= 2
                (j - 1) + (this->nx(1))*
            # endif
                (
                    (i - 1)
                )
            )
        )
    );
}


//------------------------------------------------------------------------------------------------------------
//! \brief  Returns a reference to the coefficient inside the array Abstract::DensityFunction::density.
//!
//! \param[in]      i           Spatial cell index (including ghost cells) in the \f$ x_1 \f$ dimension.
//! \param[in]      j           Spatial cell index (including ghost cells) in the \f$ x_2 \f$ dimension.
//! \param[in]      l           Index of local degree of freedom in given mesh cell.
//!
//! \return     Reference to the specified coefficient.
//!
//! \see    Abstract::DensityFunction::Index()
//------------------------------------------------------------------------------------------------------------
inline double & DensityFunction::operator() (

        const int64_t i
# if SPACE_DIMS >= 2 || defined (DOXYCOMPILE)
    ,   const int64_t j
# endif
# if SPACE_DIMS == 3
    ,   const int64_t k
# endif

    ,   const int64_t l

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
                                (i,l)
                            # elif SPACE_DIMS == 2
                                (i,j,l)
                            # elif SPACE_DIMS == 3
                                (i,j,k,l)
                            # endif
                        ];
}


} // namespace Abstract


# endif // ifndef __ABSTRACT__DENSITY_FUNCTION_HPP__
