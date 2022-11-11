//------------------------------------------------------------------------------------------------------------
//! \file   objects/Abstract/OrdinateFlux.hpp
//! \brief  Header for Abstract::OrdinateFlux class.
//!
//! \author Michael Crockatt
//! \date   February 2018
//------------------------------------------------------------------------------------------------------------

# ifndef __ABSTRACT__ORDINATE_FLUX_HPP__
# define __ABSTRACT__ORDINATE_FLUX_HPP__


# include "objects/Abstract/DensityFunction.hpp"
# include "utils/Quadrule/OrdinateSet.hpp"


namespace Abstract {


//------------------------------------------------------------------------------------------------------------
//! \brief  Abstract class used as a base class for objects used to store spatially-varying scalar-valued
//!         angular flux functions with discrete ordinates angular discretizations.
//!
//------------------------------------------------------------------------------------------------------------
class OrdinateFlux : public virtual Abstract::DensityFunction, public Quadrule::OrdinateSet {

public:

    //!
    //! \brief  Used to keep track of the state of overlapping halo cells between MPI ranks _in the upwind
    //!         direction_.
    //!
    //! This value should be marked true whenever the halo cells in the upwind direction of each ordinate are
    //! changed without synchronization between MPI ranks by, e.g., calling OrdinateFlux::SynchronizeHalos().
    //!
    bool upwind_halo_cells_dirty;


    //========================================================================================================
    //=== CONSTRUCTORS AND DESTRUCTOR ========================================================================
    //========================================================================================================

    //--------------------------------------------------------------------------------------------------------
    //! \brief  Constructs an empty Abstract::OrdinateFlux object with zero-initialized parameters.
    //--------------------------------------------------------------------------------------------------------
    OrdinateFlux( void );

    //--------------------------------------------------------------------------------------------------------
    //! \brief  Delete copy constructor.
    //--------------------------------------------------------------------------------------------------------
    OrdinateFlux( const OrdinateFlux & ) = delete;

    //--------------------------------------------------------------------------------------------------------
    //! \brief  Constructs an Abstract::OrdinateFlux object with the given parameters.
    //--------------------------------------------------------------------------------------------------------
    OrdinateFlux(
        const int64_t ang_order,
        const bool symmetric_reduce = false,
        const OrdinateType ordinate_type
            # if SPACE_DIMS == 1
                = OrdinateType::GaussLegendre
            # elif SPACE_DIMS >= 2
                = OrdinateType::ChebyshevLegendre
            # endif
    );

    //
    // Acknowledge use of Reconfigure routines from base class.
    //
    using Abstract::DensityFunction::Reconfigure;

    //--------------------------------------------------------------------------------------------------------
    //! \brief  Reconfigures an Abstract::OrdinateFlux object with a different set of angular parameters.
    //--------------------------------------------------------------------------------------------------------
    OrdinateFlux & Reconfigure(
        const int64_t ang_order,
        const bool symmetric_reduce = false,
        const OrdinateType ordinate_type
            # if SPACE_DIMS == 1
                = OrdinateType::GaussLegendre
            # elif SPACE_DIMS >= 2
                = OrdinateType::ChebyshevLegendre
            # endif
    ) override;

    //--------------------------------------------------------------------------------------------------------
    //! \brief  Reconfigures an Abstract::OrdinateFlux object with a different set of angular parameters.
    //--------------------------------------------------------------------------------------------------------
    OrdinateFlux & Reconfigure( const OrdinateSet & ordinate_set ) override;


    //========================================================================================================
    //=== OPERATOR OVERLOADS =================================================================================
    //========================================================================================================

    //--------------------------------------------------------------------------------------------------------
    //! \brief  Delete copy assignment operator.
    //--------------------------------------------------------------------------------------------------------
    OrdinateFlux & operator=( const OrdinateFlux & ) = delete;

    //--------------------------------------------------------------------------------------------------------
    //! \brief  Overload equality operator with pointer equality.
    //!
    //! \see    Abstract::OrdinateFlux::operator!=()
    //--------------------------------------------------------------------------------------------------------
    bool operator==( const OrdinateFlux & that ) const {  return this == &that;  }

    //--------------------------------------------------------------------------------------------------------
    //! \brief  Overload inequality operator using negation of equality operator.
    //!
    //! \see    Abstract::OrdinateFlux::operator==()
    //--------------------------------------------------------------------------------------------------------
    bool operator!=( const OrdinateFlux & that ) const {  return !(*this == that);  }


    //========================================================================================================
    //=== STATIC MEMBER FUNCTIONS ============================================================================
    //========================================================================================================

    //
    // Acknowledge use of Abstract::DensityFunction::AreCompatible.
    //
    using Abstract::DensityFunction::AreCompatible;

    //--------------------------------------------------------------------------------------------------------
    //! \brief  Performs a deep copy of one object to another.
    //--------------------------------------------------------------------------------------------------------
    static void Copy(
        OrdinateFlux &,
        const OrdinateFlux &,
        const OpDomain = OpDomain::All
    );

    //--------------------------------------------------------------------------------------------------------
    //! \brief  \f$ \vec{y} \gets \beta \vec{y} + \alpha \vec{x} \f$.
    //--------------------------------------------------------------------------------------------------------
    static void AXPY(
        const double,
        const double,
        const OrdinateFlux &,
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
    //! \brief  Synchronizes the values of all halo cells between MPI ranks and across periodic boundaries.
    //--------------------------------------------------------------------------------------------------------
    void SynchronizeHalos( void ) override;

    //--------------------------------------------------------------------------------------------------------
    //! \brief  Returns the number of spatial degrees of freedom per mesh element. This includes spatial and
    //!         temporal, but not angular degrees of freedom; i.e., this is the number of degrees of freedom
    //!         for each mesh cell and ordinate.
    //--------------------------------------------------------------------------------------------------------
    virtual int64_t SpatialDOFPerCell( void ) const;

    //--------------------------------------------------------------------------------------------------------
    //! \brief  Updates ghost regions on reflecting boundaries.
    //--------------------------------------------------------------------------------------------------------
    virtual OrdinateFlux & ReflectBoundaries( void ) = 0;

    //--------------------------------------------------------------------------------------------------------
    //! \brief  Returns a pointer to the location within DensityFunction::density where the coefficients for
    //!         the given spatial cell begin.
    //--------------------------------------------------------------------------------------------------------
    virtual double * PointerAtOrdinateCell(
            const int64_t q
        ,   const int64_t i
    # if SPACE_DIMS >= 2
        ,   const int64_t j
    # endif
    # if SPACE_DIMS == 3
        ,   const int64_t k
    # endif
    ) const = 0;

# if SPACE_DIMS == 2

    //--------------------------------------------------------------------------------------------------------
    //! \brief  Returns a pointer to the location within DensityFunction::density where the coefficients for
    //!         the given spatial cell begin.
    //--------------------------------------------------------------------------------------------------------
    virtual double * PointerAtQuadrantCell(
        const int64_t quad,
        const int64_t i,
        const int64_t j
    ) const = 0;

# elif SPACE_DIMS == 3

    //--------------------------------------------------------------------------------------------------------
    //! \brief  Returns a pointer to the location within DensityFunction::density where the coefficients for
    //!         the given spatial cell begin.
    //--------------------------------------------------------------------------------------------------------
    virtual double * PointerAtOctantCell(
        const int64_t oct,
        const int64_t i,
        const int64_t j,
        const int64_t k
    ) const = 0;

# endif // if SPACE_DIMS == ?

    //--------------------------------------------------------------------------------------------------------
    //! \brief  Returns the state of \c Abstract::OrdinateFlux::upwind_halo_cells_dirty.
    //--------------------------------------------------------------------------------------------------------
    virtual bool UpwindHalosAreDirty( void ) const;

    //--------------------------------------------------------------------------------------------------------
    //! \brief  Marks all halo state variables as dirty.
    //--------------------------------------------------------------------------------------------------------
    void MarkAllHalosDirty( void ) override;

    //--------------------------------------------------------------------------------------------------------
    //! \brief  Marks all halo state variables as clean.
    //--------------------------------------------------------------------------------------------------------
    void MarkAllHalosClean( void ) override;

    //--------------------------------------------------------------------------------------------------------
    //! \brief  Performs a deep copy of one object to another.
    //--------------------------------------------------------------------------------------------------------
    void Copy(
        const OrdinateFlux &,
        const OpDomain = OpDomain::All
    );


    //========================================================================================================
    //=== ROUTINES FOR INTERFACING WITH EXTERNAL LIBRARIES ===================================================
    //========================================================================================================

    //--------------------------------------------------------------------------------------------------------
    //! \brief  Unpacks the coefficients stored in an external vector object into the given object.
    //--------------------------------------------------------------------------------------------------------
    template< typename VectorType >
    void UnpackExternalVector(
        const VectorType &,
        const size_t = 0
    );

};


//========================================================================================================
//=== TEMPLATE DEFINITIONS ===============================================================================
//========================================================================================================


//--------------------------------------------------------------------------------------------------------
//! \brief  Unpacks the coefficients stored in an external vector object into the given object.
//!
//! This function template is used to override DensityFunction::UnpackExternalVector to ensure that
//! OrdinateFlux::ReflectBoundaries is called on the unpacked object.
//!
//! \param[in]      source  External vector object to unpack coefficients from.
//! \param[in]      column  (optional) <br>
//!                         Column of external vector object to unpack. Defaults to 0 (first column).
//!
//! \see    DensityFunction::PackExternalVector()
//--------------------------------------------------------------------------------------------------------
template< typename VectorType >
void OrdinateFlux::UnpackExternalVector (

    const VectorType & source,
    const size_t column         // = 0
) {

    PRINT_STATUS( "Executing Abstract::OrdinateFlux::%s.\n", __func__ )

    this->Abstract::DensityFunction::UnpackExternalVector( source, column );
    this->ReflectBoundaries();
}


} // namespace Abstract


# endif // ifndef __ABSTRACT__ORDINATE_FLUX_HPP__
