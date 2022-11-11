//------------------------------------------------------------------------------------------------------------
//! \file   objects/SpatialDomain.hpp
//! \brief  Header for SpatialDomain class.
//!
//! \author Michael M. Crockatt
//! \date   November 2017
//------------------------------------------------------------------------------------------------------------

# ifndef __SPATIAL_DOMAIN_HPP__
# define __SPATIAL_DOMAIN_HPP__


# include <cstddef>
# include <cstdint>
# include <cstdio>
# include <initializer_list>
# include <string>

# include "utils/CLog.hpp"
# include "utils/global.hpp"


//------------------------------------------------------------------------------------------------------------
//! \brief  Stores mesh information for a spatial domain.
//------------------------------------------------------------------------------------------------------------
class SpatialDomain {

public:

    //========================================================================================================
    //=== CONSTRUCTORS AND DESTRUCTOR ========================================================================
    //========================================================================================================

    //--------------------------------------------------------------------------------------------------------
    //! \brief  Constructs a SpatialDomain with default (and generally physically invalid) values.
    //--------------------------------------------------------------------------------------------------------
    SpatialDomain( void );

    //--------------------------------------------------------------------------------------------------------
    //! \brief  Creates a SpatialDomain object using the provided values.
    //--------------------------------------------------------------------------------------------------------
    SpatialDomain(
        const int64_t (&) [SPACE_DIMS],
        const double (&) [SPACE_DIMS],
        const double (&) [SPACE_DIMS]
    );

    //--------------------------------------------------------------------------------------------------------
    //! \brief  Creates a SpatialDomain object using the provided values.
    //--------------------------------------------------------------------------------------------------------
    SpatialDomain(
        std::initializer_list<int64_t>,
        std::initializer_list<double>,
        std::initializer_list<double>
    );

    //--------------------------------------------------------------------------------------------------------
    //! \brief  Attempts to construct a SpatialDomain object using values read from a ParameterList.
    //--------------------------------------------------------------------------------------------------------
    SpatialDomain( const ParameterList & input_list );

    //--------------------------------------------------------------------------------------------------------
    //! \brief  Default POD copy constructor.
    //--------------------------------------------------------------------------------------------------------
    SpatialDomain( const SpatialDomain & ) = default;

    //--------------------------------------------------------------------------------------------------------
    //! \brief  Default POD destructor.
    //--------------------------------------------------------------------------------------------------------
    ~SpatialDomain( void ) = default;


    //========================================================================================================
    //=== OPERATOR OVERLOADS =================================================================================
    //========================================================================================================

    //--------------------------------------------------------------------------------------------------------
    //! \brief  Default POD copy assignment operator.
    //--------------------------------------------------------------------------------------------------------
    SpatialDomain & operator=( const SpatialDomain & ) = default;


    //========================================================================================================
    //=== STATIC MEMBER FUNCTIONS ============================================================================
    //========================================================================================================

    //--------------------------------------------------------------------------------------------------------
    //! \brief  Two SpatialDomain objects "are matching" if all of their parameters are equal.
    //--------------------------------------------------------------------------------------------------------
    static bool AreMatching(
        const SpatialDomain &,
        const SpatialDomain &
    );


    //========================================================================================================
    //=== ADDITIONAL MEMBER FUNCTIONS ========================================================================
    //========================================================================================================

    //--------------------------------------------------------------------------------------------------------
    //! \brief  Returns a const reference to an element of SpatialDomain::nx_var.
    //!
    //! \param[in]  dim     Index of spatial dimension to return value for.
    //!
    //! \return     Returns the value of SpatialDomain::nx_var[dim].
    //--------------------------------------------------------------------------------------------------------
    inline const int64_t & nx (

        const int64_t dim

    ) const {

    # if defined (STRICT_CHECK)
        CHECK_DIM_VALUE(dim)
    # endif

        return this->nx_var[dim];
    }

    //--------------------------------------------------------------------------------------------------------
    //! \brief  Returns a const reference to the array SpatialDomain::nx_var.
    //--------------------------------------------------------------------------------------------------------
    inline auto nx( void ) const -> const int64_t (&) [SPACE_DIMS] {  return this->nx_var;  }

    //--------------------------------------------------------------------------------------------------------
    //! \brief  Returns a const reference to an element of SpatialDomain::ax_var.
    //!
    //! \param[in]  dim     Index of spatial dimension to return value for.
    //!
    //! \return     Returns the value of SpatialDomain::ax_var[dim].
    //--------------------------------------------------------------------------------------------------------
    inline const double & ax (

        const int64_t dim

    ) const {

    # if defined (STRICT_CHECK)
        CHECK_DIM_VALUE(dim)
    # endif

        return this->ax_var[dim];
    }

    //--------------------------------------------------------------------------------------------------------
    //! \brief  Returns a const reference to the array SpatialDomain::ax_var.
    //--------------------------------------------------------------------------------------------------------
    inline auto ax( void ) const -> const double (&) [SPACE_DIMS] {  return this->ax_var;  }

    //--------------------------------------------------------------------------------------------------------
    //! \brief  Returns a const reference to an element of SpatialDomain::bx_var.
    //!
    //! \param[in]  dim     Index of spatial dimension to return value for.
    //!
    //! \return     Returns the value of SpatialDomain::bx_var[dim].
    //--------------------------------------------------------------------------------------------------------
    inline const double & bx (

        const int64_t dim

    ) const {

    # if defined (STRICT_CHECK)
        CHECK_DIM_VALUE(dim)
    # endif

        return this->bx_var[dim];
    }

    //--------------------------------------------------------------------------------------------------------
    //! \brief  Returns a const reference to the array SpatialDomain::bx_var.
    //--------------------------------------------------------------------------------------------------------
    inline auto bx( void ) const -> const double (&) [SPACE_DIMS] {  return this->bx_var;  }

    //--------------------------------------------------------------------------------------------------------
    //! \brief  Returns a const reference to an element of SpatialDomain::dx_var.
    //!
    //! \param[in]  dim     Index of spatial dimension to return value for.
    //!
    //! \return     Returns the value of SpatialDomain::dx_var[dim].
    //--------------------------------------------------------------------------------------------------------
    inline const double & dx (

        const int64_t dim

    ) const {

    # if defined (STRICT_CHECK)
        CHECK_DIM_VALUE(dim)
    # endif

        return this->dx_var[dim];
    }

    //--------------------------------------------------------------------------------------------------------
    //! \brief  Returns a const reference to the array SpatialDomain::dx_var.
    //--------------------------------------------------------------------------------------------------------
    inline auto dx( void ) const -> const double (&) [SPACE_DIMS] {  return this->dx_var;  }

    //--------------------------------------------------------------------------------------------------------
    //! \brief  Returns a string describing the values stored in the object for printing.
    //--------------------------------------------------------------------------------------------------------
    std::string PrintString( void ) const;

    //--------------------------------------------------------------------------------------------------------
    //! \brief  Reads parameters for a SpatialDomain object from the given file pointer.
    //--------------------------------------------------------------------------------------------------------
    void ReadFromDisk(
    # if defined (ENABLE_MPI)
        MPI_File &,
    # else
        std::FILE * &,
    # endif
        const int64_t version = 0
    );

    //--------------------------------------------------------------------------------------------------------
    //! \brief  Writes the parameters in a SpatialDomain object to the given file pointer.
    //--------------------------------------------------------------------------------------------------------
    void WriteToDisk(
    # if defined (ENABLE_MPI)
        MPI_File &
    # else
        std::FILE * &
    # endif
    ) const;


protected:

    //========================================================================================================
    //=== INTERNAL VARIABLES =================================================================================
    //========================================================================================================

    int64_t nx_var [SPACE_DIMS];    //!< Number of spatial cells in each dimension of mesh.
    double  ax_var [SPACE_DIMS];    //!< Lower limit of spatial domain in each dimension.
    double  bx_var [SPACE_DIMS];    //!< Upper limit of spatial domain in each dimension.
    double  dx_var [SPACE_DIMS];    //!< Mesh size in each spatial dimension.

};


# endif // ifndef __SPATIAL_DOMAIN_HPP__
