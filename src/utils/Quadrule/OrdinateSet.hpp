//------------------------------------------------------------------------------------------------------------
//! \file   utils/Quadrule/OrdinateSet.hpp
//! \brief  Header for Quadrule::OrdinateSet class.
//!
//! \author Michael M. Crockatt
//! \date   December 2017
//------------------------------------------------------------------------------------------------------------

# ifndef __ORDINATE_SET_HPP__
# define __ORDINATE_SET_HPP__


# include <cstdint>
# include <map>
# include <string>
# include <tuple>

# include "utils/CLog.hpp"


namespace Quadrule {


//------------------------------------------------------------------------------------------------------------
//! \brief  Class used for storing discrete ordinates quadrature sets.
//!
//! \todo   Use mu for 1D.
//------------------------------------------------------------------------------------------------------------
class OrdinateSet {

public:

    //--------------------------------------------------------------------------------------------------------
    //! \brief  Specifies the type or family of a discrete ordinate quadrature set.
    //--------------------------------------------------------------------------------------------------------
    enum class OrdinateType {

        //!
        //! String descriptor: \e "none"
        //!
        //! Null value. Not a valid type of ordinate set.
        //!
        None,

        // 1D only.
    # if SPACE_DIMS == 1 || defined (DOXYCOMPILE)

        //!
        //! String descriptor: \e "gauss-legendre"
        //!
        //! Standard Gauss-Legendre quadrature over \f$ [-1,1] \f$.
        //!
        GaussLegendre,

        //!
        //! String descriptor: \e "gauss-lobatto"
        //!
        //! Standard Gauss-Lobatto quadrature over \f$ [-1,1] \f$.
        //!
        GaussLobatto,

        //!
        //! String descriptor: \e "double-gauss"
        //!
        //! Double Gauss-Legendre quadrature over \f$ [-1,1] \f$.
        //!
        DoubleGauss,

        //!
        //! String descriptor: \e "double-radau"
        //!
        //! Double Gauss-Radau quadrature over \f$ [-1,1] \f$.
        //!
        DoubleRadau

    # endif // if SPACE_DIMS == 1

        // Trick to get all ordinate sets to compile into documentation
    # if defined (DOXYCOMPILE)
        ,
    # endif

        // 2D & 3D only.
    # if SPACE_DIMS >= 2 || defined (DOXYCOMPILE)

        //!
        //! String descriptor: \e "chebyshev-legendre"
        //!
        //! Chebyshev-Legendre quadrature over \f$ \mathbb{S}^2 \f$.
        //!
        ChebyshevLegendre,

        //!
        //! String descriptor: \e "spherical-triangle"
        //!
        //! \f$ T_N \f$ quadrature based on a tessellation of \f$ \mathbb{S}^2 \f$ into spherical triangles
        //! \cite Thurgood1992 \cite Thurgood1995 .
        //!
        SphericalTriangle,

    # if defined (ENABLE_LEBEDEV)

        //!
        //! String descriptor: \e "lebedev"
        //!
        //! Lebedev quadrature over \f$ \mathbb{S}^2 \f$.
        //!
        //! \see    Quadrule::ComputeLebedevQuadrature
        //!
        Lebedev

    # endif // if defined (ENABLE_LEBEDEV)

    # endif // if SPACE_DIMS >= 2
    };


    //--------------------------------------------------------------------------------------------------------
    //! \brief  Used to convert OrdinateSet::OrdinateType values to their corresponding string descriptors.
    //--------------------------------------------------------------------------------------------------------
    static std::map< OrdinateSet::OrdinateType, std::string > OrdinateType_to_String;

    //--------------------------------------------------------------------------------------------------------
    //! \brief  Used to convert string descriptors to OrdinateSet::OrdinateType values.
    //--------------------------------------------------------------------------------------------------------
    static std::map< std::string, OrdinateSet::OrdinateType > String_to_OrdinateType;


    //========================================================================================================
    //=== CONSTRUCTORS, DESTRUCTOR, AND ASSOCIATED HELPER ROUTINES ===========================================
    //========================================================================================================

    //--------------------------------------------------------------------------------------------------------
    //! \brief  Constructs an empty OrdinateSet object with zero-initialized parameters.
    //--------------------------------------------------------------------------------------------------------
    OrdinateSet( void );

    //--------------------------------------------------------------------------------------------------------
    //! \brief  Constructs an OrdinateSet object given an OrdinateType and order.
    //--------------------------------------------------------------------------------------------------------
    OrdinateSet( const int64_t ang_order,
                 const bool symmetric_reduce = false,

            # if SPACE_DIMS == 1
                 const OrdinateType ordinate_type = OrdinateType::GaussLegendre
            # elif SPACE_DIMS >= 2
                 const OrdinateType ordinate_type = OrdinateType::ChebyshevLegendre
            # endif
                  );

    //--------------------------------------------------------------------------------------------------------
    //! \brief  Constructs an OrdinateSet object using parameters of a given OrdinateSet object.
    //--------------------------------------------------------------------------------------------------------
    OrdinateSet( const OrdinateSet & );

    //--------------------------------------------------------------------------------------------------------
    //! \brief  Default destructor.
    //--------------------------------------------------------------------------------------------------------
    ~OrdinateSet( void ) = default;

    //--------------------------------------------------------------------------------------------------------
    //! \brief  Reconfigures an OrdinateSet object with the parameters of a given OrdinateSet object.
    //--------------------------------------------------------------------------------------------------------
    virtual OrdinateSet & Reconfigure( const OrdinateSet & );

    //--------------------------------------------------------------------------------------------------------
    //! \brief  Reconfigures an OrdinateSet object with the given set of parameters.
    //--------------------------------------------------------------------------------------------------------
    virtual OrdinateSet & Reconfigure( const int64_t ang_order,
                                       const bool symmetric_reduce = false,

                                    # if SPACE_DIMS == 1
                                       const OrdinateType ordinate_type = OrdinateType::GaussLegendre
                                    # elif SPACE_DIMS >= 2 || defined (DOXYCOMPILE)
                                       const OrdinateType ordinate_type = OrdinateType::ChebyshevLegendre
                                    # endif
                                        );


    //========================================================================================================
    //=== OPERATOR OVERLOADS =================================================================================
    //========================================================================================================

    //--------------------------------------------------------------------------------------------------------
    //! \brief  Reconfigures an OrdinateSet object by copying another.
    //--------------------------------------------------------------------------------------------------------
    OrdinateSet & operator=( const OrdinateSet & ) = default;

    //--------------------------------------------------------------------------------------------------------
    //! \brief  Two OrdinateSet objects are said to be "equal" if they contain quadrature sets of the same
    //!         type and order.
    //--------------------------------------------------------------------------------------------------------
    inline bool operator==( const OrdinateSet & that ) const {

        return      ( this->ordinate_data.type == that.ordinate_data.type )
                &&  ( this->ordinate_data.order == that.ordinate_data.order )
                &&  ( this->ordinate_data.symmetric_reduce == that.ordinate_data.symmetric_reduce );
    }

    //--------------------------------------------------------------------------------------------------------
    //! \brief  Two OrdinateSet objects are said to be "not equal" if they contain quadrature sets of
    //!         different types and/or orders.
    //--------------------------------------------------------------------------------------------------------
    inline bool operator!=( const OrdinateSet & that ) const {  return !( *this == that );  }


    //========================================================================================================
    //=== STATIC MEMBER FUNCTIONS ============================================================================
    //========================================================================================================

    //--------------------------------------------------------------------------------------------------------
    //! \brief  Two OrdinateSet objects "are matching" if they are equal.
    //--------------------------------------------------------------------------------------------------------
    static bool AreMatching( const OrdinateSet & thing1, const OrdinateSet & thing2 ) {

        return thing1 == thing2;
    }

# if SPACE_DIMS >= 2

    //--------------------------------------------------------------------------------------------------------
    //! \brief  Sets an array of boolean values specifying the signs on the ordinates in a given quadrant.
    //--------------------------------------------------------------------------------------------------------
    static void SetSigns( const int64_t idx, bool (& signs) [SPACE_DIMS] );

# endif // if SPACE_DIMS >= 2


    //========================================================================================================
    //=== ADDITIONAL MEMBER FUNCTIONS ========================================================================
    //========================================================================================================

    //--------------------------------------------------------------------------------------------------------
    //! \brief  Returns a const reference to the value of \f$ \xi_q \f$.
    //--------------------------------------------------------------------------------------------------------
    inline const double & xi( const int64_t q ) const {

    # if defined (STRICT_CHECK)

        if ( q < 0 || q >= this->nq() ) {

            std::string error_message = "Index " + std::to_string(q)
                                        + " not in range of valid ordinate indices [0,"
                                        + std::to_string( this->nq() ) + ").\n";

            PRINT_ERROR( error_message.c_str() )
            throw std::out_of_range( error_message );
        }

    # endif // if defined (STRICT_CHECK)

        return this->ordinate_data.xi_ptr[q];
    }


# if SPACE_DIMS >= 2 || defined (DOXYCOMPILE)

    //--------------------------------------------------------------------------------------------------------
    //! \brief  Returns a const reference to the value of \f$ \eta_q \f$.
    //--------------------------------------------------------------------------------------------------------
    inline const double & eta( const int64_t q ) const {

    # if defined (STRICT_CHECK)

        if ( q < 0 || q >= this->nq() ) {

            std::string error_message = "Index " + std::to_string(q)
                                        + " not in range of valid ordinate indices [0,"
                                        + std::to_string( this->nq() ) + ").\n";

            PRINT_ERROR( error_message.c_str() )
            throw std::out_of_range( error_message );
        }

    # endif // if defined (STRICT_CHECK)

        return this->ordinate_data.eta_ptr[q];
    }

    //--------------------------------------------------------------------------------------------------------
    //! \brief  Returns a const reference to the value of \f$ \mu_q \f$.
    //--------------------------------------------------------------------------------------------------------
    inline const double & mu( const int64_t q ) const {

    # if defined (STRICT_CHECK)

        if ( q < 0 || q >= this->nq() ) {

            std::string error_message = "Index " + std::to_string(q)
                                        + " not in range of valid ordinate indices [0,"
                                        + std::to_string( this->nq() ) + ").\n";

            PRINT_ERROR( error_message.c_str() )
            throw std::out_of_range( error_message );
        }

    # endif // if defined (STRICT_CHECK)

        return this->ordinate_data.mu_ptr[q];
    }

    //--------------------------------------------------------------------------------------------------------
    //! \brief  Returns a const reference to the value of \f$ \phi_q \f$.
    //--------------------------------------------------------------------------------------------------------
    inline const double & OrdinatePhi( const int64_t q ) const {

    # if defined (STRICT_CHECK)

        if ( q < 0 || q >= this->nq() ) {

            std::string error_message = "Index " + std::to_string(q)
                                        + " not in range of valid ordinate indices [0,"
                                        + std::to_string( this->nq() ) + ").\n";

            PRINT_ERROR( error_message.c_str() )
            throw std::out_of_range( error_message );
        }

    # endif // if defined (STRICT_CHECK)

        return this->ordinate_data.phi_ptr[q];
    }

    //--------------------------------------------------------------------------------------------------------
    //! \brief  Returns a const reference to the value of \c ordinate_data.octant[i].
    //!
    //! This yields the ordinate indices which begin each octant group.
    //--------------------------------------------------------------------------------------------------------
    inline const int64_t & Octants( const int64_t i ) const {

    # if defined (STRICT_CHECK)

        if ( i < 0 || i > 8 ) {

            std::string error_message = "Index " + std::to_string(i)
                                        + " not in valid range of [0,8] for ordinate octant boundaries.\n";

            PRINT_ERROR( error_message.c_str() )
            throw std::out_of_range( error_message );
        }

    # endif // if defined (STRICT_CHECK)

        return this->ordinate_data.octants[i];
    }

    //--------------------------------------------------------------------------------------------------------
    //! \brief  Returns a const reference to the value of \c ordinate_data.octant[2*i].
    //!
    //! This yields the ordinate indices which being each quadrant group.
    //--------------------------------------------------------------------------------------------------------
    inline const int64_t & Quadrants( const int64_t i ) const {

    # if defined (STRICT_CHECK)

        if ( i < 0 || i > 4 ) {

            std::string error_message = "Index " + std::to_string(i)
                                        + " not in valid range of [0,4] for ordinate quadrant boundaries.\n";

            PRINT_ERROR( error_message.c_str() )
            throw std::out_of_range( error_message );
        }

    # endif // if defined (STRICT_CHECK)

        return this->ordinate_data.octants[2*i];
    }

    //--------------------------------------------------------------------------------------------------------
    //! \brief  Returns a const reference to the corresponding element of the reflection map for the
    //!         quadrature rule.
    //!
    //! \param[in]  dim     Spatial dimension in which reflection occurs.
    //! \param[in]  q       Index of ordinate to reflect.
    //!
    //! \return     Returns the index of the ordinate into which the ordinate with index \pp{q} reflects
    //!             with respect to the spatial dimension \pp{dim}.
    //--------------------------------------------------------------------------------------------------------
    inline const int64_t & ReflectMap (

        const int64_t dim,
        const int64_t q

    ) const {

        return this->ordinate_data.Reflect(dim,q);
    }

# endif // if SPACE_DIMS >= 2


    //--------------------------------------------------------------------------------------------------------
    //! \brief  Returns a const reference to the quadrature weight \f$ \omega_q \f$.
    //--------------------------------------------------------------------------------------------------------
    inline const double & w( const int64_t q ) const {

    # if defined (STRICT_CHECK)

        if ( q < 0 || q >= this->nq() ) {

            std::string error_message = "Index " + std::to_string(q)
                                        + " not in range of valid ordinate indices [0,"
                                        + std::to_string( this->nq() ) + ").\n";

            PRINT_ERROR( error_message.c_str() )
            throw std::out_of_range( error_message );
        }

    # endif // if defined (STRICT_CHECK)

        return this->ordinate_data.weights[q];
    }

    //--------------------------------------------------------------------------------------------------------
    //! \brief  Returns a const reference to the number of quadrature nodes in the quadrature rule.
    //--------------------------------------------------------------------------------------------------------
    inline const int64_t & nq( void ) const {  return this->ordinate_data.num_nodes;  };

    //--------------------------------------------------------------------------------------------------------
    //! \brief  Returns a const reference to the order of the quadrature rule.
    //--------------------------------------------------------------------------------------------------------
    inline const int64_t & GetAngOrder( void ) const {  return this->ordinate_data.order;  };

    //--------------------------------------------------------------------------------------------------------
    //! \brief  Returns a const reference to the OrdinateType of the quadrature rule.
    //--------------------------------------------------------------------------------------------------------
    inline const OrdinateType & GetOrdinateType( void ) const {  return this->ordinate_data.type;  };

    //--------------------------------------------------------------------------------------------------------
    //! \brief  Returns a boolean value indicating whether or not the number of ordinates has been reduced
    //!         using symmetry. Only applicable for 2D problems.
    //--------------------------------------------------------------------------------------------------------
    inline bool GetOrdinateSymmetry( void ) const {  return this->ordinate_data.symmetric_reduce;  };

    //--------------------------------------------------------------------------------------------------------
    //! \brief  Prints the quadrature rule to the logging interface.
    //--------------------------------------------------------------------------------------------------------
    void Print( const std::string & = "  " ) const;


private:

    //========================================================================================================
    //=== HELPER ROUTINES FOR INITIALIZATION =================================================================
    //========================================================================================================

# if SPACE_DIMS == 1 || defined (DOXYCOMPILE)

    //--------------------------------------------------------------------------------------------------------
    //! \brief  Adds the Gauss-Legendre quadrature of order \pp{order} to OrdinateSet::cache.
    //--------------------------------------------------------------------------------------------------------
    static void SetupGaussLegendre( const int64_t order );

    //--------------------------------------------------------------------------------------------------------
    //! \brief  Adds the Gauss-Lobatto quadrature of order \pp{order} to OrdinateSet::cache.
    //--------------------------------------------------------------------------------------------------------
    static void SetupGaussLobatto( const int64_t order );

    //--------------------------------------------------------------------------------------------------------
    //! \brief  Adds the double-Gauss quadrature of order \pp{order} to OrdinateSet::cache.
    //--------------------------------------------------------------------------------------------------------
    static void SetupDoubleGauss( const int64_t order );

    //--------------------------------------------------------------------------------------------------------
    //! \brief  Adds the double-Radau quadrature of order \pp{order} to OrdinateSet::cache.
    //--------------------------------------------------------------------------------------------------------
    static void SetupDoubleRadau( const int64_t order );

# endif // if SPACE_DIMS == 1

# if SPACE_DIMS >= 2 || defined (DOXYCOMPILE)

    //--------------------------------------------------------------------------------------------------------
    //! \brief  Initializes \pp{this} with Chebyshev-Legendre quadrature.
    //--------------------------------------------------------------------------------------------------------
    static void SetupChebyshevLegendre( const int64_t order, const bool symmetric_reduce );

    //--------------------------------------------------------------------------------------------------------
    //! \brief  Initializes \pp{this} with spherical triangle (\f$ T_N \f$) quadrature.
    //--------------------------------------------------------------------------------------------------------
    static void SetupSphericalTriangle( const int64_t order, const bool symmetric_reduce );

# if defined (ENABLE_LEBEDEV)

    //--------------------------------------------------------------------------------------------------------
    //! \brief  Initializes \pp{this} with Lebedev quadrature.
    //--------------------------------------------------------------------------------------------------------
    static void SetupLebedev( const int64_t order, const bool symmetric_reduce );

# endif // if defined (ENABLE_LEBEDEV)

# endif // if SPACE_DIMS >= 2

    //========================================================================================================
    //=== PRIVATE DATA MEMBERS ===============================================================================
    //========================================================================================================

    struct OrdinateData {

        OrdinateType type;  //!< Specifies the type of the discrete ordinate quadrature.
        int64_t order;      //!< Specifies the order of the discrete ordinate quadrature.
        int64_t num_nodes;  //!< Specifies the number of nodes in the discrete ordinate quadrature for the given type and order.

        //!
        //! \brief  Specifies whether or not the number of ordinates is reduced using symmetry.
        //!
        //! Only used for 2D codes.
        //!
        bool symmetric_reduce;

        int64_t octants [9];    //!< Stores the ordinate indices at which each spherical octant begins.

        double * xi_ptr;    //!< Pointer to array of coordinates for each ordinate in the \f$ x_1 \f$ direction.
        double * eta_ptr;   //!< Pointer to array of coordinates for each ordinate in the \f$ x_2 \f$ direction.
        double * mu_ptr;    //!< Pointer to array of coordinates for each ordinate in the \f$ x_3 \f$ direction.
        double * phi_ptr;   //!< Pointer to array of azimuthal angles for each ordinate.
        double * weights;   //!< Pointer to array of quadrature weights for each ordinate.

        int64_t * reflect_map;  //!< Contains array for mapping between ordinates for reflecting boundary conditions.

        //
        // NOTE:    The default, implicitly-defined copy and move constructors and overloads are used; i.e.,
        //          all copy and assignment operations of these pointers are INTENTIONALLY SHALLOW.
        //

    # if SPACE_DIMS >= 2

        //----------------------------------------------------------------------------------------------------
        //! \brief  Generates nodes and weights for octants II through VIII from those of octant I by
        //!         reflections.
        //----------------------------------------------------------------------------------------------------
        void FillByOctahedralReflection( void ) const;

        //----------------------------------------------------------------------------------------------------
        //! \brief  Computes azimuthal angles for each ordinate from the Cartesian coordinates.
        //----------------------------------------------------------------------------------------------------
        void ComputePhi( void ) const;

        //----------------------------------------------------------------------------------------------------
        //! \brief  Returns an unsigned integer encoding a binary description of the octant in which the
        //!         ordinate lies.
        //----------------------------------------------------------------------------------------------------
        uint64_t DetermineOctant( const int64_t q ) const;

        //----------------------------------------------------------------------------------------------------
        //! \brief  Locates starting indices of each octant and stores the values in OrdinateData::octants.
        //----------------------------------------------------------------------------------------------------
        void LocateOctants( void );

        //----------------------------------------------------------------------------------------------------
        //! \brief  Compares the coordinates of two ordinates with given indices to define an ordering on the
        //!         ordinates suitable to sort them by octant.
        //----------------------------------------------------------------------------------------------------
        int64_t Compare( const int64_t i, const int64_t j ) const;

        //----------------------------------------------------------------------------------------------------
        //! \brief  Swaps the position of two ordinates with given indices.
        //----------------------------------------------------------------------------------------------------
        void Swap( const int64_t i, const int64_t j ) const;

        //----------------------------------------------------------------------------------------------------
        //! \brief  Enforces the min-heap property on the subtree with root index \pp{i}, considering only the
        //!         sub-array of length \pp{n}.
        //----------------------------------------------------------------------------------------------------
        void Heapify( const int64_t i, const int64_t n ) const;

        //----------------------------------------------------------------------------------------------------
        //! \brief  Sorts ordinates by octant.
        //----------------------------------------------------------------------------------------------------
        void Sort( void ) const;

        //----------------------------------------------------------------------------------------------------
        //! \brief  Returns a reference to the coordinate of the ordinate with index \pp{q} in the dimension
        //!         \pp{dim}.
        //!
        //! \param[in]  dim     Spatial dimension to get coordinate for.
        //! \param[in]  q       Index of ordinate to get coordinate for.
        //!
        //! \return     Returns a reference to the coordinate of the ordinate with index \pp{q} in the
        //!             dimension \pp{dim}.
        //----------------------------------------------------------------------------------------------------
        inline double & Coordinate (

            const int64_t dim,
            const int64_t q
        ) {

        # if defined (STRICT_CHECK)

            if ( dim < 0 || dim >= 3 ) {

                std::string error_message =   "Spatial dimension "
                                            + std::to_string(dim)
                                            + " not in valid range [0,"
                                            + std::to_string(3)
                                            + ") for "
                                            + std::string(__func__)
                                            + ".\n";

                PRINT_ERROR( error_message.c_str() )
                throw std::out_of_range( error_message );
            }

            if ( q < 0 || q >= this->num_nodes ) {

                std::string error_message =   "Index "
                                            + std::to_string(q)
                                            + " not in range of valid ordinate indices [0,"
                                            + std::to_string( this->num_nodes )
                                            + ").\n";

                PRINT_ERROR( error_message.c_str() )
                throw std::out_of_range( error_message );
            }

        # endif // if defined (STRICT_CHECK)

            switch ( dim ) {

                case 0: return this->xi_ptr [q];
                case 1: return this->eta_ptr[q];
                case 2: return this->mu_ptr [q];

                default:
                {   std::string error_message =   "Invalid spatial dimension "
                                                + std::to_string(dim)
                                                + " in "
                                                + std::string(__func__)
                                                + ".\n";

                    PRINT_ERROR( error_message.c_str() )
                    throw std::out_of_range( error_message );
                }
            }
        }

        //----------------------------------------------------------------------------------------------------
        //! \brief  Accessor for element of the reflection map stored in OrdinateData::reflect_map.
        //!
        //! \param[in]  dim     Spatial dimension in which reflection occurs.
        //! \param[in]  q       Index of ordinate to reflect.
        //!
        //! \return     Returns the index of the ordinate into which the ordinate with index \pp{q} reflects
        //!             with respect to the spatial dimension \pp{dim}.
        //----------------------------------------------------------------------------------------------------
        inline int64_t & Reflect (

            const int64_t dim,
            const int64_t q

        ) const {

        # if defined (STRICT_CHECK)

            if ( this->reflect_map == nullptr ) {

                std::string error_message =   "Pointer OrdinateData:reflect_map is NULL in '"
                                            + std::string(__func__)
                                            + "'.\n";

                PRINT_ERROR( error_message.c_str() )
                throw std::runtime_error( error_message );
            }

            CHECK_DIM_VALUE(dim)

            if ( q < 0 || q >= this->num_nodes ) {

                std::string error_message =   "Index "
                                            + std::to_string(q)
                                            + " not in range of valid ordinate indices [0,"
                                            + std::to_string( this->num_nodes )
                                            + ").\n";

                PRINT_ERROR( error_message.c_str() )
                throw std::out_of_range( error_message );
            }

        # endif // if defined (STRICT_CHECK)

            return this->reflect_map[ q + this->num_nodes * dim ];
        }

        //----------------------------------------------------------------------------------------------------
        //! \brief  Allocates memory for and fills the array specifying mappings between ordinate indices at
        //!         reflecting boundary conditions.
        //----------------------------------------------------------------------------------------------------
        void ComputeReflectMap( void );

    # if SPACE_DIMS == 2

        //----------------------------------------------------------------------------------------------------
        //! \brief  Reduces the number of ordinates for 2D problems using symmetry in the \f$ z \f$ dimension.
        //----------------------------------------------------------------------------------------------------
        void SymmetricReduce( void );

        //----------------------------------------------------------------------------------------------------
        //! \brief  Removes the ordinate with specified index.
        //----------------------------------------------------------------------------------------------------
        void Remove( const int64_t i );

    # endif // if SPACE_DIMS == 2
    # endif // if SPACE_DIMS >= 2
    };


    //--------------------------------------------------------------------------------------------------------
    //! \brief  Caches previously computed OrdinateData structs so that they may be reused with new objects.
    //--------------------------------------------------------------------------------------------------------
    static std::map< std::tuple< OrdinateType, int64_t, bool >, OrdinateData > ordinate_cache;

    //!
    //! \brief  Contains the specifications and values defining the discrete ordinate quadrature rule for the
    //!         current object.
    //!
    //! \attention  The contents of this struct are \e not exclusively owned by any given object, but are
    //!             shared between all objects that use the same OrdinateType and order pair.
    //!
    OrdinateData ordinate_data;
};


} // namespace Quadrule


# endif // ifndef __ORDINATE_SET_HPP__
