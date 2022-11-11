//------------------------------------------------------------------------------------------------------------
//! \file   utils/Mollifier.hpp
//! \brief  Header file for Mollifier class.
//!
//! \author Michael Crockatt
//! \date   May 2017
//------------------------------------------------------------------------------------------------------------

# ifndef __MOLLIFIER_HPP__
# define __MOLLIFIER_HPP__


# include <cstdint>

# include "utils/global.hpp"


//------------------------------------------------------------------------------------------------------------
//! \brief  Contains a discontinuous Galerkin (DG) finite element approximation of the exponential mollifier
//!
//! \todo   Description of mollifier.
//------------------------------------------------------------------------------------------------------------
class Mollifier {

public:

    //========================================================================================================
    //=== CONSTRUCTORS AND DESTRUCTOR ========================================================================
    //========================================================================================================

    //--------------------------------------------------------------------------------------------------------
    //! \brief  Delete default constructor.
    //--------------------------------------------------------------------------------------------------------
    Mollifier( void ) = delete;

    //--------------------------------------------------------------------------------------------------------
    //! \brief  Constructs a Mollifier object with the given parameters.
    //--------------------------------------------------------------------------------------------------------
    Mollifier(
        const double radius_in,
        const int64_t DG_degree_in = 2,
        const double tol_in = 1e-15
    );

    //--------------------------------------------------------------------------------------------------------
    //! \brief  Destructor for Mollifier class.
    //--------------------------------------------------------------------------------------------------------
    ~Mollifier( void );


    //========================================================================================================
    //=== MOLLIFIER EVALUATION ROUTINES ======================================================================
    //========================================================================================================

# if SPACE_DIMS == 1

    //--------------------------------------------------------------------------------------------------------
    //! \brief  Evaluates the mollifier at the given point.
    //--------------------------------------------------------------------------------------------------------
    double operator()( const double x ) const;

# elif SPACE_DIMS == 2

    //--------------------------------------------------------------------------------------------------------
    //! \brief  Evaluates the mollifier at the given point.
    //--------------------------------------------------------------------------------------------------------
    double operator()( const double x, const double y ) const;

# elif SPACE_DIMS == 3

    //--------------------------------------------------------------------------------------------------------
    //! \brief  Evaluates the mollifier at the given point.
    //--------------------------------------------------------------------------------------------------------
    double operator()( const double x, const double y, const double z ) const;

# endif // if SPACE_DIMS == ?

    //--------------------------------------------------------------------------------------------------------
    //! \brief  Evaluates the exact value of the one-dimensional mollifier at the given point.
    //--------------------------------------------------------------------------------------------------------
    double Exact( const double x ) const;

    //--------------------------------------------------------------------------------------------------------
    //! \brief  Evaluates the exact value of the two-dimensional mollifier at the given point.
    //--------------------------------------------------------------------------------------------------------
    double Exact( const double x, const double y ) const;


    //========================================================================================================
    //=== ADDITIONAL PUBLIC MEMBER FUNCTIONS =================================================================
    //========================================================================================================

    //--------------------------------------------------------------------------------------------------------
    //! \brief  Prints the parameters of the data stored in a Mollifier object to the global logging
    //!         interface.
    //--------------------------------------------------------------------------------------------------------
    void Print( void ) const;

    //--------------------------------------------------------------------------------------------------------
    //! \brief  Returns the radius of the mollifier.
    //--------------------------------------------------------------------------------------------------------
    double Radius( void ) const {  return this->radius;  }


private:

    //========================================================================================================
    //=== PRIVATE MEMBER VARIABLES ===========================================================================
    //========================================================================================================

    int64_t nx;                 //!< Number of cells to use in DG approximation of mollifier.

    double radius;              //!< Radius of the support of the mollifier.
    double dx;                  //!< Cell width in DG approximation of mollifier.
    double cell_factor;         //!< Scaling factor used to compute the cell in which a given value lies.

    int64_t DG_degree;          //!< Maximum degree of DG basis polynomials.

    double * density;           //!< Pointer at which DG coefficients of mollifier are stored.

    size_t sizeof_density;      //!< Number of bytes allocated at Mollifier::density.
    int64_t dimof_density;      //!< Number of doubles allocated at Mollifier::density.

    double tol;                 //!< Tolerance used when constructing mollifier object.
    double normalization_const; //!< Normalization constant of exact mollifier function.


    //========================================================================================================
    //=== PRIVATE MEMBER FUNCTIONS ===========================================================================
    //========================================================================================================

    //--------------------------------------------------------------------------------------------------------
    //! \brief  Indexing function for coefficient array stored at Mollifier::density.
    //!
    //! This function converts multidimensional indices which describe a given density coefficient into a
    //! one-dimensional index for Mollifier::density. The coefficients are indexed by spatial cell and degree
    //! of basis polynomial. The spatial cells are indexed radially from the center of the mollifier (without
    //! ghost cells).
    //!
    //! \param[in]      i           Spatial cell index.
    //! \param[in]      d           Degree of spatial basis polynomial.
    //!
    //! \return     Appropriate index within Mollifier::density array for specified coefficient.
    //!
    //! \see    Mollifier::operator()()
    //--------------------------------------------------------------------------------------------------------
    inline size_t Index (

        const int64_t i,
        const int64_t d

    ) const {

    # if defined (STRICT_CHECK)

        if ( i < 0 || i >= this->nx ) {

            std::string error_message
                = "Value " + std::to_string(i) + " for index i in '" + std::string(__func__)
                  + "' outsize permissible range [ 0, " + std::to_string( this->nx ) + " ).\n";

            PRINT_ERROR( error_message.c_str() )
            throw std::out_of_range( error_message );
        }

        if ( d < 0 || d > this->DG_degree ) {

            std::string error_message
                = "Value " + std::to_string(d) + " for degree d in '" + std::string(__func__)
                  + "' outsize permissible range [ 0, " + std::to_string( this->DG_degree ) + " ].\n";

            PRINT_ERROR( error_message.c_str() )
            throw std::out_of_range( error_message );
        }

    # endif // if defined (STRICT_CHECK)

        return ( (d) + (this->DG_degree + 1)*(i) );
    }

    //--------------------------------------------------------------------------------------------------------
    //! \brief  Returns a reference to the specified DG coefficient inside the array Mollifier::density.
    //!
    //! \param[in]      i           Spatial cell index.
    //! \param[in]      d           Degree of spatial basis polynomial.
    //!
    //! \return     Reference to the specified coefficient.
    //!
    //! \see    Mollifier::Index()
    //--------------------------------------------------------------------------------------------------------
    inline double & operator() (

        const int64_t i,
        const int64_t d

    ) const {

    # if defined (STRICT_CHECK)

        if ( this->density == nullptr ) {

            std::string error_message = "Pointer this->density is NULL in '" + std::string(__func__) + "'.\n";

            PRINT_ERROR( error_message.c_str() )
            throw std::runtime_error( error_message );
        }

    # endif // if defined (STRICT_CHECK)

        return this->density[ Index(i,d) ];
    }

    //--------------------------------------------------------------------------------------------------------
    //! \brief  Evaluates the radially-defined mollifier at the specified radius from the origin.
    //--------------------------------------------------------------------------------------------------------
    double EvalAtRadius( const double r ) const;

    //--------------------------------------------------------------------------------------------------------
    //! \brief  One-dimensional \f$ C^{\infty} \f$ exponential function with compact support in
    //!         \f$ |x| \leq r \f$.
    //--------------------------------------------------------------------------------------------------------
    double ExponentialBump( const double x ) const;

    //--------------------------------------------------------------------------------------------------------
    //! \brief  Computes a numerical approximation to the \f$ L^2 \f$ norm of the exact mollifier.
    //--------------------------------------------------------------------------------------------------------
    double ComputeNorm( const double ax, const double bx, const int64_t GK_order,
                        const double * const GK_nodes, const double * const G_weights,
                        const double * const K_weights ) const;

};


# endif // ifndef __MOLLIFIER_HPP__
