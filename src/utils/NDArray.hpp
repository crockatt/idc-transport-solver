//------------------------------------------------------------------------------------------------------------
//! \file   utils/NDArray.hpp
//! \brief  Header for \f$ N \f$-dimensional array class template.
//!
//! \author Michael Crockatt
//! \date   February 2018
//------------------------------------------------------------------------------------------------------------

# ifndef __ND_ARRAY_HPP__
# define __ND_ARRAY_HPP__

# include <initializer_list>
# include <string>


//------------------------------------------------------------------------------------------------------------
//! \brief  Declaration of simple \f$ N \f$-dimensional array class template.
//!
//! Elements are stored in row-major order.
//------------------------------------------------------------------------------------------------------------
template< typename ScalarType >
class NDArray {

public:

    //========================================================================================================
    //=== CONSTRUCTORS AND DESTRUCTOR ========================================================================
    //========================================================================================================

    //--------------------------------------------------------------------------------------------------------
    //! \brief  Delete default constructor.
    //--------------------------------------------------------------------------------------------------------
    NDArray( void ) = delete;

    //--------------------------------------------------------------------------------------------------------
    //! \brief  Constructs an NDArray object using the specified parameters.
    //--------------------------------------------------------------------------------------------------------
    NDArray( std::initializer_list<size_t> );

    //--------------------------------------------------------------------------------------------------------
    //! \brief  Delete copy constructor.
    //--------------------------------------------------------------------------------------------------------
    NDArray( const NDArray<ScalarType> & ) = delete;

    //--------------------------------------------------------------------------------------------------------
    //! \brief  Destructor for NDArray class.
    //--------------------------------------------------------------------------------------------------------
    virtual ~NDArray( void );


    //========================================================================================================
    //=== OPERATOR OVERLOADS =================================================================================
    //========================================================================================================

    //--------------------------------------------------------------------------------------------------------
    //! \brief  Delete copy assignment operation.
    //--------------------------------------------------------------------------------------------------------
    NDArray<ScalarType> & operator=( const NDArray<ScalarType> & );

    //--------------------------------------------------------------------------------------------------------
    //! \brief  Accessor overload without bounds checking.
    //--------------------------------------------------------------------------------------------------------
    ScalarType & operator[]( std::initializer_list<size_t> ) const;


    //========================================================================================================
    //=== ADDITIONAL MEMBER FUNCTIONS ========================================================================
    //========================================================================================================

    //--------------------------------------------------------------------------------------------------------
    //! \brief  Accessor routine with bounds checking.
    //--------------------------------------------------------------------------------------------------------
    ScalarType & at( std::initializer_list<size_t> ) const;


private:

    //========================================================================================================
    //=== PRIVATE HELPER FUNCTIONS ===========================================================================
    //========================================================================================================

    //--------------------------------------------------------------------------------------------------------
    //! \brief  Indexing function for linearized array.
    //--------------------------------------------------------------------------------------------------------
    size_t Index( std::initializer_list<size_t> ) const;


    //========================================================================================================
    //=== INTERNAL VARIABLES =================================================================================
    //========================================================================================================

    size_t N;                       //!< Number of dimensions in array.
    std::vector<size_t> dims;       //!< Holds the number of elements in each dimension.
    ScalarType * array;             //!< Pointer to the actual array.

};


//============================================================================================================
//=== CONSTRUCTORS AND DESTRUCTOR ============================================================================
//============================================================================================================


//------------------------------------------------------------------------------------------------------------
//! \brief  Constructs an NDArray object using the specified parameters.
//!
//! Elements of the array are value-initialized using the default constructor for \c ScalarType.
//------------------------------------------------------------------------------------------------------------
template< typename ScalarType >
NDArray<ScalarType>::NDArray (

    std::initializer_list<size_t> dim_sizes
) :
    N{ dim_sizes.size() },
    dims{ dim_sizes }
{

    // Check that number of dimensions is positive.
    if ( N <= 0 )
        throw std::invalid_argument( "Number of dimensions for NDArray object must be positive.\n" );

    // Check that the size of each dimension is positive.
    for ( size_t dim = 0; dim < this->N; ++dim )
        if ( this->dims.at(dim) <= 0 )
            throw std::invalid_argument(   "Size of dimensions for NDArray object must be positive: Size "
                                         + std::to_string( this->dims.at(dim) )
                                         + " given for dimension "
                                         + std::to_string(dim)
                                         + ".\n"
                                       );
    // Allocate memory.
    size_t num_elements = 1;
    for ( size_t size_of_dim : dims ) {  num_elements *= size_of_dim;  }

    // Elements of array are value-initialized using the default constructor for ScalarType.
    this->array = new ScalarType[ num_elements ]();
}


//------------------------------------------------------------------------------------------------------------
//! \brief  Destructor for NDArray class.
//------------------------------------------------------------------------------------------------------------
template< typename ScalarType >
NDArray<ScalarType>::~NDArray ( void ) {

    delete [] this->array;
}


//============================================================================================================
//=== ACCESSOR ROUTINES ======================================================================================
//============================================================================================================


//------------------------------------------------------------------------------------------------------------
//! \brief  Accessor overload without bounds checking.
//------------------------------------------------------------------------------------------------------------
template< typename ScalarType >
ScalarType & NDArray<ScalarType>::operator[] (

    std::initializer_list<size_t> idx

) const {

    return this->array[ this->Index(idx) ];
}


//------------------------------------------------------------------------------------------------------------
//! \brief  Accessor routine with bounds checking.
//------------------------------------------------------------------------------------------------------------
template< typename ScalarType >
ScalarType & NDArray<ScalarType>::at (

    std::initializer_list<size_t> idx

) const {

    // Check number of indices.
    if ( idx.size() != this->N )
        throw std::invalid_argument(   "Number of indices provided ("
                                     + std::to_string( idx.size() )
                                     + ") not equal to number of dimensions ("
                                     + std::to_string( this->N )
                                     + ") in NDArray object.\n"
                                   );

    // Check bounds on each index.
    size_t dim = 0;

    for ( size_t i : idx ) {

        if ( i >= this->dims.at(dim) )
            throw std::invalid_argument(   "Value ("
                                         + std::to_string(i)
                                         + ") for index in dimension "
                                         + std::to_string(dim)
                                         + " outside valid range [0, "
                                         + std::to_string( this->dims.at(dim) )
                                         + ").\n"
                                       );

        ++dim;
    }

    return (*this)[ idx ];
}


//============================================================================================================
//=== PRIVATE HELPER FUNCTIONS ===============================================================================
//============================================================================================================


//------------------------------------------------------------------------------------------------------------
//! \brief  Indexing function for linearized array.
//------------------------------------------------------------------------------------------------------------
template< typename ScalarType >
size_t NDArray<ScalarType>::Index (

    std::initializer_list<size_t> idx

) const {

    size_t linear_index = 0;

    size_t dim = 0;
    for ( size_t i : idx ) {

        linear_index *= this->dims.at(dim);
        linear_index += i;

        ++dim;
    }

    return linear_index;
}


# endif // ifndef __ND_ARRAY_HPP__
