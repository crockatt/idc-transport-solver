//------------------------------------------------------------------------------------------------------------
//! \file   utils/MVWrapper.hpp
//! \brief  Header for MVWrapper class template.
//!
//! \author Michael M. Crockatt
//! \date   February 2018
//------------------------------------------------------------------------------------------------------------

# ifndef __MVWRAPPER_HPP__
# define __MVWRAPPER_HPP__


//------------------------------------------------------------------------------------------------------------
//! \brief  Class template used to wrap different vector and multivector objects. This class defines a
//!         uniform interface into the underlying data structure, regardless of the lower-level library
//!         implementation.
//!
//! This class template is configured to provide access to only one vector (column) of the multivector. The
//! vector to which access is provided is set during object construction and cannot be changed.
//!
//------------------------------------------------------------------------------------------------------------
template< typename ViewType >
class MVWrapper {

public:

    //========================================================================================================
    //=== CONSTRUCTORS AND DESTRUCTOR ========================================================================
    //========================================================================================================

    //--------------------------------------------------------------------------------------------------------
    //! \brief  Delete default constructor.
    //--------------------------------------------------------------------------------------------------------
    MVWrapper( void ) = delete;

    //--------------------------------------------------------------------------------------------------------
    //! \brief  Delete copy constructor.
    //--------------------------------------------------------------------------------------------------------
    MVWrapper( const MVWrapper & ) = delete;

    //--------------------------------------------------------------------------------------------------------
    //! \brief  Construct a wrapper using the provided values.
    //--------------------------------------------------------------------------------------------------------
    MVWrapper(
        ViewType & data_in,
        const int64_t LD_data_in = 0,
        const int64_t column_in = 0
    );

    //--------------------------------------------------------------------------------------------------------
    //! \brief  Default destructor.
    //--------------------------------------------------------------------------------------------------------
    ~MVWrapper( void ) = default;


    //========================================================================================================
    //=== OPERATOR OVERLOADS =================================================================================
    //========================================================================================================

    //--------------------------------------------------------------------------------------------------------
    //! \brief  Delete copy assignment operator.
    //--------------------------------------------------------------------------------------------------------
    MVWrapper & operator=( const MVWrapper & ) = delete;

    //--------------------------------------------------------------------------------------------------------
    //! \brief  Provides access to the specified element of the multivector.
    //--------------------------------------------------------------------------------------------------------
    inline double & operator()( const int64_t ) const;


private:

    //========================================================================================================
    //=== INTERNAL VARIABLES =================================================================================
    //========================================================================================================

    ViewType & data;                //!< Reference to underlying data representation.
    const int64_t LD_data;          //!< Leading dimension for the array data.
    const int64_t column;           //!< Column of multivector to provide access to.

};


//------------------------------------------------------------------------------------------------------------
//! \brief  Construct a wrapper using the provided values.
//!
//! \param[in]  data_in     Reference to underlying data representation.
//! \param[in]  LD_data_in  Dimension of each vector.
//! \param[in]  column_in   Column of multivector to provide access to.
//------------------------------------------------------------------------------------------------------------
template< typename ViewType >
MVWrapper<ViewType>::MVWrapper (

    ViewType & data_in,
    const int64_t LD_data_in,   // = 0
    const int64_t column_in     // = 0
) :
    data{ data_in },
    LD_data{ LD_data_in },
    column{ column_in }
{
    PRINT_STATUS( "Executing MVWrapper<%s>::s.\n", typeid(ViewType).name(), __func__ )
}


//------------------------------------------------------------------------------------------------------------
//! \brief  Provides access to the specified element of the multivector.
//!
//! Specialization for the case where MVWrapper::data is a pointer to type \c double. Assumes a column-major
//! ordering of the multivector stored at MVWrapper::data.
//!
//! \param[in]  i   Index of element of vector to access.
//------------------------------------------------------------------------------------------------------------
template<>
inline double & MVWrapper<double *>::operator() (

    const int64_t i

) const {

    PRINT_STATUS( "Executing MVWrapper<double *>::s.\n", __func__ )

    return this->data[ i + this->LD_data * this->column ];
}


//------------------------------------------------------------------------------------------------------------
//! \brief  Provides access to the specified element of the multivector.
//!
//! This generic interface assumes that the object containing the underlying data (\pp{data}) has overloaded
//! the function call operator to accept two indices: the second index specifies the column (vector) to access
//! and the first index specifies the element of that column to access.
//!
//! \param[in]  i   Index of element of vector to access.
//------------------------------------------------------------------------------------------------------------
template< typename ViewType >
inline double & MVWrapper<ViewType>::operator() (

    const int64_t i

) const {

    PRINT_STATUS( "Executing MVWrapper<%s>::s.\n", typeid(ViewType).name(), __func__ )

    return this->data( i, column );
}


# endif // ifndef __MVWRAPPER_HPP__
