//------------------------------------------------------------------------------------------------------------
//! \file   utils/ParameterList.hpp
//! \brief  Header for ParameterList class.
//!
//! \author Michael Crockatt
//! \date   May 2017
//------------------------------------------------------------------------------------------------------------

# ifndef __INPUT_LIST_HPP__
# define __INPUT_LIST_HPP__


# include <cstdint>
# include <map>
# include <sstream>
# include <string>
# include <typeinfo>


//------------------------------------------------------------------------------------------------------------
//! \brief  Manages reading, writing, and processing of lists of parameters stored as key-value pairs.
//!
//! \attention  This class is compatible _only_ with ASCII plain text files. You are hereby warned that for
//!             any characters not in the ASCII set the behavior of this library is _undefined_.
//!
//! \note   The handling of operator overloads (operator+ and operator+=) is setup such that values in the list
//!         on the right side of the operator should _always_ override values that may already be in the list
//!         on the left side of the operator.
//!
//------------------------------------------------------------------------------------------------------------
class ParameterList {

public:

    //========================================================================================================
    //=== CONSTRUCTORS AND DESTRUCTOR ========================================================================
    //========================================================================================================


    //--------------------------------------------------------------------------------------------------------
    //! \brief  Default constructor.
    //--------------------------------------------------------------------------------------------------------
    ParameterList( void ) = default;


    //--------------------------------------------------------------------------------------------------------
    //! \brief  Default copy constructor.
    //--------------------------------------------------------------------------------------------------------
    ParameterList( const ParameterList & ) = default;


    //--------------------------------------------------------------------------------------------------------
    //! \brief  Default move constructor.
    //--------------------------------------------------------------------------------------------------------
    ParameterList( ParameterList && ) = default;


    //--------------------------------------------------------------------------------------------------------
    //! \brief  Constructor for enabling bracket initialization from an std::map object.
    //--------------------------------------------------------------------------------------------------------
    ParameterList( std::map< std::string, std::string > && );


    //--------------------------------------------------------------------------------------------------------
    //! \brief  Default destructor.
    //--------------------------------------------------------------------------------------------------------
    ~ParameterList( void ) = default;


    //========================================================================================================
    //=== OPERATOR OVERLOADS =================================================================================
    //========================================================================================================


    //--------------------------------------------------------------------------------------------------------
    //! \brief  Default copy assignment operator.
    //--------------------------------------------------------------------------------------------------------
    ParameterList & operator=( const ParameterList & ) = default;


    //--------------------------------------------------------------------------------------------------------
    //! \brief  Default move assignment operator.
    //--------------------------------------------------------------------------------------------------------
    ParameterList & operator=( ParameterList && ) = default;


    //--------------------------------------------------------------------------------------------------------
    //! \brief  Adds the items of one list to another.
    //--------------------------------------------------------------------------------------------------------
    ParameterList & operator+=( const ParameterList & );


    //========================================================================================================
    //=== PUBLIC MEMBER FUNCTIONS ============================================================================
    //========================================================================================================


    //--------------------------------------------------------------------------------------------------------
    //! \brief  Clears all items from the object.
    //--------------------------------------------------------------------------------------------------------
    ParameterList & Clear( void );


    //--------------------------------------------------------------------------------------------------------
    //! \brief  Returns the number of key-value pairs in the list.
    //--------------------------------------------------------------------------------------------------------
    std::size_t Size( void );


    //--------------------------------------------------------------------------------------------------------
    //! \brief  Returns an iterator to the beginning of the list of key-value pairs.
    //--------------------------------------------------------------------------------------------------------
    std::map< std::string, std::string >::const_iterator begin( void ) const;


    //--------------------------------------------------------------------------------------------------------
    //! \brief  Returns an iterator to the end of the list of key-value pairs.
    //--------------------------------------------------------------------------------------------------------
    std::map< std::string, std::string >::const_iterator end( void ) const;


    //--------------------------------------------------------------------------------------------------------
    //! \brief  Prints the contents of the list to the logging interface.
    //--------------------------------------------------------------------------------------------------------
    void Print( void ) const;


    //--------------------------------------------------------------------------------------------------------
    //! \brief  Adds the given key-value pair to the list. Values of existing keys are overwritten.
    //--------------------------------------------------------------------------------------------------------
    void SetValue( const std::string &, const std::string & );


    //--------------------------------------------------------------------------------------------------------
    //! \brief  Adds the given key-value pair to the list. Values of existing keys are overwritten.
    //--------------------------------------------------------------------------------------------------------
    void SetValue( const std::pair< std::string, std::string > & );


    //--------------------------------------------------------------------------------------------------------
    //! \brief  Removes a given key from the list.
    //--------------------------------------------------------------------------------------------------------
    void RemoveValue( const std::string &, const bool throw_if_not_found = true );


    //--------------------------------------------------------------------------------------------------------
    //! \brief  Read a list of key-value pairs from a file, adding those found to any already in the list.
    //--------------------------------------------------------------------------------------------------------
    ParameterList & ReadFromFile( const std::string & filename );


    //--------------------------------------------------------------------------------------------------------
    //! \brief  Writes the contents of the list to a file with the given name.
    //--------------------------------------------------------------------------------------------------------
    ParameterList & WriteToFile( const std::string & filename );


    //========================================================================================================
    //=== TEMPLATE MEMBER FUNCTIONS ==========================================================================
    //========================================================================================================


    //--------------------------------------------------------------------------------------------------------
    //! \brief  Retrieves a value from the list and converts it into the given type.
    //--------------------------------------------------------------------------------------------------------
    template< typename T >
    void GetValue( const std::string & key, T & value ) const;


    //--------------------------------------------------------------------------------------------------------
    //! \brief  Retrieves a value from the list and converts it into the given type.
    //--------------------------------------------------------------------------------------------------------
    template< typename T >
    T GetValue( const std::string & key ) const;


    //--------------------------------------------------------------------------------------------------------
    //! \brief  Retrieves a value from the list and converts it into an enum class value using the given
    //!         conversion map.
    //--------------------------------------------------------------------------------------------------------
    template< typename EnumClass >
    void GetValue( const std::string & key, EnumClass & value,
                   std::map< std::string, EnumClass > String_to_EnumClass ) const;


    //--------------------------------------------------------------------------------------------------------
    //! \brief  Retrieves a value from the list and converts it into an enum class value using the given
    //!         conversion map.
    //--------------------------------------------------------------------------------------------------------
    template< typename EnumClass >
    EnumClass GetValue( const std::string & key,
                        std::map< std::string, EnumClass > String_to_EnumClass ) const;


private:

    //========================================================================================================
    //=== PRIVATE MEMBER VARIABLES ===========================================================================
    //========================================================================================================

    //!
    //! \brief  Contains the key-value pairs of strings in the list.
    //!
    std::map< std::string, std::string > the_list;

};


//============================================================================================================
//=== NON-MEMBER OPERATOR OVERLOADS ==========================================================================
//============================================================================================================


//------------------------------------------------------------------------------------------------------------
//! \brief  Adds the contents of two ParameterList objects together.
//------------------------------------------------------------------------------------------------------------
inline ParameterList operator+ (

    ParameterList left,
    const ParameterList & right
) {

    left += right;
    return left;
}


//============================================================================================================
//=== TEMPLATE DEFINITIONS ===================================================================================
//============================================================================================================


//------------------------------------------------------------------------------------------------------------
//! \brief  Retrieves a value from the list and converts it into the given type.
//!
//! The variable \pp{value} is unmodified if this routine throws an exception.
//!
//! \param[in]      key     Key with which to search for value with.
//! \param[out]     value   Reference to object in which to store converted value with key \pp{key}.
//------------------------------------------------------------------------------------------------------------
template< typename T >
void ParameterList::GetValue (

    const std::string & key,
    T & value

) const {

    // Get value as string.
    std::string value_str;

    try {
        value_str = this->the_list.at( key );
    }
    catch ( std::out_of_range & e ) {

        std::string error_message = "Key '" + key + "' not found in input list.\n";
        throw std::out_of_range( error_message );
    }

    // Convert value string to given type.
    if ( std::stringstream( value_str ) >> value ) {/* empty */}
    else {
        std::string error_message = "Error converting value '" + value_str + "' for key '" + key +
                                    "' to type '" + std::string( typeid(T).name() ) + "'.";
        throw std::runtime_error( error_message );
    }
}

//!
//! \brief  Template specialization for type bool.
//!
template<> void ParameterList::GetValue<bool>( const std::string & key, bool & value ) const;


//------------------------------------------------------------------------------------------------------------
//! \brief  Retrieves a value from the list and converts it into the given type.
//!
//! \param[in]      key     Key with which to search for value with.
//!
//! \return     Returns the converted value with key \pp{key} from ParameterList::the_list.
//------------------------------------------------------------------------------------------------------------
template< typename T >
T ParameterList::GetValue (

    const std::string & key

) const {

    // Set default null value.
    T value = {0};

    // Get value as string.
    std::string value_str;

    try {
        value_str = this->the_list.at( key );
    }
    catch ( std::out_of_range & e ) {

        std::string error_message = "Key '" + key + "' not found in input list.\n";
        throw std::out_of_range( error_message );
    }

    // Convert value string to given type.
    if ( std::stringstream( value_str ) >> value ) {/* empty */}
    else {
        std::string error_message = "Error converting value '" + value_str + "' for key '" + key +
                                    "' to type '" + std::string( typeid(T).name() ) + "'.";
        throw std::runtime_error( error_message );
    }

    return value;
}

//!
//! \brief  Template specialization for type bool.
//!
template<> bool ParameterList::GetValue( const std::string & ) const;


//------------------------------------------------------------------------------------------------------------
//! \brief  Retrieves a value from the list and converts it into an enum class value using the given
//!         conversion map.
//!
//! The variable \pp{value} is unmodified if this routine throws an exception.
//!
//! \param[in]      key                     Key with which to search for value with.
//! \param[out]     value                   Reference to object in which to store converted value with key
//!                                         \pp{key}.
//! \param[in]      String_to_EnumClass     Map used to convert string descriptors to values of the given
//!                                         enum class.
//------------------------------------------------------------------------------------------------------------
template< typename EnumClass >
void ParameterList::GetValue (

    const std::string & key,
    EnumClass & value,
    std::map< std::string, EnumClass > String_to_EnumClass

) const {

    // Get value as string.
    std::string value_str;

    try {
        value_str = this->the_list.at( key );
    }
    catch ( std::out_of_range & e ) {

        std::string error_message = "Key '" + key + "' not found in input list.\n";
        throw std::out_of_range( error_message );
    }

    // Convert value string to given type.
    try {
        value = String_to_EnumClass.at( value_str );
    }
    catch (...) {

        std::string error_message = "Error converting value '" + value_str + "' for key '" + key +
                                    "' to type '" + std::string( typeid(EnumClass).name() ) + "'.";
        throw std::runtime_error( error_message );
    }
}


//------------------------------------------------------------------------------------------------------------
//! \brief  Retrieves a value from the list and converts it into an enum class value using the given
//!         conversion map.
//!
//! \param[in]      key                     Key with which to search for value with.
//! \param[in]      String_to_EnumClass     Map used to convert string descriptors to values of the given
//!                                         enum class.
//------------------------------------------------------------------------------------------------------------
template< typename EnumClass >
EnumClass ParameterList::GetValue (

    const std::string & key,
    std::map< std::string, EnumClass > String_to_EnumClass

) const {

    // Set default null value.
    EnumClass value = EnumClass::None;

    // Get value as string.
    std::string value_str;

    try {
        value_str = this->the_list.at( key );
    }
    catch ( std::out_of_range & e ) {

        std::string error_message = "Key '" + key + "' not found in input list.\n";
        throw std::out_of_range( error_message );
    }

    // Convert value string to given type.
    try {
        value = String_to_EnumClass.at( value_str );
    }
    catch (...) {

        std::string error_message = "Error converting value '" + value_str + "' for key '" + key +
                                    "' to type '" + std::string( typeid(EnumClass).name() ) + "'.";
        throw std::runtime_error( error_message );
    }

    return value;
}


# endif // ifndef __INPUT_LIST_HPP__
