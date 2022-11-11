//------------------------------------------------------------------------------------------------------------
//! \file   ParameterList.cpp
//! \brief  Implementation of ParameterList class.
//!
//! \author Michael Crockatt
//! \date   May 2017
//------------------------------------------------------------------------------------------------------------


# include <algorithm>
# include <fstream>

# include "utils/CLog.hpp"
# include "utils/ParameterList.hpp"
# include "utils/global.hpp"


//============================================================================================================
//=== STATIC DEFINITIONS =====================================================================================
//============================================================================================================


//------------------------------------------------------------------------------------------------------------
//! \brief  Wrapper to convert string to all lower case.
//------------------------------------------------------------------------------------------------------------
static void StrLower (

    std::string & str
) {
    std::transform( str.begin(), str.end(), str.begin(), ::tolower );
}


//============================================================================================================
//=== CONSTRUCTORS AND DESTRUCTOR ============================================================================
//============================================================================================================


//--------------------------------------------------------------------------------------------------------
//! \brief  Constructor for enabling bracket initialization from an std::map object.
//--------------------------------------------------------------------------------------------------------
ParameterList::ParameterList (

    std::map< std::string, std::string > && input_map
) {
    this->the_list = input_map;
}



//============================================================================================================
//=== OPERATOR OVERLOADS =====================================================================================
//============================================================================================================


//------------------------------------------------------------------------------------------------------------
//! \brief  Adds the items of one list to another.
//!
//! \param[in]  that    ParameterList containing key-value pairs to add to this.
//------------------------------------------------------------------------------------------------------------
ParameterList & ParameterList::operator+= (

    const ParameterList & that
) {

    // We want to do insert_or_assign here, but not necessarily require C++17 support for now.
    for ( const auto & item : that.the_list )
        try         {  this->the_list.at( item.first ) = item.second;  }
        catch (...) {  this->the_list.insert( item );  }

    return *this;
}


//============================================================================================================
//=== PUBLIC MEMBER FUNCTIONS ================================================================================
//============================================================================================================


//------------------------------------------------------------------------------------------------------------
//! \brief  Clears all items from the object.
//------------------------------------------------------------------------------------------------------------
ParameterList & ParameterList::Clear ( void ) {

    the_list.clear();

    return *this;
}


//------------------------------------------------------------------------------------------------------------
//! \brief  Returns the number of key-value pairs in the list.
//------------------------------------------------------------------------------------------------------------
std::size_t ParameterList::Size ( void ) {

    return the_list.size();
}


//------------------------------------------------------------------------------------------------------------
//! \brief  Returns an iterator to the beginning of the list of key-value pairs.
//------------------------------------------------------------------------------------------------------------
std::map< std::string, std::string >::const_iterator ParameterList::begin ( void ) const {

    return this->the_list.begin();
}


//------------------------------------------------------------------------------------------------------------
//! \brief  Returns an iterator to the end of the list of key-value pairs.
//------------------------------------------------------------------------------------------------------------
std::map< std::string, std::string >::const_iterator ParameterList::end ( void ) const {

    return this->the_list.end();
}


//------------------------------------------------------------------------------------------------------------
//! \brief  Adds the given key-value pair to the list.  Values of existing keys are overwritten.
//!
//! \param[in]  new_pair    Key-value pair to add to the list.
//------------------------------------------------------------------------------------------------------------
void ParameterList::SetValue (

    const std::pair< std::string, std::string > & new_pair
) {

    this->SetValue( new_pair.first, new_pair.second );
}


//------------------------------------------------------------------------------------------------------------
//! \brief  Adds the given key-value pair to the list.  Values of existing keys are overwritten.
//!
//! \param[in]  key     Key of pair to add to list.
//! \param[in]  value   Value of pair to add to list.
//------------------------------------------------------------------------------------------------------------
void ParameterList::SetValue (

    const std::string & key,
    const std::string & value
) {

    if ( this->the_list.count( key ) > 0 ) {

        this->the_list.at( key ) = value;

    } else {

        this->the_list.insert( std::pair< std::string, std::string >( key, value ) );
    }
}


//------------------------------------------------------------------------------------------------------------
//! \brief  Removes items from the list by key.
//!
//! \param[in]  key                 Key of item to remove from the list.
//! \param[in]  throw_if_not_found  (optional) Throw an error if the key to be removed is not present in the
//!                                 list. Defaults to true.
//------------------------------------------------------------------------------------------------------------
void ParameterList::RemoveValue (

    const std::string & key,
    const bool throw_if_not_found // = true
) {

    const auto key_it = this->the_list.find( key );

    if (    (key_it == this->the_list.end())
        &&  throw_if_not_found
    ) {
        std::string error_message = "Key value '" + key + "' passed to " + std::string(__func__)
                                  + " does not exist in the list!";

        PRINT_ERROR( error_message.c_str() )
        throw std::runtime_error( error_message );
    }

    if ( key_it != this->the_list.end() )
        this->the_list.erase( key_it );
}


//------------------------------------------------------------------------------------------------------------
//! \brief  Prints the contents of the list to the logging interface.
//------------------------------------------------------------------------------------------------------------
void ParameterList::Print ( void ) const {

    for ( auto item : this->the_list ) {

        PRINT_LOG( "%s %s\n", item.first.c_str(), item.second.c_str() )
    }
}


//------------------------------------------------------------------------------------------------------------
//! \brief  Read a list of key-value pairs from a file, adding those found to any already in the list.
//!
//! \param[in]      filename    Name of file to read key-value pairs from.
//------------------------------------------------------------------------------------------------------------
ParameterList & ParameterList::ReadFromFile (

    const std::string & filename
) {

    // Open file.
    std::fstream file( filename, std::ios::in );

    if ( !file.is_open() ) {

        std::string error_message = "Failed to open file '" + filename + "' in '" + std::string(__func__)
                                    + "'.\n";

        PRINT_ERROR( error_message.c_str() )
        throw std::runtime_error( error_message );
    }

    // Read and process lines.
    std::string line, key, value;

    while ( std::getline( file, line ) ) {

        // Strip comments from line.
        line = line.substr( 0, line.find( '#' ) );

        // Convert line to lower case.
        StrLower( line );

        // Extract key and value from line.
        std::stringstream( line ) >> key >> value;

        // Insert key and value into list.
        the_list.insert( std::make_pair( key, value ) );
    }

    // Close file and return.
    file.close();

    return *this;
}


//------------------------------------------------------------------------------------------------------------
//! \brief  Writes the contents of the list to a file with the given name.
//!
//! \attention  If a file with name \pp{filename} already exists, it is overwritten without warning.
//!
//! \param[in]      filename    Name of the file to write key-value pairs to.
//------------------------------------------------------------------------------------------------------------
ParameterList & ParameterList::WriteToFile (

    const std::string & filename
) {

    // Open file.
    std::fstream file( filename, std::ios::out );

    if ( !file.is_open() ) {

        std::string error_message = "Failed to open file '" + filename + "' in '" + std::string(__func__)
                                    + "'.\n";

        PRINT_ERROR( error_message.c_str() )
        throw std::runtime_error( error_message );
    }

    // Write key-value pairs.
    for ( auto item : this->the_list ) {

        file << item.first << " " << item.second << std::endl;
    }

    // Close file and return.
    file.close();

    return *this;
}


//============================================================================================================
//=== TEMPLATE SPECIALIZATIONS ===============================================================================
//============================================================================================================


//------------------------------------------------------------------------------------------------------------
//! \brief  Retrieves a value from the list and converts it into boolean type.
//!
//! \param[in]      key     Key with which to search for value with.
//! \param[out]     value   Reference to object in which to store converted value with key \pp{key}.
//------------------------------------------------------------------------------------------------------------
template<>
void ParameterList::GetValue (

    const std::string & key,
    bool & value

) const {

    // Set default null value.
    value = false;

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
    if ( std::stringstream( value_str ) >> std::boolalpha >> value ) {/* empty */}
    else {
        std::string error_message = "Error converting value '" + value_str + "' for key '" + key +
                                    "' to type 'bool'";
        throw std::runtime_error( error_message );
    }
}


//------------------------------------------------------------------------------------------------------------
//! \brief  Retrieves a value from the list and converts it into boolean type.
//!
//! \param[in]      key     Key with which to search for value with.
//!
//! \return     Returns the converted value with key \pp{key} from ParameterList::the_list.
//------------------------------------------------------------------------------------------------------------
template<>
bool ParameterList::GetValue (

    const std::string & key

) const {

    // Set default null value.
    bool value = false;

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
    if ( std::stringstream( value_str ) >> std::boolalpha >> value ) {/* empty */}
    else {
        std::string error_message = "Error converting value '" + value_str + "' for key '" + key +
                                    "' to type 'bool'.";
        throw std::runtime_error( error_message );
    }

    return value;
}
