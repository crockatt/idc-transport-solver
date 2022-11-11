//------------------------------------------------------------------------------------------------------------
//! \file   utils/Factory/Factory.hpp
//! \brief  Header file for generic factory interface.
//!
//! \author Michael M. Crockatt
//! \date   April 2018
//------------------------------------------------------------------------------------------------------------

# ifndef __FACTORY_HPP__
# define __FACTORY_HPP__


# include "utils/CLog.hpp"
# include "utils/ParameterList.hpp"


//============================================================================================================
//=== CLASS DECLARATION ======================================================================================
//============================================================================================================


//------------------------------------------------------------------------------------------------------------
//! \brief  Factory class template used to create factories for class hierarchies that are constructed from
//!         values stored in a ParameterList.
//!
//! This factory interface is intended to be used to construct objects of concrete types derived from some
//! abstract type \c BaseType. To use this template:
//! -# Define a static member function \c GetInputKey in \c BaseType that returns a value of type std::string.
//!    This value will serve as the key used to search the ParameterList provided to Factory::CreateObject for
//!    the string descriptor of the derived type that is desired.
//! -# In each concrete type derived from \c BaseType that should be constructable through this factory
//!    interface, define a static member function \c Descriptor that returns a value of type std::string.
//!    This value will serve as the string descriptor for the derived type.
//! -# Provide an explicit instantiation definition of DerivedFactory for each derived type of \c BaseType
//!    in the implementation file of each derived type that should be constructable through this factory
//!    interface.
//!
//------------------------------------------------------------------------------------------------------------
template< class BaseType >
class Factory {

    //========================================================================================================
    //=== CONSTRUCTORS AND DESTRUCTOR ========================================================================
    //========================================================================================================

protected:

    Factory( void ) = default;

    Factory( const Factory & ) = delete;
    Factory & operator=( const Factory & ) = delete;

public:

    ~Factory( void ) = default;


    //========================================================================================================
    //=== PUBLIC INTERFACE ROUTINES ==========================================================================
    //========================================================================================================

    //--------------------------------------------------------------------------------------------------------
    //! \brief  Creates and returns an object of a type derived from BaseType.
    //--------------------------------------------------------------------------------------------------------
    static BaseType * CreateObject( const ParameterList & input_list );


protected:

    //========================================================================================================
    //=== PROTECTED INTERFACE ROUTINES =======================================================================
    //========================================================================================================

    //--------------------------------------------------------------------------------------------------------
    //! \brief  Returns a reference to the map stored at Factory::factory_map_ptr.
    //--------------------------------------------------------------------------------------------------------
    static std::map< std::string, Factory<BaseType> * > & GetFactoryMap( void );

    //--------------------------------------------------------------------------------------------------------
    //! \brief  Adds a reference to a factory class to the map Factory::factory_map_ptr.
    //--------------------------------------------------------------------------------------------------------
    void AddFactory (
        const std::string & descriptor,
        Factory<BaseType> & factory
    );

    //--------------------------------------------------------------------------------------------------------
    //! \brief  Virtual interface for derived factory classes to create objects of their respective types.
    //--------------------------------------------------------------------------------------------------------
    virtual BaseType * Create( const ParameterList & input_list ) const = 0;


    //========================================================================================================
    //=== INTERNAL VARIABLES =================================================================================
    //========================================================================================================

    //!
    //! \brief  Pointer to map for mapping between string descriptors and factories for classes deriving from
    //!         BaseType.
    //!
    //! Factories are implemented as singletons and are responsible for adding themselves to this map
    //! through the Factory::AddFactory method.
    //!
    static std::map< std::string, Factory<BaseType> * > * factory_map_ptr;
};


//============================================================================================================
//=== MEMBER DEFINITIONS FOR Factory CLASS ===================================================================
//============================================================================================================


//!
//! \brief  Definition of Factory::factory_map_ptr.
//!
template< class BaseType >
std::map< std::string, Factory<BaseType> * > *
Factory<BaseType>::factory_map_ptr{ nullptr };


//------------------------------------------------------------------------------------------------------------
//! \brief  Returns a reference to the map stored at Factory::factory_map_ptr.
//!
//! Responsible for allocating the map if Factory::factory_map_ptr is null.
//------------------------------------------------------------------------------------------------------------
template< class BaseType >
std::map< std::string, Factory<BaseType> * > &
Factory<BaseType>::GetFactoryMap ( void ) {

    if ( factory_map_ptr == nullptr )
        factory_map_ptr = new std::map< std::string, Factory<BaseType> * >();

    return *factory_map_ptr;
}


//------------------------------------------------------------------------------------------------------------
//! \brief  Adds a reference to a factory class to the map Factory::factory_map.
//!
//! Factory classes deriving from Factory (i.e., using DerivedFactory) call this routine to add a mapping to
//! their singleton instance from the string descriptor of the subclass which they are intended to construct.
//!
//! \note   Because this function uses the \c insert routine of std::map each factory may only be added once.
//!
//! \param[in]  descriptor  String descriptor of the derived type associated with the factory \pp{factory}.
//! \param[in]  factory     Factory object for creating objects of types deriving from BaseType.
//------------------------------------------------------------------------------------------------------------
template< class BaseType >
void Factory<BaseType>::AddFactory (

    const std::string & descriptor,
    Factory<BaseType> & factory
) {

    PRINT_STATUS( "Executing Factory<%s>::%s.\n", typeid(BaseType).name(), __func__ )

    auto & factory_map = this->GetFactoryMap();

    if ( factory_map.count( descriptor ) ) {

        std::string error_message =   "Factory for derived type '"
                                    + descriptor
                                    + "' already present in Factory<"
                                    + typeid(BaseType).name()
                                    + ">::factory_map.\n";

        PRINT_ERROR( error_message.c_str() )
        throw std::runtime_error( error_message );
    }

    factory_map.insert( std::pair< std::string, Factory<BaseType> * >( descriptor, &factory ) );
}


//------------------------------------------------------------------------------------------------------------
//! \brief  Calls creator for desired base type implementation.
//!
//! Uses the map Factory::factory_map to determine which factory to call to create object.
//!
//! \param[in]  input_list      Contains parameters for object construction.
//!
//! \return     Returns a pointer to the created object.
//------------------------------------------------------------------------------------------------------------
template< class BaseType >
BaseType * Factory<BaseType>::CreateObject (

    const ParameterList & input_list
) {

    PRINT_STATUS( "Executing Factory<%s>::%s.\n", typeid(BaseType).name(), __func__ )

    std::string derived_type_str;

    // Get string descriptor of derived type from input list.
    try {

        derived_type_str = input_list.GetValue<std::string>( BaseType::GetInputKey() );

    } catch (...) {

        std::string error_message =   "Failed to find key '"
                                    + BaseType::GetInputKey()
                                    + "' for constructing object of type "
                                    + typeid(BaseType).name()
                                    + " in ParameterList provided to Factory<"
                                    + typeid(BaseType).name()
                                    + ">::"
                                    + std::string(__func__)
                                    + ".\n";

        throw std::out_of_range( error_message );
    }

    Factory * derived_factory = nullptr;

    // Retrieve factory pointer with given string descriptor.
    try {

        auto & factory_map = GetFactoryMap();
        derived_factory = factory_map.at( derived_type_str );

    } catch (...) {

        std::string error_message =   "Failed to access factory for derived class with descriptor '"
                                    + derived_type_str
                                    + "' in Factory<"
                                    + typeid(BaseType).name()
                                    + ">::"
                                    + std::string(__func__)
                                    + ".\n";

        throw std::out_of_range( error_message );
    }

    // Create desired object using factory.
    return derived_factory->Create( input_list );
}


# endif // ifndef __FACTORY_HPP__
