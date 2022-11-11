//------------------------------------------------------------------------------------------------------------
//! \file   utils/Factory/DerivedFactory.hpp
//! \brief  Header file for derived component of factory interface.
//!
//! \author Michael M. Crockatt
//! \date   April 2018
//------------------------------------------------------------------------------------------------------------

# ifndef __DERIVED_FACTORY_HPP__
# define __DERIVED_FACTORY_HPP__


# include "utils/Factory/Factory.hpp"


//============================================================================================================
//=== CLASS DECLARATION ======================================================================================
//============================================================================================================


//------------------------------------------------------------------------------------------------------------
//! \brief  Class template for factories of classes derived from \c BaseType.
//!
//! Used to enable each derived type of \c BaseType to be added to the factory of \c BaseType.
//------------------------------------------------------------------------------------------------------------
template< class BaseType, class DerivedType >
class DerivedFactory :
    public Factory<BaseType>
{

    //========================================================================================================
    //=== CONSTRUCTORS AND DESTRUCTOR ========================================================================
    //========================================================================================================

protected:

    DerivedFactory( void );

    DerivedFactory( const DerivedFactory< BaseType, DerivedType > & ) = delete;

    DerivedFactory< BaseType, DerivedType > &
        operator=( const DerivedFactory< BaseType, DerivedType > & ) = delete;

public:

    ~DerivedFactory( void ) = default;


protected:

    //========================================================================================================
    //=== PROTECTED INTERFACE ROUTINES =======================================================================
    //========================================================================================================

    //!
    //! \brief  Inherit Factory::AddFactory.
    //!
    using Factory<BaseType>::AddFactory;

    //--------------------------------------------------------------------------------------------------------
    //! \brief  Creates and returns an object of a type derived from BaseType.
    //--------------------------------------------------------------------------------------------------------
    BaseType * Create( const ParameterList & input_list ) const override;


    //========================================================================================================
    //=== INTERNAL VARIABLES =================================================================================
    //========================================================================================================

    //!
    //! \brief  Singleton instance.
    //!
    static DerivedFactory< BaseType, DerivedType > instance;
};


//============================================================================================================
//=== DEFINITION OF MEMBER FUNCTIONS =========================================================================
//============================================================================================================


//!
//! \brief  Define static singleton instance DerivedFactory::instance.
//!
template< class BaseType, class DerivedType >
DerivedFactory< BaseType, DerivedType >
DerivedFactory< BaseType, DerivedType >::instance;


//------------------------------------------------------------------------------------------------------------
//! \brief  Default constructor.
//!
//! Adds singleton object to Factory::factory_map using DerivedType::Descriptor().
//!
//------------------------------------------------------------------------------------------------------------
template< class BaseType, class DerivedType >
DerivedFactory< BaseType, DerivedType >::DerivedFactory ( void ) {

    PRINT_STATUS( "Executing DerivedFactory<%s,%s>::%s.\n",
                  typeid(BaseType).name(), typeid(DerivedType).name(), __func__ )

    AddFactory( DerivedType::Descriptor(), this->instance );
}


//------------------------------------------------------------------------------------------------------------
//! \brief  Creates and returns an object of type derived from ImplicitSolver.
//!
//! \param[in]  domain_decomposition    Contains parameters of spatial discretization.
//! \param[in]  ordinate_set            Contains parameters of angular discretization.
//! \param[in]  input_list              Contains additional parameters.
//! \param[in]  dt                      Initial timestep size.
//!
//! \return     Returns a pointer to the created object.
//------------------------------------------------------------------------------------------------------------
template< class BaseType, class DerivedType >
BaseType *
DerivedFactory< BaseType, DerivedType >::Create (

    const ParameterList & input_list

) const {

    PRINT_STATUS( "Executing DerivedFactory<%s,%s>::%s.\n",
                  typeid(BaseType).name(), typeid(DerivedType).name(), __func__ )

    return new DerivedType( input_list );
}


# endif // ifndef __DERIVED_FACTORY_HPP__
