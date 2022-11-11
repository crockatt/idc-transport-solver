//------------------------------------------------------------------------------------------------------------
//! \file   linear_solvers/Abstract/SweepPattern/SweepPatternDerivedFactory.hpp
//! \brief  Header file containing declaration of SweepPatternDerivedFactory class template.
//!
//! \author Michael M. Crockatt
//! \date   February 2018
//------------------------------------------------------------------------------------------------------------

# ifndef __ABSTRACT__SWEEP_PATTERN_DERIVED_FACTORY_HPP__
# define __ABSTRACT__SWEEP_PATTERN_DERIVED_FACTORY_HPP__


//============================================================================================================
//=== CLASS DECLARATION ======================================================================================
//============================================================================================================


//------------------------------------------------------------------------------------------------------------
//! \brief  Class template for factories of classes derived from SweepPattern.
//!
//! Implements the interface defined by SweepPatternFactory.
//!
//! Let \c DerivedPattern denote a class that has been derived from SweepPattern. Use of this template
//! to register \c DerivedPattern in SweepPatternFactory::factory_map so that objects of type
//! \c DerivedPattern may be created requires two steps:
//!     1. Declare a specialization of this template for \c DerivedPattern.
//!     2. Perform an explicit instantiation of this specialization in the implementation file for
//!        \c DerivedPattern.
//!
//------------------------------------------------------------------------------------------------------------
template< class OrdinateFlux >
template< class DerivedPattern >
class SweepOperator<OrdinateFlux>::SweepPattern::SweepPatternDerivedFactory :
    public SweepOperator<OrdinateFlux>::SweepPattern::SweepPatternFactory
{

    //========================================================================================================
    //=== CONSTRUCTORS AND DESTRUCTOR ========================================================================
    //========================================================================================================

protected:

    SweepPatternDerivedFactory( void );

    SweepPatternDerivedFactory( const SweepPatternDerivedFactory<DerivedPattern> & ) = delete;

    SweepPatternDerivedFactory<DerivedPattern> &
        operator=( const SweepPatternDerivedFactory<DerivedPattern> & ) = delete;

public:

    ~SweepPatternDerivedFactory( void ) = default;


protected:

    //========================================================================================================
    //=== PROTECTED INTERFACE ROUTINES =======================================================================
    //========================================================================================================

    //!
    //! \brief  Inherit SweepPatternFactory::AddFactory.
    //!
    using SweepPattern::SweepPatternFactory::AddFactory;

    //--------------------------------------------------------------------------------------------------------
    //! \brief  Creates and returns an object of a type derived from SweepPattern.
    //--------------------------------------------------------------------------------------------------------
    SweepPattern * Create (
        const SweepOperator<OrdinateFlux> & enclosing,
        const ParameterList & input_list
    ) const override;


    //========================================================================================================
    //=== INTERNAL VARIABLES =================================================================================
    //========================================================================================================

    //!
    //! \brief  Singleton instance.
    //!
    static SweepPatternDerivedFactory<DerivedPattern> instance;
};


//============================================================================================================
//=== DEFINITION OF MEMBER FUNCTIONS =========================================================================
//============================================================================================================


//------------------------------------------------------------------------------------------------------------
//! \brief  Default constructor.
//!
//! Adds singleton object to SweepPatternFactory::factory_map using
//! SweepPatternDerivedFactory::Descriptor().
//!
//------------------------------------------------------------------------------------------------------------
template< class OrdinateFlux >
template< class DerivedPattern >
SweepOperator<OrdinateFlux>::SweepPattern::SweepPatternDerivedFactory<DerivedPattern>::
SweepPatternDerivedFactory ( void ) {

    PRINT_STATUS( "Executing SweepOperator<%s>::SweepPattern::SweepPatternDerivedFactory<%s>::%s.\n",
                  typeid(OrdinateFlux).name(), typeid(DerivedPattern).name(), __func__ )

    AddFactory( DerivedPattern::Descriptor(), this->instance );
}


//------------------------------------------------------------------------------------------------------------
//! \brief  Creates and returns an object of type derived from SweepPattern.
//!
//! \param[in]  enclosing   Reference to instance of enclosing SweepOperator class.
//! \param[in]  input_list  Contains parameters for initialization of object.
//!
//! \return     Returns a pointer to the created object.
//------------------------------------------------------------------------------------------------------------
template< class OrdinateFlux >
template< class DerivedPattern >
typename SweepOperator<OrdinateFlux>::SweepPattern *
SweepOperator<OrdinateFlux>::SweepPattern::SweepPatternDerivedFactory<DerivedPattern>::Create (

    const SweepOperator<OrdinateFlux> & enclosing,
    const ParameterList & input_list

) const {

    PRINT_STATUS( "Executing SweepOperator<%s>::SweepPattern::SweepPatternDerivedFactory<%s>::%s.\n",
                  typeid(OrdinateFlux).name(), typeid(DerivedPattern).name(), __func__ )

    return new DerivedPattern( enclosing, input_list );
}


//============================================================================================================
//=== INSTANTIATION OF STATIC MEMBERS ========================================================================
//============================================================================================================


//!
//! \brief  Define static singleton instance SweepPatternDerivedFactory::instance.
//!
template< class OrdinateFlux >
template< class DerivedPattern >
typename SweepOperator<OrdinateFlux>::SweepPattern::
    template SweepPatternDerivedFactory<DerivedPattern>
SweepOperator<OrdinateFlux>::SweepPattern::SweepPatternDerivedFactory<DerivedPattern>::instance;


# endif // ifndef __ABSTRACT__SWEEP_PATTERN_DERIVED_FACTORY_HPP__
