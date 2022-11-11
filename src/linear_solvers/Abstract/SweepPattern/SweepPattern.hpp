//------------------------------------------------------------------------------------------------------------
//! \file   linear_solvers/Abstract/SweepPattern/SweepPattern.hpp
//! \brief  Header file containing declaration of SweepPattern abstract class and associated factory class.
//!
//! \note   The declarations in this header are used for the internal mechanics of the SweepPattern classes.
//!         This header is therefore not intended to be included into user implementations: Instead, include
//!         SweepOperator.hpp and interface with the SweepOperator class defined there.
//!
//! \author Michael M. Crockatt
//! \date   March 2018
//------------------------------------------------------------------------------------------------------------

# ifndef __ABSTRACT__SWEEP_PATTERN_HPP__
# define __ABSTRACT__SWEEP_PATTERN_HPP__


# include <map>
# include <memory>
# include <string>
# include <vector>

# include "linear_solvers/Abstract/SweepOperator.hpp"
# include "utils/ParameterList.hpp"
# include "utils/SIMD.hpp"


//------------------------------------------------------------------------------------------------------------
//! \brief  Abstract class outlining objects implementing mesh-traversal algorithms used for discrete
//!         ordinates transport sweeps.
//!
//! Constructors are declared protected and nested class SweepPatternFactory is declared public. This is done
//! to force construction of SweepPattern objects through the factory interface.
//!
//------------------------------------------------------------------------------------------------------------
template< class OrdinateFlux >
class SweepOperator<OrdinateFlux>::SweepPattern {

public:

    //!
    //! \brief  Declare SweepPatternFactory as public to expose factory interface.
    //!
    class SweepPatternFactory;


protected:

    //!
    //! \brief  Declare protected template factory class for derived pattern classes.
    //!
    template< class DerivedPattern >
    class SweepPatternDerivedFactory;


    //========================================================================================================
    //=== CONSTRUCTORS AND DESTRUCTOR ========================================================================
    //========================================================================================================

    //--------------------------------------------------------------------------------------------------------
    //! \brief  Delete default constructor.
    //--------------------------------------------------------------------------------------------------------
    SweepPattern( void ) = delete;

    //--------------------------------------------------------------------------------------------------------
    //! \brief  Delete copy constructor.
    //--------------------------------------------------------------------------------------------------------
    SweepPattern( const SweepPattern & ) = delete;

    //--------------------------------------------------------------------------------------------------------
    //! \brief  Initializes a SweepPattern object.
    //--------------------------------------------------------------------------------------------------------
    SweepPattern (
        const SweepOperator<OrdinateFlux> & enclosing,
        const ParameterList & // input_list
    );

public:

    //--------------------------------------------------------------------------------------------------------
    //! \brief  Virtual destructor.
    //--------------------------------------------------------------------------------------------------------
    virtual ~SweepPattern( void );


    //========================================================================================================
    //=== OPERATOR OVERLOADS =================================================================================
    //========================================================================================================

    //--------------------------------------------------------------------------------------------------------
    //! \brief  Delete copy assignment operator.
    //--------------------------------------------------------------------------------------------------------
    SweepPattern & operator=( const SweepPattern & ) = delete;


    //========================================================================================================
    //=== PUBLIC INTERFACE ROUTINES ==========================================================================
    //========================================================================================================

    //--------------------------------------------------------------------------------------------------------
    //! \brief  Returns the string descriptor for an objects type.
    //--------------------------------------------------------------------------------------------------------
    virtual const std::string GetDescriptor( void ) const = 0;

    //--------------------------------------------------------------------------------------------------------
    //! \brief  Prints a summary of the object's configuration to the logging interface.
    //--------------------------------------------------------------------------------------------------------
    virtual void Print( const std::string & = "  " ) const;

    //--------------------------------------------------------------------------------------------------------
    //! \brief  Sets parameters given in a ParameterList.
    //--------------------------------------------------------------------------------------------------------
    virtual void SetParameters( const ParameterList & );

    //--------------------------------------------------------------------------------------------------------
    //! \brief  Performs a transport sweep.
    //--------------------------------------------------------------------------------------------------------
    virtual void Sweep (
        const OrdinateFlux & source,
        OrdinateFlux & result,
        const RKDG::CrossSection & sigma,
        const RKDG::OrdinateFlux * const initial_condition = nullptr
    ) const = 0;


protected:

    //========================================================================================================
    //=== INTERNAL VARIABLES =================================================================================
    //========================================================================================================

    //!
    //! \brief  Reference to instance of enclosing SweepOperator class.
    //!
    const SweepOperator<OrdinateFlux> & sw_op;


    //========================================================================================================
    //=== PROTECTED MEMBER FUNCTIONS =========================================================================
    //========================================================================================================

    //--------------------------------------------------------------------------------------------------------
    //! \brief  Construct lists for SIMD blocking across angular ordinates.
    //--------------------------------------------------------------------------------------------------------
    # if SPACE_DIMS == 1
        std::vector< SIMD_BlkIdx<0> >
    # else
        std::vector< std::vector< SIMD_BlkIdx<0> > >
    # endif
    ConstructSIMDAngleBlocks( void ) const;

};


//------------------------------------------------------------------------------------------------------------
//! \brief  Abstract class defining interface for factories used to construct objects derived from
//!         SweepPattern class.
//!
//! Concrete classes derived from SweepPattern should use the class template SweepPatternDerivedFactory (which
//! derives from this class) to enable their construction.
//!
//------------------------------------------------------------------------------------------------------------
template< class OrdinateFlux >
class SweepOperator<OrdinateFlux>::SweepPattern::SweepPatternFactory {

    //========================================================================================================
    //=== CONSTRUCTORS AND DESTRUCTOR ========================================================================
    //========================================================================================================

protected:

    SweepPatternFactory( void ) = default;

    SweepPatternFactory( const SweepPatternFactory & ) = delete;
    SweepPatternFactory & operator=( const SweepPatternFactory & ) = delete;

public:

    ~SweepPatternFactory( void ) = default;


    //========================================================================================================
    //=== PUBLIC INTERFACE ROUTINES ==========================================================================
    //========================================================================================================

    //--------------------------------------------------------------------------------------------------------
    //! \brief  Creates and returns an object of a type derived from SweepPattern.
    //--------------------------------------------------------------------------------------------------------
    static SweepPattern * CreatePattern (
        const SweepOperator<OrdinateFlux> & enclosing,
        const ParameterList & input_list
    );


protected:

    //========================================================================================================
    //=== PROTECTED INTERFACE ROUTINES =======================================================================
    //========================================================================================================

    //--------------------------------------------------------------------------------------------------------
    //! \brief  Returns a reference to the map stored at SweepPatternFactory::factory_map_ptr.
    //--------------------------------------------------------------------------------------------------------
    static std::map< std::string, SweepPatternFactory * > & GetFactoryMap( void );

    //--------------------------------------------------------------------------------------------------------
    //! \brief  Adds a reference to a factory class to the map SweepPatternFactory::factory_map_ptr.
    //--------------------------------------------------------------------------------------------------------
    void AddFactory (
        const std::string & descriptor,
        SweepPatternFactory & factory
    );

    //--------------------------------------------------------------------------------------------------------
    //! \brief  Virtual interface for derived factory classes to create objects of their respective types.
    //--------------------------------------------------------------------------------------------------------
    virtual SweepPattern * Create (
        const SweepOperator<OrdinateFlux> & enclosing,
        const ParameterList & input_list
    ) const = 0;


    //========================================================================================================
    //=== INTERNAL VARIABLES =================================================================================
    //========================================================================================================

    //!
    //! \brief  Pointer to map for mapping between string descriptors and factories for classes deriving from
    //!         SweepPattern.
    //!
    //! Factories are implemented as singletons and are responsible for adding themselves to this map
    //! through the SweepPattern::AddFactory method.
    //!
    static std::map< std::string, SweepPatternFactory * > * factory_map_ptr;

};


# endif // ifndef __ABSTRACT__SWEEP_PATTERN_HPP__

