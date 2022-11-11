//------------------------------------------------------------------------------------------------------------
//! \file   linear_solvers/Abstract/ImplicitSolver.hpp
//! \brief  Header file containing declaration of ImplicitSolver class template.
//!
//! \author Michael M. Crockatt
//! \date   March 2018
//------------------------------------------------------------------------------------------------------------

# ifndef __ABSTRACT__IMPLICIT_SOLVER_HPP__
# define __ABSTRACT__IMPLICIT_SOLVER_HPP__


# include "linear_solvers/Abstract/SweepOperator.hpp"
# include "objects/DomainDecomposition.hpp"
# include "utils/Quadrule/OrdinateSet.hpp"


//------------------------------------------------------------------------------------------------------------
//! \brief  Abstract class template defining interface for implicit solvers.
//!
//! This class defines the interface used by objects of the ImplicitSolver class hierarchy for solving
//! globally-coupled systems of equations resulting from application of implicit time discretizations to the
//! transport system.
//------------------------------------------------------------------------------------------------------------
template< class AngularFlux >
class ImplicitSolver :
    public DomainDecomposition,
    public Quadrule::OrdinateSet
{

public:

    //!
    //! \brief  Declare ImplicitSolverFactory as public to expose factory interface.
    //!
    class ImplicitSolverFactory;


protected:

    //!
    //! \brief  Declare protected template factory class for derived solver classes.
    //!
    template< class DerivedSolver >
    class ImplicitSolverDerivedFactory;


    //========================================================================================================
    //=== CONSTRUCTORS AND DESTRUCTOR ========================================================================
    //========================================================================================================

    //--------------------------------------------------------------------------------------------------------
    //! \brief  Delete default constructor.
    //--------------------------------------------------------------------------------------------------------
    ImplicitSolver( void ) = delete;

    //--------------------------------------------------------------------------------------------------------
    //! \brief  Delete copy constructor.
    //--------------------------------------------------------------------------------------------------------
    ImplicitSolver( const ImplicitSolver<AngularFlux> & ) = delete;

    //--------------------------------------------------------------------------------------------------------
    //! \brief  Initializes an ImplicitSolver object.
    //--------------------------------------------------------------------------------------------------------
    ImplicitSolver (
        const DomainDecomposition & domain_decomposition,
        const Quadrule::OrdinateSet & ordinate_set,
        const ParameterList & input_list,
        const double dt = std::numeric_limits<double>::infinity()
    );

public:

    //--------------------------------------------------------------------------------------------------------
    //! \brief  Virtual destructor.
    //--------------------------------------------------------------------------------------------------------
    virtual ~ImplicitSolver( void );


    //========================================================================================================
    //=== OPERATOR OVERLOADS =================================================================================
    //========================================================================================================

    //--------------------------------------------------------------------------------------------------------
    //! \brief  Delete copy assignment operator.
    //--------------------------------------------------------------------------------------------------------
    ImplicitSolver<AngularFlux> & operator=( const ImplicitSolver<AngularFlux> & ) = delete;


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
    void Print( const std::string & = "  " ) const override;

    //--------------------------------------------------------------------------------------------------------
    //! \brief  Sets parameters given in a ParameterList.
    //--------------------------------------------------------------------------------------------------------
    virtual void SetParameters( const ParameterList & );

    //--------------------------------------------------------------------------------------------------------
    //! \brief  Sets the initial guess for the next solve.
    //--------------------------------------------------------------------------------------------------------
    virtual void SetInitialGuess( const AngularFlux * const = nullptr );

    //--------------------------------------------------------------------------------------------------------
    //! \brief  Solves a system of equations using the specified parameters.
    //--------------------------------------------------------------------------------------------------------
    virtual void Solve (
        AngularFlux & source,
        const RKDG::CrossSection & sigma_t,
        const RKDG::CrossSection & sigma_s,
        const double dt = std::numeric_limits<double>::infinity(),
        const RKDG::OrdinateFlux * const initial = nullptr
    ) = 0;

};


//------------------------------------------------------------------------------------------------------------
//! \brief  Abstract class defining interface for factories used to construct objects derived from
//!         ImplicitSolver class.
//!
//! Concrete classes derived from ImplicitSolver should use the class template ImplicitSolverDerivedFactory
//! (which derives from this class) to enable their construction.
//!
//------------------------------------------------------------------------------------------------------------
template< class AngularFlux >
class ImplicitSolver<AngularFlux>::ImplicitSolverFactory {

    //========================================================================================================
    //=== CONSTRUCTORS AND DESTRUCTOR ========================================================================
    //========================================================================================================

protected:

    ImplicitSolverFactory( void ) = default;

    ImplicitSolverFactory( const ImplicitSolverFactory & ) = delete;
    ImplicitSolverFactory & operator=( const ImplicitSolverFactory & ) = delete;

public:

    ~ImplicitSolverFactory( void ) = default;


    //========================================================================================================
    //=== PUBLIC INTERFACE ROUTINES ==========================================================================
    //========================================================================================================

    //--------------------------------------------------------------------------------------------------------
    //! \brief  Creates and returns an object of a type derived from ImplicitSolver.
    //--------------------------------------------------------------------------------------------------------
    static ImplicitSolver<AngularFlux> * CreateSolver (
        const DomainDecomposition & domain_decomposition,
        const Quadrule::OrdinateSet & ordinate_set,
        const ParameterList & input_list,
        const double dt = std::numeric_limits<double>::infinity()
    );


protected:

    //========================================================================================================
    //=== PROTECTED INTERFACE ROUTINES =======================================================================
    //========================================================================================================

    //--------------------------------------------------------------------------------------------------------
    //! \brief  Returns a reference to the map stored at ImplicitSolverFactory::factory_map_ptr.
    //--------------------------------------------------------------------------------------------------------
    static std::map< std::string, ImplicitSolverFactory * > & GetFactoryMap( void );

    //--------------------------------------------------------------------------------------------------------
    //! \brief  Adds a reference to a factory class to the map ImplicitSolverFactory::factory_map_ptr.
    //--------------------------------------------------------------------------------------------------------
    void AddFactory (
        const std::string & descriptor,
        ImplicitSolverFactory & factory
    );

    //--------------------------------------------------------------------------------------------------------
    //! \brief  Virtual interface for derived factory classes to create objects of their respective types.
    //--------------------------------------------------------------------------------------------------------
    virtual ImplicitSolver<AngularFlux> * Create (
        const DomainDecomposition & domain_decomposition,
        const Quadrule::OrdinateSet & ordinate_set,
        const ParameterList & input_list,
        const double dt = std::numeric_limits<double>::infinity()
    ) const = 0;


    //========================================================================================================
    //=== INTERNAL VARIABLES =================================================================================
    //========================================================================================================

    //!
    //! \brief  Pointer to map for mapping between string descriptors and factories for classes deriving from
    //!         ImplicitSolver.
    //!
    //! Factories are implemented as singletons and are responsible for adding themselves to this map
    //! through the ImplicitSolverFactory::AddFactory method.
    //!
    static std::map< std::string, ImplicitSolverFactory * > * factory_map_ptr;

};


# endif // ifndef __ABSTRACT__IMPLICIT_SOLVER_HPP__

