//------------------------------------------------------------------------------------------------------------
//! \file   linear_solvers/Abstract/ImplicitSolver/SweepSolver.hpp
//! \brief  Header file containing declaration of SweepSolver class template.
//!
//! \author Michael M. Crockatt
//! \date   March 2018
//------------------------------------------------------------------------------------------------------------

# ifndef __ABSTRACT__SWEEP_SOLVER_HPP__
# define __ABSTRACT__SWEEP_SOLVER_HPP__


# include "linear_solvers/Abstract/ImplicitSolver.hpp"
# include "linear_solvers/Abstract/SweepOperator.hpp"


//------------------------------------------------------------------------------------------------------------
//! \brief  Class template for implicit solver composed of a single transport sweep.
//!
//!
//------------------------------------------------------------------------------------------------------------
template< class AngularFlux >
class SweepSolver :
    public ImplicitSolver<AngularFlux>
{

    //!
    //! \brief  Typedef for objects holding initial conditions.
    //!
    typedef RKDG::OrdinateFlux AngularFluxInit;

    //!
    //! \brief  Typedef for cross section objects.
    //!
    typedef RKDG::CrossSection CrossSection;

    //!
    //! \brief  Friend instantiation of derived factory class.
    //!
    friend typename ImplicitSolver<AngularFlux>::template ImplicitSolverDerivedFactory<SweepSolver>;


    //========================================================================================================
    //=== CONSTRUCTORS AND DESTRUCTOR ========================================================================
    //========================================================================================================

protected:

    //--------------------------------------------------------------------------------------------------------
    //! \brief  Delete default constructor.
    //--------------------------------------------------------------------------------------------------------
    SweepSolver( void ) = delete;

    //--------------------------------------------------------------------------------------------------------
    //! \brief  Delete copy constructor.
    //--------------------------------------------------------------------------------------------------------
    SweepSolver( const SweepSolver<AngularFlux> & ) = delete;

    //--------------------------------------------------------------------------------------------------------
    //! \brief  Initializes a SweepSolver object.
    //--------------------------------------------------------------------------------------------------------
    SweepSolver (
        const DomainDecomposition & domain_decomposition,
        const Quadrule::OrdinateSet & ordinate_set,
        const ParameterList & input_list,
        const double dt = std::numeric_limits<double>::infinity()
    );

public:

    //--------------------------------------------------------------------------------------------------------
    //! \brief  Virtual destructor.
    //--------------------------------------------------------------------------------------------------------
    virtual ~SweepSolver( void );


    //========================================================================================================
    //=== OPERATOR OVERLOADS =================================================================================
    //========================================================================================================

    //--------------------------------------------------------------------------------------------------------
    //! \brief  Delete copy assignment operator.
    //--------------------------------------------------------------------------------------------------------
    SweepSolver<AngularFlux> & operator=( const SweepSolver<AngularFlux> & ) = delete;


    //========================================================================================================
    //=== PUBLIC INTERFACE ===================================================================================
    //========================================================================================================

    //--------------------------------------------------------------------------------------------------------
    //! \brief  Returns the string descriptor of the class type.
    //--------------------------------------------------------------------------------------------------------
    static const std::string Descriptor( void );

    //--------------------------------------------------------------------------------------------------------
    //! \brief  Returns the string descriptor for an objects type.
    //--------------------------------------------------------------------------------------------------------
    const std::string GetDescriptor( void ) const override;

    //--------------------------------------------------------------------------------------------------------
    //! \brief  Prints a summary of the object's configuration to the logging interface.
    //--------------------------------------------------------------------------------------------------------
    void Print( const std::string & = "  " ) const override;

    //--------------------------------------------------------------------------------------------------------
    //! \brief  Sets parameters given in a ParameterList.
    //--------------------------------------------------------------------------------------------------------
    void SetParameters( const ParameterList & ) override;

    //--------------------------------------------------------------------------------------------------------
    //! \brief  Solves a system of equations using the specified parameters.
    //--------------------------------------------------------------------------------------------------------
    void Solve (
        AngularFlux & source,
        const CrossSection & sigma_t,
        const CrossSection &, // sigma_s,
        const double dt = std::numeric_limits<double>::infinity(),
        const AngularFluxInit * const initial = nullptr
    ) override;


protected:

    //========================================================================================================
    //=== INTERNAL VARIABLES =================================================================================
    //========================================================================================================

    //!
    //! \brief  Pointer to SweepOperator object for performing transport sweeps.
    //!
    std::shared_ptr<SweepOperator<AngularFlux>> sweep_operator;

};


# endif // ifndef __ABSTRACT__SWEEP_SOLVER_HPP__
