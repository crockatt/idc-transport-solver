//------------------------------------------------------------------------------------------------------------
//! \file   linear_solvers/Abstract/ImplicitSolver/BelosGMRESSolver.hpp
//! \brief  Header file containing declaration of BelosGMRESSolver class template.
//!
//! \author Michael M. Crockatt
//! \date   March 2018
//------------------------------------------------------------------------------------------------------------

# ifndef __ABSTRACT__BELOS_GMRES_SOLVER_HPP__
# define __ABSTRACT__BELOS_GMRES_SOLVER_HPP__

# if defined (ENABLE_BELOS_TPETRA) || defined (DOXYCOMPILE)


# include <type_traits>
# include <vector>


# include "linear_solvers/Abstract/ImplicitSolver.hpp"
# include "linear_solvers/Abstract/SweepOperator.hpp"
# include "objects/RKDG/OrdinateFlux.hpp"
# include "objects/RKDG/ScalarFlux.hpp"
# include "objects/STDG/OrdinateFlux.hpp"
# include "objects/STDG/ScalarFlux.hpp"
# include "utils/BelosTpetraInterface.hpp"


//------------------------------------------------------------------------------------------------------------
//! \brief  Class template for solver based on inverting the operator
//!         \f$ (\mathcal{I} + \mathcal{L}_0^{-1}( \mathcal{L}_1 - \mathcal{SP})) \f$
//!         using Belos' GMRES implementation.
//!
//! \tparam     AngularFlux     Type of angular flux vector. The appropriate type for scalar flux vectors is
//!                             deduced from this type.
//!
//------------------------------------------------------------------------------------------------------------
template< class AngularFlux >
class BelosGMRESSolver :
    public ImplicitSolver<AngularFlux>
{

    //!
    //! \brief  Determine type of scalar flux vector.
    //!
    typedef
        typename std::conditional<
                std::is_same< AngularFlux, RKDG::OrdinateFlux >::value,
                RKDG::ScalarFlux,
                STDG::ScalarFlux
            >::type
        ScalarFlux;

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
    friend typename ImplicitSolver<AngularFlux>::template ImplicitSolverDerivedFactory<BelosGMRESSolver>;

    //!
    //! \brief  Declare operator for solve as nested class.
    //!
    class ILLSP;


    //========================================================================================================
    //=== CONSTRUCTORS AND DESTRUCTOR ========================================================================
    //========================================================================================================

protected:

    //--------------------------------------------------------------------------------------------------------
    //! \brief  Delete default constructor.
    //--------------------------------------------------------------------------------------------------------
    BelosGMRESSolver( void ) = delete;

    //--------------------------------------------------------------------------------------------------------
    //! \brief  Delete copy constructor.
    //--------------------------------------------------------------------------------------------------------
    BelosGMRESSolver( const BelosGMRESSolver<AngularFlux> & ) = delete;

    //--------------------------------------------------------------------------------------------------------
    //! \brief  Initializes a BelosGMRESSolver object.
    //--------------------------------------------------------------------------------------------------------
    BelosGMRESSolver (
        const DomainDecomposition & domain_decomposition,
        const Quadrule::OrdinateSet & ordinate_set,
        const ParameterList & input_list,
        const double dt = std::numeric_limits<double>::infinity()
    );

public:

    //--------------------------------------------------------------------------------------------------------
    //! \brief  Virtual destructor.
    //--------------------------------------------------------------------------------------------------------
    virtual ~BelosGMRESSolver( void );


    //========================================================================================================
    //=== OPERATOR OVERLOADS =================================================================================
    //========================================================================================================

    //--------------------------------------------------------------------------------------------------------
    //! \brief  Delete copy assignment operator.
    //--------------------------------------------------------------------------------------------------------
    BelosGMRESSolver<AngularFlux> & operator=( const BelosGMRESSolver<AngularFlux> & ) = delete;


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
    //! \brief  Sets the initial guess for the next solve.
    //--------------------------------------------------------------------------------------------------------
    void SetInitialGuess( const AngularFlux * const = nullptr ) override;

    //--------------------------------------------------------------------------------------------------------
    //! \brief  Solves a system of equations using the specified parameters.
    //--------------------------------------------------------------------------------------------------------
    void Solve (
        AngularFlux & source,
        const CrossSection & sigma_t,
        const CrossSection & sigma_s,
        const double dt = std::numeric_limits<double>::infinity(),
        const AngularFluxInit * const initial_condition = nullptr
    ) override;


protected:

    //========================================================================================================
    //=== INTERNAL VARIABLES =================================================================================
    //========================================================================================================

    double dt;                              //!< Timestep size for current implicit solve.
    const CrossSection * sigma_t;           //!< Total cross-section \f$ \sigma_{\mathrm{t}} \f$.
    const CrossSection * sigma_s;           //!< Scattering cross-section \f$ \sigma_{\mathrm{s}} \f$.

    std::shared_ptr<AngularFlux> psi;       //!< Pointer to temporary angular flux vector.
    std::shared_ptr<ScalarFlux> phi;        //!< Pointer to temporary scalar flux vector.

    Teuchos::RCP<Teuchos::ParameterList> transp_params; //!< Parameters for Belos GMRES solver.
    Teuchos::RCP<BelosProblem> transp_problem;          //!< Belos::LinearProblem context.
    Teuchos::RCP<BelosSolver> transp_solver;            //!< Belos::SolverManager context.
    Teuchos::RCP<TpetraOperator> transp_op;             //!< Tpetra::Operator defining matrix to invert.
    Teuchos::RCP<TpetraVec> transp_vec_sol;             //!< Tpetra::MultiVector for the solution.
    Teuchos::RCP<TpetraVec> transp_vec_rhs;             //!< Tpetra::MultiVector for right-hand vector of system.

    //!
    //! \brief  Pointer to SweepOperator object for performing transport sweeps.
    //!
    std::shared_ptr<SweepOperator<AngularFlux>> sweep_operator;

};


//------------------------------------------------------------------------------------------------------------
//! \brief  Operator class used to define interface for applying the operator
//!         \f$ (\mathcal{I} + \mathcal{L}_0^{-1}( \mathcal{L}_1 - \mathcal{SP})) \f$ to a given vector to
//!         obtain vectors for the Krylov basis used by the GMRES solver.
//------------------------------------------------------------------------------------------------------------
template< class AngularFlux >
class BelosGMRESSolver<AngularFlux>::ILLSP :
    public TpetraOperator
{

public:

    //========================================================================================================
    //=== CONSTRUCTORS AND DESTRUCTOR ========================================================================
    //========================================================================================================

    //--------------------------------------------------------------------------------------------------------
    //! \brief  Delete default constructor.
    //--------------------------------------------------------------------------------------------------------
    ILLSP( void ) = delete;

    //--------------------------------------------------------------------------------------------------------
    //! \brief  Constructs an ILLSP object from the given parameters.
    //--------------------------------------------------------------------------------------------------------
    ILLSP( const BelosGMRESSolver<AngularFlux> & );

    //--------------------------------------------------------------------------------------------------------
    //! \brief  Destructor for class Abstract::ILLSP.
    //--------------------------------------------------------------------------------------------------------
    virtual ~ILLSP( void ) {/* empty */};


    //========================================================================================================
    //=== INTERFACE ROUTINES =================================================================================
    //========================================================================================================

    //--------------------------------------------------------------------------------------------------------
    //! \brief  Returns the domain Map of the operator.
    //--------------------------------------------------------------------------------------------------------
    Teuchos::RCP<const TpetraMap> getDomainMap( void ) const {  return domain_decomp_map;  }

    //--------------------------------------------------------------------------------------------------------
    //! \brief  Returns the range Map of the operator.
    //--------------------------------------------------------------------------------------------------------
    Teuchos::RCP<const TpetraMap> getRangeMap( void ) const {  return domain_decomp_map;  }

    //--------------------------------------------------------------------------------------------------------
    //! \brief  Computes
    //!         \f$ \vec{y} \gets \beta \vec{y} + \alpha (\mathcal{I} + \mathcal{L}_0^{-1}( \mathcal{L}_1 - \mathcal{SP})) \vec{x} \f$.
    //--------------------------------------------------------------------------------------------------------
    virtual void apply (
        const TpetraVec & X,
        TpetraVec & Y,
        Teuchos::ETransp mode = Teuchos::NO_TRANS,
        TpetraScalar alpha = Teuchos::ScalarTraits<TpetraScalar>::one(),
        TpetraScalar beta = Teuchos::ScalarTraits<TpetraScalar>::zero()
    ) const;


protected:

    //========================================================================================================
    //=== INTERNAL VARIABLES =================================================================================
    //========================================================================================================

    //!
    //! \brief  Map over MPI ranks for operator.
    //!
    Teuchos::RCP<const TpetraMap> domain_decomp_map;

    //!
    //! \brief  Reference to enclosing solver object.
    //!
    const BelosGMRESSolver<AngularFlux> & enc;
};


# endif // if defined (ENABLE_BELOS_TPETRA)
# endif // ifndef __ABSTRACT__BELOS_GMRES_SOLVER_HPP__
