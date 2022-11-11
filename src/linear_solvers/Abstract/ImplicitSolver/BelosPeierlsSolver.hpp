//------------------------------------------------------------------------------------------------------------
//! \file   linear_solvers/Abstract/ImplicitSolver/BelosPeierlsSolver.hpp
//! \brief  Header file containing declaration of BelosPeierlsSolver class template.
//!
//! \author Michael M. Crockatt
//! \date   March 2018
//------------------------------------------------------------------------------------------------------------

# ifndef __ABSTRACT__BELOS_PEIERLS_SOLVER_HPP__
# define __ABSTRACT__BELOS_PEIERLS_SOLVER_HPP__

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
//!         \f$ (\mathcal{I} - \mathcal{PL}^{-1}\mathcal{S}) \f$ using Belos' GMRES implementation.
//!
//! \tparam     AngularFlux     Type of angular flux vector. The appropriate type for scalar flux vectors is
//!                             deduced from this type.
//!
//------------------------------------------------------------------------------------------------------------
template< class AngularFlux >
class BelosPeierlsSolver :
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
    friend typename ImplicitSolver<AngularFlux>::template ImplicitSolverDerivedFactory<BelosPeierlsSolver>;

    //!
    //! \brief  Declare operator for solve as nested class.
    //!
    class IPLS;


    //========================================================================================================
    //=== CONSTRUCTORS AND DESTRUCTOR ========================================================================
    //========================================================================================================

protected:

    //--------------------------------------------------------------------------------------------------------
    //! \brief  Delete default constructor.
    //--------------------------------------------------------------------------------------------------------
    BelosPeierlsSolver( void ) = delete;

    //--------------------------------------------------------------------------------------------------------
    //! \brief  Delete copy constructor.
    //--------------------------------------------------------------------------------------------------------
    BelosPeierlsSolver( const BelosPeierlsSolver<AngularFlux> & ) = delete;

    //--------------------------------------------------------------------------------------------------------
    //! \brief  Initializes a BelosPeierlsSolver object.
    //--------------------------------------------------------------------------------------------------------
    BelosPeierlsSolver (
        const DomainDecomposition & domain_decomposition,
        const Quadrule::OrdinateSet & ordinate_set,
        const ParameterList & input_list,
        const double dt = std::numeric_limits<double>::infinity()
    );

public:

    //--------------------------------------------------------------------------------------------------------
    //! \brief  Virtual destructor.
    //--------------------------------------------------------------------------------------------------------
    virtual ~BelosPeierlsSolver( void );


    //========================================================================================================
    //=== OPERATOR OVERLOADS =================================================================================
    //========================================================================================================

    //--------------------------------------------------------------------------------------------------------
    //! \brief  Delete copy assignment operator.
    //--------------------------------------------------------------------------------------------------------
    BelosPeierlsSolver<AngularFlux> & operator=( const BelosPeierlsSolver<AngularFlux> & ) = delete;


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

    std::shared_ptr<AngularFlux> psi;               //!< Pointer to temporary angular flux vector.
    std::vector< std::shared_ptr<ScalarFlux> > phi; //!< Pointer to temporary scalar flux vector(s).

    Teuchos::RCP<Teuchos::ParameterList> solver_params;     //!< Parameters for Belos GMRES solver.
    Teuchos::RCP<BelosProblem> belos_problem;               //!< Belos::LinearProblem context.
    Teuchos::RCP<BelosSolver> belos_solver;                 //!< Belos::SolverManager context.
    Teuchos::RCP<TpetraOperator> trilinos_mat;              //!< Tpetra::Operator defining matrix to invert.
    Teuchos::RCP<TpetraVec> trilinos_vec_x;                 //!< Tpetra::MultiVector for the solution.
    Teuchos::RCP<TpetraVec> trilinos_vec_b;                 //!< Tpetra::MultiVector for right-hand vector of system.

    //!
    //! \brief  Pointer to SweepOperator object for performing transport sweeps.
    //!
    std::shared_ptr<SweepOperator<AngularFlux>> sweep_operator;

};


//------------------------------------------------------------------------------------------------------------
//! \brief  Operator class used to define interface for applying the operator
//!         \f$ (\mathcal{I} - \mathcal{PL}^{-1}\mathcal{S}) \f$ to a given vector
//!         \f$ \mathcal{P} \Psi = \Phi \f$ to obtain vectors for the Krylov basis used by the GMRES solver.
//------------------------------------------------------------------------------------------------------------
template< class AngularFlux >
class BelosPeierlsSolver<AngularFlux>::IPLS :
    public TpetraOperator
{

public:

    //========================================================================================================
    //=== CONSTRUCTORS AND DESTRUCTOR ========================================================================
    //========================================================================================================

    //--------------------------------------------------------------------------------------------------------
    //! \brief  Delete default constructor.
    //--------------------------------------------------------------------------------------------------------
    IPLS( void ) = delete;

    //--------------------------------------------------------------------------------------------------------
    //! \brief  Constructs an IPLS object from the given parameters.
    //--------------------------------------------------------------------------------------------------------
    IPLS( const BelosPeierlsSolver<AngularFlux> & );

    //--------------------------------------------------------------------------------------------------------
    //! \brief  Destructor for class Abstract::IPLS.
    //--------------------------------------------------------------------------------------------------------
    virtual ~IPLS( void ) {/* empty */};


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
    //!         \f$ \vec{y} \gets \beta \vec{y} + \alpha (\mathcal{I} - \mathcal{PL}^{-1}\mathcal{S}) \vec{x} \f$.
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
    const BelosPeierlsSolver<AngularFlux> & enc;
};


# endif // if defined (ENABLE_BELOS_TPETRA)
# endif // ifndef __ABSTRACT__BELOS_PEIERLS_SOLVER_HPP__
