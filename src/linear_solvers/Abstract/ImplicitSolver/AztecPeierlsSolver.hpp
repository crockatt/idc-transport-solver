//------------------------------------------------------------------------------------------------------------
//! \file   linear_solvers/Abstract/ImplicitSolver/AztecPeierlsSolver.hpp
//! \brief  Header file containing declaration of AztecPeierlsSolver class template.
//!
//! \author Michael M. Crockatt
//! \date   March 2018
//------------------------------------------------------------------------------------------------------------

# ifndef __ABSTRACT__AZTEC_PEIERLS_SOLVER_HPP__
# define __ABSTRACT__AZTEC_PEIERLS_SOLVER_HPP__

# if defined (ENABLE_AZTEC_EPETRA) || defined (DOXYCOMPILE)


# include <memory>
# include <type_traits>
# include <vector>

# include <AztecOO.h>
# include <Epetra_LinearProblem.h>
# include <Epetra_Operator.h>
# include <Epetra_MultiVector.h>

# include "linear_solvers/Abstract/ImplicitSolver.hpp"
# include "linear_solvers/Abstract/SweepOperator.hpp"
# include "objects/RKDG/OrdinateFlux.hpp"
# include "objects/RKDG/ScalarFlux.hpp"
# include "objects/STDG/OrdinateFlux.hpp"
# include "objects/STDG/ScalarFlux.hpp"
# include "utils/BelosTpetraInterface.hpp"


//------------------------------------------------------------------------------------------------------------
//! \brief  Class template for solver based on inverting the operator
//!         \f$ (\mathcal{I} - \mathcal{PL}^{-1}\mathcal{S}) \f$ using AztecOO's GMRES implementation.
//!
//! \tparam     AngularFlux     Type of angular flux vector. The appropriate type for scalar flux vectors is
//!                             deduced from this type.
//!
//------------------------------------------------------------------------------------------------------------
template< class AngularFlux >
class AztecPeierlsSolver :
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
    friend typename ImplicitSolver<AngularFlux>::template ImplicitSolverDerivedFactory<AztecPeierlsSolver>;

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
    AztecPeierlsSolver( void ) = delete;

    //--------------------------------------------------------------------------------------------------------
    //! \brief  Delete copy constructor.
    //--------------------------------------------------------------------------------------------------------
    AztecPeierlsSolver( const AztecPeierlsSolver<AngularFlux> & ) = delete;

    //--------------------------------------------------------------------------------------------------------
    //! \brief  Initializes a AztecPeierlsSolver object.
    //--------------------------------------------------------------------------------------------------------
    AztecPeierlsSolver (
        const DomainDecomposition & domain_decomposition,
        const Quadrule::OrdinateSet & ordinate_set,
        const ParameterList & input_list,
        const double dt = std::numeric_limits<double>::infinity()
    );

public:

    //--------------------------------------------------------------------------------------------------------
    //! \brief  Virtual destructor.
    //--------------------------------------------------------------------------------------------------------
    virtual ~AztecPeierlsSolver( void );


    //========================================================================================================
    //=== OPERATOR OVERLOADS =================================================================================
    //========================================================================================================

    //--------------------------------------------------------------------------------------------------------
    //! \brief  Delete copy assignment operator.
    //--------------------------------------------------------------------------------------------------------
    AztecPeierlsSolver<AngularFlux> & operator=( const AztecPeierlsSolver<AngularFlux> & ) = delete;


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

    //!
    //! /\brief List of parameters for AztecOO GMRES transport solve.
    //!
    std::shared_ptr<Teuchos::ParameterList> transport_params;

    std::shared_ptr<Epetra_LinearProblem> transport_prob;   //!< Epetra_LinearProblem object for solve.
    std::shared_ptr<AztecOO> transport_solver;              //!< AztecOO solver object.
    std::shared_ptr<Epetra_Operator> ipls;                  //!< Epetra_Operator object for solve.
    std::shared_ptr<Epetra_MultiVector> epetra_vec_x;       //!< Epetra_MultiVector for the solution.
    std::shared_ptr<Epetra_MultiVector> epetra_vec_b;       //!< Epetra_MultiVector for right-hand vector.

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
class AztecPeierlsSolver<AngularFlux>::IPLS :
    public virtual Epetra_Operator
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
    IPLS( const AztecPeierlsSolver<AngularFlux> & );

    //--------------------------------------------------------------------------------------------------------
    //! \brief  Destructor for class Abstract::IPLS.
    //--------------------------------------------------------------------------------------------------------
    virtual ~IPLS( void ) {/* empty */};


    //========================================================================================================
    //=== ATTRIBUTE SET METHODS ==============================================================================
    //========================================================================================================

    //--------------------------------------------------------------------------------------------------------
    //! \brief  Returns -1 because this implementation of the Epetra_Operator interface does not support
    //!         transpose use.
    //--------------------------------------------------------------------------------------------------------
    int SetUseTranspose( bool ) override;


    //========================================================================================================
    //=== MATHEMATICAL FUNCTIONS =============================================================================
    //========================================================================================================

    //--------------------------------------------------------------------------------------------------------
    //! \brief  Upon return, \pp{Y} contains the result of applying the operator to \pp{X}.
    //--------------------------------------------------------------------------------------------------------
    int Apply( const Epetra_MultiVector &, Epetra_MultiVector & ) const override;

    //--------------------------------------------------------------------------------------------------------
    //! \brief  Returns -1 because this implementation does not support inverse use.
    //--------------------------------------------------------------------------------------------------------
    int ApplyInverse( const Epetra_MultiVector &, Epetra_MultiVector & ) const override;

    //--------------------------------------------------------------------------------------------------------
    //! \brief  Returns -1 because this implementation does not support use of the infinity norm.
    //--------------------------------------------------------------------------------------------------------
    double NormInf( void ) const override;


    //========================================================================================================
    //=== ATTRIBUTE ACCESS FUNCTIONS =========================================================================
    //========================================================================================================

    //--------------------------------------------------------------------------------------------------------
    //! \brief  Returns a character string describing the operator.
    //--------------------------------------------------------------------------------------------------------
    const char * Label( void ) const override;

    //--------------------------------------------------------------------------------------------------------
    //! \brief  Returns false because this implementation does not support transpose use.
    //--------------------------------------------------------------------------------------------------------
    bool UseTranspose( void ) const override;

    //--------------------------------------------------------------------------------------------------------
    //! \brief  Returns false because this object cannot provide an approximate Inf-norm.
    //--------------------------------------------------------------------------------------------------------
    bool HasNormInf( void ) const override;

    //--------------------------------------------------------------------------------------------------------
    //! \brief  Returns a reference to the Epetra_Comm communicator associated with this operator.
    //--------------------------------------------------------------------------------------------------------
    const Epetra_Comm & Comm( void ) const override;

    //--------------------------------------------------------------------------------------------------------
    //! \brief  Returns a reference to the Epetra_Map object associated with the domain of this operator.
    //--------------------------------------------------------------------------------------------------------
    const Epetra_Map & OperatorDomainMap( void ) const override;

    //--------------------------------------------------------------------------------------------------------
    //! \brief  Returns a reference to the Epetra_Map object associated with the range of this operator.
    //--------------------------------------------------------------------------------------------------------
    const Epetra_Map & OperatorRangeMap( void ) const override;


protected:

    //========================================================================================================
    //=== INTERNAL VARIABLES =================================================================================
    //========================================================================================================

    //!
    //! \brief  Character string describing the operator.
    //!
    const char label [10] = "\0";

    //!
    //! \brief  Epetra communicator for operator.
    //!
    std::shared_ptr<Epetra_Comm> comm;

    //!
    //! \brief  Epetra map associated with domain and range of operator.
    //!
    std::shared_ptr<Epetra_Map> domain_decomp_map;

    //!
    //! \brief  Reference to enclosing solver object.
    //!
    const AztecPeierlsSolver<AngularFlux> & enc;
};


# endif // if defined (ENABLE_AZTEC_EPETRA)
# endif // ifndef __ABSTRACT__AZTEC_PEIERLS_SOLVER_HPP__
