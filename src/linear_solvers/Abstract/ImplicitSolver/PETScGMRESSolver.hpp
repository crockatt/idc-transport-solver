//------------------------------------------------------------------------------------------------------------
//! \file   linear_solvers/Abstract/ImplicitSolver/PETScGMRESSolver.hpp
//! \brief  Header file containing declaration of PETScGMRESSolver class template.
//!
//! \author Michael M. Crockatt
//! \date   March 2018
//------------------------------------------------------------------------------------------------------------

# ifndef __ABSTRACT__PETSC_GMRES_SOLVER_HPP__
# define __ABSTRACT__PETSC_GMRES_SOLVER_HPP__

# if defined (ENABLE_PETSC) || defined (DOXYCOMPILE)


# include <type_traits>
# include <vector>

# include <petscksp.h>

# include "linear_solvers/Abstract/ImplicitSolver.hpp"
# include "linear_solvers/Abstract/SweepOperator.hpp"
# include "objects/RKDG/OrdinateFlux.hpp"
# include "objects/RKDG/ScalarFlux.hpp"
# include "objects/STDG/OrdinateFlux.hpp"
# include "objects/STDG/ScalarFlux.hpp"
# include "utils/PETScInterface.hpp"


//------------------------------------------------------------------------------------------------------------
//! \brief  Class template for solver based on inverting the operator
//!         \f$ (\mathcal{I} - \mathcal{PL}^{-1}\mathcal{S}) \f$ using PETSc's GMRES implementation.
//!
//! \tparam     AngularFlux     Type of angular flux vector. The appropriate type for scalar flux vectors is
//!                             deduced from this type.
//!
//------------------------------------------------------------------------------------------------------------
template< class AngularFlux >
class PETScGMRESSolver :
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
    friend typename ImplicitSolver<AngularFlux>::template ImplicitSolverDerivedFactory<PETScGMRESSolver>;


    //========================================================================================================
    //=== CONSTRUCTORS AND DESTRUCTOR ========================================================================
    //========================================================================================================

protected:

    //--------------------------------------------------------------------------------------------------------
    //! \brief  Delete default constructor.
    //--------------------------------------------------------------------------------------------------------
    PETScGMRESSolver( void ) = delete;

    //--------------------------------------------------------------------------------------------------------
    //! \brief  Delete copy constructor.
    //--------------------------------------------------------------------------------------------------------
    PETScGMRESSolver( const PETScGMRESSolver<AngularFlux> & ) = delete;

    //--------------------------------------------------------------------------------------------------------
    //! \brief  Initializes an PETScGMRESSolver object.
    //--------------------------------------------------------------------------------------------------------
    PETScGMRESSolver (
        const DomainDecomposition & domain_decomposition,
        const Quadrule::OrdinateSet & ordinate_set,
        const ParameterList & input_list,
        const double dt = std::numeric_limits<double>::infinity()
    );

public:

    //--------------------------------------------------------------------------------------------------------
    //! \brief  Virtual destructor.
    //--------------------------------------------------------------------------------------------------------
    virtual ~PETScGMRESSolver( void );


    //========================================================================================================
    //=== OPERATOR OVERLOADS =================================================================================
    //========================================================================================================

    //--------------------------------------------------------------------------------------------------------
    //! \brief  Delete copy assignment operator.
    //--------------------------------------------------------------------------------------------------------
    PETScGMRESSolver<AngularFlux> & operator=( const PETScGMRESSolver<AngularFlux> & ) = delete;


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


    //========================================================================================================
    //=== ADDITIONAL HELPER FUNCTIONS ========================================================================
    //========================================================================================================

    //--------------------------------------------------------------------------------------------------------
    //! \brief  Applies the operator \f$ (\mathcal{I} - \mathcal{PL}^{-1}\mathcal{S}) \f$ to a given vector.
    //--------------------------------------------------------------------------------------------------------
    static PetscErrorCode ILLSP( PETScMat, PETScVec, PETScVec );


protected:

    //========================================================================================================
    //=== INTERNAL VARIABLES =================================================================================
    //========================================================================================================

    double dt_save;                         //!< Timestep size for current implicit solve.
    const CrossSection * sigma_t_save;      //!< Total cross-section \f$ \sigma_{\mathrm{t}} \f$.
    const CrossSection * sigma_s_save;      //!< Scattering cross-section \f$ \sigma_{\mathrm{s}} \f$.

    //!
    //! \brief  Pointer to currently executing solver object.
    //!
    //! This is required in order to obtain the current solver context containing solve-specific
    //! parameters (PETScGMRESSolver::dt, PETScGMRESSolver::sigma_t_save, and 
    //! PETScGMRESSolver::sigma_s_save) due to the requirement that a PETSc shell matrix is defined with a
    //! function pointer.
    //!
    //! \attention  Because this declared static <em>only one object of each instantiated class type can
    //!             execute a solve operation at any one time.</em>
    //!
    static PETScGMRESSolver<AngularFlux> * instance;

    std::shared_ptr<AngularFlux> psi;               //!< Pointer to temporary angular flux vector.
    std::shared_ptr<ScalarFlux> phi;                //!< Pointer to temporary scalar flux vector.

    // Variables for transport solve.
    std::shared_ptr<ParameterList> transp_params;   //!< Contains current parameters for transport solver.
    PETScKSP transp_solver;                         //!< PETSc iterative solver context.
    PETScPC transp_precond;                         //!< PETSc preconditioner context.
    PETScMat transp_mat;                            //!< PETSc matrix object to be inverted.
    PETScVec transp_vec_sol;                        //!< PETSc vector for storing solution of system.
    PETScVec transp_vec_rhs;                        //!< PETSc vector for holding the right side of the system.

    //!
    //! \brief  Pointer to SweepOperator object for performing transport sweeps.
    //!
    std::shared_ptr<SweepOperator<AngularFlux>> sweep_operator;

};


# endif // if defined (ENABLE_PETSC)
# endif // ifndef __ABSTRACT__PETSC_GMRES_SOLVER_HPP__
