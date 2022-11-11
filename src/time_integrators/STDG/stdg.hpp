//------------------------------------------------------------------------------------------------------------
//! \file   time_integrators/STDG/stdg.hpp
//! \brief  Header file for space-time DG (STDG) implementation.
//!
//! \author Michael M. Crockatt
//! \date   February 2018
//------------------------------------------------------------------------------------------------------------

# ifndef __STDG_INTEGRATOR_HPP__
# define __STDG_INTEGRATOR_HPP__


# include <memory>
# include <cstdint>

# include "linear_solvers/Abstract/ImplicitSolver.hpp"
# include "utils/global.hpp"
# include "objects/RKDG/CrossSection.hpp"
# include "objects/RKDG/OrdinateFlux.hpp"
# include "objects/STDG/OrdinateFlux.hpp"
# include "objects/STDG/ScalarFlux.hpp"
# include "operators/RelabelOperator.hpp"
# include "time_integrators/Abstract/TimeIntegrator.hpp"


//------------------------------------------------------------------------------------------------------------
//! \brief  Class for STDG time integration scheme.
//!
//!
//------------------------------------------------------------------------------------------------------------
class STDGIntegrator : public Abstract::TimeIntegrator {

public:


    //--------------------------------------------------------------------------------------------------------
    //! \brief  Specifies the type of hybrid splitting (if any) used by the STDG integrator.
    //--------------------------------------------------------------------------------------------------------
    enum class HybridMethod {

        //!
        //! String descriptor: \e "none"
        //!
        //! No hybrid splitting used. The transport system is solved using a standard discrete ordinates
        //! approximation.
        //!
        //! \see    STDGIntegrator::Step_Nonhybrid()
        //!
        None,

        //!
        //! String descriptor: \e "hybrid-a"
        //!
        //! Applies a hybrid method that, in the case of DG time integration, can be equivalently derived by
        //! splitting the transport equation at either the continuum level or the discrete level.
        //! Relabeling of the collided flux is (optionally) performed at the end of each timestep via an
        //! interpolatory reconstruction procedure.
        //!
        //! The string descriptors \e "hybrid-ia" and \e "hybrid-iia" will be converted to hybrid methods of
        //! this type.
        //!
        //! \see    STDGIntegrator::Step_Hybrid()
        //!
        Hybrid_a,

        //!
        //! String descriptor: \e "hybrid-c"
        //!
        //! Applies a hybrid method that, in the case of DG time integration, can be equivalently derived by
        //! splitting the transport equation at either the continuum level or the discrete level.
        //! Relabeling of the collided flux is (optionally) performed at the end of each timestep via a
        //! reconstruction procedure motivated by the Nyström interpolation technique.
        //!
        //! The string descriptors \e "hybrid-ic" and \e "hybrid-iic" will be converted to hybrid methods of
        //! this type.
        //!
        //! \see    STDGIntegrator::Step_Hybrid()
        //!
        Hybrid_c
    };


    //--------------------------------------------------------------------------------------------------------
    //! \brief  Used to convert STDGIntegrator::HybridMethod values to their corresponding string descriptors.
    //--------------------------------------------------------------------------------------------------------
    static std::map< STDGIntegrator::HybridMethod, std::string > HybridMethod_to_String;


    //--------------------------------------------------------------------------------------------------------
    //! \brief  Used to convert string descriptors to STDGIntegrator::HybridMethod values.
    //--------------------------------------------------------------------------------------------------------
    static std::map< std::string, STDGIntegrator::HybridMethod > String_to_HybridMethod;


    //========================================================================================================
    //=== CONSTRUCTORS AND DESTRUCTOR ========================================================================
    //========================================================================================================


    //--------------------------------------------------------------------------------------------------------
    //! \brief  Delete default constructor.
    //--------------------------------------------------------------------------------------------------------
    STDGIntegrator( void ) = delete;


    //--------------------------------------------------------------------------------------------------------
    //! \brief  Delete copy constructor.
    //--------------------------------------------------------------------------------------------------------
    STDGIntegrator( const STDGIntegrator & ) = delete;


    //--------------------------------------------------------------------------------------------------------
    //! \brief  Constructs an STDGIntegrator object using the provided values.
    //--------------------------------------------------------------------------------------------------------
    STDGIntegrator( const ParameterList & );


    //--------------------------------------------------------------------------------------------------------
    //! \brief  Destructor for class STDGIntegrator.
    //--------------------------------------------------------------------------------------------------------
    ~STDGIntegrator( void ) override;


    //========================================================================================================
    //=== OPERATOR OVERLOADS =================================================================================
    //========================================================================================================


    //--------------------------------------------------------------------------------------------------------
    //! \brief  Delete copy assignment operator.
    //--------------------------------------------------------------------------------------------------------
    STDGIntegrator & operator=( const STDGIntegrator & ) = delete;


    //========================================================================================================
    //=== INTERFACE ROUTINES =================================================================================
    //========================================================================================================


    //--------------------------------------------------------------------------------------------------------
    //! \brief  Interface routine for computing a single timestep with step size \pp{step_size}.
    //--------------------------------------------------------------------------------------------------------
    STDGIntegrator & Step (

        RKDG::OrdinateFlux & initial_condition,
        const RKDG::OrdinateFlux & source,
        const RKDG::CrossSection & sigma_t,
        const RKDG::CrossSection & sigma_s,
        const double step_size

    ) override;


    //--------------------------------------------------------------------------------------------------------
    //! \brief  Interface routine for computing a single timestep with step size \pp{step_size}.
    //--------------------------------------------------------------------------------------------------------
    STDGIntegrator & Step (

        RKDG::OrdinateFlux & u_initial_condition,
        RKDG::OrdinateFlux * c_initial_condition,
        const RKDG::OrdinateFlux & source,
        const RKDG::CrossSection & sigma_t,
        const RKDG::CrossSection & sigma_s,
        const double step_size

    ) override;


    //--------------------------------------------------------------------------------------------------------
    //! \brief  Prints a summary of the integrator configuration to the logging interface.
    //--------------------------------------------------------------------------------------------------------
    void Print( const std::string & = "  " ) const override;


    //--------------------------------------------------------------------------------------------------------
    //! \brief  Returns true if the integrator has been constructed with any form of hybrid splitting and
    //!         false otherwise.
    //--------------------------------------------------------------------------------------------------------
    bool IsHybrid( void ) const override;


    //--------------------------------------------------------------------------------------------------------
    //! \brief  Returns the TimeIntegratorType of the object.
    //--------------------------------------------------------------------------------------------------------
    TimeIntegratorType GetIntegratorType( void ) const override {  return TimeIntegratorType::STDG;  }


private:

    //========================================================================================================
    //=== PRIVATE MEMBER VARIABLES ===========================================================================
    //========================================================================================================


    STDG::ScalarFlux phi;

    // Nonhybrid fluxes and solvers.
    STDG::OrdinateFlux * psi;           //!< Stores the STDG coefficients of the angular flux.

    std::shared_ptr<ImplicitSolver<STDG::OrdinateFlux>> solver;

    // Variables for hybrid methods.
    HybridMethod hybrid_method;         //!< Specifies the type of hybrid splitting (if any) to use.

    STDG::OrdinateFlux * uncollided;    //!< Stores the STDG coefficients of the uncollided flux.
    STDG::OrdinateFlux * collided;      //!< Stores the STDG coefficients of the collided flux.

    std::shared_ptr<ImplicitSolver<STDG::OrdinateFlux>> u_solver;   //!< Solver for uncollided systems.
    std::shared_ptr<ImplicitSolver<STDG::OrdinateFlux>> c_solver;   //!< Solver for collided systems.

    bool relabel = true;                //!< Specifies whether relabeling should be performed or not. Defaults to true.
    RelabelOperator * relabel_operator; //!< RelabelOperator for mapping collided fluxes to uncollided.


    //========================================================================================================
    //=== PRIVATE TIME INTEGRATION ROUTINES ==================================================================
    //========================================================================================================


    //--------------------------------------------------------------------------------------------------------
    //! \brief  Computes a one timestep update without hybrid splitting.
    //--------------------------------------------------------------------------------------------------------
    void Step_Nonhybrid (

        RKDG::OrdinateFlux & initial_condition,
        const RKDG::OrdinateFlux & source,
        const RKDG::CrossSection & sigma_t,
        const RKDG::CrossSection & sigma_s,
        const double h
    );


    //--------------------------------------------------------------------------------------------------------
    //! \brief  Computes a one timestep update using hybrid splitting.
    //--------------------------------------------------------------------------------------------------------
    void Step_Hybrid (

        RKDG::OrdinateFlux & u_initial_condition,
        RKDG::OrdinateFlux * c_initial_condition,
        const RKDG::OrdinateFlux & source,
        const RKDG::CrossSection & sigma_t,
        const RKDG::CrossSection & sigma_s,
        const double h
    );
};


# endif // ifndef __STDG_INTEGRATOR_HPP__
