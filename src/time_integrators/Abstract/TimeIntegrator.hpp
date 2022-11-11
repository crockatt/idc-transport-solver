//------------------------------------------------------------------------------------------------------------
//! \file   time_integrators/Abstract/TimeIntegrator.hpp
//! \brief  Header for Abstract::TimeIntegrator class.
//!
//! \author Michael Crockatt
//! \date   February 2018
//------------------------------------------------------------------------------------------------------------

# ifndef __ABSTRACT__TIME_INTEGRATOR_HPP__
# define __ABSTRACT__TIME_INTEGRATOR_HPP__


# include <map>
# include <string>

# include "objects/DomainDecomposition.hpp"
# include "objects/RKDG/OrdinateFlux.hpp"
# include "objects/RKDG/CrossSection.hpp"


//------------------------------------------------------------------------------------------------------------
//! \brief  Specifies the type of a time integrator.
//!
//! Used to determine which class derived from Abstract::TimeIntegrator should be used when initializing the
//! time integrator for a simulation.
//------------------------------------------------------------------------------------------------------------
enum class TimeIntegratorType {

    //!
    //! String descriptor: \e "none"
    //!
    //! Default null value. Not a valid time integrator. Cannot be used to initialize a time integrator class.
    //!
    None,

    //!
    //! String descriptor: \e "idc"
    //!
    //! Integral deferred correction time integration methods.
    //!
    //! \see    IDCIntegrator
    //!
    IDC,

    //!
    //! String descriptor: \e "ls-idc"
    //!
    //! Low-storage integral deferred correction time integration methods.
    //!
    //! \see    LSIDCIntegrator
    //!
    LSIDC,

    //!
    //! String descriptor: \e "dirk"
    //!
    //! Diagonally implicit Runge-Kutta (DIRK) time integration methods.
    //!
    //! \see    DIRKIntegrator
    //!
    DIRK,

    //!
    //! String descriptor: \e "stdg"
    //!
    //! Space-time discontinuous Galerkin (STDG) methods.
    //!
    //! \see    STDGIntegrator
    //!
    STDG,

    //!
    //! String descriptor: \e "erk"
    //!
    //! Explicit Runge-Kutta (ERK) time integration methods.
    //!
    //! \see    ERKIntegrator
    //!
    ERK
};


//------------------------------------------------------------------------------------------------------------
//! \brief  Used to convert TimeIntegratorType values to their corresponding string descriptors.
//------------------------------------------------------------------------------------------------------------
const std::map< TimeIntegratorType, std::string > TimeIntegratorType_to_String = {

    { TimeIntegratorType::None,             "none"              },
    { TimeIntegratorType::IDC,              "idc"               },
    { TimeIntegratorType::LSIDC,            "ls-idc"            },
    { TimeIntegratorType::DIRK,             "dirk"              },
    { TimeIntegratorType::STDG,             "stdg"              },
    { TimeIntegratorType::ERK,              "erk"               }
};


//------------------------------------------------------------------------------------------------------------
//! \brief  Used to convert string descriptors to TimeIntegratorType values.
//------------------------------------------------------------------------------------------------------------
const std::map< std::string, TimeIntegratorType > String_to_TimeIntegratorType = {

    { "none",               TimeIntegratorType::None            },
    { "idc",                TimeIntegratorType::IDC             },
    { "ls-idc",             TimeIntegratorType::LSIDC           },
    { "dirk",               TimeIntegratorType::DIRK            },
    { "stdg",               TimeIntegratorType::STDG            },
    { "erk",                TimeIntegratorType::ERK             }
};


//------------------------------------------------------------------------------------------------------------
//! \brief  Declaration of signal-handling function for checkpoint.
//------------------------------------------------------------------------------------------------------------
extern "C"
void CheckpointSignalHandler( int );


namespace Abstract {


//------------------------------------------------------------------------------------------------------------
//! \brief  Abstract time integrator class defining interface for time integration schemes.
//!
//!
//------------------------------------------------------------------------------------------------------------
class TimeIntegrator : public DomainDecomposition {

    //!
    //! \brief  Friend declaration for signal handler.
    //!
    friend void ::CheckpointSignalHandler( int );

public:

    //========================================================================================================
    //=== CONSTRUCTORS, DESTRUCTOR, AND RECONFIGURATION ROUTINES =============================================
    //========================================================================================================

    //--------------------------------------------------------------------------------------------------------
    //! \brief  Delete default constructor.
    //--------------------------------------------------------------------------------------------------------
    TimeIntegrator( void ) = delete;

    //--------------------------------------------------------------------------------------------------------
    //! \brief  Delete copy constructor.
    //--------------------------------------------------------------------------------------------------------
    TimeIntegrator( const TimeIntegrator & ) = delete;

    //--------------------------------------------------------------------------------------------------------
    //! \brief  Constructs an Abstract::TimeIntegrator object with the given set of parameters.
    //--------------------------------------------------------------------------------------------------------
    TimeIntegrator( const ParameterList & );

    //--------------------------------------------------------------------------------------------------------
    //! \brief  Destructor for class Abstract::TimeIntegrator.
    //--------------------------------------------------------------------------------------------------------
    virtual ~TimeIntegrator( void ) {};

    //--------------------------------------------------------------------------------------------------------
    //! \brief  Creates an integrator of the specified type.
    //--------------------------------------------------------------------------------------------------------
    static TimeIntegrator * CreateIntegrator( const ParameterList & );


    //========================================================================================================
    //=== OPERATOR OVERLOADS =================================================================================
    //========================================================================================================

    //--------------------------------------------------------------------------------------------------------
    //! \brief  Delete copy assignment operator.
    //--------------------------------------------------------------------------------------------------------
    TimeIntegrator & operator=( const TimeIntegrator & ) = delete;


    //========================================================================================================
    //=== INTERFACE ROUTINES FOR TIME MARCHING ===============================================================
    //========================================================================================================

    //--------------------------------------------------------------------------------------------------------
    //! \brief  Computes a single timestep with step size \pp{step_size}.
    //--------------------------------------------------------------------------------------------------------
    virtual TimeIntegrator & Step (
        RKDG::OrdinateFlux & initial_condition,
        const RKDG::OrdinateFlux & source,
        const RKDG::CrossSection & sigma_t,
        const RKDG::CrossSection & sigma_s,
        const double step_size
    ) = 0;


    //--------------------------------------------------------------------------------------------------------
    //! \brief  Computes a single timestep with step size \pp{step_size}.
    //--------------------------------------------------------------------------------------------------------
    virtual TimeIntegrator & Step (
        RKDG::OrdinateFlux & u_initial_condition,
        RKDG::OrdinateFlux * c_initial_condition,
        const RKDG::OrdinateFlux & source,
        const RKDG::CrossSection & sigma_t,
        const RKDG::CrossSection & sigma_s,
        const double step_size
    ) = 0;


    //--------------------------------------------------------------------------------------------------------
    //! \brief  Integrates from time zero to time \pp{final_time} using a step size \pp{max_step_size}.
    //--------------------------------------------------------------------------------------------------------
    TimeIntegrator & Integrate (
        RKDG::OrdinateFlux & initial_condition,
        const RKDG::OrdinateFlux & source,
        const RKDG::CrossSection & sigma_t,
        const RKDG::CrossSection & sigma_s,
        const double final_time
    );

    //--------------------------------------------------------------------------------------------------------
    //! \brief  Integrates from time zero to time \pp{final_time} using a step size \pp{max_step_size}.
    //--------------------------------------------------------------------------------------------------------
    TimeIntegrator & Integrate (
        RKDG::OrdinateFlux & u_initial_condition,
        RKDG::OrdinateFlux * c_initial_condition,
        const RKDG::OrdinateFlux & source,
        const RKDG::CrossSection & sigma_t,
        const RKDG::CrossSection & sigma_s,
        const double final_time
    );


    //========================================================================================================
    //=== STATIC MEMBER FUNCTIONS ============================================================================
    //========================================================================================================

    //--------------------------------------------------------------------------------------------------------
    //! \brief  Prepares a list for use when constructing objects for the collided component.
    //--------------------------------------------------------------------------------------------------------
    static ParameterList MakeCollidedList( const ParameterList & );

    //--------------------------------------------------------------------------------------------------------
    //! \brief  Prepares a list for use when constructing objects for the uncollided component.
    //--------------------------------------------------------------------------------------------------------
    static ParameterList MakeUncollidedList( const ParameterList & );


    //========================================================================================================
    //=== ADDITIONAL MEMBER FUNCTIONS ========================================================================
    //========================================================================================================

    //--------------------------------------------------------------------------------------------------------
    //! \brief  Returns the key used to search input lists for determining the time integrator type to
    //!         construct.
    //--------------------------------------------------------------------------------------------------------
    static std::string GetInputKey( void );

    //--------------------------------------------------------------------------------------------------------
    //! \brief  Prints a summary of the integrator configuration to the logging interface.
    //--------------------------------------------------------------------------------------------------------
    void Print( const std::string & = "  " ) const override;

    //--------------------------------------------------------------------------------------------------------
    //! \brief  Returns true if the integrator has been configured with some form of hybrid splitting.
    //--------------------------------------------------------------------------------------------------------
    virtual bool IsHybrid( void ) const = 0;

    //--------------------------------------------------------------------------------------------------------
    //! \brief  Returns the TimeIntegratorType of the object.
    //--------------------------------------------------------------------------------------------------------
    virtual TimeIntegratorType GetIntegratorType( void ) const = 0;

    //--------------------------------------------------------------------------------------------------------
    //! \brief  Outputs the solution at the current time.
    //--------------------------------------------------------------------------------------------------------
    void OutputSolution(
        const RKDG::OrdinateFlux &,
        const RKDG::OrdinateFlux * const = nullptr,
        const bool = false
    ) const;

    //--------------------------------------------------------------------------------------------------------
    //! \brief  Returns the total walltime elapsed while executing inside the time-dependent solver.
    //--------------------------------------------------------------------------------------------------------
    double GetSolveTime( void ) {  return this->solve_time.GetElapsedTime();  }


protected:

    //========================================================================================================
    //=== PROTECTED MEMBER VARIABLES =========================================================================
    //========================================================================================================

    bool output_steps;          //!< Specifies whether or not output should happen after each timestep.
    bool output_sf;             //!< Specifies whether or not output should include the scalar flux distribution.
    bool output_af;             //!< Specifies whether or not output should include the angular flux distribution.

    int64_t current_step;       //!< Number of current step.
    double current_time;        //!< Current time.

    const double CFL;           //!< CFL number used to determine default step size.
    const double max_step_size; //!< Default step size.

    Walltimer solve_time;       //!< Timer for total solve time.

    int64_t checkpoint_steps;   //!< Number of steps between checkpoints.

    double checkpoint_wtime;    //!< Walltime between checkpoints.
    double last_wtime;          //!< Walltime at the end of the last timestep.
    double est_wtime;           //!< Current walltime estimate for step.
    const double A_wtime;       //!< Weight for updating walltime estimate.
    bool first_step;            //!< Used to keep track of whether it is the first step or not, for walltime estimate.


    //========================================================================================================
    //=== PROTECTED MEMBER FUNCTIONS =========================================================================
    //========================================================================================================

    //--------------------------------------------------------------------------------------------------------
    //! \brief  Writes a checkpoint file for the current time.
    //--------------------------------------------------------------------------------------------------------
    void WriteCheckpoint(
        const RKDG::OrdinateFlux &,
        const RKDG::OrdinateFlux * const = nullptr
    ) const;

    //--------------------------------------------------------------------------------------------------------
    //! \brief  Loads program state from a checkpoint file.
    //--------------------------------------------------------------------------------------------------------
    void ResumeCheckpoint(
        RKDG::OrdinateFlux &,
        RKDG::OrdinateFlux * const = nullptr
    );

    //--------------------------------------------------------------------------------------------------------
    //! \brief  Sets a signal handler to checkpoint the program state using the given objects.
    //--------------------------------------------------------------------------------------------------------
    void SetSignalHandler(
        RKDG::OrdinateFlux &,
        RKDG::OrdinateFlux * const = nullptr
    );

    //--------------------------------------------------------------------------------------------------------
    //! \brief  Checks walltime estimate to see if we need to checkpoint.
    //--------------------------------------------------------------------------------------------------------
    bool CheckWalltimeCheckpoint( void );

};


} // namespace Abstract


# endif // ifndef __ABSTRACT__TIME_INTEGRATOR_HPP__
