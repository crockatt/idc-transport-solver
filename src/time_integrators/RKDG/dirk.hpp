//------------------------------------------------------------------------------------------------------------
//! \file   time_integrators/RKDG/dirk.hpp
//! \brief  Header file for diagonally implicit Runge-Kutta (DIRK) time integrator implementation.
//!
//! \author Michael M. Crockatt
//! \date   February 2018
//------------------------------------------------------------------------------------------------------------

# ifndef __RKDG__DIRK_INTEGRATOR_HPP__
# define __RKDG__DIRK_INTEGRATOR_HPP__


# include <cstdint>
# include <limits>
# include <map>
# include <memory>
# include <string>

# include "linear_solvers/Abstract/ImplicitSolver.hpp"
# include "objects/RKDG/CrossSection.hpp"
# include "objects/RKDG/OrdinateFlux.hpp"
# include "objects/RKDG/ScalarFlux.hpp"
# include "operators/RelabelOperator.hpp"
# include "time_integrators/Abstract/TimeIntegrator.hpp"
# include "utils/global.hpp"


//------------------------------------------------------------------------------------------------------------
//! \brief  Class for DIRK time integration schemes.
//!
//!
//------------------------------------------------------------------------------------------------------------
class DIRKIntegrator : public Abstract::TimeIntegrator {

public:

    //--------------------------------------------------------------------------------------------------------
    //! \brief  Specifies the DIRK method.
    //!
    //! \todo   Include references for DIRK methods.
    //--------------------------------------------------------------------------------------------------------
    enum class DIRKMethod {

        //!
        //! String descriptor: \e "none"
        //!
        //! Null value. Not a valid DIRK method.
        //!
        None,

        //!
        //! String descriptor: \e "euler"
        //!
        //! Stages: 1               <br>
        //! Order: 1                <br>
        //! Stability: L-stable     <br>
        //! Stiffly accurate: yes
        //!
        Euler,

        //!
        //! String descriptor: \e "sdirk2"
        //!
        //! Stages: 2               <br>
        //! Order: 2                <br>
        //! Stability: L-stable     <br>
        //! Stiffly accurate: yes
        //!
        SDIRK2,

        //!
        //! String descriptor: \e "sdirk3"
        //!
        //! Stages: 3               <br>
        //! Order: 3                <br>
        //! Stability: L-stable     <br>
        //! Stiffly accurate: yes
        //!
        SDIRK3,

        //!
        //! String descriptor: \e "sdirk5"
        //!
        //! Stages: 5               <br>
        //! Order: 5                <br>
        //! Stability: L-stable     <br>
        //! Stiffly accurate: yes
        //!
        SDIRK5,

        //!
        //! String descriptor: \e "sdirk-5-3-4"
        //!
        //! Stages: 5               <br>
        //! Order: 4                <br>
        //! Stability: L-stable     <br>
        //! Stiffly accurate: yes
        //!
        SDIRK_5_3_4,

        //!
        //! String descriptor: \e "kvaerno-4-2-3"
        //!
        //! Stages: 4               <br>
        //! Order: 3                <br>
        //! Stability: L-stable     <br>
        //! Stiffly accurate: yes
        //!
        KVAERNO_4_2_3,

        //!
        //! String descriptor: \e "kvaerno-7-4-5"
        //!
        //! Stages: 7               <br>
        //! Order: 5                <br>
        //! Stability: L-stable     <br>
        //! Stiffly accurate: yes
        //!
        KVAERNO_7_4_5,

        //!
        //! String descriptor: \e "ark-8-4-5"
        //!
        //! Stages: 8               <br>
        //! Order: 5                <br>
        //! Stability: L-stable     <br>
        //! Stiffly accurate: yes
        //!
        ARK_8_4_5,

        //!
        //! String descriptor: \e "gsbp-dirk-3"
        //!
        //! Stages: 3               <br>
        //! Order: 3                <br>
        //! Stability: L-stable     <br>
        //! Stiffly accurate: no
        //!
        GSBP_DIRK_3,

        //!
        //! String descriptor: \e "gsbp-dirk-4"
        //!
        //! Stages: 4               <br>
        //! Order: 4                <br>
        //! Stability: L-stable     <br>
        //! Stiffly accurate: no
        //!
        GSBP_DIRK_4
    };


    //--------------------------------------------------------------------------------------------------------
    //! \brief  Specifies the type of hybrid splitting (if any) used by the DIRK integrator.
    //--------------------------------------------------------------------------------------------------------
    enum class HybridMethod {

        //!
        //! String descriptor: \e "none"
        //!
        //! No hybrid splitting used. The transport system is solved using a standard discrete ordinates
        //! approximation.
        //!
        //! \see    DIRKIntegrator::Step_Nonhybrid()
        //!
        None,

        //!
        //! String descriptor: \e "hybrid-ia"
        //!
        //! Applies a hybrid method derived by splitting the transport equation at the continuum level.
        //! Relabeling of the collided flux is (optionally) performed at the end of each timestep via an
        //! interpolatory reconstruction procedure.
        //!
        //! \see    DIRKIntegrator::Step_HybridI()
        //!
        HybridIa,

        //!
        //! String descriptor: \e "hybrid-ic"
        //!
        //! Applies a hybrid method derived by splitting the transport equation at the continuum level.
        //! Relabeling of the collided flux is (optionally) performed at the end of each timestep via a
        //! reconstruction procedure motivated by the Nyström interpolation technique.
        //!
        //! \see    DIRKIntegrator::Step_HybridI()
        //!
        HybridIc,

        //!
        //! String descriptor: \e "hybrid-iia"
        //!
        //! Applies a hybrid method derived by splitting the transport equation at the discrete level.
        //! Relabeling of the collided flux is always performed after the implicit solve performed at each
        //! stage of the DIRK method via an interpolatory reconstruction procedure.
        //!
        //! \see    DIRKIntegrator::Step_HybridII()
        //!
        HybridIIa,

        //!
        //! String descriptor: \e "hybrid-iic"
        //!
        //! Applies a hybrid method derived by splitting the transport equation at the discrete level.
        //! Relabeling of the collided flux is always performed after the implicit solve performed at each
        //! stage of the DIRK method via a reconstruction procedure motivated by the Nyström interpolation
        //! technique.
        //!
        //! \see    DIRKIntegrator::Step_HybridII()
        //!
        HybridIIc
    };


    //--------------------------------------------------------------------------------------------------------
    //! \brief  Used to convert DIRKIntegrator::DIRKMethod values to their corresponding string descriptors.
    //--------------------------------------------------------------------------------------------------------
    static std::map< DIRKIntegrator::DIRKMethod, std::string > DIRKMethod_to_String;


    //--------------------------------------------------------------------------------------------------------
    //! \brief  Used to convert string descriptors to DIRKIntegrator::DIRKMethod values.
    //--------------------------------------------------------------------------------------------------------
    static std::map< std::string, DIRKIntegrator::DIRKMethod > String_to_DIRKMethod;


    //--------------------------------------------------------------------------------------------------------
    //! \brief  Used to convert DIRKIntegrator::HybridMethod values to their corresponding string descriptors.
    //--------------------------------------------------------------------------------------------------------
    static std::map< DIRKIntegrator::HybridMethod, std::string > HybridMethod_to_String;


    //--------------------------------------------------------------------------------------------------------
    //! \brief  Used to convert string descriptors to DIRKIntegrator::HybridMethod values.
    //--------------------------------------------------------------------------------------------------------
    static std::map< std::string, DIRKIntegrator::HybridMethod > String_to_HybridMethod;


    //========================================================================================================
    //=== CONSTRUCTORS AND DESTRUCTOR ========================================================================
    //========================================================================================================


    //--------------------------------------------------------------------------------------------------------
    //! \brief  Delete default constructor.
    //--------------------------------------------------------------------------------------------------------
    DIRKIntegrator( void ) = delete;


    //--------------------------------------------------------------------------------------------------------
    //! \brief  Delete copy constructor.
    //--------------------------------------------------------------------------------------------------------
    DIRKIntegrator( const DIRKIntegrator & ) = delete;


    //--------------------------------------------------------------------------------------------------------
    //! \brief  Constructs a DIRKIntegrator object using the provided values.
    //--------------------------------------------------------------------------------------------------------
    DIRKIntegrator( const ParameterList & );


    //--------------------------------------------------------------------------------------------------------
    //! \brief  Destructor for class DIRKIntegrator.
    //--------------------------------------------------------------------------------------------------------
    ~DIRKIntegrator( void ) override;


    //========================================================================================================
    //=== OPERATOR OVERLOADS =================================================================================
    //========================================================================================================


    //--------------------------------------------------------------------------------------------------------
    //! \brief  Delete copy assignment operator.
    //--------------------------------------------------------------------------------------------------------
    DIRKIntegrator & operator=( const DIRKIntegrator & ) = delete;


    //========================================================================================================
    //=== INTERFACE ROUTINES =================================================================================
    //========================================================================================================


    //--------------------------------------------------------------------------------------------------------
    //! \brief  Interface routine for computing a single timestep with step size \pp{step_size}.
    //--------------------------------------------------------------------------------------------------------
    DIRKIntegrator & Step (

        RKDG::OrdinateFlux & initial_condition,
        const RKDG::OrdinateFlux & source,
        const RKDG::CrossSection & sigma_t,
        const RKDG::CrossSection & sigma_s,
        const double step_size

    ) override;


    //--------------------------------------------------------------------------------------------------------
    //! \brief  Interface routine for computing a single timestep with step size \pp{step_size}.
    //--------------------------------------------------------------------------------------------------------
    DIRKIntegrator & Step (

        RKDG::OrdinateFlux & u_initial_condition,
        RKDG::OrdinateFlux * c_initial_condition,
        const RKDG::OrdinateFlux & source,
        const RKDG::CrossSection & sigma_t,
        const RKDG::CrossSection & sigma_s,
        const double step_size

    ) override;

    //--------------------------------------------------------------------------------------------------------
    //! \brief  Returns the key used to determine the type of DIRK method.
    //--------------------------------------------------------------------------------------------------------
    static std::string GetInputKey( void );

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
    TimeIntegratorType GetIntegratorType( void ) const override {  return TimeIntegratorType::DIRK;  }


    //--------------------------------------------------------------------------------------------------------
    //! \brief  Prints the Butcher tableau of the DIRK method to the logging interface.
    //--------------------------------------------------------------------------------------------------------
    void PrintTableau( const std::string prefix = "  " ) const;


private:

    //========================================================================================================
    //=== PRIVATE MEMBER VARIABLES ===========================================================================
    //========================================================================================================


    DIRKMethod dirk_method;             //!< Specifies the DIRK method used.

    double * A;                         //!< \f$ A \f$ matrix of Butcher tableau.
    double * b;                         //!< \f$ b \f$ vector of Butcher tableau.
    double * c;                         //!< \f$ c \f$ vector of Butcher tableau.

    bool stiffly_accurate;              //!< Specifies if the DIRK method is stiffly accurate or not.
    int64_t num_stages;                 //!< Number of interior stages of the DIRK method.
    int64_t stages_allocd;              //!< Number of stage vectors allocated by the DIRK method.

    //!
    //! \brief  Defines a map between stage vectors to allow re-use of allocated objects for multiple stages
    //!         when possible.
    //!
    int64_t * stage_alloc_map;

    RKDG::ScalarFlux * phi;             //!< Temporary storage for scalar fluxes.

    // Nonhybrid fluxes and solvers.
    RKDG::OrdinateFlux ** stage_flux;   //!< Angular flux stage values.

    //! Solver for nonhybrid systems.
    std::shared_ptr<ImplicitSolver<RKDG::OrdinateFlux>> solver;

    // Variables for hybrid methods.
    HybridMethod hybrid_method;         //!< Specifies the type of hybrid splitting (if any) to use.

    RKDG::OrdinateFlux ** u_stage_flux; //!< Uncollided flux stage values.

    //! Solver for uncollided systems.
    std::shared_ptr<ImplicitSolver<RKDG::OrdinateFlux>> u_solver;

    RKDG::OrdinateFlux ** c_stage_flux; //!< Collided flux stage values.

    //! Solver for collided systems.
    std::shared_ptr<ImplicitSolver<RKDG::OrdinateFlux>> c_solver;

    RKDG::OrdinateFlux * c_temp;        //!< Collided flux for hybrid-II DIRK methods.

    bool relabel = true;                //!< Specifies whether relabeling should be performed or not. Defaults to true.
    RelabelOperator * relabel_operator; //!< RelabelOperator for mapping collided fluxes to uncollided.


    //========================================================================================================
    //=== PRIVATE HELPER ROUTINES ============================================================================
    //========================================================================================================


    //--------------------------------------------------------------------------------------------------------
    //! \brief  Indexing function for the Butcher matrix \f$ A \f$ stored at DIRKIntegrator::A.
    //!
    //! \param[in]      i       Row index of element.
    //! \param[in]      j       Column index of element.
    //!
    //! \return     Returns a reference to the (\pp{i},\pp{j})th element of the Butcher matrix \f$ A \f$ for
    //!             the DIRK method.
    //--------------------------------------------------------------------------------------------------------
    inline double & IA( const int64_t i, const int64_t j ) const {

    # if defined (STRICT_CHECK)

        if ( i < 0 || i >= this->num_stages ) {

            std::string error_message
                = "Value " + std::to_string(i) + " for index i in '" + std::string(__func__)
                  + "' outsize permissible range [ 0, " + std::to_string( this->num_stages ) + " ).\n";

            PRINT_ERROR( error_message.c_str() )
            throw std::out_of_range( error_message );
        }

        if ( j < 0 || j >= this->num_stages ) {

            std::string error_message
                = "Value " + std::to_string(j) + " for index j in '" + std::string(__func__)
                  + "' outsize permissible range [ 0, " + std::to_string( this->num_stages ) + " ).\n";

            PRINT_ERROR( error_message.c_str() )
            throw std::out_of_range( error_message );
        }

        if ( j > i ) {

            std::string warning_message
                = "Value " + std::to_string(j) + " for index j in '" + std::string(__func__)
                  + "' is larger than value " + std::to_string(i) + " for index i.\n";

            PRINT_WARNING( warning_message.c_str() )
        }

    # endif // if defined (STRICT_CHECK)

        return this->A[ j + this->num_stages*i ];
    }


    //--------------------------------------------------------------------------------------------------------
    //! \brief  Computes the final output value of the DIRK method using the intermediate stage values of the
    //!         non-hybrid or reconstructed hybrid approximation.
    //--------------------------------------------------------------------------------------------------------
    void ComputeOutput_Nonhybrid (

        const RKDG::OrdinateFlux & initial_condition,
        const RKDG::OrdinateFlux & source,
        const RKDG::CrossSection & sigma_t,
        const RKDG::CrossSection & sigma_s,
        const double h
    );


    //--------------------------------------------------------------------------------------------------------
    //! \brief  Computes the final output values of the DIRK method using the intermediate stage values of the
    //!         non-reconstructed hybrid approximation.
    //--------------------------------------------------------------------------------------------------------
    void ComputeOutput_Hybrid (

        const RKDG::OrdinateFlux & u_initial_condition,
        const RKDG::OrdinateFlux * c_initial_condition,
        const RKDG::OrdinateFlux & source,
        const RKDG::CrossSection & sigma_t,
        const RKDG::CrossSection & sigma_s,
        const double h
    );


    //--------------------------------------------------------------------------------------------------------
    //! \brief  Computes reconstructed approximations of the intermediate stage values of a set of hybrid
    //!         solutions using a Nyström interpolation technique.
    //--------------------------------------------------------------------------------------------------------
    void ComputeNystromInterpolants (

        const RKDG::OrdinateFlux & u_initial_condition,
        const RKDG::OrdinateFlux & source,
        const RKDG::CrossSection & sigma_t,
        const RKDG::CrossSection & sigma_s,
        const double h
    );


    //--------------------------------------------------------------------------------------------------------
    //! \brief  Constructs the Butcher tableau for the DIRK method set by DIRKIntegrator::dirk_method.
    //--------------------------------------------------------------------------------------------------------
    void SetTableau( void );


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
    //! \brief  Computes a one timestep update using hybrid I type splitting.
    //--------------------------------------------------------------------------------------------------------
    void Step_HybridI (

        RKDG::OrdinateFlux & u_initial_condition,
        RKDG::OrdinateFlux * c_initial_condition,
        const RKDG::OrdinateFlux & source,
        const RKDG::CrossSection & sigma_t,
        const RKDG::CrossSection & sigma_s,
        const double h
    );


    //--------------------------------------------------------------------------------------------------------
    //! \brief  Computes a one timestep update using hybrid II type splitting.
    //--------------------------------------------------------------------------------------------------------
    void Step_HybridII (

        RKDG::OrdinateFlux & initial_condition,
        const RKDG::OrdinateFlux & source,
        const RKDG::CrossSection & sigma_t,
        const RKDG::CrossSection & sigma_s,
        const double h
    );
};


# endif // ifndef __DIRK_INTEGRATOR_HPP__
