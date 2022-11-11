//------------------------------------------------------------------------------------------------------------
//! \file   time_integrators/RKDG/erk.hpp
//! \brief  Header file for explicit Runge-Kutta (ERK) time integrator implementation.
//!
//! \author Michael M. Crockatt
//! \date   February 2018
//------------------------------------------------------------------------------------------------------------

# ifndef __ERK_INTEGRATOR_HPP__
# define __ERK_INTEGRATOR_HPP__


# include <cstdint>
# include <map>
# include <string>

# include "objects/RKDG/CrossSection.hpp"
# include "objects/RKDG/OrdinateFlux.hpp"
# include "time_integrators/Abstract/TimeIntegrator.hpp"
# include "utils/global.hpp"


//------------------------------------------------------------------------------------------------------------
//! \brief  Class for ERK time integration schemes.
//!
//!
//------------------------------------------------------------------------------------------------------------
class ERKIntegrator : public Abstract::TimeIntegrator {

public:

    //--------------------------------------------------------------------------------------------------------
    //! \brief  Specifies the EKR method.
    //!
    //! \todo   Descriptions of ERK methods.
    //--------------------------------------------------------------------------------------------------------
    enum class ERKMethod {

        None,
        RK4
    };


    //--------------------------------------------------------------------------------------------------------
    //! \brief  Specifies the type of hybrid splitting (if any) used by the ERK integrator.
    //!
    //! \todo   Include hybrid descriptions.
    //--------------------------------------------------------------------------------------------------------
    enum class HybridMethod {

        None
    };


    //--------------------------------------------------------------------------------------------------------
    //! \brief  Used to convert ERKIntegrator::ERKMethod values to their corresponding string descriptors.
    //--------------------------------------------------------------------------------------------------------
    static std::map< ERKIntegrator::ERKMethod, std::string > ERKMethod_to_String;


    //--------------------------------------------------------------------------------------------------------
    //! \brief  Used to convert string descriptors to ERKIntegrator::ERKMethod values.
    //--------------------------------------------------------------------------------------------------------
    static std::map< std::string, ERKIntegrator::ERKMethod > String_to_ERKMethod;


    //--------------------------------------------------------------------------------------------------------
    //! \brief  Used to convert ERKIntegrator::HybridMethod values to their corresponding string descriptors.
    //--------------------------------------------------------------------------------------------------------
    static std::map< ERKIntegrator::HybridMethod, std::string > HybridMethod_to_String;


    //--------------------------------------------------------------------------------------------------------
    //! \brief  Used to convert string descriptors to ERKIntegrator::HybridMethod values.
    //--------------------------------------------------------------------------------------------------------
    static std::map< std::string, ERKIntegrator::HybridMethod > String_to_HybridMethod;


    //========================================================================================================
    //=== CONSTRUCTORS AND DESTRUCTOR ========================================================================
    //========================================================================================================


    //--------------------------------------------------------------------------------------------------------
    //! \brief  Delete default constructor.
    //--------------------------------------------------------------------------------------------------------
    ERKIntegrator( void ) = delete;


    //--------------------------------------------------------------------------------------------------------
    //! \brief  Delete copy constructor.
    //--------------------------------------------------------------------------------------------------------
    ERKIntegrator( const ERKIntegrator & ) = delete;


    //--------------------------------------------------------------------------------------------------------
    //! \brief  Constructs an ERKIntegrator object using the provided values.
    //--------------------------------------------------------------------------------------------------------
    ERKIntegrator( const ParameterList & );


    //--------------------------------------------------------------------------------------------------------
    //! \brief  Destructor for class ERKIntegrator.
    //--------------------------------------------------------------------------------------------------------
    ~ERKIntegrator( void ) override;


    //========================================================================================================
    //=== OPERATOR OVERLOADS =================================================================================
    //========================================================================================================


    //--------------------------------------------------------------------------------------------------------
    //! \brief  Delete copy assignment operator.
    //--------------------------------------------------------------------------------------------------------
    ERKIntegrator & operator=( const ERKIntegrator & ) = delete;


    //========================================================================================================
    //=== INTERFACE ROUTINES =================================================================================
    //========================================================================================================


    //--------------------------------------------------------------------------------------------------------
    //! \brief  Interface routine for computing a single timestep with step size \pp{step_size}.
    //--------------------------------------------------------------------------------------------------------
    ERKIntegrator & Step (

        RKDG::OrdinateFlux & initial_condition,
        const RKDG::OrdinateFlux & source,
        const RKDG::CrossSection & sigma_t,
        const RKDG::CrossSection & sigma_s,
        const double step_size

    ) override;


    //--------------------------------------------------------------------------------------------------------
    //! \brief  Interface routine for computing a single timestep with step size \pp{step_size}.
    //--------------------------------------------------------------------------------------------------------
    ERKIntegrator & Step (

        RKDG::OrdinateFlux & u_initial_condition,
        RKDG::OrdinateFlux *,
        const RKDG::OrdinateFlux & source,
        const RKDG::CrossSection & sigma_t,
        const RKDG::CrossSection & sigma_s,
        const double step_size

    ) override;


    //--------------------------------------------------------------------------------------------------------
    //! \brief  Prints a summary of the integrator configuration to the logging interface.
    //--------------------------------------------------------------------------------------------------------
    void Print( const std::string & = "  " ) const override {};


    //--------------------------------------------------------------------------------------------------------
    //! \brief  Returns true if the integrator has been constructed with any form of hybrid splitting and
    //!         false otherwise.
    //--------------------------------------------------------------------------------------------------------
    bool IsHybrid( void ) const override;


    //--------------------------------------------------------------------------------------------------------
    //! \brief  Returns the TimeIntegratorType of the object.
    //--------------------------------------------------------------------------------------------------------
    TimeIntegratorType GetIntegratorType( void ) const override {  return TimeIntegratorType::ERK;  }


    //--------------------------------------------------------------------------------------------------------
    //! \brief  Prints the Butcher tableau of the ERK method to the logging interface.
    //--------------------------------------------------------------------------------------------------------
    void PrintTableau( void ) const;


private:

    //========================================================================================================
    //=== PRIVATE MEMBER VARIABLES ===========================================================================
    //========================================================================================================

    ERKMethod erk_method;               //!< Specifies the ERK method used.

    double * A;                         //!< \f$ A \f$ matrix of Butcher tableau.
    double * b;                         //!< \f$ b \f$ vector of Butcher tableau.
    double * c;                         //!< \f$ c \f$ vector of Butcher tableau.

    int num_stages;                     //!< Number of interior stages of the DIRK method.

    RKDG::ScalarFlux phi;               //!< Temporary storage for scalar fluxes.

    // Nonhybrid fluxes.
    RKDG::OrdinateFlux ** stage_flux;   //!< Angular flux stage values.

    // Variables for hybrid methods.
    HybridMethod hybrid_method;         //!< Specifies the type of hybrid splitting (if any) to use.


    //========================================================================================================
    //=== PRIVATE HELPER ROUTINES ============================================================================
    //========================================================================================================

    //--------------------------------------------------------------------------------------------------------
    //! \brief  Indexing function for the Butcher matrix \f$ A \f$ stored at ERKIntegrator::A.
    //!
    //! \param[in]      i       Row index of element.
    //! \param[in]      j       Column index of element.
    //!
    //! \return     Returns a reference to the (\pp{i},\pp{j})th element of the Butcher matrix \f$ A \f$ for
    //!             the ERK method.
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

        if ( j >= i ) {

            std::string warning_message
                = "Value " + std::to_string(j) + " for index j in '" + std::string(__func__)
                  + "' is at least as large as the value " + std::to_string(i) + " for index i.\n";

            PRINT_WARNING( warning_message.c_str() )
        }

    # endif // if defined (STRICT_CHECK)

        return this->A[ j + this->num_stages*i ];
    }


    //--------------------------------------------------------------------------------------------------------
    //! \brief  Constructs the Butcher tableau for the ERK method set by ERKIntegrator::erk_method.
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
};


# endif // ifndef __ERK_INTEGRATOR_HPP__
