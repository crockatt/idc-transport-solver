//------------------------------------------------------------------------------------------------------------
//! \file   time_integrators/RKDG/erk.cpp
//! \brief  Implementation of explicit Runge-Kutta (ERK) time integration schemes.
//!
//! \author Michael M. Crockatt
//! \date   February 2018
//------------------------------------------------------------------------------------------------------------

# include <cinttypes>
# include <cstdint>
# include <cstdlib>
# include <cstring>

# include "objects/RKDG/CrossSection.hpp"
# include "objects/RKDG/OrdinateFlux.hpp"
# include "operators/TransportOperator.hpp"
# include "time_integrators/RKDG/erk.hpp"
# include "utils/CLog.hpp"
# include "utils/global.hpp"


using namespace RKDG;
using namespace Quadrule;


//============================================================================================================
//=== STATIC DEFINITIONS =====================================================================================
//============================================================================================================


//------------------------------------------------------------------------------------------------------------
//! \brief  Used to convert ERKIntegrator::ERKMethod values to their corresponding string descriptors.
//------------------------------------------------------------------------------------------------------------
std::map< ERKIntegrator::ERKMethod, std::string > ERKIntegrator::ERKMethod_to_String = {

    { ERKIntegrator::ERKMethod::None,       "none"      },
    { ERKIntegrator::ERKMethod::RK4,        "rk4"       }
};


//------------------------------------------------------------------------------------------------------------
//! \brief  Used to convert string descriptors to ERKIntegrator::ERKMethod values.
//------------------------------------------------------------------------------------------------------------
std::map< std::string, ERKIntegrator::ERKMethod > ERKIntegrator::String_to_ERKMethod = {

    { "none",       ERKIntegrator::ERKMethod::None      },
    { "rk4",        ERKIntegrator::ERKMethod::RK4       }
};


//------------------------------------------------------------------------------------------------------------
//! \brief  Used to convert ERKIntegrator::HybridMethod values to their corresponding string descriptors.
//------------------------------------------------------------------------------------------------------------
std::map< ERKIntegrator::HybridMethod, std::string > ERKIntegrator::HybridMethod_to_String = {

    { ERKIntegrator::HybridMethod::None,        "none"      }
};


//------------------------------------------------------------------------------------------------------------
//! \brief  Used to convert string descriptors to ERKIntegrator::HybridMethod values.
//------------------------------------------------------------------------------------------------------------
std::map< std::string, ERKIntegrator::HybridMethod > ERKIntegrator::String_to_HybridMethod = {

    { "none",       ERKIntegrator::HybridMethod::None       }
};


//============================================================================================================
//=== CONSTRUCTORS AND DESTRUCTOR ============================================================================
//============================================================================================================


//------------------------------------------------------------------------------------------------------------
//! \brief  Constructs an ERKIntegrator object using the provided values.
//!
//! \param[in]  input_list      List of parameters to use for initialization.
//------------------------------------------------------------------------------------------------------------
ERKIntegrator::ERKIntegrator (

    const ParameterList & input_list
) :
    TimeIntegrator( input_list ),

    erk_method(),
    A{ nullptr },
    b{ nullptr },
    c{ nullptr },
    num_stages(),
    phi(),
    stage_flux{ nullptr },
    hybrid_method()
{

    bool symmetric_reduce = false;

    ParameterList   params = input_list;
    ParameterList u_params = MakeUncollidedList( input_list );
    ParameterList c_params = MakeCollidedList( input_list );

    // Read polynomial degree for DG approximation.
    int64_t DG_degree = GetDGDegreeX( params );

    // Determine type of hybrid splitting (default to none).
    try         {  params.GetValue( "hybrid_method", this->hybrid_method, String_to_HybridMethod );  }
    catch (...) {  this->hybrid_method = HybridMethod::None;  }

    // Determine angular discretization(s).
# if SPACE_DIMS == 2

    try         {  params.GetValue( "ordinate_sym_reduce", symmetric_reduce );  }
    catch (...) {  symmetric_reduce = false;  }

# endif // if SPACE_DIMS == 2

    OrdinateSet ordinate_set, u_ordinate_set, c_ordinate_set;

    // Angular discretizations for hybrid methods.
    if ( IsHybrid() ) {

    // Angular discretizations for nonhybrid methods.
    } else {

        // Set default ordinate type.
        OrdinateSet::OrdinateType ordinate_type =
            # if SPACE_DIMS == 1
                OrdinateSet::OrdinateType::GaussLegendre;
            # elif SPACE_DIMS >= 2
                OrdinateSet::OrdinateType::ChebyshevLegendre;
            # endif

        // Try reading alternative specifications.
        try {  ordinate_type = params.GetValue( "ordinate_type",
                                                OrdinateSet::String_to_OrdinateType );
        } catch (...) {/* empty */}

        // Configure ordinate specification object.
        ordinate_set.Reconfigure( params.GetValue<int64_t>( "ang_order" ),
                                  symmetric_reduce,
                                  ordinate_type );
    }


    // --- Setup ERK method. -------------------------------------------------------------------------- //

    params.GetValue( "erk_method", this->erk_method, ERKIntegrator::String_to_ERKMethod );

    SetTableau();


    // --- Setup objects used by ERK integrator. ----------------------------------------------------- //

    this->phi.Reconfigure( *this, DG_degree );

    // Hybrid/non-hybrid setup.
    switch ( this->hybrid_method ) {

        case ERKIntegrator::HybridMethod::None:
        {
            this->stage_flux = new OrdinateFlux*[ this->num_stages ];

            for ( int64_t i = 0; i < this->num_stages; ++i )
                this->stage_flux[i] = new OrdinateFlux( *this, ordinate_set, DG_degree );
        } break;

        default:
        {   std::string error_message = "Invalid hybrid splitting '"
                                        + HybridMethod_to_String.at( this->hybrid_method )
                                        + "' for initializing ERK integrator in '" + std::string(__func__)
                                        + "'.\n";

            PRINT_ERROR( error_message.c_str() )
            throw std::invalid_argument( error_message );
        }
    }
}


//------------------------------------------------------------------------------------------------------------
//! \brief  Destructor for class ERKIntegrator.
//------------------------------------------------------------------------------------------------------------
ERKIntegrator::~ERKIntegrator ( void ) {

    delete [] this->A;
    delete [] this->b;
    delete [] this->c;

    delete [] this->stage_flux;
}


//============================================================================================================
//=== INTERFACE ROUTINES =====================================================================================
//============================================================================================================


//------------------------------------------------------------------------------------------------------------
//! \brief  Calls the appropriate \c Step_* function corresponding to the desired ERK hybrid splitting as set
//!         by ERKIntegrator::hybrid_method.
//!
//! \param[in,out]  initial_condition   Initially contains the initial condition \f$ \Psi^n \f$ for the
//!                                     current timestep. Upon return contains the result \f$ \Psi^{n+1} \f$
//!                                     of the timestep update.
//! \param[in]      sigma_t             RKDG::CrossSection object containing the total cross section
//!                                     \f$ \sigma_{\mathrm{t}} \f$.
//! \param[in]      sigma_s             RKDG::CrossSection object containing the scattering cross section
//!                                     \f$ \sigma_{\mathrm{s}} \f$.
//! \param[in]      source              Source term of the transport system. <br>
//!                                     The source term is assumed to be constant across the timestep interval.
//! \param[in]      step_size           Timestep size.
//!
//! \see    EKIntegrator::Step_Nonhybrid()
//------------------------------------------------------------------------------------------------------------
ERKIntegrator & ERKIntegrator::Step (

    OrdinateFlux & initial_condition,
    const OrdinateFlux & source,
    const CrossSection & sigma_t,
    const CrossSection & sigma_s,
    const double step_size
) {

    switch ( this->hybrid_method ) {

        case HybridMethod::None:

            Step_Nonhybrid( initial_condition, source, sigma_t, sigma_s, step_size );
            break;

        default:
        {   std::string error_message = "Invalid hybrid splitting '"
                                        + HybridMethod_to_String.at( this->hybrid_method )
                                        + "' for ERK step.\n";

            PRINT_ERROR( error_message.c_str() )
            throw std::invalid_argument( error_message );
        }
    }

    return *this;
}


//------------------------------------------------------------------------------------------------------------
//! \brief  Calls the appropriate \c Step_* function corresponding to the desired ERK hybrid splitting as set
//!         by ERKIntegrator::hybrid_method.
//!
//! \param[in,out]  u_initial_condition Initially contains the initial condition \f$ \Psi^n \f$ for the
//!                                     uncollided flux in the current timestep. Upon return contains the
//!                                     result \f$ \Psi^{n+1} \f$ of the timestep update if relabeling is used
//!                                     or only the uncollided component of the result if relabeling is not
//!                                     used.
// \param[in,out]  c_initial_condition Initially contains the initial condition \f$ \Psi^n \f$ for the
//                                     collided flux in the current timestep. Upon return contains either
//                                     a zero flux if relabeling is used or the collided component of the
//                                     result \f$ \Psi^{n+1} \f$ if relabeling is not used. <br/>
//                                     If this is \c null, then relabeling is forced.
//! \param[in]      sigma_t             RKDG::CrossSection object containing the total cross section
//!                                     \f$ \sigma_{\mathrm{t}} \f$.
//! \param[in]      sigma_s             RKDG::CrossSection object containing the scattering cross section
//!                                     \f$ \sigma_{\mathrm{s}} \f$.
//! \param[in]      source              Source term of the transport system. <br>
//!                                     The source term is assumed to be constant across the timestep interval.
//! \param[in]      step_size           Timestep size.
//!
//! \see    EKIntegrator::Step_Nonhybrid()
//------------------------------------------------------------------------------------------------------------
ERKIntegrator & ERKIntegrator::Step (

    OrdinateFlux & u_initial_condition,
    OrdinateFlux *,                     // c_initial_condition,
    const OrdinateFlux & source,
    const CrossSection & sigma_t,
    const CrossSection & sigma_s,
    const double step_size
) {

    switch ( this->hybrid_method ) {

        case HybridMethod::None:

        # if defined (STRICT_CHECK) && 0

            if ( c_initial_condition != nullptr ) {

                PRINT_WARNING( "ERK integrator with hybrid splitting '%s' ignores pointer to collided "
                               "initial condition.\n",
                               HybridMethod_to_String.at( this->hybrid_method ).c_str() )
            }

        # endif // if defined (STRICT_CHECK)

            Step_Nonhybrid( u_initial_condition, source, sigma_t, sigma_s, step_size );
            break;

        default:
        {   std::string error_message = "Invalid hybrid splitting '"
                                        + HybridMethod_to_String.at( this->hybrid_method )
                                        + "' for ERK step.\n";

            PRINT_ERROR( error_message.c_str() )
            throw std::invalid_argument( error_message );
        }
    }

    return *this;
}


//------------------------------------------------------------------------------------------------------------
//! \brief  Returns true if the integrator has been constructed with hybrid splitting and false otherwise.
//------------------------------------------------------------------------------------------------------------
bool ERKIntegrator::IsHybrid ( void ) const {

    switch ( this->hybrid_method ) {

        default:
            return false;
    }
}


//------------------------------------------------------------------------------------------------------------
//! \brief  Prints the Butcher tableau of the ERK method to the logging interface.
//------------------------------------------------------------------------------------------------------------
void ERKIntegrator::PrintTableau ( void ) const {

    // --- Check to make sure integrator is properly initialized. ------------------------------------- //

    if ( this->erk_method == ERKIntegrator::ERKMethod::None ) {

        std::string error_message = "No Butcher tableau to print for ERKMethod '"
                                    + ERKMethod_to_String.at( this->erk_method )
                                    + "'.\n";

        PRINT_ERROR( error_message.c_str() )
        throw std::invalid_argument( error_message );
    }

    // --- Print tableau. ----------------------------------------------------------------------------- //

    char buffer[100];

    buffer[0] = '\0';
    strncat( buffer, "-------------------------------------------------", Global::col_width );

    PRINT_LOG( "\n" )

    for ( int64_t i = 0; i < this->num_stages; ++i ) {

        // Print c.
        PRINT_LOG( "  % -*.4e |", Global::col_width, this->c[i] )

        // Print row of A.
        for ( int64_t j = 0; j < i; ++j )
            PRINT_LOG( " % -*.4e ", Global::col_width, IA(i,j) )

        PRINT_LOG( "\n" )
    }

    // Print horizontal line.
    PRINT_LOG( " -%*s-|", Global::col_width, buffer )

    for ( int64_t i = 0; i < this->num_stages; ++i ) {  PRINT_LOG( "-%*s-", Global::col_width, buffer )  }

    PRINT_LOG( "\n" )

    // Print last row for b.
    buffer[0] = '\0';

    PRINT_LOG( "  %*s |", Global::col_width, buffer )

    for ( int64_t j = 0; j < this->num_stages; ++j )
        PRINT_LOG( " % -*.4e ", Global::col_width, this->b[j] )

    PRINT_LOG( "\n" )
    PRINT_LOG( "\n" )
}


//============================================================================================================
//=== PRIVATE TIME INTEGRATION ROUTINES ======================================================================
//============================================================================================================


//------------------------------------------------------------------------------------------------------------
//! \brief  Computes a one timestep update without hybrid splitting.
//!
//! \param[in,out]  initial_condition   Initially contains the initial condition \f$ \Psi^n \f$ for the
//!                                     current timestep. Upon return contains the result \f$ \Psi^{n+1} \f$
//!                                     of the timestep update.
//! \param[in]      sigma_t             RKDG::CrossSection object containing the total cross section
//!                                     \f$ \sigma_{\mathrm{t}} \f$.
//! \param[in]      sigma_s             RKDG::CrossSection object containing the scattering cross section
//!                                     \f$ \sigma_{\mathrm{s}} \f$.
//! \param[in]      source              Source term of the transport system. <br>
//!                                     The source term is assumed to be constant across the timestep interval.
//! \param[in]      h                   Timestep size.
//!
//! \see    ERKIntegrator::Step()
//------------------------------------------------------------------------------------------------------------
void ERKIntegrator::Step_Nonhybrid (

    OrdinateFlux & initial_condition,
    const OrdinateFlux & source,
    const CrossSection & sigma_t,
    const CrossSection & sigma_s,
    const double h
) {

    double stage_coeff, source_coeff;

    // --- Construct stage values in "Y" formulation. ------------------------------------------------- //

    // Compute each internal stage.
    for ( int64_t i = 0; i < this->num_stages; ++i ) {

        stage_coeff = h;
        OrdinateFlux::Copy( *this->stage_flux[i], initial_condition );
        source_coeff = 0.0;

        this->phi.ZeroDensity();

        // Compute weighted sum of transport operator applied to previous stages.
        for ( int64_t j = 0; j < i; ++j ) {

            if ( IA(i,j) == 0.0 ) {  continue;  }

            source_coeff += IA(i,j);

            TransportOperator::Lmv( -stage_coeff * IA(i,j), 1.0, sigma_t, *this->stage_flux[j],
                                    *this->stage_flux[i] );
            TransportOperator::Pmv( stage_coeff * IA(i,j), 1.0, *this->stage_flux[j], this->phi );
        }

        TransportOperator::Smv( 1.0, 1.0, sigma_s, this->phi, *this->stage_flux[i] );

        // Include source from previous stages.
        OrdinateFlux::AXPY( stage_coeff * source_coeff, 1.0, source, *this->stage_flux[i] );

        this->stage_flux[i]->Copy( source, OpDomain::Boundary );
    }

    // Compute result from internal stages.
    this->phi.ZeroDensity();

    for ( int64_t i = 0; i < this->num_stages; ++i ) {

        if ( this->b[i] == 0.0 ) {  continue;  }

        TransportOperator::Lmv( -h * this->b[i], 1.0, sigma_t, *this->stage_flux[i], initial_condition );
        TransportOperator::Pmv( h * this->b[i], 1.0, *this->stage_flux[i], this->phi );
        OrdinateFlux::AXPY( h * this->b[i], 1.0, source, initial_condition );
    }

    TransportOperator::Smv( 1.0, 1.0, sigma_s, this->phi, initial_condition );

    initial_condition.ZeroDensity( OpDomain::Boundary );
}


//============================================================================================================
//=== PRIVATE HELPER ROUTINES ================================================================================
//============================================================================================================


//------------------------------------------------------------------------------------------------------------
//! \brief  Constructs the Butcher tableau for the ERK method set by ERKIntegrator::erk_method.
//!
//! Sets ERKIntegrator::num_stages. Allocates memory for ERKIntegrator::A, ERKIntegrator::b, and
//! ERKIntegrator::c, and initializes allocated memory with the proper values from the Butcher tableau of the
//! method.
//------------------------------------------------------------------------------------------------------------
void ERKIntegrator::SetTableau ( void ) {

    switch( this->erk_method ) {

        case ERKIntegrator::ERKMethod::RK4:
        {
            //
            // RK4
            //
            // Stages:      four
            // Order:       fourth
            //

            this->num_stages = 4;

            this->A = new double[ this->num_stages * this->num_stages ];
            this->b = new double[ this->num_stages ];
            this->c = new double[ this->num_stages ];

            IA(0,0) = 0.0;
            IA(0,1) = 0.0;
            IA(0,2) = 0.0;
            IA(0,3) = 0.0;

            IA(1,0) = 0.5;
            IA(1,1) = 0.0;
            IA(1,2) = 0.0;
            IA(1,3) = 0.0;

            IA(2,0) = 0.0;
            IA(2,1) = 0.5;
            IA(2,2) = 0.0;
            IA(2,3) = 0.0;

            IA(3,0) = 0.0;
            IA(3,1) = 0.0;
            IA(3,2) = 1.0;
            IA(3,3) = 0.0;

            this->b[0] = 1.0 / 6.0;
            this->b[1] = 1.0 / 3.0;
            this->b[2] = 1.0 / 3.0;
            this->b[3] = 1.0 / 6.0;

            this->c[0] = 0.0;
            this->c[1] = 0.5;
            this->c[2] = 0.5;
            this->c[3] = 1.0;

            break;
        }

        default:
        {   std::string error_message = "Invalid ERKMethod '"
                                        + ERKMethod_to_String.at( this->erk_method )
                                        + "' for initializing ERK integrator this in '"
                                        + std::string(__func__) + "'.\n";

            PRINT_ERROR( error_message.c_str() )
            throw std::invalid_argument( error_message );
        }
    }
}
