//------------------------------------------------------------------------------------------------------------
//! \file   time_integrators/STDG/stdg.cpp
//! \brief  Implementation file for space-time DG (STDG) time-integration schemes.
//!
//! \author Michael M. Crockatt
//! \date   January 2018
//------------------------------------------------------------------------------------------------------------


# include <cinttypes>

# include "operators/TransportOperator.hpp"
# include "time_integrators/STDG/stdg.hpp"
# include "utils/CLog.hpp"


using namespace Quadrule;


//============================================================================================================
//=== STATIC DEFINITIONS =====================================================================================
//============================================================================================================


//------------------------------------------------------------------------------------------------------------
//! \brief  Used to convert STDGIntegrator::HybridMethod values to their corresponding string descriptors.
//------------------------------------------------------------------------------------------------------------
std::map< STDGIntegrator::HybridMethod, std::string > STDGIntegrator::HybridMethod_to_String = {

    { STDGIntegrator::HybridMethod::None,       "none"      },
    { STDGIntegrator::HybridMethod::Hybrid_a,   "hybrid-a"  },
    { STDGIntegrator::HybridMethod::Hybrid_c,   "hybrid-c"  }
};


//------------------------------------------------------------------------------------------------------------
//! \brief  Used to convert string descriptors to STDGIntegrator::HybridMethod values.
//------------------------------------------------------------------------------------------------------------
std::map< std::string, STDGIntegrator::HybridMethod > STDGIntegrator::String_to_HybridMethod = {

    { "none",       STDGIntegrator::HybridMethod::None      },
    { "hybrid-a",   STDGIntegrator::HybridMethod::Hybrid_a  },
    { "hybrid-ia",  STDGIntegrator::HybridMethod::Hybrid_a  },
    { "hybrid-iia", STDGIntegrator::HybridMethod::Hybrid_a  },
    { "hybrid-c",   STDGIntegrator::HybridMethod::Hybrid_c  },
    { "hybrid-ic",  STDGIntegrator::HybridMethod::Hybrid_c  },
    { "hybrid-iic", STDGIntegrator::HybridMethod::Hybrid_c  }
};


//============================================================================================================
//=== CONSTRUCTORS AND DESTRUCTOR ============================================================================
//============================================================================================================


//------------------------------------------------------------------------------------------------------------
//! \brief  Constructs an STDGIntegrator object using the provided values.
//!
//! \param[in]  input_list      List of parameters to use for initialization.
//------------------------------------------------------------------------------------------------------------
STDGIntegrator::STDGIntegrator (

    const ParameterList & input_list
) :
    TimeIntegrator( input_list ),
    phi(),
    psi{ nullptr },
    solver{ nullptr },
    hybrid_method(),
    uncollided{ nullptr },
    collided{ nullptr },
    u_solver{ nullptr },
    c_solver{ nullptr },
    relabel_operator{ nullptr }
{

    bool symmetric_reduce = false;

    ParameterList   params = input_list;
    ParameterList u_params = MakeUncollidedList( input_list );
    ParameterList c_params = MakeCollidedList( input_list );

    // Read polynomial degrees for DG approximations.
    int64_t DG_degree_x = GetDGDegreeX( params ),
            DG_degree_t = GetDGDegreeT( params );

    // Determine type of hybrid splitting (default to none).
    try         {  params.GetValue( "hybrid_method", this->hybrid_method,
                                                String_to_HybridMethod );              }
    catch (...) {  this->hybrid_method = HybridMethod::None;                           }

    // Determine angular discretization(s).
# if SPACE_DIMS == 2

    try         {  params.GetValue( "ordinate_sym_reduce", symmetric_reduce );  }
    catch (...) {  symmetric_reduce = false;                                                }

# endif // if SPACE_DIMS == 2

    OrdinateSet ordinate_set, u_ordinate_set, c_ordinate_set;

    // Angular discretizations for hybrid methods.
    if ( IsHybrid() ) {

        // Set default ordinate types.
        OrdinateSet::OrdinateType u_ordinate_type =
            # if SPACE_DIMS == 1
                OrdinateSet::OrdinateType::GaussLegendre;
            # elif SPACE_DIMS >= 2
                OrdinateSet::OrdinateType::ChebyshevLegendre;
            # endif

        OrdinateSet::OrdinateType c_ordinate_type =
            # if SPACE_DIMS == 1
                OrdinateSet::OrdinateType::GaussLegendre;
            # elif SPACE_DIMS >= 2
                OrdinateSet::OrdinateType::ChebyshevLegendre;
            # endif

        // Try reading alternative specifications.
        try {  u_ordinate_type = u_params.GetValue( "ordinate_type",
                                                    OrdinateSet::String_to_OrdinateType );
        } catch (...) {/* empty */}

        try {  c_ordinate_type = c_params.GetValue( "ordinate_type",
                                                    OrdinateSet::String_to_OrdinateType );
        } catch (...) {/* empty */}

        // Configure ordinate specification objects.
        u_ordinate_set.Reconfigure( u_params.GetValue<int64_t>( "ang_order" ),
                                    symmetric_reduce,
                                    u_ordinate_type );

        c_ordinate_set.Reconfigure( c_params.GetValue<int64_t>( "ang_order" ),
                                    symmetric_reduce,
                                    c_ordinate_type );

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


    // --- Setup STDG method. ------------------------------------------------------------------------- //

    switch ( this->hybrid_method ) {

        case HybridMethod::Hybrid_a:
        case HybridMethod::Hybrid_c:
        {
            try         {  params.GetValue( "relabel", this->relabel );  }
            catch (...) {  this->relabel = true;                                     }

        } break;

        case HybridMethod::None:
        default:
        break;
    }


    // --- Setup objects used by STDG integrator. ----------------------------------------------------- //

    this->phi.Reconfigure( *this, DG_degree_x, DG_degree_t );

    switch ( this->hybrid_method ) {

        case HybridMethod::None:
        {
            this->psi = new STDG::OrdinateFlux( *this, ordinate_set, DG_degree_x, DG_degree_t );

            this->solver = std::shared_ptr<ImplicitSolver<STDG::OrdinateFlux>>(
                    ImplicitSolver<STDG::OrdinateFlux>::ImplicitSolverFactory::CreateSolver(
                        *this, ordinate_set, params
                ) );
        } break;

        case HybridMethod::Hybrid_a:
        case HybridMethod::Hybrid_c:
        {
            this->uncollided = new STDG::OrdinateFlux( *this, u_ordinate_set, DG_degree_x, DG_degree_t );
            this->collided   = new STDG::OrdinateFlux( *this, c_ordinate_set, DG_degree_x, DG_degree_t );

            // Uncollided solver.
            this->u_solver = std::shared_ptr<ImplicitSolver<STDG::OrdinateFlux>>(
                    ImplicitSolver<STDG::OrdinateFlux>::ImplicitSolverFactory::CreateSolver(
                        *this, u_ordinate_set, u_params
                ) );

            // Collided solver.
            this->c_solver = std::shared_ptr<ImplicitSolver<STDG::OrdinateFlux>>(
                    ImplicitSolver<STDG::OrdinateFlux>::ImplicitSolverFactory::CreateSolver(
                        *this, c_ordinate_set, c_params
                ) );

            // Relabel operator.
            if ( this->hybrid_method == HybridMethod::Hybrid_a )
                this->relabel_operator = new RelabelOperator( c_ordinate_set, u_ordinate_set );
        } break;

        default:
        {   std::string error_message = "Invalid hybrid splitting '"
                                        + HybridMethod_to_String.at( this->hybrid_method )
                                        + "' for initializing STDG integrator in '" + std::string(__func__)
                                        + "'.\n";

            PRINT_ERROR( error_message.c_str() )
            throw std::invalid_argument( error_message );
        }
    }
}


//------------------------------------------------------------------------------------------------------------
//! \brief  Destructor for class STDGIntegrator.
//------------------------------------------------------------------------------------------------------------
STDGIntegrator::~STDGIntegrator( void ) {

    delete psi;

    delete uncollided;
    delete collided;

    delete relabel_operator;
}


//============================================================================================================
//=== INTERFACE ROUTINES =====================================================================================
//============================================================================================================


//------------------------------------------------------------------------------------------------------------
//! \brief  Calls the appropriate \c Step_* function corresponding to the desired STDG hybrid splitting as set
//!         by STDGIntegrator::hybrid_method.
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
//! \see    STDGIntegrator::Step_Nonhybrid()
//! \see    STDGIntegrator::Step_Hybrid()
//------------------------------------------------------------------------------------------------------------
STDGIntegrator & STDGIntegrator::Step (

    RKDG::OrdinateFlux & initial_condition,
    const RKDG::OrdinateFlux & source,
    const RKDG::CrossSection & sigma_t,
    const RKDG::CrossSection & sigma_s,
    const double step_size
) {

# if LOGLEVEL >= 1
    PRINT_LOG( "\n" )
# endif

    switch ( this->hybrid_method ) {

        case HybridMethod::None:

            Step_Nonhybrid( initial_condition, source, sigma_t, sigma_s, step_size );
            break;

        case HybridMethod::Hybrid_a:
        case HybridMethod::Hybrid_c:

            Step_Hybrid( initial_condition, nullptr, source, sigma_t, sigma_s, step_size );
            break;

        default:
        {   std::string error_message = "Invalid hybrid splitting '"
                                        + HybridMethod_to_String.at( this->hybrid_method )
                                        + "' for STDG step.\n";

            PRINT_ERROR( error_message.c_str() )
            throw std::invalid_argument( error_message );
        }
    }

# if LOGLEVEL >= 1
    PRINT_LOG( "\n" )
# endif

    return *this;
}


//------------------------------------------------------------------------------------------------------------
//! \brief  Calls the appropriate \c Step_* function corresponding to the desired STDG hybrid splitting as set
//!         by STDGIntegrator::hybrid_method.
//!
//! \param[in,out]  u_initial_condition Initially contains the initial condition \f$ \Psi^n \f$ for the
//!                                     uncollided flux in the current timestep. Upon return contains the
//!                                     result \f$ \Psi^{n+1} \f$ of the timestep update if relabeling is used
//!                                     or only the uncollided component of the result if relabeling is not
//!                                     used.
//! \param[in,out]  c_initial_condition Initially contains the initial condition \f$ \Psi^n \f$ for the
//!                                     collided flux in the current timestep. Upon return contains either
//!                                     a zero flux if relabeling is used or the collided component of the
//!                                     result \f$ \Psi^{n+1} \f$ if relabeling is not used. <br/>
//!                                     If this is \c null, then relabeling is forced.
//! \param[in]      sigma_t             RKDG::CrossSection object containing the total cross section
//!                                     \f$ \sigma_{\mathrm{t}} \f$.
//! \param[in]      sigma_s             RKDG::CrossSection object containing the scattering cross section
//!                                     \f$ \sigma_{\mathrm{s}} \f$.
//! \param[in]      source              Source term of the transport system. <br>
//!                                     The source term is assumed to be constant across the timestep interval.
//! \param[in]      step_size           Timestep size.
//!
//! \see    STDGIntegrator::Step_Nonhybrid()
//! \see    STDGIntegrator::Step_Hybrid()
//------------------------------------------------------------------------------------------------------------
STDGIntegrator & STDGIntegrator::Step (

    RKDG::OrdinateFlux & u_initial_condition,
    RKDG::OrdinateFlux * c_initial_condition,
    const RKDG::OrdinateFlux & source,
    const RKDG::CrossSection & sigma_t,
    const RKDG::CrossSection & sigma_s,
    const double step_size
) {

# if LOGLEVEL >= 1
    PRINT_LOG( "\n" )
# endif

    switch ( this->hybrid_method ) {

        case HybridMethod::None:

        # if defined (STRICT_CHECK)

            if ( c_initial_condition != nullptr ) {

                PRINT_WARNING( "STDG integrator with hybrid splitting '%s' ignores pointer to collided "
                               "initial condition.\n",
                               HybridMethod_to_String.at( this->hybrid_method ).c_str() )
            }

        # endif // if defined (STRICT_CHECK)

            Step_Nonhybrid( u_initial_condition, source, sigma_t, sigma_s, step_size );
            break;

        case HybridMethod::Hybrid_a:
        case HybridMethod::Hybrid_c:

            Step_Hybrid( u_initial_condition, c_initial_condition, source, sigma_t, sigma_s, step_size );
            break;

        default:
        {   std::string error_message = "Invalid hybrid splitting '"
                                        + HybridMethod_to_String.at( this->hybrid_method )
                                        + "' for STDG step.\n";

            PRINT_ERROR( error_message.c_str() )
            throw std::invalid_argument( error_message );
        }
    }

# if LOGLEVEL >= 1
    PRINT_LOG( "\n" )
# endif

    return *this;
}


//------------------------------------------------------------------------------------------------------------
//! \brief  Prints a summary of the integrator configuration to the logging interface.
//!
//! \param[in]  prefix  String prepended to each line of output, intended to be used for controlling
//!                     indentation of output.
//------------------------------------------------------------------------------------------------------------
void STDGIntegrator::Print (

    const std::string & prefix // = "  "

) const {

    int64_t DG_degree_t = -1;

    if ( IsHybrid() )
        DG_degree_t = this->uncollided->DG_degree_t;
    else
        DG_degree_t = this->psi->DG_degree_t;

    PRINT_LOG( "%s%-*s  %s\n", prefix.c_str(),
               Global::col_width, "Time integrator:",
               TimeIntegratorType_to_String.at( GetIntegratorType() ).c_str() )

    PRINT_LOG( "%s%-*s % " PRId64 "\n", prefix.c_str(), Global::col_width, "Max DG Degree (t):", DG_degree_t )

    this->Abstract::TimeIntegrator::Print( prefix );

    PRINT_LOG( "%s%-*s  %s\n", prefix.c_str(),
               Global::col_width, "Hybrid:",
               STDGIntegrator::HybridMethod_to_String.at( this->hybrid_method ).c_str() )

    if ( IsHybrid() ) {

        PRINT_LOG( "\n" )
        PRINT_LOG( "%s%-*s  %s\n", prefix.c_str(),
                   Global::col_width, "Relabel:", ( this->relabel ? "yes" : "no" ) )

        if (    this->relabel
             && this->hybrid_method == HybridMethod::Hybrid_a
        ) {
            this->relabel_operator->Print();
        }

        PRINT_LOG( "\n" )
        PRINT_LOG( "%sUncollided:\n", prefix.c_str() )

        this->u_solver->Print( prefix + "    " );

        PRINT_LOG( "\n" )
        PRINT_LOG( "%sCollided:\n", prefix.c_str() )

        this->c_solver->Print( prefix + "    " );

    } else {

        PRINT_LOG( "\n" )

        this->solver->Print( prefix );
    }
}


//------------------------------------------------------------------------------------------------------------
//! \brief  Returns true if the integrator has been constructed with any form of hybrid splitting and false
//!         otherwise.
//------------------------------------------------------------------------------------------------------------
bool STDGIntegrator::IsHybrid ( void ) const {

    switch ( this->hybrid_method ) {

        case HybridMethod::Hybrid_a:
        case HybridMethod::Hybrid_c:

            return true;

        default:
            return false;
    }
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
//! \see    STDGIntegrator::Step()
//------------------------------------------------------------------------------------------------------------
void STDGIntegrator::Step_Nonhybrid (

    RKDG::OrdinateFlux & initial_condition,
    const RKDG::OrdinateFlux & source,
    const RKDG::CrossSection & sigma_t,
    const RKDG::CrossSection & sigma_s,
    const double h
) {

    // Set initial guess for implicit solve from initial condition.
    STDG::OrdinateFlux::AXPY( 1.0, 0.0, initial_condition, *this->psi );
    this->solver->SetInitialGuess( this->psi );

    // Initialize source for implicit solve.
    STDG::OrdinateFlux::AXPY( 1.0, 0.0, source, *this->psi );

    // Perform STDG implicit solve.
    this->solver->Solve( *this->psi, sigma_t, sigma_s, h, &initial_condition );

    // Construct output value from DG temporal expansion.
    STDG::OrdinateFlux::AXPY( 1.0, 0.0, *this->psi, initial_condition );

    initial_condition.ZeroDensity( OpDomain::Boundary );
}


//------------------------------------------------------------------------------------------------------------
//! \brief  Computes a one timestep update with hybrid splitting.
//!
//! \param[in,out]  u_initial_condition Initially contains the initial condition \f$ \Psi^n \f$ for the
//!                                     uncollided flux in the current timestep. Upon return contains the
//!                                     result \f$ \Psi^{n+1} \f$ of the timestep update if relabeling is used
//!                                     or only the uncollided component of the result if relabeling is not
//!                                     used.
//! \param[in,out]  c_initial_condition Initially contains the initial condition \f$ \Psi^n \f$ for the
//!                                     collided flux in the current timestep. Upon return contains either
//!                                     a zero flux if relabeling is used or the collided component of the
//!                                     result \f$ \Psi^{n+1} \f$ if relabeling is not used. <br/>
//!                                     If this is \c null, then relabeling is forced.
//! \param[in]      sigma_t             RKDG::CrossSection object containing the total cross section
//!                                     \f$ \sigma_{\mathrm{t}} \f$.
//! \param[in]      sigma_s             RKDG::CrossSection object containing the scattering cross section
//!                                     \f$ \sigma_{\mathrm{s}} \f$.
//! \param[in]      source              Source term of the transport system. <br>
//!                                     The source term is assumed to be constant across the timestep interval.
//! \param[in]      h                   Timestep size.
//!
//! \see    STDGIntegrator::Step()
//------------------------------------------------------------------------------------------------------------
void STDGIntegrator::Step_Hybrid (

    RKDG::OrdinateFlux & u_initial_condition,
    RKDG::OrdinateFlux * c_initial_condition,
    const RKDG::OrdinateFlux & source,
    const RKDG::CrossSection & sigma_t,
    const RKDG::CrossSection & sigma_s,
    const double h
) {

    // --- Uncollided --------------------------------------------------------------------------------- //

    Global::TMR_uncollided.Start();

    // Initialize source for implicit solve.
    STDG::OrdinateFlux::AXPY( 1.0, 0.0, source, *this->uncollided );

    // Perform uncollided solve.
    this->u_solver->Solve( *this->uncollided, sigma_t, sigma_s, h, &u_initial_condition );

    Global::TMR_uncollided.Stop();

    // --- Collided ----------------------------------------------------------------------------------- //

    Global::TMR_collided.Start();

    // Set initial guess for collided solve from initial condition (if present).
    if ( c_initial_condition != nullptr ) {

        STDG::OrdinateFlux::AXPY( 1.0, 0.0, *c_initial_condition, *this->collided );
        this->c_solver->SetInitialGuess( this->collided );

    } else {

        this->c_solver->SetInitialGuess( nullptr );
    }

    // Construct scattering source from uncollided flux.
    TransportOperator::Pmv( 1.0, 0.0, *this->uncollided, this->phi );
    this->phi.ZeroDensity( OpDomain::Boundary );
    TransportOperator::Smv( 1.0, 0.0, sigma_s, this->phi, *this->collided );

    // Perform collided solve.
    this->c_solver->Solve( *this->collided, sigma_t, sigma_s, h, c_initial_condition );

    Global::TMR_collided.Stop();

    // --- Relabel ------------------------------------------------------------------------------------ //

    if (    this->relabel == false
         && c_initial_condition != nullptr
    ) {

        Global::TMR_uncollided.Start();
        STDG::OrdinateFlux::AXPY( 1.0, 0.0, *this->uncollided, u_initial_condition );
        Global::TMR_uncollided.Stop();

        Global::TMR_collided.Start();
        STDG::OrdinateFlux::AXPY( 1.0, 0.0, *this->collided, *c_initial_condition );
        Global::TMR_collided.Stop();

        c_initial_condition->ZeroDensity( OpDomain::Boundary );

    } else {

        switch ( this->hybrid_method ) {

            case HybridMethod::Hybrid_a:
            {
                Global::TMR_uncollided.Start();
                STDG::OrdinateFlux::AXPY( 1.0, 0.0, *this->uncollided, u_initial_condition );
                Global::TMR_uncollided.Stop();

                this->relabel_operator->Relabel( 1.0, *this->collided, u_initial_condition );

                break;
            }

            case HybridMethod::Hybrid_c:
            {
                Global::TMR_uncollided.Start();

                TransportOperator::Pmv( 1.0, 0.0, *this->uncollided, this->phi );
                TransportOperator::Pmv( 1.0, 1.0, *this->collided,   this->phi );
                this->phi.ZeroDensity( OpDomain::Boundary );

                STDG::OrdinateFlux::AXPY( 1.0, 0.0, source, *this->uncollided );
                TransportOperator::Smv( 1.0, 1.0, sigma_s, this->phi, *this->uncollided );

                this->u_solver->Solve( *this->uncollided, sigma_t, sigma_s, h, &u_initial_condition );

                STDG::OrdinateFlux::AXPY( 1.0, 0.0, *this->uncollided, u_initial_condition );

                Global::TMR_uncollided.Stop();

            } break;

            default:
            {   std::string error_message = "Invalid hybrid splitting '"
                                            + HybridMethod_to_String.at( this->hybrid_method )
                                            + "' in '" + std::string(__func__) + "'.\n";

                PRINT_ERROR( error_message.c_str() )
                throw std::invalid_argument( error_message );
            }
        }
    }

    u_initial_condition.ZeroDensity( OpDomain::Boundary );
}
