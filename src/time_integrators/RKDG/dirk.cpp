//------------------------------------------------------------------------------------------------------------
//! \file   time_integrators/RKDG/dirk.cpp
//! \brief  Implementation of diagonally implicit Runge-Kutta (DIRK) time integration schemes.
//!
//! \author Michael Crockatt
//! \date   January 2018
//------------------------------------------------------------------------------------------------------------

# include <cinttypes>
# include <cmath>
# include <cstdint>
# include <cstdlib>
# include <cstring>

# include "objects/RKDG/OrdinateFlux.hpp"
# include "objects/RKDG/CrossSection.hpp"
# include "operators/TransportOperator.hpp"
# include "time_integrators/RKDG/dirk.hpp"
# include "utils/CLog.hpp"
# include "utils/global.hpp"


using namespace RKDG;
using namespace Quadrule;


//============================================================================================================
//=== STATIC DEFINITIONS =====================================================================================
//============================================================================================================


//------------------------------------------------------------------------------------------------------------
//! \brief  Used to convert DIRKIntegrator::DIRKMethod values to their corresponding string descriptors.
//------------------------------------------------------------------------------------------------------------
std::map< DIRKIntegrator::DIRKMethod, std::string > DIRKIntegrator::DIRKMethod_to_String = {

    { DIRKIntegrator::DIRKMethod::None,             "none"          },
    { DIRKIntegrator::DIRKMethod::Euler,            "euler"         },
    { DIRKIntegrator::DIRKMethod::SDIRK2,           "sdirk2"        },
    { DIRKIntegrator::DIRKMethod::SDIRK3,           "sdirk3"        },
    { DIRKIntegrator::DIRKMethod::SDIRK5,           "sdirk5"        },
    { DIRKIntegrator::DIRKMethod::SDIRK_5_3_4,      "sdirk-5-3-4"   },
    { DIRKIntegrator::DIRKMethod::KVAERNO_4_2_3,    "kvaerno-4-2-3" },
    { DIRKIntegrator::DIRKMethod::KVAERNO_7_4_5,    "kvaerno-7-4-5" },
    { DIRKIntegrator::DIRKMethod::ARK_8_4_5,        "ark-8-4-5"     },
    { DIRKIntegrator::DIRKMethod::GSBP_DIRK_3,      "gsbp-dirk-3"   },
    { DIRKIntegrator::DIRKMethod::GSBP_DIRK_4,      "gsbp-dirk-4"   }
};


//------------------------------------------------------------------------------------------------------------
//! \brief  Used to convert string descriptors to DIRKIntegrator::DIRKMethod values.
//------------------------------------------------------------------------------------------------------------
std::map< std::string, DIRKIntegrator::DIRKMethod > DIRKIntegrator::String_to_DIRKMethod = {

    { "none",           DIRKIntegrator::DIRKMethod::None            },
    { "euler",          DIRKIntegrator::DIRKMethod::Euler           },
    { "sdirk2",         DIRKIntegrator::DIRKMethod::SDIRK2          },
    { "sdirk3",         DIRKIntegrator::DIRKMethod::SDIRK3          },
    { "sdirk5",         DIRKIntegrator::DIRKMethod::SDIRK5          },
    { "sdirk-5-3-4",    DIRKIntegrator::DIRKMethod::SDIRK_5_3_4     },
    { "kvaerno-4-2-3",  DIRKIntegrator::DIRKMethod::KVAERNO_4_2_3   },
    { "kvaerno-7-4-5",  DIRKIntegrator::DIRKMethod::KVAERNO_7_4_5   },
    { "ark-8-4-5",      DIRKIntegrator::DIRKMethod::ARK_8_4_5       },
    { "gsbp-dirk-3",    DIRKIntegrator::DIRKMethod::GSBP_DIRK_3     },
    { "gsbp-dirk-4",    DIRKIntegrator::DIRKMethod::GSBP_DIRK_4     }
};


//------------------------------------------------------------------------------------------------------------
//! \brief  Used to convert DIRKIntegrator::HybridMethod values to their corresponding string descriptors.
//------------------------------------------------------------------------------------------------------------
std::map< DIRKIntegrator::HybridMethod, std::string > DIRKIntegrator::HybridMethod_to_String = {

    { DIRKIntegrator::HybridMethod::None,       "none"          },
    { DIRKIntegrator::HybridMethod::HybridIa,   "hybrid-ia"     },
    { DIRKIntegrator::HybridMethod::HybridIc,   "hybrid-ic"     },
    { DIRKIntegrator::HybridMethod::HybridIIa,  "hybrid-iia"    },
    { DIRKIntegrator::HybridMethod::HybridIIc,  "hybrid-iic"    }
};


//------------------------------------------------------------------------------------------------------------
//! \brief  Used to convert string descriptors to DIRKIntegrator::HybridMethod values.
//------------------------------------------------------------------------------------------------------------
std::map< std::string, DIRKIntegrator::HybridMethod > DIRKIntegrator::String_to_HybridMethod = {

    { "none",       DIRKIntegrator::HybridMethod::None      },
    { "hybrid-ia",  DIRKIntegrator::HybridMethod::HybridIa  },
    { "hybrid-ic",  DIRKIntegrator::HybridMethod::HybridIc  },
    { "hybrid-iia", DIRKIntegrator::HybridMethod::HybridIIa },
    { "hybrid-iic", DIRKIntegrator::HybridMethod::HybridIIc }
};


//============================================================================================================
//=== CONSTRUCTORS AND DESTRUCTOR ============================================================================
//============================================================================================================


//------------------------------------------------------------------------------------------------------------
//! \brief  Constructs a DIRKIntegrator object using the provided values.
//!
//! \param[in]  input_list      List of parameters to use for initialization.
//------------------------------------------------------------------------------------------------------------
DIRKIntegrator::DIRKIntegrator (

    const ParameterList & input_list
) :
    TimeIntegrator( input_list ),

    dirk_method(),

    A{ nullptr },
    b{ nullptr },
    c{ nullptr },

    stiffly_accurate{},
    num_stages{},
    stages_allocd{},
    stage_alloc_map{ nullptr },

    phi{ nullptr },

    stage_flux{ nullptr },
    solver{ nullptr },

    hybrid_method(),

    u_stage_flux{ nullptr },
    u_solver{ nullptr },
    c_stage_flux{ nullptr },
    c_solver{ nullptr },
    c_temp{ nullptr },

    relabel_operator{ nullptr }
{

    bool symmetric_reduce = false;

    ParameterList   params = input_list;
    ParameterList u_params = MakeUncollidedList( input_list );
    ParameterList c_params = MakeCollidedList( input_list );

    // Read polynomial degree for DG approximation.
    int64_t DG_degree = GetDGDegreeX( params );

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


    // --- Setup DIRK method. ------------------------------------------------------------------------- //

    params.GetValue( GetInputKey(), this->dirk_method, DIRKIntegrator::String_to_DIRKMethod );

    switch ( this->hybrid_method ) {

        case HybridMethod::HybridIa:
        case HybridMethod::HybridIc:
        {
            try         {  params.GetValue( "relabel", this->relabel );  }
            catch (...) {  this->relabel = true;                                     }

        } /* fall through */
        case HybridMethod::HybridIIa:
        case HybridMethod::HybridIIc:
        break;

        case HybridMethod::None:
        default:
        break;
    }

    SetTableau();


    // --- Setup objects used by DIRK integrator. ----------------------------------------------------- //

    // Hybrid/non-hybrid setup.
    switch ( this->hybrid_method ) {

        case HybridMethod::None:
        {
            this->stage_flux = new OrdinateFlux*[ this->num_stages ];

            for ( int64_t i = 0; i < this->num_stages; ++i ) {

                if ( this->stage_alloc_map[i] == i )
                    this->stage_flux[i] = new OrdinateFlux( *this, ordinate_set, DG_degree );
                else
                    this->stage_flux[i] = this->stage_flux[ this->stage_alloc_map[i] ];
            }

            this->phi = new ScalarFlux( *this, DG_degree );

            this->solver = std::shared_ptr<ImplicitSolver<OrdinateFlux>>(
                    ImplicitSolver<OrdinateFlux>::ImplicitSolverFactory::CreateSolver(
                        *this, ordinate_set, params
                ) );
        } break;

        case HybridMethod::HybridIa:
        case HybridMethod::HybridIc:
        case HybridMethod::HybridIIa:
        case HybridMethod::HybridIIc:
        {
            // --- Options common to all hybrid methods. --- //

            this->u_stage_flux = new OrdinateFlux*[ this->num_stages ];

            for ( int64_t i = 0; i < this->num_stages; ++i ) {

                if ( this->stage_alloc_map[i] == i )
                    this->u_stage_flux[i] = new OrdinateFlux( *this, u_ordinate_set, DG_degree );
                else
                    this->u_stage_flux[i] = this->u_stage_flux[ this->stage_alloc_map[i] ];
            }

            // Pointer aliasing for ComputeOutput_Nonhybrid.
            this->stage_flux = this->u_stage_flux;

            // Uncollided solver.
            this->u_solver = std::shared_ptr<ImplicitSolver<OrdinateFlux>>(
                    ImplicitSolver<OrdinateFlux>::ImplicitSolverFactory::CreateSolver(
                        *this, u_ordinate_set, u_params
                ) );

            // Collided solver.
            this->c_solver = std::shared_ptr<ImplicitSolver<OrdinateFlux>>(
                    ImplicitSolver<OrdinateFlux>::ImplicitSolverFactory::CreateSolver(
                        *this, c_ordinate_set, c_params
                ) );

            // --- Options differing between hybrid methods. --- //

            // Allocate collided flux vectors.
            if (    this->hybrid_method == HybridMethod::HybridIa
                 || this->hybrid_method == HybridMethod::HybridIc
            ) {

                this->c_stage_flux = new OrdinateFlux*[ this->num_stages ];

                for ( int64_t i = 0; i < this->num_stages; ++i ) {

                    if ( this->stage_alloc_map[i] == i )
                        this->c_stage_flux[i] = new OrdinateFlux( *this, c_ordinate_set, DG_degree );
                    else
                        this->c_stage_flux[i] = this->c_stage_flux[ this->stage_alloc_map[i] ];
                }

            } else if (    this->hybrid_method == HybridMethod::HybridIIa
                        || this->hybrid_method == HybridMethod::HybridIIc
            ) {
                this->c_temp = new OrdinateFlux( *this, c_ordinate_set, DG_degree );
            }

            // Allocate scalar flux vectors.
            if ( this->hybrid_method == HybridMethod::HybridIc ) {

                this->phi = new ScalarFlux[ this->num_stages ];

                for ( int64_t i = 0; i < this->num_stages; ++i )
                    this->phi[i].Reconfigure( *this, DG_degree );

            } else {

                this->phi = new ScalarFlux( *this, DG_degree );
            }

            // Setup relabel operator.
            if (    this->hybrid_method == HybridMethod::HybridIa
                 || this->hybrid_method == HybridMethod::HybridIIa
            ) {
                this->relabel_operator = new RelabelOperator( c_ordinate_set, u_ordinate_set );
            }
        } break;

        default:
        {   std::string error_message =   "Invalid hybrid splitting '"
                                        + HybridMethod_to_String.at( this->hybrid_method )
                                        + "' for initializing DIRK integrator in '"
                                        + std::string(__func__)
                                        + "'.\n";

            PRINT_ERROR( error_message.c_str() )
            throw std::invalid_argument( error_message );
        }
    }
}


//------------------------------------------------------------------------------------------------------------
//! \brief  Destructor for class DIRKIntegrator.
//------------------------------------------------------------------------------------------------------------
DIRKIntegrator::~DIRKIntegrator ( void ) {

    delete [] this->A;
    delete [] this->b;
    delete [] this->c;

    delete this->c_temp;
    delete this->relabel_operator;

    if ( this->hybrid_method == HybridMethod::HybridIc )
        delete [] this->phi;
    else
        delete this->phi;

    if ( this->hybrid_method == HybridMethod::None ) {

        if ( this->stage_flux != nullptr ) {

            for ( int64_t i = 0; i < this->num_stages; ++i ) {

                if ( this->stage_alloc_map[i] == i )
                    delete this->stage_flux[i];
            }

            delete [] this->stage_flux;
        }

    } else {

        if ( this->u_stage_flux != nullptr ) {

            for ( int64_t i = 0; i < this->num_stages; ++i ) {

                if ( this->stage_alloc_map[i] == i )
                    delete this->u_stage_flux[i];
            }

            delete [] this->u_stage_flux;
        }

        if ( this->c_stage_flux != nullptr ) {

            for ( int64_t i = 0; i < this->num_stages; ++i ) {

                if ( this->stage_alloc_map[i] == i )
                    delete this->c_stage_flux[i];
            }

            delete [] this->c_stage_flux;
        }
    }

    delete [] this->stage_alloc_map;
}


//============================================================================================================
//=== INTERFACE ROUTINES =====================================================================================
//============================================================================================================


//------------------------------------------------------------------------------------------------------------
//! \brief  Calls the appropriate \c Step_* function corresponding to the desired DIRK hybrid splitting as set
//!         by DIRKIntegrator::hybrid_method.
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
//! \see    DIRKIntegrator::Step_Nonhybrid()
//! \see    DIRKIntegrator::Step_HybridI()
//! \see    DIRKIntegrator::Step_HybridII()
//------------------------------------------------------------------------------------------------------------
DIRKIntegrator & DIRKIntegrator::Step (

    OrdinateFlux & initial_condition,
    const OrdinateFlux & source,
    const CrossSection & sigma_t,
    const CrossSection & sigma_s,
    const double step_size
) {

# if LOGLEVEL >= 1
    PRINT_LOG( "\n" )
# endif

    switch ( this->hybrid_method ) {

        case HybridMethod::None:

            Step_Nonhybrid( initial_condition, source, sigma_t, sigma_s, step_size );
            break;

        case HybridMethod::HybridIa:
        case HybridMethod::HybridIc:

            Step_HybridI( initial_condition, nullptr, source, sigma_t, sigma_s, step_size );
            break;

        case HybridMethod::HybridIIa:
        case HybridMethod::HybridIIc:

            Step_HybridII( initial_condition, source, sigma_t, sigma_s, step_size );
            break;

        default:
        {   std::string error_message = "Invalid hybrid splitting '"
                                        + HybridMethod_to_String.at( this->hybrid_method )
                                        + "' for DIRK step.\n";

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
//! \brief  Calls the appropriate \c Step_* function corresponding to the desired DIRK hybrid splitting as set
//!         by DIRKIntegrator::hybrid_method.
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
//! \see    DIRKIntegrator::Step_Nonhybrid()
//! \see    DIRKIntegrator::Step_HybridI()
//! \see    DIRKIntegrator::Step_HybridII()
//------------------------------------------------------------------------------------------------------------
DIRKIntegrator & DIRKIntegrator::Step (

    OrdinateFlux & u_initial_condition,
    OrdinateFlux * c_initial_condition,
    const OrdinateFlux & source,
    const CrossSection & sigma_t,
    const CrossSection & sigma_s,
    const double step_size
) {

# if LOGLEVEL >= 1
    PRINT_LOG( "\n" )
# endif

    switch ( this->hybrid_method ) {

        case HybridMethod::None:

        # if defined (STRICT_CHECK)

            if ( c_initial_condition != nullptr ) {

                PRINT_WARNING( "DIRK integrator with hybrid splitting '%s' ignores pointer to collided "
                               "initial condition.\n",
                               HybridMethod_to_String.at( this->hybrid_method ).c_str() )
            }

        # endif // if defined (STRICT_CHECK)

            Step_Nonhybrid( u_initial_condition, source, sigma_t, sigma_s, step_size );
            break;

        case HybridMethod::HybridIa:
        case HybridMethod::HybridIc:

            Step_HybridI( u_initial_condition, c_initial_condition, source, sigma_t, sigma_s, step_size );
            break;

        case HybridMethod::HybridIIa:
        case HybridMethod::HybridIIc:

        # if defined (STRICT_CHECK)

            if ( c_initial_condition != nullptr ) {

                PRINT_WARNING( "DIRK integrator with hybrid splitting '%s' ignores pointer to collided "
                               "initial condition.\n",
                               HybridMethod_to_String.at( this->hybrid_method ).c_str() )
            }

        # endif // if defined (STRICT_CHECK)

            Step_HybridII( u_initial_condition, source, sigma_t, sigma_s, step_size );
            break;

        default:
        {   std::string error_message = "Invalid hybrid splitting '"
                                        + HybridMethod_to_String.at( this->hybrid_method )
                                        + "' for DIRK step.\n";

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
//! \brief  Returns the key used to determine the type of DIRK method.
//------------------------------------------------------------------------------------------------------------
std::string DIRKIntegrator::GetInputKey( void ) {

    PRINT_STATUS( "Executing DIRKIntegrator::%s.\n", __func__ )

    return "dirk_method";
}


//------------------------------------------------------------------------------------------------------------
//! \brief  Prints a summary of the integrator configuration to the logging interface.
//!
//! \param[in]  prefix  String prepended to each line of output, intended to be used for controlling
//!                     indentation of output.
//------------------------------------------------------------------------------------------------------------
void DIRKIntegrator::Print (

    const std::string & prefix // = "  "

) const {

    PRINT_LOG( "%s%-*s  %s\n", prefix.c_str(),
               Global::col_width, "Time integrator:",
               TimeIntegratorType_to_String.at( GetIntegratorType() ).c_str() )

    PRINT_LOG( "%s%-*s  %s\n", prefix.c_str(),
               Global::col_width, "DIRK method:",
               DIRKIntegrator::DIRKMethod_to_String.at( this->dirk_method ).c_str() )

    this->Abstract::TimeIntegrator::Print( prefix );

# if LOGLEVEL >= 2
    this->PrintTableau( prefix );
# endif

    PRINT_LOG( "%s%-*s  %s\n", prefix.c_str(),
               Global::col_width, "Hybrid:",
               HybridMethod_to_String.at( this->hybrid_method ).c_str() )

    if ( IsHybrid() ) {

        PRINT_LOG( "\n" )
        PRINT_LOG( "%s%-*s  %s\n", prefix.c_str(),
                   Global::col_width, "Relabel:", ( this->relabel ? "yes" : "no" ) )

        if (    this->relabel
             && (    this->hybrid_method == HybridMethod::HybridIa
                  || this->hybrid_method == HybridMethod::HybridIIa
                )
        ) {
            this->relabel_operator->Print( prefix );
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
//! \brief  Returns true if the integrator has been constructed with hybrid splitting and false otherwise.
//------------------------------------------------------------------------------------------------------------
bool DIRKIntegrator::IsHybrid ( void ) const {

    switch ( this->hybrid_method ) {

        case HybridMethod::HybridIa:
        case HybridMethod::HybridIc:
        case HybridMethod::HybridIIa:
        case HybridMethod::HybridIIc:

            return true;

        default:
            return false;
    }
}


//------------------------------------------------------------------------------------------------------------
//! \brief  Prints the Butcher tableau of the DIRK method to the logging interface.
//!
//! \param[in]  prefix  String prepended to each line of output, intended to be used for controlling
//!                     indentation of output.
//------------------------------------------------------------------------------------------------------------
void DIRKIntegrator::PrintTableau (

    const std::string prefix    // = "  "

) const {

    // --- Check to make sure integrator is properly initialized. ------------------------------------- //

    if ( this->dirk_method == DIRKIntegrator::DIRKMethod::None ) {

        std::string error_message =   "No Butcher tableau to print for DIRKMethod '"
                                    + DIRKMethod_to_String.at( this->dirk_method )
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
        PRINT_LOG( "%s % -*.4e |", prefix.c_str(), Global::col_width, this->c[i] )

        // Print row of A.
        for ( int64_t j = 0; j <= i; ++j )
            PRINT_LOG( " % -*.4e ", Global::col_width, IA(i,j) )

        PRINT_LOG( "\n" )
    }

    // Print horizontal line.
    PRINT_LOG( "%s-%*s-|", prefix.c_str(), Global::col_width, buffer )

    for ( int64_t i = 0; i < this->num_stages; ++i ) {  PRINT_LOG( "-%*s-", Global::col_width, buffer )  }

    PRINT_LOG( "\n" )

    // Print last row for b.
    buffer[0] = '\0';

    PRINT_LOG( "%s %*s |", prefix.c_str(), Global::col_width, buffer )

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
//! \see    DIRKIntegrator::Step()
//! \see    DIRKIntegrator::ComputeOutput_Nonhybrid()
//------------------------------------------------------------------------------------------------------------
void DIRKIntegrator::Step_Nonhybrid (

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

        // Explicit stage.
        if ( IA(i,i) == 0.0 ) {

            stage_coeff = h;
            OrdinateFlux::Copy( *this->stage_flux[i], initial_condition );

        // Implicit stage.
        } else {

            stage_coeff = 1.0 / IA(i,i);
            OrdinateFlux::AXPY( 1.0 / (h * IA(i,i)), 0.0, initial_condition, *this->stage_flux[i] );
        }

        source_coeff = 0.0;
        this->phi->ZeroDensity();

        // Compute weighted sum of transport operator applied to previous stages.
        for ( int64_t j = 0; j < i; ++j ) {

            if ( IA(i,j) == 0.0 ) {  continue;  }

            source_coeff += IA(i,j);

            TransportOperator::Lmv( - stage_coeff * IA(i,j), 1.0, sigma_t, *this->stage_flux[j],
                                    *this->stage_flux[i] );
            TransportOperator::Pmv( stage_coeff * IA(i,j), 1.0, *this->stage_flux[j], *this->phi );
        }

        this->phi->ZeroDensity( OpDomain::Boundary );
        TransportOperator::Smv( 1.0, 1.0, sigma_s, *this->phi, *this->stage_flux[i] );

        // Include source from previous stages.
        OrdinateFlux::AXPY( stage_coeff * source_coeff, 1.0, source, *this->stage_flux[i] );
        this->stage_flux[i]->ZeroDensity( OpDomain::Boundary );

        // Explicit stage: Just set boundary conditions on stage vector.
        if ( IA(i,i) == 0.0 ) {

            this->stage_flux[i]->Copy( source, OpDomain::Boundary );

            if ( BitmaskHasAny( Global::reflecting_boundaries, BoundaryEdge::All ) )
                this->stage_flux[i]->ReflectBoundaries();

        // Implicit stage: Add source term and solve system.
        } else {

            // Include source term for current stage.
            OrdinateFlux::AXPY( 1.0, 1.0, source, *this->stage_flux[i] );

            // Set initial guess for solve.
            if ( i > 0 )
                this->solver->SetInitialGuess( this->stage_flux[i-1] );
            else
                this->solver->SetInitialGuess( &initial_condition );

            // Perform implicit solve.
            this->solver->Solve( *this->stage_flux[i], sigma_t, sigma_s, h * IA(i,i) );
        }
    }

    // --- Compute return value from intermediate stages. --------------------------------------------- //

    if ( !this->stiffly_accurate )
        ComputeOutput_Nonhybrid( initial_condition, source, sigma_t, sigma_s, h );

    this->stage_flux[this->num_stages - 1]->ZeroDensity( OpDomain::Boundary );
    OrdinateFlux::swap( initial_condition, *this->stage_flux[this->num_stages - 1] );
}


//------------------------------------------------------------------------------------------------------------
//! \brief  Computes a one timestep update using hybrid-I type splitting.
//!
//! \attention  If DIRKIntegrator::relabel is "true" then the initial collided distribution is assumed to be
//!             zero and \pp{c_initial_condition} is not dereferenced.
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
//! \see    DIRKIntegrator::HybridMethod::HybridI
//! \see    DIRKIntegrator::Step()
//------------------------------------------------------------------------------------------------------------
void DIRKIntegrator::Step_HybridI (

    OrdinateFlux & u_initial_condition,
    OrdinateFlux * const c_initial_condition,
    const OrdinateFlux & source,
    const CrossSection & sigma_t,
    const CrossSection & sigma_s,
    const double h
) {

    double stage_coeff, source_coeff;

    // --- Construct stage values in "Y" formulation. ------------------------------------------------- //

    // Compute each internal stage.
    for ( int64_t i = 0; i < this->num_stages; ++i ) {

        // Explicit stage.
        if ( IA(i,i) == 0.0 ) {

            stage_coeff = h;

            // --- High-resolution (uncollided). --- //

            Global::TMR_uncollided.Start();
            OrdinateFlux::Copy( *this->u_stage_flux[i], u_initial_condition );
            Global::TMR_uncollided.Stop();

            // --- Low-resolution (collided), if not relabeling. --- //

            Global::TMR_collided.Start();

            if (    this->relabel == false
                 && c_initial_condition != nullptr
            ) {
                OrdinateFlux::Copy( *this->c_stage_flux[i], *c_initial_condition );

            } else {

                this->c_stage_flux[i]->ZeroDensity();
            }

            Global::TMR_collided.Stop();

        // Implicit stage.
        } else {

            stage_coeff = 1.0 / IA(i,i);

            // --- High-resolution (uncollided). --- //

            Global::TMR_uncollided.Start();
            OrdinateFlux::AXPY( 1.0 / (h * IA(i,i)), 0.0, u_initial_condition, *this->u_stage_flux[i] );
            Global::TMR_uncollided.Stop();

            // --- Low-resolution (collided), if not relabeling. --- //

            Global::TMR_collided.Start();

            if (    this->relabel == false
                 && c_initial_condition != nullptr
            ) {
                OrdinateFlux::AXPY( 1.0 / (h * IA(i,i)), 0.0, *c_initial_condition, *this->c_stage_flux[i] );

            } else {

                this->c_stage_flux[i]->ZeroDensity();
            }

            Global::TMR_collided.Stop();
        }

        source_coeff = 0.0;
        ScalarFlux & current_phi = ( this->hybrid_method == HybridMethod::HybridIc ? this->phi[i] : *this->phi );
        current_phi.ZeroDensity();

        // Compute weighted sum of hybrid transport operator applied to previous stages.
        for ( int64_t j = 0; j < i; ++j ) {

            if ( IA(i,j) == 0.0 ) {  continue;  }

            source_coeff += IA(i,j);

            // --- High-resolution (uncollided). --- //

            Global::TMR_uncollided.Start();

            TransportOperator::Lmv( - stage_coeff * IA(i,j), 1.0, sigma_t, *this->u_stage_flux[j],
                                    *this->u_stage_flux[i] );

            Global::TMR_uncollided.Stop();

            // --- Low-resolution (collided). --- //

            Global::TMR_collided.Start();

            TransportOperator::Lmv( - stage_coeff * IA(i,j), 1.0, sigma_t, *this->c_stage_flux[j],
                                    *this->c_stage_flux[i] );
            TransportOperator::Pmv( stage_coeff * IA(i,j), 1.0, *this->u_stage_flux[j], current_phi );
            TransportOperator::Pmv( stage_coeff * IA(i,j), 1.0, *this->c_stage_flux[j], current_phi );

            Global::TMR_collided.Stop();
        }

        Global::TMR_collided.Start();

        current_phi.ZeroDensity( OpDomain::Boundary );
        TransportOperator::Smv( 1.0, 1.0, sigma_s, current_phi, *this->c_stage_flux[i] );

        Global::TMR_collided.Stop();
        Global::TMR_uncollided.Start();

        // Include source from previous stages.
        OrdinateFlux::AXPY( stage_coeff * source_coeff, 1.0, source, *this->u_stage_flux[i] );
        this->u_stage_flux[i]->ZeroDensity( OpDomain::Boundary );

        Global::TMR_uncollided.Stop();

        // Explicit stage: Just set boundary conditions on uncollided flux.
        if ( IA(i,i) == 0.0 ) {

            this->u_stage_flux[i]->Copy( source, OpDomain::Boundary );

            if ( BitmaskHasAny( Global::reflecting_boundaries, BoundaryEdge::All ) )
                this->u_stage_flux[i]->ReflectBoundaries();

        // Implicit stage: Add source term and solve systems.
        } else {

            Global::TMR_uncollided.Start();

            // Include source term for current stage.
            OrdinateFlux::AXPY( 1.0, 1.0, source, *this->u_stage_flux[i] );

            // Perform high-resolution (uncollided) implicit solve.
            this->u_solver->Solve( *this->u_stage_flux[i], sigma_t, sigma_s, h * IA(i,i) );

            Global::TMR_uncollided.Stop();
            Global::TMR_collided.Start();

            // Compute source for low-resolution (collided) implicit solve.
            TransportOperator::Pmv( 1.0, 0.0, *this->u_stage_flux[i], current_phi );
            current_phi.ZeroDensity( OpDomain::Boundary );
            TransportOperator::Smv( 1.0, 1.0, sigma_s, current_phi, *this->c_stage_flux[i] );

            // Set initial guess for solve.
            if ( i > 0 )
                this->c_solver->SetInitialGuess( this->c_stage_flux[i-1] );
            else
                this->c_solver->SetInitialGuess( c_initial_condition );

            // Perform low-resolution (collided) implicit solve.
            this->c_solver->Solve( *this->c_stage_flux[i], sigma_t, sigma_s, h * IA(i,i) );

            Global::TMR_collided.Stop();
        }

        // For hybrid-Ic methods, construct scattering source of current stage for Nyström reconstruction.
        if ( this->hybrid_method == HybridMethod::HybridIc ) {

            this->phi[i].ZeroDensity();

            if ( IA(i,i) != 0.0 ) {

                TransportOperator::Pmv( 1.0, 1.0, *this->u_stage_flux[i], this->phi[i] );
                TransportOperator::Pmv( 1.0, 1.0, *this->c_stage_flux[i], this->phi[i] );
            }
        }
    }

    // --- Compute return value(s) from intermediate stages. ------------------------------------------ //

    switch ( this->hybrid_method ) {

        case HybridMethod::HybridIa:
        {
            if ( !this->stiffly_accurate )
                ComputeOutput_Hybrid( u_initial_condition, c_initial_condition, source, sigma_t, sigma_s, h );

            OrdinateFlux::swap( u_initial_condition, *this->u_stage_flux[this->num_stages - 1] );

            if (    this->relabel == true
                 || c_initial_condition == nullptr
            ) {
                this->relabel_operator->Relabel( 1.0, *this->c_stage_flux[this->num_stages - 1],
                                                 u_initial_condition );
            } else {

                OrdinateFlux::swap( *c_initial_condition, *this->c_stage_flux[this->num_stages - 1] );
                c_initial_condition->ZeroDensity( OpDomain::Boundary );
            }
        } break;

        case HybridMethod::HybridIc:
        {
            // Do not perform reconstruction, and return both collided and uncollided fluxes.
            if (    this->relabel == false
                 && c_initial_condition != nullptr
            ) {

                if ( !this->stiffly_accurate ) {

                    ComputeOutput_Hybrid( u_initial_condition, c_initial_condition,
                                          source, sigma_t, sigma_s, h );
                }

                this->u_stage_flux[this->num_stages - 1]->ZeroDensity( OpDomain::Boundary );
                this->c_stage_flux[this->num_stages - 1]->ZeroDensity( OpDomain::Boundary );

                OrdinateFlux::swap(  u_initial_condition, *this->u_stage_flux[this->num_stages - 1] );
                OrdinateFlux::swap( *c_initial_condition, *this->c_stage_flux[this->num_stages - 1] );

            // Compute and return reconstructed solution.
            } else {

                Global::TMR_uncollided.Start();

                ComputeNystromInterpolants( u_initial_condition, source, sigma_t, sigma_s, h );

                if ( !this->stiffly_accurate )
                    ComputeOutput_Nonhybrid( u_initial_condition, source, sigma_t, sigma_s, h );

                this->u_stage_flux[this->num_stages - 1]->ZeroDensity( OpDomain::Boundary );
                OrdinateFlux::swap( u_initial_condition, *this->u_stage_flux[this->num_stages - 1] );

                Global::TMR_uncollided.Stop();
            }
        } break;

        default:
        {   std::string error_message = "Invalid hybrid splitting '"
                                        + HybridMethod_to_String.at( this->hybrid_method )
                                        + "' in '" + std::string(__func__) + "'.\n";

            PRINT_ERROR( error_message.c_str() )
            throw std::invalid_argument( error_message );
        }
    }

    // Zero ghost cells before return.
    u_initial_condition.ZeroDensity( OpDomain::Boundary );
}


//------------------------------------------------------------------------------------------------------------
//! \brief  Computes a one timestep update using hybrid-II type splitting.
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
//! \see    DIRKIntegrator::HybridMethod::HybridII
//! \see    DIRKIntegrator::Step()
//------------------------------------------------------------------------------------------------------------
void DIRKIntegrator::Step_HybridII (

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

        Global::TMR_uncollided.Start();

        // Explicit stage.
        if ( IA(i,i) == 0.0 ) {

            stage_coeff = h;
            OrdinateFlux::Copy( *this->u_stage_flux[i], initial_condition );

        // Implicit stage.
        } else {

            stage_coeff = 1.0 / IA(i,i);
            OrdinateFlux::AXPY( 1.0 / (h * IA(i,i)), 0.0, initial_condition, *this->u_stage_flux[i] );
        }

        source_coeff = 0.0;
        this->phi->ZeroDensity();

        // Compute weighted sum of transport operator applied to previous stages.
        for ( int64_t j = 0; j < i; ++j ) {

            if ( IA(i,j) == 0.0 ) {  continue;  }

            source_coeff += IA(i,j);

            TransportOperator::Lmv( - stage_coeff * IA(i,j), 1.0, sigma_t, *this->u_stage_flux[j],
                                    *this->u_stage_flux[i] );
            TransportOperator::Pmv( stage_coeff * IA(i,j), 1.0, *this->u_stage_flux[j], *this->phi );
        }

        this->phi->ZeroDensity( OpDomain::Boundary );
        TransportOperator::Smv( 1.0, 1.0, sigma_s, *this->phi, *this->u_stage_flux[i] );

        // Include source from previous stages.
        OrdinateFlux::AXPY( stage_coeff * source_coeff, 1.0, source, *this->u_stage_flux[i] );
        this->u_stage_flux[i]->ZeroDensity( OpDomain::Boundary );

        Global::TMR_uncollided.Stop();

        // Explicit stage: Just set boundary conditions on uncollided flux.
        if ( IA(i,i) == 0.0 ) {

            this->u_stage_flux[i]->Copy( source, OpDomain::Boundary );

            if ( BitmaskHasAny( Global::reflecting_boundaries, BoundaryEdge::All ) )
                this->u_stage_flux[i]->ReflectBoundaries();

        // Implicit stage: Add source term and solve systems.
        } else {

            Global::TMR_uncollided.Start();

            OrdinateFlux::AXPY( 1.0, 1.0, source, *this->u_stage_flux[i] );

            // Perform high-resolution (uncollided) implicit solve.
            this->u_solver->Solve( *this->u_stage_flux[i], sigma_t, sigma_s, h * IA(i,i) );

            Global::TMR_uncollided.Stop();
            Global::TMR_collided.Start();

            // Compute source for low-resolution (collided) implicit solve.
            TransportOperator::Pmv( 1.0, 0.0, *this->u_stage_flux[i], *this->phi );
            this->phi->ZeroDensity( OpDomain::Boundary );
            TransportOperator::Smv( 1.0, 0.0, sigma_s, *this->phi, *this->c_temp );

            // Set initial guess for solve.
            this->c_solver->SetInitialGuess( nullptr );

            // Perform low-resolution (collided) implicit solve.
            this->c_solver->Solve( *this->c_temp, sigma_t, sigma_s, h * IA(i,i) );

            Global::TMR_collided.Stop();

            // Compute reconstructed stage value.
            switch ( this->hybrid_method ) {

                case HybridMethod::HybridIIa:
                {
                    this->relabel_operator->Relabel( 1.0, *this->c_temp, *this->u_stage_flux[i] );
                    break;
                }

                case HybridMethod::HybridIIc:
                {
                    Global::TMR_uncollided.Start();

                    // Compute scalar flux distribution for current stage from hybrid components.
                    TransportOperator::Pmv( 1.0, 0.0, *this->u_stage_flux[i], *this->phi );
                    TransportOperator::Pmv( 1.0, 1.0, *this->c_temp, *this->phi );

                    OrdinateFlux::AXPY( 1.0 / (h * IA(i,i)), 0.0, initial_condition, *this->u_stage_flux[i] );

                    // Compute weighted sum of transport operator applied to previous stages.
                    for ( int64_t j = 0; j < i; ++j ) {

                        const double stage_weight = IA(i,j) / IA(i,i);

                        TransportOperator::Pmv( stage_weight, 1.0, *this->u_stage_flux[j], *this->phi );
                        TransportOperator::Lmv( -stage_weight, 1.0, sigma_t, *this->u_stage_flux[j],
                                                *this->u_stage_flux[i] );
                    }

                    this->phi->ZeroDensity( OpDomain::Boundary );
                    TransportOperator::Smv( 1.0, 1.0, sigma_s, *this->phi, *this->u_stage_flux[i] );

                    // Include source term from previous stages.
                    OrdinateFlux::AXPY( stage_coeff * source_coeff, 1.0, source, *this->u_stage_flux[i] );
                    this->u_stage_flux[i]->ZeroDensity( OpDomain::Boundary );

                    // Include source and boundary conditions for current stage.
                    OrdinateFlux::AXPY( 1.0, 1.0, source, *this->u_stage_flux[i] );

                    this->u_solver->Solve( *this->u_stage_flux[i], sigma_t, sigma_s, h * IA(i,i) );

                    Global::TMR_uncollided.Stop();
                    break;
                }

                default:
                {   std::string error_message = "Invalid hybrid splitting '"
                                                + HybridMethod_to_String.at( this->hybrid_method )
                                                + "' in '" + std::string(__func__) + "'.\n";

                    PRINT_ERROR( error_message.c_str() )
                    throw std::invalid_argument( error_message );
                }
            }
        }
    }

    // --- Compute return value from intermediate stages. --------------------------------------------- //

    if ( !this->stiffly_accurate ) {

        Global::TMR_uncollided.Start();
        ComputeOutput_Nonhybrid( initial_condition, source, sigma_t, sigma_s, h );
        Global::TMR_uncollided.Stop();
    }

    OrdinateFlux::swap( initial_condition, *this->u_stage_flux[this->num_stages - 1] );
    initial_condition.ZeroDensity( OpDomain::Boundary );
}


//============================================================================================================
//=== PRIVATE HELPER ROUTINES ================================================================================
//============================================================================================================


//------------------------------------------------------------------------------------------------------------
//! \brief  Computes the final output value of the DIRK method using the intermediate stage values of the
//!         non-hybrid or reconstructed hybrid approximation.
//!
//! The output value is stored in <code> this->stage_flux[ this->num_stages - 1 ] </code>.
//!
//! \param[in]      initial_condition   The initial condition \f$ \Psi^n \f$ for the current timestep.
//! \param[in]      sigma_t             RKDG::CrossSection object containing the total cross section
//!                                     \f$ \sigma_{\mathrm{t}} \f$.
//! \param[in]      sigma_s             RKDG::CrossSection object containing the scattering cross section
//!                                     \f$ \sigma_{\mathrm{s}} \f$.
//! \param[in]      source              Source term of the transport system. <br>
//!                                     The source term is assumed to be constant across the timestep interval.
//! \param[in]      h                   Timestep size.
//!
//! \see    DIRKIntegrator::Step_Nonhybrid()
//! \see    DIRKIntegrator::Step_HybridI()
//! \see    DIRKIntegrator::Step_HybridII()
//------------------------------------------------------------------------------------------------------------
void DIRKIntegrator::ComputeOutput_Nonhybrid (

    const OrdinateFlux & initial_condition,
    const OrdinateFlux & source,
    const CrossSection & sigma_t,
    const CrossSection & sigma_s,
    const double h
) {

    OrdinateFlux & result = *this->stage_flux[this->num_stages - 1];
    this->phi->ZeroDensity();

    // Apply transport operators to stage vectors.
    for ( int64_t i = this->num_stages - 1; i >= 0; --i ) {

        if ( this->b[i] == 0.0 ) {  continue;  }

        TransportOperator::Pmv( h * this->b[i], 1.0, *this->stage_flux[i], *this->phi );
        TransportOperator::Lmv( -h * this->b[i], ( i == this->num_stages - 1 ? 0.0 : 1.0 ),  sigma_t,
                                *this->stage_flux[i], result );
    }

    TransportOperator::Smv( 1.0, 1.0, sigma_s, *this->phi, result );

    // Include initial condition and source for timestep.
    OrdinateFlux::AXPY( 1.0, 1.0, initial_condition, result );
    OrdinateFlux::AXPY( h, 1.0, source, result );
}


//------------------------------------------------------------------------------------------------------------
//! \brief  Computes the final output values of the DIRK method using the intermediate stage values of the
//!         hybrid approximation.
//!
//! The output values are stored in <code> this->u_stage_flux[this->num_stages - 1] </code> and
//! <code> this->c_stage_flux[this->num_stages - 1] </code>.
//!
//! \attention  If \c DIRKIntegrator::relabel is "true" then the initial collided distribution is assumed to
//!             be zero and \pp{c_initial_condition} is not dereferenced.
//!
//! \param[in]      u_initial_condition Contains the initial condition \f$ \Psi^n \f$ for the uncollided flux
//!                                     in the current timestep.
//! \param[in]      c_initial_condition Contains the initial condition \f$ \Psi^n \f$ for the collided flux in
//!                                     the current timestep. <br>
//!                                     If this is \c null, then the initial value for the collided flux is
//!                                     assumed to be zero.
//! \param[in]      sigma_t             RKDG::CrossSection object containing the total cross section
//!                                     \f$ \sigma_{\mathrm{t}} \f$.
//! \param[in]      sigma_s             RKDG::CrossSection object containing the scattering cross section
//!                                     \f$ \sigma_{\mathrm{s}} \f$.
//! \param[in]      source              Source term of the transport system. <br>
//!                                     The source term is assumed to be constant across the timestep interval.
//! \param[in]      h                   Timestep size.
//!
//! \see    DIRKIntegrator::Step_HybridI()
//------------------------------------------------------------------------------------------------------------
void DIRKIntegrator::ComputeOutput_Hybrid (

    const OrdinateFlux & u_initial_condition,
    const OrdinateFlux * const c_initial_condition,
    const OrdinateFlux & source,
    const CrossSection & sigma_t,
    const CrossSection & sigma_s,
    const double h
) {

    OrdinateFlux & u_result = *this->u_stage_flux[this->num_stages - 1];
    OrdinateFlux & c_result = *this->c_stage_flux[this->num_stages - 1];

    this->phi->ZeroDensity();

    // --- Construct uncollided result. --------------------------------------------------------------- //

    Global::TMR_uncollided.Start();

    // Apply transport operators to stage vectors.
    for ( int64_t i = this->num_stages - 1; i >= 0; --i ) {

        if ( this->b[i] == 0.0 ) {  continue;  }

        TransportOperator::Pmv( h * this->b[i], 1.0, *this->u_stage_flux[i], *this->phi );
        TransportOperator::Lmv( -h * this->b[i], ( i == this->num_stages - 1 ? 0.0 : 1.0 ), sigma_t,
                                *this->u_stage_flux[i], u_result );
    }

    // Include initial condition and source for timestep.
    OrdinateFlux::AXPY( 1.0, 1.0, u_initial_condition, u_result );
    OrdinateFlux::AXPY( h, 1.0, source, u_result );

    Global::TMR_uncollided.Stop();

    // --- Construct collided result. ----------------------------------------------------------------- //

    Global::TMR_collided.Start();

    // Apply transport operators to stage vectors.
    for ( int64_t i = this->num_stages - 1; i >= 0; --i ) {

        if ( this->b[i] == 0.0 ) {  continue;  }

        TransportOperator::Pmv( h * this->b[i], 1.0, *this->c_stage_flux[i], *this->phi );
        TransportOperator::Lmv( -h * this->b[i], ( i == this->num_stages - 1 ? 0.0 : 1.0 ), sigma_t,
                                *this->c_stage_flux[i], c_result );
    }

    this->phi->ZeroDensity( OpDomain::Boundary );
    TransportOperator::Smv( 1.0, 1.0, sigma_s, *this->phi, c_result );

    // Include collided initial condition if necessary (collided source is always zero).
    if ( c_initial_condition != nullptr )
        OrdinateFlux::AXPY( 1.0, 1.0, *c_initial_condition, c_result );

    Global::TMR_collided.Stop();
}


//------------------------------------------------------------------------------------------------------------
//! \brief  Computes reconstructed approximations of the intermediate stage values of a set of hybrid
//!         solutions using a Nyström interpolation technique.
//!
//! The reconstructed approximations are stored in the elements of <code> this->u_stage_flux </code>.
//!
//! This routine is only called by HybridIc methods.
//!
//! \param[in,out]  u_initial_condition Contains the initial condition \f$ \Psi^n \f$ for the uncollided flux
//!                                     in the current timestep. The initial condition for the collided flux
//!                                     is assumed to be zero.
//! \param[in]      sigma_t             RKDG::CrossSection object containing the total cross section
//!                                     \f$ \sigma_{\mathrm{t}} \f$.
//! \param[in]      sigma_s             RKDG::CrossSection object containing the scattering cross section
//!                                     \f$ \sigma_{\mathrm{s}} \f$.
//! \param[in]      source              Source term of the transport system. <br>
//!                                     The source term is assumed to be constant across the timestep interval.
//! \param[in]      h                   Timestep size.
//!
//! \see    DIRKIntegrator::Step_HybridI()
//------------------------------------------------------------------------------------------------------------
void DIRKIntegrator::ComputeNystromInterpolants (

    const OrdinateFlux & u_initial_condition,
    const OrdinateFlux & source,
    const CrossSection & sigma_t,
    const CrossSection & sigma_s,
    const double h
) {

    double stage_coeff, source_coeff;

    for ( int64_t i = 0; i < this->num_stages; ++i ) {

        // Explicit stage.
        if ( IA(i,i) == 0.0 ) {

            stage_coeff = h;
            OrdinateFlux::Copy( *this->u_stage_flux[i], u_initial_condition );

        // Implicit stage.
        } else {

            stage_coeff = 1.0 / IA(i,i);
            OrdinateFlux::AXPY( 1.0 / (h * IA(i,i)), 0.0, u_initial_condition, *this->u_stage_flux[i] );
        }

        source_coeff = 0.0;

        // Compute weighted sum of transport operator applied to previous stages.
        for ( int64_t j = 0; j < i; ++j ) {

            if ( IA(i,j) == 0.0 ) {  continue;  }

            source_coeff += IA(i,j);

            TransportOperator::Pmv( stage_coeff * IA(i,j), 1.0, *this->u_stage_flux[j], this->phi[i] );
            TransportOperator::Lmv( -stage_coeff * IA(i,j), 1.0, sigma_t, *this->u_stage_flux[j],
                                    *this->u_stage_flux[i] );
        }

        this->phi[i].ZeroDensity( OpDomain::Boundary );
        TransportOperator::Smv( 1.0, 1.0, sigma_s, this->phi[i], *this->u_stage_flux[i] );

        // Include source from previous stages
        OrdinateFlux::AXPY( stage_coeff * source_coeff, 1.0, source, *this->u_stage_flux[i] );
        this->u_stage_flux[i]->ZeroDensity( OpDomain::Boundary );

        // If the stage is implicit, include source terms for current stage and do sweep.
        if ( IA(i,i) != 0.0 ) {

            OrdinateFlux::AXPY( 1.0, 1.0, source, *this->u_stage_flux[i] );
            this->u_solver->Solve( *this->u_stage_flux[i], sigma_t, sigma_s, h * IA(i,i) );

        } else {

            this->u_stage_flux[i]->Copy( source, OpDomain::Boundary );
        }
    }
}


//------------------------------------------------------------------------------------------------------------
//! \brief  Constructs the Butcher tableau for the DIRK method set by DIRKIntegrator::dirk_method.
//!
//! Sets DIRKIntegrator::num_stages and DIRKIntegrator::stiffly_accurate. Allocates memory for
//! DIRKIntegrator::A, DIRKIntegrator::b, and DIRKIntegrator::c, and initializes allocated memory with the
//! proper values from the Butcher tableau of the method.
//------------------------------------------------------------------------------------------------------------
void DIRKIntegrator::SetTableau ( void ) {

    switch( this->dirk_method ) {

        case DIRKIntegrator::DIRKMethod::Euler:
        {
            //
            // Euler
            //
            // Stages:      one
            // Order:       first
            // Stability:   L-stable
            // Stiffly Acc: yes
            //

            this->stiffly_accurate = true;
            this->num_stages = 1;
            this->stages_allocd = 1;

            this->A = new double[ this->num_stages * this->num_stages ];
            this->b = new double[ this->num_stages ];
            this->c = new double[ this->num_stages ];

            this->stage_alloc_map = new int64_t[ this->num_stages ];

            IA(0,0) = 1.0;

            this->b[0] = 1.0;

            this->c[0] = 1.0;

            this->stage_alloc_map[0] = 0;

        } break;

        case DIRKIntegrator::DIRKMethod::SDIRK2:
        {
            //
            // SDIRK2
            //
            // Stages:      two
            // Order:       second
            // Stability:   L-stable
            // Stiffly Acc: yes
            //

            this->stiffly_accurate = true;
            this->num_stages = 2;
            this->stages_allocd = 2;

            this->A = new double[ this->num_stages * this->num_stages ];
            this->b = new double[ this->num_stages ];
            this->c = new double[ this->num_stages ];

            this->stage_alloc_map = new int64_t[ this->num_stages ];

            const double alpha = 1.0 - sqrt(2.0) / 2.0;

            IA(0,0) = alpha;
            IA(0,1) = 0.0;
            IA(1,0) = 1.0 - alpha;
            IA(1,1) = alpha;

            this->b[0] = 1.0 - alpha;
            this->b[1] = alpha;

            this->c[0] = alpha;
            this->c[1] = 1.0;

            this->stage_alloc_map[0] = 0;
            this->stage_alloc_map[1] = 1;

        } break;

        case DIRKIntegrator::DIRKMethod::SDIRK3:
        {
            //
            // SDIRK3
            //
            // Stages:      three
            // Order:       third
            // Stability:   L-stable
            // Stiffly Acc: yes
            //

            this->stiffly_accurate = true;
            this->num_stages = 3;
            this->stages_allocd = 3;

            this->A = new double[ this->num_stages * this->num_stages ];
            this->b = new double[ this->num_stages ];
            this->c = new double[ this->num_stages ];

            this->stage_alloc_map = new int64_t[ this->num_stages ];

            const double alpha = 0.435866521508459;
            const double tau   = ( 1.0 + alpha ) / 2.0;
            const double b1    = - ( 6.0 * alpha*alpha - 16.0 * alpha + 1.0 ) / 4.0;
            const double b2    = ( 6.0 * alpha*alpha - 20.0 * alpha + 5.0 ) / 4.0;

            IA(0,0) = alpha;
            IA(0,1) = 0.0;
            IA(0,2) = 0.0;

            IA(1,0) = tau - alpha;
            IA(1,1) = alpha;
            IA(1,2) = 0.0;

            IA(2,0) = b1;
            IA(2,1) = b2;
            IA(2,2) = alpha;

            this->b[0] = b1;
            this->b[1] = b2;
            this->b[2] = alpha;

            this->c[0] = alpha;
            this->c[1] = tau;
            this->c[2] = 1.0;

            this->stage_alloc_map[0] = 0;
            this->stage_alloc_map[1] = 1;
            this->stage_alloc_map[2] = 2;

        } break;

        case DIRKIntegrator::DIRKMethod::SDIRK5:
        {
            //
            // SDIRK5
            //
            // Stages:      five
            // Order:       fifth
            // Stability:   L-stable
            // Stiffly Acc: yes
            //

            this->stiffly_accurate = true;
            this->num_stages = 5;
            this->stages_allocd = 4;

            this->A = new double[ this->num_stages * this->num_stages ];
            this->b = new double[ this->num_stages ];
            this->c = new double[ this->num_stages ];

            this->stage_alloc_map = new int64_t[ this->num_stages ];

            IA(0,0) = 0.25;
            IA(0,1) = 0.0;
            IA(0,2) = 0.0;
            IA(0,3) = 0.0;
            IA(0,4) = 0.0;

            IA(1,0) = -1.0 / 12.0;
            IA(1,1) = 0.25;
            IA(1,2) = 0.0;
            IA(1,3) = 0.0;
            IA(1,4) = 0.0;

            IA(2,0) = ( 73.0 + 12.0 * sqrt(41.0) ) / 150.0;
            IA(2,1) = ( 24.0 - 19.0 * sqrt(41.0) ) / 300.0;
            IA(2,2) = 0.25;
            IA(2,3) = 0.0;
            IA(2,4) = 0.0;

            IA(3,0) = 0.16336450515843079546;   // ( 59765462671.0 - 2469071899.0 * sqrt(41.0) ) / 269065110000.0;
            IA(3,1) = 0.31739121890801644520;   // ( 26775007261.0 + 244199891.0 * sqrt(41.0) ) / 89286180000.0;
            IA(3,2) = -0.020807794690328052104; // ( 889326089143.0 - 203592224167.0 * sqrt(41.0) ) / 19910818140000.0;
            IA(3,3) = 0.25;
            IA(3,4) = 0.0;

            IA(4,0) = 0.0;
            IA(4,1) = 15.0 / 37.0;
            IA(4,2) = ( 2091.0 - 879.0 * sqrt(41.0) ) / 12136.0;
            IA(4,3) = ( 2091.0 + 879.0 * sqrt(41.0) ) / 12136.0;
            IA(4,4) = 0.25;

            this->b[0] = 0.0;
            this->b[1] = 15.0 / 37.0;
            this->b[2] = ( 2091.0 - 879.0 * sqrt(41.0) ) / 12136.0;
            this->b[3] = ( 2091.0 + 879.0 * sqrt(41.0) ) / 12136.0;
            this->b[4] = 0.25;

            this->c[0] = 0.25;
            this->c[1] = 1.0 / 6.0;
            this->c[2] = ( 49.0 + sqrt(41.0) ) / 60.0;
            this->c[3] = ( 49.0 + sqrt(41.0) ) / 60.0;
            this->c[4] = 1.0;

            this->stage_alloc_map[0] = 0;
            this->stage_alloc_map[1] = 1;
            this->stage_alloc_map[2] = 2;
            this->stage_alloc_map[3] = 3;
            this->stage_alloc_map[4] = 0;

        } break;

        case DIRKIntegrator::DIRKMethod::SDIRK_5_3_4:
        {
            //
            // SDIRK-5-3-4
            //
            // Stages:      five
            // Order:       fourth
            // Stability:   L-stable
            // Stiffly Acc: yes
            //

            this->stiffly_accurate = true;
            this->num_stages = 5;
            this->stages_allocd = 5;

            this->A = new double[ this->num_stages * this->num_stages ];
            this->b = new double[ this->num_stages ];
            this->c = new double[ this->num_stages ];

            this->stage_alloc_map = new int64_t[ this->num_stages ];

            IA(0,0) = 0.25;
            IA(0,1) = 0.0;
            IA(0,2) = 0.0;
            IA(0,3) = 0.0;
            IA(0,4) = 0.0;

            IA(1,0) = 0.5;
            IA(1,1) = 0.25;
            IA(1,2) = 0.0;
            IA(1,3) = 0.0;
            IA(1,4) = 0.0;

            IA(2,0) = 17.0 / 50.0;
            IA(2,1) = - 1.0 / 25.0;
            IA(2,2) = 0.25;
            IA(2,3) = 0.0;
            IA(2,4) = 0.0;

            IA(3,0) = 371.0 / 1360.0;
            IA(3,1) = - 137.0 / 2720.0;
            IA(3,2) = 15.0 / 544.0;
            IA(3,3) = 0.25;
            IA(3,4) = 0.0;

            IA(4,0) = 25.0 / 24.0;
            IA(4,1) = - 49.0 / 48.0;
            IA(4,2) = 125.0 / 16.0;
            IA(4,3) = - 85.0 / 12.0;
            IA(4,4) = 0.25;

            this->b[0] = 25.0 / 24.0;
            this->b[1] = - 49.0 / 48.0;
            this->b[2] = 125.0 / 16.0;
            this->b[3] = - 85.0 / 12.0;
            this->b[4] = 0.25;

            this->c[0] = 0.25;
            this->c[1] = 0.75;
            this->c[2] = 11.0 / 20.0;
            this->c[3] = 0.5;
            this->c[4] = 1.0;

            this->stage_alloc_map[0] = 0;
            this->stage_alloc_map[1] = 1;
            this->stage_alloc_map[2] = 2;
            this->stage_alloc_map[3] = 3;
            this->stage_alloc_map[4] = 4;

        } break;

        case DIRKIntegrator::DIRKMethod::KVAERNO_4_2_3:
        {
            //
            // KVAERNO-4-2-3
            //
            // Stages:      four
            // Order:       third
            // Stability:   L-stable
            // Stiffly Acc: yes
            //

            this->stiffly_accurate = true;
            this->num_stages = 4;
            this->stages_allocd = 4;

            this->A = new double[ this->num_stages * this->num_stages ];
            this->b = new double[ this->num_stages ];
            this->c = new double[ this->num_stages ];

            this->stage_alloc_map = new int64_t[ this->num_stages ];

            IA(0,0) = 0.0;
            IA(0,1) = 0.0;
            IA(0,2) = 0.0;
            IA(0,3) = 0.0;

            IA(1,0) = 0.4358665215;
            IA(1,1) = 0.4358665215;
            IA(1,2) = 0.0;
            IA(1,3) = 0.0;

            IA(2,0) = 0.490563388419108;
            IA(2,1) = 0.073570090080892;
            IA(2,2) = 0.4358665215;
            IA(2,3) = 0.0;

            IA(3,0) = 0.308809969973036;
            IA(3,1) = 1.490563388254106;
            IA(3,2) = -1.235239879727145;
            IA(3,3) = 0.4358665215;

            this->b[0] = 0.308809969973036;
            this->b[1] = 1.490563388254106;
            this->b[2] = -1.235239879727145;
            this->b[3] = 0.4358665215;

            this->c[0] = 0.0;
            this->c[1] = 0.871733043;
            this->c[2] = 1.0;
            this->c[3] = 1.0;

            this->stage_alloc_map[0] = 0;
            this->stage_alloc_map[1] = 1;
            this->stage_alloc_map[2] = 2;
            this->stage_alloc_map[3] = 3;

        } break;

        case DIRKIntegrator::DIRKMethod::KVAERNO_7_4_5:
        {
            //
            // KVAERNO-7-4-5
            //
            // Stages:      seven
            // Order:       fifth
            // Stability:   L-stable
            // Stiffly Acc: yes
            //

            this->stiffly_accurate = true;
            this->num_stages = 7;
            this->stages_allocd = 6;

            this->A = new double[ this->num_stages * this->num_stages ];
            this->b = new double[ this->num_stages ];
            this->c = new double[ this->num_stages ];

            this->stage_alloc_map = new int64_t[ this->num_stages ];

            IA(0,0) = 0.0;
            IA(0,1) = 0.0;
            IA(0,2) = 0.0;
            IA(0,3) = 0.0;
            IA(0,4) = 0.0;
            IA(0,5) = 0.0;
            IA(0,6) = 0.0;

            IA(1,0) = 0.26;
            IA(1,1) = 0.26;
            IA(1,2) = 0.0;
            IA(1,3) = 0.0;
            IA(1,4) = 0.0;
            IA(1,5) = 0.0;
            IA(1,6) = 0.0;

            IA(2,0) = 0.13;
            IA(2,1) = 0.84033320996790809;
            IA(2,2) = 0.26;
            IA(2,3) = 0.0;
            IA(2,4) = 0.0;
            IA(2,5) = 0.0;
            IA(2,6) = 0.0;

            IA(3,0) = 0.22371961478320505;
            IA(3,1) = 0.47675532319799699;
            IA(3,2) = -0.06470895363112615;
            IA(3,3) = 0.26;
            IA(3,4) = 0.0;
            IA(3,5) = 0.0;
            IA(3,6) = 0.0;

            IA(4,0) = 0.16648564323248321;
            IA(4,1) = 0.10450018841591720;
            IA(4,2) = 0.03631482272098715;
            IA(4,3) = -0.13090704451073998;
            IA(4,4) = 0.26;
            IA(4,5) = 0.0;
            IA(4,6) = 0.0;

            IA(5,0) = 0.13855640231268224;
            IA(5,1) = 0.0;
            IA(5,2) = -0.04245337201752043;
            IA(5,3) = 0.02446657898003141;
            IA(5,4) = 0.61943039072480676;
            IA(5,5) = 0.26;
            IA(5,6) = 0.0;

            IA(6,0) = 0.13659751177640291;
            IA(6,1) = 0.0;
            IA(6,2) = -0.05496908796538376;
            IA(6,3) = -0.04118626728321046;
            IA(6,4) = 0.62993304899016403;
            IA(6,5) = 0.06962479448202728;
            IA(6,6) = 0.26;

            this->b[0] = 0.13659751177640291;
            this->b[1] = 0.0;
            this->b[2] = -0.05496908796538376;
            this->b[3] = -0.04118626728321046;
            this->b[4] = 0.62993304899016403;
            this->b[5] = 0.06962479448202728;
            this->b[6] = 0.26;

            this->c[0] = 0.0;
            this->c[1] = 0.52;
            this->c[2] = 1.230333209967908;
            this->c[3] = 0.895765984350076;
            this->c[4] = 0.436393609858648;
            this->c[5] = 1.0;
            this->c[6] = 1.0;

            this->stage_alloc_map[0] = 0;
            this->stage_alloc_map[1] = 1;
            this->stage_alloc_map[2] = 2;
            this->stage_alloc_map[3] = 3;
            this->stage_alloc_map[4] = 4;
            this->stage_alloc_map[5] = 5;
            this->stage_alloc_map[6] = 1;

        } break;

        case DIRKIntegrator::DIRKMethod::ARK_8_4_5:
        {
            //
            // ARK-8-4-5
            //
            // Stages:      eight
            // Order:       fifth
            // Stability:   L-stable
            // Stiffly Acc: yes
            //

            this->stiffly_accurate = true;
            this->num_stages = 8;
            this->stages_allocd = 6;

            this->A = new double[ this->num_stages * this->num_stages ];
            this->b = new double[ this->num_stages ];
            this->c = new double[ this->num_stages ];

            this->stage_alloc_map = new int64_t[ this->num_stages ];

            IA(0,0) = 0.0;
            IA(0,1) = 0.0;
            IA(0,2) = 0.0;
            IA(0,3) = 0.0;
            IA(0,4) = 0.0;
            IA(0,5) = 0.0;
            IA(0,6) = 0.0;
            IA(0,7) = 0.0;

            IA(1,0) = 41.0 / 200.0;
            IA(1,1) = 41.0 / 200.0;
            IA(1,2) = 0.0;
            IA(1,3) = 0.0;
            IA(1,4) = 0.0;
            IA(1,5) = 0.0;
            IA(1,6) = 0.0;
            IA(1,7) = 0.0;

            IA(2,0) = 41.0 / 400.0;
            IA(2,1) = -567603406766.0 / 11931857230679.0;
            IA(2,2) = 41.0 / 200.0;
            IA(2,3) = 0.0;
            IA(2,4) = 0.0;
            IA(2,5) = 0.0;
            IA(2,6) = 0.0;
            IA(2,7) = 0.0;

            IA(3,0) = 683785636431.0 / 9252920307686.0;
            IA(3,1) = 0.0;
            IA(3,2) = -110385047103.0 / 1367015193373.0;
            IA(3,3) = 41.0 / 200.0;
            IA(3,4) = 0.0;
            IA(3,5) = 0.0;
            IA(3,6) = 0.0;
            IA(3,7) = 0.0;

            IA(4,0) = 3016520224154.0 / 10081342136671.0;
            IA(4,1) = 0.0;
            IA(4,2) = 30586259806659.0 / 12414158314087.0;
            IA(4,3) = -22760509404356.0 / 11113319521817.0;
            IA(4,4) = 41.0 / 200.0;
            IA(4,5) = 0.0;
            IA(4,6) = 0.0;
            IA(4,7) = 0.0;

            IA(5,0) = 218866479029.0 / 1489978393911.0;
            IA(5,1) = 0.0;
            IA(5,2) = 638256894668.0 / 5436446318841.0;
            IA(5,3) = -1179710474555.0 / 5321154724896.0;
            IA(5,4) = -60928119172.0 / 8023461067671.0;
            IA(5,5) = 41.0 / 200.0;
            IA(5,6) = 0.0;
            IA(5,7) = 0.0;

            IA(6,0) = 1020004230633.0 / 5715676835656.0;
            IA(6,1) = 0.0;
            IA(6,2) = 25762820946817.0 / 25263940353407.0;
            IA(6,3) = -2161375909145.0 / 9755907335909.0;
            IA(6,4) = -211217309593.0 / 5846859502534.0;
            IA(6,5) = -4269925059573.0 / 7827059040749.0;
            IA(6,6) = 41.0 / 200.0;
            IA(6,7) = 0.0;

            IA(7,0) = -872700587467.0 / 9133579230613.0;
            IA(7,1) = 0.0;
            IA(7,2) = 0.0;
            IA(7,3) = 22348218063261.0 / 9555858737531.0;
            IA(7,4) = -1143369518992.0 / 8141816002931.0;
            IA(7,5) = -39379526789629.0 / 19018526304540.0;
            IA(7,6) = 32727382324388.0 / 42900044865799.0;
            IA(7,7) = 41.0 / 200.0;

            this->b[0] = -872700587467.0 / 9133579230613.0;
            this->b[1] = 0.0;
            this->b[2] = 0.0;
            this->b[3] = 22348218063261.0 / 9555858737531.0;
            this->b[4] = -1143369518992.0 / 8141816002931.0;
            this->b[5] = -39379526789629.0 / 19018526304540.0;
            this->b[6] = 32727382324388.0 / 42900044865799.0;
            this->b[7] = 41.0 / 200.0;

            this->c[0] = 0.0;
            this->c[1] = 41.0 / 100.0;
            this->c[2] = 2935347310677.0 / 11292855782101.0;
            this->c[3] = 1426016391358.0 / 7196633302097.0;
            this->c[4] = 92.0 / 100.0;
            this->c[5] = 24.0 / 100.0;
            this->c[6] = 3.0 / 5.0;
            this->c[7] = 1.0;

            this->stage_alloc_map[0] = 0;
            this->stage_alloc_map[1] = 1;
            this->stage_alloc_map[2] = 2;
            this->stage_alloc_map[3] = 3;
            this->stage_alloc_map[4] = 4;
            this->stage_alloc_map[5] = 5;
            this->stage_alloc_map[6] = 1;
            this->stage_alloc_map[7] = 2;

        } break;

        case DIRKIntegrator::DIRKMethod::GSBP_DIRK_3:
        {
            //
            // GSBP-DIRK-3
            //
            // Stages:      three
            // Order:       third
            // Stability:   L-stable
            // Stiffly Acc: no
            //

            this->num_stages = 3;
            this->stiffly_accurate = false;

            this->A = new double[ this->num_stages * this->num_stages ];
            this->b = new double[ this->num_stages ];
            this->c = new double[ this->num_stages ];

            this->stage_alloc_map = new int64_t[ this->num_stages ];

            IA(0,0) = 0.0585104413426586;
            IA(0,1) = 0.0;
            IA(0,2) = 0.0;

            IA(1,0) = 0.0389225469556698;
            IA(1,1) = 0.7675348853239251;
            IA(1,2) = 0.0;

            IA(2,0) = 0.1613387070350185;
            IA(2,1) = -0.5944302919004032;
            IA(2,2) = 0.7165457925008468;

            this->b[0] = 0.1008717264855379;
            this->b[1] = 0.4574278841698629;
            this->b[2] = 0.4417003893445992;

            this->c[0] = 0.0585104413419415;
            this->c[1] = 0.8064574322792799;
            this->c[2] = 0.2834542075672883;

            this->stage_alloc_map[0] = 0;
            this->stage_alloc_map[1] = 1;
            this->stage_alloc_map[2] = 2;

        } break;

        case DIRKIntegrator::DIRKMethod::GSBP_DIRK_4:
        {
            //
            // GSBP-DIRK-4
            //
            // Stages:      four
            // Order:       fourth
            // Stability:   L-stable
            // Stiffly Acc: no
            //

            this->num_stages = 4;
            this->stiffly_accurate = false;

            this->A = new double[ this->num_stages * this->num_stages ];
            this->b = new double[ this->num_stages ];
            this->c = new double[ this->num_stages ];

            this->stage_alloc_map = new int64_t[ this->num_stages ];

            IA(0,0) = 0.5975501145870646;
            IA(0,1) = 0.0;
            IA(0,2) = 0.0;
            IA(0,3) = 0.0;

            IA(1,0) = -0.3662683378362842;
            IA(1,1) = 0.4899631271029300;
            IA(1,2) = 0.0;
            IA(1,3) = 0.0;

            IA(2,0) = -0.9122346095222909;
            IA(2,1) = 1.395636663278596;
            IA(2,2) = 0.4979628247281717;
            IA(2,3) = 0.0;

            IA(3,0) = 4.870201094711127;
            IA(3,1) = -3.007233691002447;
            IA(3,2) = -2.425297972138512;
            IA(3,3) = 0.7811652842149162;

            this->b[0] = 0.5263633266867775;
            this->b[1] = 0.3002573924935185;
            this->b[2] = 0.1447678514141155;
            this->b[3] = 0.02861142940558849;

            this->c[0] = 0.5975501145870646;
            this->c[1] = 0.1236947892666459;
            this->c[2] = 0.9813648784844768;
            this->c[3] = 0.2188347157850838;

            this->stage_alloc_map[0] = 0;
            this->stage_alloc_map[1] = 1;
            this->stage_alloc_map[2] = 2;
            this->stage_alloc_map[3] = 3;

        } break;

        default:
        {   std::string error_message = "Invalid DIRKMethod '"
                                        + DIRKMethod_to_String.at( this->dirk_method )
                                        + "' for initializing DIRK integrator this in '"
                                        + std::string(__func__) + "'.\n";

            PRINT_ERROR( error_message.c_str() )
            throw std::invalid_argument( error_message );
        }
    }
}
