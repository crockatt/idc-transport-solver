//------------------------------------------------------------------------------------------------------------
//! \file   time_integrators/RKDG/idc.cpp
//! \brief  Implementation of IDC integrators using implicit Euler steps.
//!
//! \author Michael M. Crockatt
//! \date   February 2018
//------------------------------------------------------------------------------------------------------------

# include <cinttypes>
# include <cmath>
# include <cstdint>
# include <cstdlib>
# include <cstring>

# include "operators/TransportOperator.hpp"
# include "time_integrators/RKDG/idc.hpp"
# include "utils/CLog.hpp"
# include "utils/global.hpp"
# include "utils/Quadrule/Quadrule.hpp"


using namespace RKDG;
using namespace Quadrule;


//============================================================================================================
//=== STATIC DEFINITIONS =====================================================================================
//============================================================================================================


//------------------------------------------------------------------------------------------------------------
//! \brief  Used to convert IDCIntegrator::IDCType values to their corresponding string descriptors.
//------------------------------------------------------------------------------------------------------------
std::map< IDCIntegrator::IDCType, std::string > IDCIntegrator::IDCType_to_String = {

    { IDCIntegrator::IDCType::None,     "none"      },
    { IDCIntegrator::IDCType::Error,    "error"     },
    { IDCIntegrator::IDCType::Update,   "update"    }
};


//------------------------------------------------------------------------------------------------------------
//! \brief  Used to convert string descriptors to IDCIntegrator::IDCType values.
//------------------------------------------------------------------------------------------------------------
std::map< std::string, IDCIntegrator::IDCType > IDCIntegrator::String_to_IDCType = {

    { "none",       IDCIntegrator::IDCType::None        },
    { "error",      IDCIntegrator::IDCType::Error       },
    { "update",     IDCIntegrator::IDCType::Update      }
};


//------------------------------------------------------------------------------------------------------------
//! \brief  Used to convert IDCIntegrator::HybridMethod values to their corresponding string descriptors.
//------------------------------------------------------------------------------------------------------------
std::map< IDCIntegrator::HybridMethod, std::string > IDCIntegrator::HybridMethod_to_String = {

    { IDCIntegrator::HybridMethod::None,        "none"          },
    { IDCIntegrator::HybridMethod::HybridIa,    "hybrid-ia"     },
    { IDCIntegrator::HybridMethod::HybridIb,    "hybrid-ib"     },
    { IDCIntegrator::HybridMethod::HybridIIa,   "hybrid-iia"    },
    { IDCIntegrator::HybridMethod::HybridIIb,   "hybrid-iib"    },
    { IDCIntegrator::HybridMethod::HybridIIc,   "hybrid-iic"    }
};


//------------------------------------------------------------------------------------------------------------
//! \brief  Used to convert string descriptors to IDCIntegrator::HybridMethod values.
//------------------------------------------------------------------------------------------------------------
std::map< std::string, IDCIntegrator::HybridMethod > IDCIntegrator::String_to_HybridMethod = {

    { "none",           IDCIntegrator::HybridMethod::None       },
    { "hybrid-ia",      IDCIntegrator::HybridMethod::HybridIa   },
    { "hybrid-ib",      IDCIntegrator::HybridMethod::HybridIb   },
    { "hybrid-iia",     IDCIntegrator::HybridMethod::HybridIIa  },
    { "hybrid-iib",     IDCIntegrator::HybridMethod::HybridIIb  },
    { "hybrid-iic",     IDCIntegrator::HybridMethod::HybridIIc  }
};


//============================================================================================================
//=== CONSTRUCTORS AND DESTRUCTOR ============================================================================
//============================================================================================================


//------------------------------------------------------------------------------------------------------------
//! \brief  Creates an IDCIntegrator object using specified parameters and returns a pointer to the object.
//!
//! This is a wrapper function to ensure that IDCIntegrator::Allocate is called immediately after the object
//! is constructed.
//!
//! \param[in]  input_list      List of parameters to use for initialization.
//------------------------------------------------------------------------------------------------------------
IDCIntegrator * IDCIntegrator::Create (

    const ParameterList & input_list
) {

    IDCIntegrator * result = new IDCIntegrator( input_list );
    result->Allocate();

    return result;
}


//------------------------------------------------------------------------------------------------------------
//! \brief  Constructs an IDCIntegrator object using the provided values.
//!
//! \param[in]  input_list      List of parameters to use for initialization.
//------------------------------------------------------------------------------------------------------------
IDCIntegrator::IDCIntegrator (

    const ParameterList & input_list
) :
    TimeIntegrator( input_list ),

    nodes_type{},
    nodes{ nullptr },
    weights{ nullptr },

    left_endpoint{ false },
    right_endpoint{ false },
    use_collocation{ true },

    num_stages{},
    num_corrections{},

    phi{ nullptr },

    stages{ nullptr },      residuals{ nullptr },       error{ nullptr },       solver{ nullptr },

    hybrid_method{},

    u_stages{ nullptr },    u_residuals{ nullptr },     u_error{ nullptr },     u_solver{ nullptr },
    c_stages{ nullptr },    c_residuals{ nullptr },     c_error{ nullptr },     c_solver{ nullptr },

    relabel_operator{ nullptr }
{

    bool symmetric_reduce = false;

    ParameterList   params = input_list;
    ParameterList u_params = MakeUncollidedList( input_list );
    ParameterList c_params = MakeCollidedList( input_list );

    // Read polynomial degree for DG approximation.
    DG_degree = GetDGDegreeX( params );

    // Determine type of hybrid splitting (default to none).
    try         {  params.GetValue( "hybrid_method", this->hybrid_method,
                                                String_to_HybridMethod );              }
    catch (...) {  this->hybrid_method = HybridMethod::None;                           }

# if SPACE_DIMS == 2

    try         {  params.GetValue( "ordinate_sym_reduce", symmetric_reduce );  }
    catch (...) {  symmetric_reduce = false;                                                }

# endif // if SPACE_DIMS == 2

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
        this->u_ordinate_set.Reconfigure( u_params.GetValue<int64_t>( "ang_order" ),
                                          symmetric_reduce,
                                          u_ordinate_type );

        this->c_ordinate_set.Reconfigure( c_params.GetValue<int64_t>( "ang_order" ),
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
        this->ordinate_set.Reconfigure( params.GetValue<int64_t>( "ang_order" ),
                                        symmetric_reduce,
                                        ordinate_type );
    }


    // --- Setup IDC method. -------------------------------------------------------------------------- //

    params.GetValue( "idc_type", this->idc_type, IDCIntegrator::String_to_IDCType );
    params.GetValue( "idc_num_corr", this->num_corrections );

    SetNodes();

    if ( this->right_endpoint ) {

        try         {  params.GetValue( "idc_use_collocation", this->use_collocation );  }
        catch (...) {  this->use_collocation = false;                                                }
    }

    switch ( this->hybrid_method ) {

        case HybridMethod::HybridIa:
        {
            try         {  params.GetValue( "relabel", this->relabel );  }
            catch (...) {  this->relabel = true;                                     }

        } /* fall through */
        case HybridMethod::HybridIb:
        case HybridMethod::HybridIIa:
        case HybridMethod::HybridIIb:
        case HybridMethod::HybridIIc:
        break;

        case HybridMethod::None:
        default:
        break;
    }

    // --- Setup solvers used by IDC integrator. ------------------------------------------------------ //

    // Hybrid/non-hybrid setup.
    switch ( this->hybrid_method ) {

        case HybridMethod::None:
        {
            this->solver = std::shared_ptr<ImplicitSolver<OrdinateFlux>>(
                    ImplicitSolver<OrdinateFlux>::ImplicitSolverFactory::CreateSolver(
                        *this, this->ordinate_set, params
                ) );
        } break;

        case HybridMethod::HybridIa:
        case HybridMethod::HybridIb:
        case HybridMethod::HybridIIa:
        case HybridMethod::HybridIIb:
        {
            // Relabel operator.
            this->relabel_operator = new RelabelOperator( this->c_ordinate_set, this->u_ordinate_set );

        } /* fall through */
        case HybridMethod::HybridIIc:
        {
            // Uncollided solver.
            this->u_solver = std::shared_ptr<ImplicitSolver<OrdinateFlux>>(
                    ImplicitSolver<OrdinateFlux>::ImplicitSolverFactory::CreateSolver(
                        *this, this->u_ordinate_set, u_params
                ) );

            // Collided solver.
            this->c_solver = std::shared_ptr<ImplicitSolver<OrdinateFlux>>(
                    ImplicitSolver<OrdinateFlux>::ImplicitSolverFactory::CreateSolver(
                        *this, this->c_ordinate_set, c_params
                ) );
        } break;

        default:
        {   std::string error_message = "Invalid hybrid splitting '"
                                        + HybridMethod_to_String.at( this->hybrid_method )
                                        + "' for initializing IDC integrator in '" + std::string(__func__)
                                        + "'.\n";

            PRINT_ERROR( error_message.c_str() )
            throw std::invalid_argument( error_message );
        }
    }
}


//------------------------------------------------------------------------------------------------------------
//! \brief  Allocates memory for internal stage values used by the IDCIntegrator object.
//------------------------------------------------------------------------------------------------------------
void IDCIntegrator::Allocate( void ) {

    this->phi = new ScalarFlux( *this, DG_degree );

    // Hybrid/non-hybrid setup.
    switch ( this->hybrid_method ) {

        case HybridMethod::None:
        {
            if ( idc_type == IDCType::Error )
                this->error = new OrdinateFlux( *this, this->ordinate_set, DG_degree );

            this->stages    = new OrdinateFlux*[ this->num_stages ];
            this->residuals = new OrdinateFlux*[ this->num_stages ];

            this->stages[0] = nullptr;

            if ( this->left_endpoint )
                this->residuals[0] = new OrdinateFlux( *this, this->ordinate_set, DG_degree );
            else
                this->residuals[0] = nullptr;

            for ( int64_t i = 1; i < this->num_stages; ++i ) {

                this->stages[i]    = new OrdinateFlux( *this, this->ordinate_set, DG_degree );
                this->residuals[i] = new OrdinateFlux( *this, this->ordinate_set, DG_degree );
            }
        } break;

        case HybridMethod::HybridIa:
        case HybridMethod::HybridIb:
        {
            if ( idc_type == IDCType::Error ) {

                this->u_error = new OrdinateFlux( *this, this->u_ordinate_set, DG_degree );
                this->c_error = new OrdinateFlux( *this, this->c_ordinate_set, DG_degree );
            }

            this->u_stages    = new OrdinateFlux*[ this->num_stages ];
            this->c_stages    = new OrdinateFlux*[ this->num_stages ];
            this->u_residuals = new OrdinateFlux*[ this->num_stages ];
            this->c_residuals = new OrdinateFlux*[ this->num_stages ];

            this->u_stages[0] = nullptr;
            this->c_stages[0] = nullptr;

            if ( this->left_endpoint ) {

                this->u_residuals[0] = new OrdinateFlux( *this, this->u_ordinate_set, DG_degree );
                this->c_residuals[0] = new OrdinateFlux( *this, this->c_ordinate_set, DG_degree );

            } else {

                this->u_residuals[0] = nullptr;
                this->c_residuals[0] = nullptr;
            }

            for ( int64_t i = 1; i < this->num_stages; ++i ) {

                this->u_stages[i] = new OrdinateFlux( *this, this->u_ordinate_set, DG_degree );
                this->c_stages[i] = new OrdinateFlux( *this, this->c_ordinate_set, DG_degree );

                this->u_residuals[i] = new OrdinateFlux( *this, this->u_ordinate_set, DG_degree );
                this->c_residuals[i] = new OrdinateFlux( *this, this->c_ordinate_set, DG_degree );
            }
        } break;

        case HybridMethod::HybridIIa:
        case HybridMethod::HybridIIb:
        case HybridMethod::HybridIIc:
        {
            if (    this->hybrid_method == HybridMethod::HybridIIc
                 || this->idc_type == IDCType::Error
            )
                this->u_error = new OrdinateFlux( *this, this->u_ordinate_set, DG_degree );

            this->c_error = new OrdinateFlux( *this, this->c_ordinate_set, DG_degree );

            this->u_stages    = new OrdinateFlux*[ this->num_stages ];
            this->u_residuals = new OrdinateFlux*[ this->num_stages ];

            this->stages    = this->u_stages;       // Pointer aliasing for ComputeResiduals and
            this->residuals = this->u_residuals;    // ComputeCollocation_Nonhybrid.

            this->u_stages[0] = nullptr;

            if ( this->left_endpoint )
                this->u_residuals[0] = new OrdinateFlux( *this, this->u_ordinate_set, DG_degree );
            else
                this->u_residuals[0] = nullptr;

            if (    this->hybrid_method != HybridMethod::HybridIIc
                 && this->idc_type == IDCType::Update
            ) {
                this->c_stages = new OrdinateFlux*[ this->num_stages ];
                this->c_stages[0] = nullptr;
            }

            for ( int64_t i = 1; i < this->num_stages; ++i ) {

                this->u_stages[i]    = new OrdinateFlux( *this, this->u_ordinate_set, DG_degree );
                this->u_residuals[i] = new OrdinateFlux( *this, this->u_ordinate_set, DG_degree );

                if (    this->hybrid_method != HybridMethod::HybridIIc
                     && this->idc_type == IDCType::Update
                )
                    this->c_stages[i] = new OrdinateFlux( *this, this->c_ordinate_set, DG_degree );
            }
        } break;

        default:
        {   std::string error_message = "Invalid hybrid splitting '"
                                        + HybridMethod_to_String.at( this->hybrid_method )
                                        + "' for initializing DIRK integrator in '" + std::string(__func__)
                                        + "'.\n";

            PRINT_ERROR( error_message.c_str() )
            throw std::invalid_argument( error_message );
        }
    }
}


//------------------------------------------------------------------------------------------------------------
//! \brief  Destructor for class IDCIntegrator.
//------------------------------------------------------------------------------------------------------------
IDCIntegrator::~IDCIntegrator( void ) {

    delete [] nodes;
    delete [] weights;

    delete error;
    delete u_error;
    delete c_error;

    if ( this->hybrid_method == HybridMethod::None ) {

        if ( this->stages != nullptr ) {

            for ( int64_t i = 0; i < this->num_stages; ++i )
                delete this->stages[i];

            delete [] this->stages;
        }

        if ( this->residuals != nullptr ) {

            for ( int64_t i = 0; i < this->num_stages; ++i )
                delete this->residuals[i];

            delete [] this->residuals;
        }

    } else {

        if ( this->u_stages != nullptr ) {

            for ( int64_t i = 0; i < this->num_stages; ++i )
                delete this->u_stages[i];

            delete [] this->u_stages;
        }

        if ( this->c_stages != nullptr ) {

            for ( int64_t i = 0; i < this->num_stages; ++i )
                delete this->c_stages[i];

            delete [] this->c_stages;
        }

        if ( this->u_residuals != nullptr ) {

            for ( int64_t i = 0; i < this->num_stages; ++i )
                delete this->u_residuals[i];

            delete [] this->u_residuals;
        }

        if ( this->c_residuals != nullptr ) {

            for ( int64_t i = 0; i < this->num_stages; ++i )
                delete this->c_residuals[i];

            delete [] this->c_residuals;
        }
    }

    delete this->relabel_operator;
}


//============================================================================================================
//=== INTERFACE ROUTINES =====================================================================================
//============================================================================================================


//------------------------------------------------------------------------------------------------------------
//! \brief  Calls the appropriate \c Step_* function corresponding to the desired IDC hybrid splitting as set
//!         by IDCIntegrator::hybrid_method.
//!
//! \param[in,out]  initial_condition   Initially contains the initial condition \f$ \Psi^n \f$ for the
//!                                     current timestep. Upon return contains the result \f$ \Psi^{n+1} \f$
//!                                     of the timestep update.
//! \param[in]      sigma_t_in          RKDG::CrossSection object containing the total cross section
//!                                     \f$ \sigma_{\mathrm{t}} \f$.
//! \param[in]      sigma_s_in          RKDG::CrossSection object containing the scattering cross section
//!                                     \f$ \sigma_{\mathrm{s}} \f$.
//! \param[in]      source              Source term of the transport system. <br>
//!                                     The source term is assumed to be constant across the timestep interval.
//! \param[in]      step_size           Timestep size.
//!
//! \see    IDCIntegrator::Step_Nonhybrid()
//! \see    IDCIntegrator::Step_HybridI()
//! \see    IDCIntegrator::Step_HybridII()
//------------------------------------------------------------------------------------------------------------
IDCIntegrator & IDCIntegrator::Step (

    OrdinateFlux & initial_condition,
    const OrdinateFlux & source,
    const CrossSection & sigma_t_in,
    const CrossSection & sigma_s_in,
    const double step_size
) {

    this->sigma_s = &sigma_s_in;
    this->sigma_t = &sigma_t_in;

# if LOGLEVEL >= 1
    PRINT_LOG( "\n" )
# endif

    switch ( this->hybrid_method ) {

        case HybridMethod::None:

            Step_Nonhybrid( initial_condition, source, step_size );
            break;

        case HybridMethod::HybridIa:
        case HybridMethod::HybridIb:

            Step_HybridI( initial_condition, nullptr, source, step_size );
            break;

        case HybridMethod::HybridIIa:
        case HybridMethod::HybridIIb:

            Step_HybridII( initial_condition, source, step_size );
            break;

        default:
        {   std::string error_message = "Invalid hybrid splitting '"
                                        + HybridMethod_to_String.at( this->hybrid_method )
                                        + "' for IDC step.\n";

            PRINT_ERROR( error_message.c_str() )
            throw std::invalid_argument( error_message );
        }
    }

# if LOGLEVEL >= 1
    PRINT_LOG( "\n" )
# endif

    this->sigma_s = nullptr;
    this->sigma_t = nullptr;

    return *this;
}


//------------------------------------------------------------------------------------------------------------
//! \brief  Calls the appropriate \c Step_* function corresponding to the desired IDC hybrid splitting as set
//!         by IDCIntegrator::hybrid_method.
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
//! \param[in]      sigma_t_in          RKDG::CrossSection object containing the total cross section
//!                                     \f$ \sigma_{\mathrm{t}} \f$.
//! \param[in]      sigma_s_in          RKDG::CrossSection object containing the scattering cross section
//!                                     \f$ \sigma_{\mathrm{s}} \f$.
//! \param[in]      source              Source term of the transport system. <br>
//!                                     The source term is assumed to be constant across the timestep interval.
//! \param[in]      step_size           Timestep size.
//!
//! \see    IDCIntegrator::Step_Nonhybrid()
//! \see    IDCIntegrator::Step_HybridI()
//! \see    IDCIntegrator::Step_HybridII()
//------------------------------------------------------------------------------------------------------------
IDCIntegrator & IDCIntegrator::Step (

    OrdinateFlux & u_initial_condition,
    OrdinateFlux * c_initial_condition,
    const OrdinateFlux & source,
    const CrossSection & sigma_t_in,
    const CrossSection & sigma_s_in,
    const double step_size
) {

    this->sigma_s = &sigma_s_in;
    this->sigma_t = &sigma_t_in;

# if LOGLEVEL >= 1
    PRINT_LOG( "\n" )
# endif

    switch ( this->hybrid_method ) {

        case HybridMethod::None:

        # if defined (STRICT_CHECK)

            if ( c_initial_condition != nullptr ) {

                PRINT_WARNING( "IDC integrator with hybrid splitting '%s' ignores pointer to collided "
                               "initial condition.\n",
                               HybridMethod_to_String.at( this->hybrid_method ).c_str() )
            }

        # endif // if defined (STRICT_CHECK)

            Step_Nonhybrid( u_initial_condition, source, step_size );
            break;

        case HybridMethod::HybridIa:
        case HybridMethod::HybridIb:

            Step_HybridI( u_initial_condition, c_initial_condition, source, step_size );
            break;

        case HybridMethod::HybridIIa:
        case HybridMethod::HybridIIb:
        case HybridMethod::HybridIIc:

        # if defined (STRICT_CHECK)

            if ( c_initial_condition != nullptr ) {

                PRINT_WARNING( "IDC integrator with hybrid splitting '%s' ignores pointer to collided "
                               "initial condition.\n",
                               HybridMethod_to_String.at( this->hybrid_method ).c_str() )
            }

        # endif // if defined (STRICT_CHECK)

            Step_HybridII( u_initial_condition, source, step_size );
            break;

        default:
        {   std::string error_message = "Invalid hybrid splitting '"
                                        + HybridMethod_to_String.at( this->hybrid_method )
                                        + "' for IDC step.\n";

            PRINT_ERROR( error_message.c_str() )
            throw std::invalid_argument( error_message );
        }
    }

# if LOGLEVEL >= 1
    PRINT_LOG( "\n" )
# endif

    this->sigma_s = nullptr;
    this->sigma_t = nullptr;

    return *this;
}


//------------------------------------------------------------------------------------------------------------
//! \brief  Prints a summary of the integrator configuration to the logging interface.
//!
//! \param[in]  prefix  String prepended to each line of output, intended to be used for controlling
//!                     indentation of output.
//------------------------------------------------------------------------------------------------------------
void IDCIntegrator::Print (

    const std::string & prefix // = "  "

) const {

    PRINT_LOG( "%s%-*s  %s\n", prefix.c_str(),
               Global::col_width, "Time integrator:",
               TimeIntegratorType_to_String.at( GetIntegratorType() ).c_str() )

    PRINT_LOG( "%s%-*s  %s\n", prefix.c_str(),
               Global::col_width, "IDC type:",
               IDCType_to_String.at( this->idc_type ).c_str() )

    PRINT_LOG( "%s%-*s  %s\n", prefix.c_str(),
               Global::col_width, "IDC nodes type:",
               NodesType_to_String.at( this->nodes_type ).c_str() )

    PRINT_LOG( "%s%-*s % " PRId64 "\n", prefix.c_str(),
               Global::col_width, "IDC stages:", this->num_stages )

    PRINT_LOG( "%s%-*s % " PRId64 "\n", prefix.c_str(),
               Global::col_width, "IDC corrections:", this->num_corrections )

    this->Abstract::TimeIntegrator::Print( prefix );

# if LOGLEVEL >= 2
    this->PrintQuadrature();
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
                  || this->hybrid_method == HybridMethod::HybridIb
                  || this->hybrid_method == HybridMethod::HybridIIa
                  || this->hybrid_method == HybridMethod::HybridIIb
                )
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
//! \brief  Returns true if the integrator has been constructed with hybrid splitting and false otherwise.
//------------------------------------------------------------------------------------------------------------
bool IDCIntegrator::IsHybrid ( void ) const {

    switch ( this->hybrid_method ) {

        case HybridMethod::HybridIa:
        case HybridMethod::HybridIb:
        case HybridMethod::HybridIIa:
        case HybridMethod::HybridIIb:
        case HybridMethod::HybridIIc:

            return true;

        default:
            return false;
    }
}


//------------------------------------------------------------------------------------------------------------
//! \brief  Prints a table of the quadrature nodes and weights used by the IDC algorithm.
//------------------------------------------------------------------------------------------------------------
void IDCIntegrator::PrintQuadrature ( void ) const {

    char buffer[100];
    double weightTotal = 0.0;

    // Table header.
    PRINT_LOG( "\n" )
    PRINT_LOG( " %-*s ", Global::col_width, "  Nodes" )

    for ( int k = 0; k < this->num_stages; ++k ) {

        sprintf( buffer, "  Weights %d", k );
        PRINT_LOG( " %-*s ", Global::col_width, buffer )
    }

    PRINT_LOG( "\n" )

    buffer[0] = '\0';
    strncat( buffer, "-------------------------------------------------", Global::col_width );

    for ( int k = 0; k <= this->num_stages; ++k )
        PRINT_LOG( " %-*s ", Global::col_width, buffer )

    // Node and weight values.
    for ( int i = 0; i < this->num_stages; ++i ) {

        PRINT_LOG( "\n" )
        PRINT_LOG( " % *.5e  ", Global::col_width-2, this->nodes[i] )

        for ( int j = 0; j < this->num_stages; ++j )  {

            PRINT_LOG( "  % *.5e  ", Global::col_width-2, Weight(j,i) )
            weightTotal += Weight(j,i);
        }
    }

    PRINT_LOG( "\n" )
    PRINT_LOG( "\n" )
    PRINT_LOG( " Weight total = % .5e\n", weightTotal )
    PRINT_LOG( " Deviation    = % .5e\n", fabs( 2.0 - weightTotal ) )
}


//============================================================================================================
//=== PRIVATE TIME INTEGRATION ROUTINES ======================================================================
//============================================================================================================


//------------------------------------------------------------------------------------------------------------
//! \brief  Computes a one timestep update without hybrid splitting.
//!
//! \param[in,out]  initial_condition   Initially contains the initial condition \f$ \Psi^n \f$ for the
//!                                     current timestep. Upon return contains the result \f$ \Psi^{n+1} \f$
//!                                     of the timestep update. Ghost cells are assumed to be zero on entry
//!                                     and are zeroed upon return.
//! \param[in]      source              Material source term and boundary conditions for the transport system.
//!                                     <br/>
//!                                     Assumed to be constant across the timestep interval.
//! \param[in]      dt                  Timestep size. Note that this is the full IDC timestep size from which
//!                                     the substep sizes are determined.
//!
//! \see    IDCIntegrator::Step_HybridI()
//! \see    IDCIntegrator::Step_HybridII()
//------------------------------------------------------------------------------------------------------------
void IDCIntegrator::Step_Nonhybrid (

    OrdinateFlux & initial_condition,
    const OrdinateFlux & source,
    const double dt
) {

    this->stages[0] = &initial_condition;

    if (    this->left_endpoint
         && BitmaskHasAny( Global::reflecting_boundaries, BoundaryEdge::All )
    ) {
        this->stages[0]->ReflectBoundaries();
    }

    // --- Prediction steps. -------------------------------------------------------------------------- //

    for ( int64_t n = 1; n < this->num_stages; ++n ) {

        // Normalized substep length.
        const double h = this->nodes[n] - this->nodes[n-1];

        // Compute source for implicit solve.
        OrdinateFlux::AXPY( 1.0/(h * dt), 0.0, *this->stages[n-1], *this->stages[n] );
        this->stages[n]->ZeroDensity( OpDomain::Boundary );
        OrdinateFlux::AXPY( 1.0, 1.0, source, *this->stages[n] );

        // Set initial guess for solve.
        this->solver->SetInitialGuess( this->stages[n-1] );

        // Perform implicit solve.
        this->solver->Solve( *this->stages[n], *this->sigma_t, *this->sigma_s, h * dt );
    }

    // --- Correction iterations. --------------------------------------------------------------------- //

    for ( int64_t p = 0; p < this->num_corrections; ++p ) {

        ComputeResiduals();

        for ( int64_t n = 1; n < this->num_stages; ++n ) {

            // Normalized substep length.
            const double h = (this->nodes[n] - this->nodes[n-1]);

            // Set initial guess for implicit solve.
            switch ( this->idc_type ) {

                case IDCType::Error:

                    this->solver->SetInitialGuess( nullptr );
                    break;

                case IDCType::Update:

                    this->solver->SetInitialGuess( this->stages[n] );
                    break;

                default:
                {   std::string error_message = "Invalid IDCType '" + IDCType_to_String.at( this->idc_type )
                                                + "' in '" + std::string(__func__) + "'.\n";

                    PRINT_ERROR( error_message.c_str() )
                    throw std::invalid_argument( error_message );
                }
            }

            // Compute source for implicit solve.
            ComputeSource( source, dt, n );

            // Perform implicit solve and apply correction as necessary.
            switch ( this->idc_type ) {

                case IDCType::Error:

                    this->solver->Solve( *this->error, *this->sigma_t, *this->sigma_s, h * dt );
                    OrdinateFlux::AXPY( 1.0, 1.0, *this->error, *this->stages[n] );
                    break;

                case IDCType::Update:

                    this->solver->Solve( *this->stages[n], *this->sigma_t, *this->sigma_s, h * dt );
                    break;

                default:
                {   std::string error_message = "Invalid IDCType '" + IDCType_to_String.at( this->idc_type )
                                                + "' in '" + std::string(__func__) + "'.\n";

                    PRINT_ERROR( error_message.c_str() )
                    throw std::invalid_argument( error_message );
                }
            }
        }
    }

    // --- Compute result of timestep. ---------------------------------------------------------------- //

    if (    !this->right_endpoint
         || this->use_collocation
    ) {
        ComputeCollocation_Nonhybrid( initial_condition, source, dt );
    }

    OrdinateFlux::swap( initial_condition, *this->stages[this->num_stages - 1] );
    initial_condition.ZeroDensity( OpDomain::Boundary );

    this->stages[0] = nullptr;
}


//------------------------------------------------------------------------------------------------------------
//! \brief  Computes a one timestep update using hybrid-I type splittings.
//!
//! \param[in,out]  u_initial_condition Initially contains the initial condition \f$ \Psi^n \f$ for the
//!                                     uncollided flux in the current timestep. Upon return contains the
//!                                     result \f$ \Psi^{n+1} \f$ of the timestep update if relabeling is used
//!                                     or only the uncollided component of the result if relabeling is not
//!                                     used. Ghost cells are assumed to be zero on entry and are zeroed upon
//!                                     return.
//! \param[in,out]  c_initial_condition Initially contains the initial condition \f$ \Psi^n \f$ for the
//!                                     collided flux in the current timestep. Upon return this value is
//!                                     unchanged if relabeling is used, or it contains the collided component
//!                                     of the result \f$ \Psi^{n+1} \f$ if relabeling is not used. Ghost
//!                                     cells are assumed to be zero on entry and are zeroed upon return.
//!                                     <br/>
//!                                     If this is \c null, then relabeling is forced.
//! \param[in]      source              Material source term and boundary conditions for the transport system.
//!                                     <br/>
//!                                     Assumed to be constant across the timestep interval.
//! \param[in]      dt                  Timestep size. Note that this is the full IDC timestep size from which
//!                                     the substep sizes are determined.
//!
//! \see    IDCIntegrator::Step_Nonhybrid()
//! \see    IDCIntegrator::Step_HybridII()
//------------------------------------------------------------------------------------------------------------
void IDCIntegrator::Step_HybridI (

    OrdinateFlux & u_initial_condition,
    OrdinateFlux * const c_initial_condition,
    const OrdinateFlux & source,
    const double dt
) {

    this->u_stages[0] = &u_initial_condition;
    this->c_stages[0] = c_initial_condition;

    if (    this->left_endpoint
         && BitmaskHasAny( Global::reflecting_boundaries, BoundaryEdge::All )
    ) {
        this->u_stages[0]->ReflectBoundaries();

        if ( this->c_stages[0] != nullptr )
            this->c_stages[0]->ReflectBoundaries();
    }

    // --- Prediction steps. -------------------------------------------------------------------------- //

    for ( int64_t n = 1; n < this->num_stages; ++n ) {

        // Normalized substep length.
        const double h = this->nodes[n] - this->nodes[n-1];

        // --- High-resolution (uncollided). --- //

        Global::TMR_uncollided.Start();

        // Compute source for implicit solve.
        OrdinateFlux::AXPY( 1.0/(h * dt), 0.0, *this->u_stages[n-1], *this->u_stages[n] );
        this->u_stages[n]->ZeroDensity( OpDomain::Boundary );
        OrdinateFlux::AXPY( 1.0, 1.0, source, *this->u_stages[n] );

        if (    this->hybrid_method == HybridMethod::HybridIb
             && this->c_stages[n-1] != nullptr
        ) {
            Global::TMR_uncollided.Stop();
            this->relabel_operator->Relabel( 1.0/(h * dt), *this->c_stages[n-1], *this->u_stages[n] );
            Global::TMR_uncollided.Start();
        }

        // Perform implicit solve.
        this->u_solver->Solve( *this->u_stages[n], *this->sigma_t, *this->sigma_s, h * dt );

        Global::TMR_uncollided.Stop();

        // --- Low-resolution (collided). --- //

        Global::TMR_collided.Start();

        // Compute source for implicit solve.
        TransportOperator::Pmv( 1.0, 0.0, *this->u_stages[n], *this->phi );
        this->phi->ZeroDensity( OpDomain::Boundary );
        TransportOperator::Smv( 1.0, 0.0, *this->sigma_s, *this->phi, *this->c_stages[n] );

        if (    this->hybrid_method == HybridMethod::HybridIa
             && this->c_stages[n-1] != nullptr
        ) {
            OrdinateFlux::AXPY( 1.0/(h * dt), 1.0, *this->c_stages[n-1], *this->c_stages[n] );
        }

        // Set initial guess for solve.
        this->c_solver->SetInitialGuess( this->c_stages[n-1] );

        // Perform implicit solve.
        this->c_solver->Solve( *this->c_stages[n], *this->sigma_t, *this->sigma_s, h * dt );

        Global::TMR_collided.Stop();
    }

    // --- Correction iterations. --------------------------------------------------------------------- //

    for ( int64_t p = 0; p < this->num_corrections; ++p ) {

        ComputeResiduals_HybridI();

        for ( int64_t n = 1; n < this->num_stages; ++n ) {

            // Normalized substep length.
            const double h = this->nodes[n] - this->nodes[n-1];

            // --- High-resolution solve (uncollided). --- //

            Global::TMR_uncollided.Start();

            // Compute source for implicit solve.
            u_ComputeSource_HybridI( source, dt, n );

            // Perform implicit solve and apply correction as necessary.
            switch ( this->idc_type ) {

                case IDCType::Error:

                    this->u_solver->Solve( *this->u_error, *this->sigma_t, *this->sigma_s, h * dt );
                    OrdinateFlux::AXPY( 1.0, 1.0, *this->u_error, *this->u_stages[n] );
                    break;

                case IDCType::Update:

                    this->u_solver->Solve( *this->u_stages[n], *this->sigma_t, *this->sigma_s, h * dt );
                    break;

                default:
                {   std::string error_message = "Invalid IDCType '" + IDCType_to_String.at( this->idc_type )
                                                + "' in '" + std::string(__func__) + "'.\n";

                    PRINT_ERROR( error_message.c_str() )
                    throw std::invalid_argument( error_message );
                }
            }

            Global::TMR_uncollided.Stop();

            // --- Low-resolution solve (collided). --- //

            Global::TMR_collided.Start();

            // Set initial guess for implicit solve.
            switch ( this->idc_type ) {

                case IDCType::Error:

                    this->c_solver->SetInitialGuess( nullptr );
                    break;

                case IDCType::Update:

                    this->c_solver->SetInitialGuess( this->c_stages[n] );
                    break;

                default:
                {   std::string error_message = "Invalid IDCType '" + IDCType_to_String.at( this->idc_type )
                                                + "' in '" + std::string(__func__) + "'.\n";

                    PRINT_ERROR( error_message.c_str() )
                    throw std::invalid_argument( error_message );
                }
            }

            // Compute source for implicit solve.
            c_ComputeSource_HybridI( dt, n );

            // Perform implicit solve.
            switch ( this->idc_type ) {

                case IDCType::Error:

                    this->c_solver->Solve( *this->c_error, *this->sigma_t, *this->sigma_s, h * dt );
                    OrdinateFlux::AXPY( 1.0, 1.0, *this->c_error, *this->c_stages[n] );
                    break;

                case IDCType::Update:

                    this->c_solver->Solve( *this->c_stages[n], *this->sigma_t, *this->sigma_s, h * dt );
                    break;

                default:
                {   std::string error_message = "Invalid IDCType '" + IDCType_to_String.at( this->idc_type )
                                                + "' in '" + std::string(__func__) + "'.\n";

                    PRINT_ERROR( error_message.c_str() )
                    throw std::invalid_argument( error_message );
                }
            }

            Global::TMR_collided.Stop();
        }
    }

    // --- Compute result of timestep. ---------------------------------------------------------------- //

    // Apply collocation formula.
    if (    !this->right_endpoint
         || this->use_collocation
    ) {
        ComputeCollocation_Hybrid( u_initial_condition, c_initial_condition, source, dt );
    }

    // Arrange uncollided output.
    OrdinateFlux::swap( u_initial_condition, *this->u_stages[this->num_stages - 1] );

    // Perform relabeling or arrange collided output.
    if (    this->hybrid_method == HybridMethod::HybridIa
         && this->relabel == false
         && c_initial_condition != nullptr
    ) {

        OrdinateFlux::swap( *c_initial_condition, *this->c_stages[this->num_stages - 1] );
        c_initial_condition->ZeroDensity( OpDomain::Boundary );

    } else {

        this->relabel_operator->Relabel( 1.0, *this->c_stages[this->num_stages - 1], u_initial_condition );
    }

    u_initial_condition.ZeroDensity( OpDomain::Boundary );

    // Cleanup.
    this->u_stages[0] = nullptr;
    this->c_stages[0] = nullptr;
}


//------------------------------------------------------------------------------------------------------------
//! \brief  Computes a one timestep update using hybrid I type splitting.
//!
//! \param[in,out]  initial_condition   Initially contains the initial condition \f$ \Psi^n \f$ for the
//!                                     current timestep. Upon return contains the result \f$ \Psi^{n+1} \f$
//!                                     of the timestep update.
//! \param[in]      source              Source term of the transport system. <br>
//!                                     The source term is assumed to be constant across the timestep interval.
//! \param[in]      dt                  Timestep size. Note that this is the full IDC timestep size from which
//!                                     the substep sizes are determined.
//!
//! \see    IDCIntegrator::Step_Nonhybrid()
//! \see    IDCIntegrator::Step_HybridI()
//------------------------------------------------------------------------------------------------------------
void IDCIntegrator::Step_HybridII (

    OrdinateFlux & initial_condition,
    const RKDG::OrdinateFlux & source,
    const double dt
) {

    this->u_stages[0] = &initial_condition;

    if (    this->left_endpoint
         && BitmaskHasAny( Global::reflecting_boundaries, BoundaryEdge::All )
    ) {
        this->u_stages[0]->ReflectBoundaries();
    }

    // --- Prediction steps. -------------------------------------------------------------------------- //

    for ( int64_t n = 1; n < this->num_stages; ++n ) {

        // Normalized substep length.
        const double h = this->nodes[n] - this->nodes[n-1];

        // --- High-resolution (uncollided). --- //

        Global::TMR_uncollided.Start();

        // Compute source for implicit solve.
        OrdinateFlux::AXPY( 1.0/(h * dt), 0.0, *this->u_stages[n-1], *this->u_stages[n] );
        this->u_stages[n]->ZeroDensity( OpDomain::Boundary );
        OrdinateFlux::AXPY( 1.0, 1.0, source, *this->u_stages[n] );

        // Perform implicit solve.
        this->u_solver->Solve( *this->u_stages[n], *this->sigma_t, *this->sigma_s, h * dt );

        Global::TMR_uncollided.Stop();

        // --- Low-resolution (collided). --- //

        Global::TMR_collided.Start();

        if (    this->hybrid_method == HybridMethod::HybridIIc
             || this->idc_type == IDCType::Error
        ) {
            // Set initial guess for solve.
            if ( n > 1 )
                this->c_solver->SetInitialGuess( this->c_error );
            else
                this->c_solver->SetInitialGuess( nullptr );

            // For hybrid-IIa methods, relabel previous substep before computing next implicit Euler step.
            if (    this->hybrid_method == HybridMethod::HybridIIa
                 && n > 1
            ) {
                // \psi_u^{n-1,[0]} + R \psi_c^{n-1,[0]}
                Global::TMR_collided.Stop();
                this->relabel_operator->Relabel( 1.0, *this->c_error, *this->u_stages[n-1] );
                Global::TMR_collided.Start();

                // \frac{1}{h_n \Delta t} \psi_c^{n-1,[0]}
                this->c_error->Scale( 1.0/(h * dt) );

            } else {

                this->c_error->ZeroDensity();
            }

            // SP psi_u^{n,[0]}
            TransportOperator::Pmv( 1.0, 0.0, *this->u_stages[n], *this->phi );
            this->phi->ZeroDensity( OpDomain::Boundary );
            TransportOperator::Smv( 1.0, 1.0, *this->sigma_s, *this->phi, *this->c_error );

            // Perform implicit solve.
            this->c_solver->Solve( *this->c_error, *this->sigma_t, *this->sigma_s, h * dt );

            Global::TMR_collided.Stop();

            // For hybrid-IIb methods, relabel step after computing it.
            if ( this->hybrid_method == HybridMethod::HybridIIb )
                this->relabel_operator->Relabel( 1.0, *this->c_error, *this->u_stages[n] );

            // For hybrid-IIc methods, perform Nyström reconstruction.
            if ( this->hybrid_method == HybridMethod::HybridIIc ) {

                Global::TMR_relabel.Start();

                // P psi_u^{n,[0]} + P psi_c^{n,[0]}
                TransportOperator::Pmv( 1.0, 0.0, *this->u_stages[n], *this->phi );
                TransportOperator::Pmv( 1.0, 1.0, *this->c_error, *this->phi );

                // Include term for previous timestep value, without boundary conditions.
                OrdinateFlux::AXPY( 1.0/(h * dt), 0.0, *this->u_stages[n-1], *this->u_stages[n] );
                this->u_stages[n]->ZeroDensity( OpDomain::Boundary );

                // Include source term, with boundary conditions.
                OrdinateFlux::AXPY( 1.0, 1.0, source, *this->u_stages[n] );

                // Include scattering source, without boundary conditions.
                this->phi->ZeroDensity( OpDomain::Boundary );
                TransportOperator::Smv( 1.0, 1.0, *this->sigma_s, *this->phi, *this->u_stages[n] );

                // Apply sweep to compute the reconstructed solution.
                this->u_solver->Solve( *this->u_stages[n], *this->sigma_t, *this->sigma_s, h * dt );

                Global::TMR_relabel.Stop();
            }

        } else if ( this->idc_type == IDCType::Update ) {

            // For hybrid-IIa methods, relabel previous substep before computing next implicit Euler step.
            this->c_stages[n]->ZeroDensity();

            if (    this->hybrid_method == HybridMethod::HybridIIa
                 && this->c_stages[n-1] != nullptr
            ) {
                // \frac{1}{h_n \Delta t} \psi_c^{n-1,[0]}
                OrdinateFlux::AXPY( 1.0/(h * dt), 1.0, *this->c_stages[n-1], *this->c_stages[n] );
            }

            // SP psi_u^{n,[0]}
            TransportOperator::Pmv( 1.0, 0.0, *this->u_stages[n], *this->phi );
            this->phi->ZeroDensity( OpDomain::Boundary );
            TransportOperator::Smv( 1.0, 1.0, *this->sigma_s, *this->phi, *this->c_stages[n] );

            // Set initial guess for solve.
            this->c_solver->SetInitialGuess( this->c_stages[n-1] );

            // Perform implicit solve.
            this->c_solver->Solve( *this->c_stages[n], *this->sigma_t, *this->sigma_s, h * dt );

            Global::TMR_collided.Stop();

            // For hybrid-IIb methods, relabel step after computing it.
            if ( this->hybrid_method == HybridMethod::HybridIIb )
                this->relabel_operator->Relabel( 1.0, *this->c_stages[n], *this->u_stages[n] );

        } else {
            std::string error_message = "Invalid IDCType '" + IDCType_to_String.at( this->idc_type )
                                        + "' in '" + std::string(__func__) + "'.\n";

            PRINT_ERROR( error_message.c_str() )
            throw std::invalid_argument( error_message );
        }
    }

    // For hybrid-IIa methods, finish relabeling collided substeps.
    if ( this->hybrid_method == HybridMethod::HybridIIa ) {

        switch( this->idc_type ) {

            case IDCType::Error:

                this->relabel_operator->Relabel( 1.0, *this->c_error, *this->u_stages[this->num_stages - 1] );
                break;

            case IDCType::Update:

                for ( int64_t n = 1; n < this->num_stages; ++n )
                    this->relabel_operator->Relabel( 1.0, *this->c_stages[n], *this->u_stages[n] );
                break;

            default:
            {   std::string error_message = "Invalid IDCType '" + IDCType_to_String.at( this->idc_type )
                                            + "' in '" + std::string(__func__) + "'.\n";

                PRINT_ERROR( error_message.c_str() )
                throw std::invalid_argument( error_message );
            }
        }
    }

    // --- Correction iterations. --------------------------------------------------------------------- //

    for ( int64_t p = 0; p < this->num_corrections; ++p ) {

        Global::TMR_uncollided.Start();
        ComputeResiduals();
        Global::TMR_uncollided.Stop();

        for ( int64_t n = 1; n < this->num_stages; ++n ) {

            // Normalized substep length.
            const double h = this->nodes[n] - this->nodes[n-1];

            // --- High-resolution (uncollided) solve. --- //

            Global::TMR_uncollided.Start();

            // Compute source for implicit solve.
            u_ComputeSource_HybridII( source, dt, n );

            // Perform implicit solve.
            if (    this->hybrid_method == HybridMethod::HybridIIc
                 || this->idc_type == IDCType::Error
            )
                this->u_solver->Solve( *this->u_error, *this->sigma_t, *this->sigma_s, h * dt );
            else
                this->u_solver->Solve( *this->u_stages[n], *this->sigma_t, *this->sigma_s, h * dt );

            // Apply correction as necessary.
            if (    this->hybrid_method != HybridMethod::HybridIIc
                 && this->idc_type == IDCType::Error
            ) {
                OrdinateFlux::AXPY( 1.0, 1.0, *this->u_error, *this->u_stages[n] );
            }

            Global::TMR_uncollided.Stop();

            // --- Low-resolution (collided) solve. --- //

            // Apply relabeled error from previous step if necessary.
            if (    this->idc_type == IDCType::Error
                 && this->hybrid_method == HybridMethod::HybridIIa
                 && n > 1
            ) {
                // \psi^{n-1,[p-1]} + e_u^{n-1,[p-1]} + R e_c^{n-1,[p-1]}
                this->relabel_operator->Relabel( 1.0, *this->c_error, *this->u_stages[n-1] );
            }

            Global::TMR_collided.Start();

            // Set initial guess for implicit solve.
            if (    this->hybrid_method == HybridMethod::HybridIIc
                 || this->idc_type == IDCType::Error
            )
                this->c_solver->SetInitialGuess( nullptr );
            else
                this->c_solver->SetInitialGuess( this->c_stages[n] );

            // Compute source for implicit solve.
            c_ComputeSource_HybridII( dt, n );

            // Perform implicit solve.
            this->c_solver->Solve( *this->c_error, *this->sigma_t, *this->sigma_s, h * dt );

            Global::TMR_collided.Stop();

            // For hybrid-IIb methods, relabel collided correction.
            if ( this->hybrid_method == HybridMethod::HybridIIb )
                this->relabel_operator->Relabel( 1.0, *this->c_error, *this->u_stages[n] );

            // For hybrid-IIc methods, perform Nyström reconstruction with collided correction.
            if ( this->hybrid_method == HybridMethod::HybridIIc ) {

                Global::TMR_relabel.Start();
                ComputeNystromSource( source, dt, n );

                switch ( this->idc_type ) {

                    case IDCType::Error:

                        this->u_solver->Solve( *this->u_error, *this->sigma_t, *this->sigma_s, h * dt );
                        OrdinateFlux::AXPY( 1.0, 1.0, *this->u_error, *this->u_stages[n] );
                        break;

                    case IDCType::Update:

                        this->u_solver->Solve( *this->u_stages[n], *this->sigma_t, *this->sigma_s, h * dt );
                        break;

                    default:
                    {   std::string error_message = "Invalid IDCType '"
                                                    + IDCType_to_String.at( this->idc_type )
                                                    + "' in '" + std::string(__func__) + "'.\n";

                        PRINT_ERROR( error_message.c_str() )
                        throw std::invalid_argument( error_message );
                    }
                }

                Global::TMR_relabel.Stop();

            // For update forms of hybrid-IIa and hybrid-IIb methods, store collided correction for later use.
            } else if ( this->idc_type == IDCType::Update ) {

                std::swap( this->c_error, this->c_stages[n] );
            }
        }

        // Finish relabeling collided corrections.
        if ( this->hybrid_method == HybridMethod::HybridIIa ) {

            switch( this->idc_type ) {

                case IDCType::Error:

                    this->relabel_operator->Relabel( 1.0, *this->c_error,
                                                     *this->u_stages[this->num_stages - 1] );
                    break;

                case IDCType::Update:

                    for ( int64_t n = 1; n < this->num_stages; ++n )
                        this->relabel_operator->Relabel( 1.0, *this->c_stages[n], *this->u_stages[n] );
                    break;

                default:
                {   std::string error_message = "Invalid IDCType '" + IDCType_to_String.at( this->idc_type )
                                                + "' in '" + std::string(__func__) + "'.\n";

                    PRINT_ERROR( error_message.c_str() )
                    throw std::invalid_argument( error_message );
                }
            }
        }
    }

    // --- Compute result of timestep. ---------------------------------------------------------------- //

    if (    !this->right_endpoint
         || this->use_collocation
    ) {
        ComputeCollocation_Nonhybrid( initial_condition, source, dt );
    }

    OrdinateFlux::swap( initial_condition, *this->u_stages[this->num_stages - 1] );
    initial_condition.ZeroDensity( OpDomain::Boundary );

    this->u_stages[0] = nullptr;
}


//============================================================================================================
//=== PRIVATE HELPER ROUTINES ================================================================================
//============================================================================================================


//------------------------------------------------------------------------------------------------------------
//! \brief  Computes the final output value of the IDC method by applying the collocation formula to the
//!         substep values of the non-hybrid or reconstructed hybrid approximation.
//!
//! The output value is stored in <code> this->stages[ this->num_stages - 1 ] </code>.
//!
//! \param[in]  initial_condition   The initial condition \f$ \Psi^n \f$ for the current timestep.
//! \param[in]  source              Source term of the transport system. <br>
//!                                 The source term is assumed to be constant across the timestep interval.
//! \param[in]  dt                  Timestep size.
//!
//! \see    IDCIntegrator::Step_Nonhybrid()
//! \see    IDCIntegrator::Step_HybridI()
//------------------------------------------------------------------------------------------------------------
void IDCIntegrator::ComputeCollocation_Nonhybrid (

    const OrdinateFlux & initial_condition,
    const OrdinateFlux & source,
    const double dt
) {

    OrdinateFlux & result = *this->stages[this->num_stages - 1];
    this->phi->ZeroDensity();

    // Apply transport operators to stage vectors.
    for ( int64_t n = this->num_stages - 1; n >= ((int)(!this->left_endpoint)); --n ) {

        TransportOperator::Pmv( dt * Weight(0,n), 1.0, *this->stages[n], *this->phi );
        TransportOperator::Lmv( -dt * Weight(0,n), ( n == this->num_stages - 1 ? 0.0 : 1.0 ), *this->sigma_t,
                                *this->stages[n], result );
    }

    this->phi->ZeroDensity( OpDomain::Boundary );
    TransportOperator::Smv( 1.0, 1.0, *this->sigma_s, *this->phi, result );

    // Include initial condition and source for timestep.
    OrdinateFlux::AXPY( 1.0, 1.0, initial_condition, result );
    OrdinateFlux::AXPY( dt, 1.0, source, result );
}


//------------------------------------------------------------------------------------------------------------
//! \brief  Computes the final output values of the IDC method by applying the collocation formula to the
//!         substep values of the hybrid approximation.
//!
//! The output values are stored in <code> this->u_stages[this->num_stages - 1] </code> and
//! <code> this->c_stages[this->num_stages - 1] </code>.
//!
//! \attention  If \c IDCIntegrator::relabel is "true" then the initial collided distribution is assumed to be
//!             zero and \pp{c_initial_condition} is not dereferenced.
//!
//! \param[in]  u_initial_condition The initial condition \f$ \Psi^n \f$ for the uncollided flux in the
//!                                 current timestep.
//! \param[in]  c_initial_condition The initial condition \f$ \Psi^n \f$ for the collided flux in the current
//!                                 timestep. <br>
//!                                 If this is \c null, then the initial value for the collided flux is
//!                                 assumed to be zero.
//! \param[in]  source              Source term of the transport system. <br>
//!                                 The source term is assumed to be constant across the timestep interval.
//! \param[in]  dt                  Timestep size.
//!
//! \see    IDCIntegrator::Step_HybridII()
//------------------------------------------------------------------------------------------------------------
void IDCIntegrator::ComputeCollocation_Hybrid (

    const OrdinateFlux & u_initial_condition,
    const OrdinateFlux * const c_initial_condition,
    const OrdinateFlux & source,
    const double dt
) {

    OrdinateFlux & u_result = *this->u_stages[this->num_stages - 1];
    OrdinateFlux & c_result = *this->c_stages[this->num_stages - 1];

    this->phi->ZeroDensity();

    // --- Construct uncollided result. -------------------------------------------------------------- //

    Global::TMR_uncollided.Start();

    // Apply transport operators.
    for ( int64_t n = this->num_stages - 1; n >= ((int)(!this->left_endpoint)); --n ) {

        TransportOperator::Pmv( dt * Weight(0,n), 1.0, *this->u_stages[n], *this->phi );
        TransportOperator::Lmv( -dt * Weight(0,n), ( n == this->num_stages - 1 ? 0.0 : 1.0 ), *this->sigma_t,
                                *this->u_stages[n], u_result );
    }

    // Include initial condition and source for timestep.
    OrdinateFlux::AXPY( 1.0, 1.0, u_initial_condition, u_result );
    OrdinateFlux::AXPY( dt, 1.0, source, u_result );

    Global::TMR_uncollided.Stop();

    // --- Construct collided result. ---------------------------------------------------------------- //

    Global::TMR_collided.Start();

    // Apply transport operators.
    for ( int64_t n = this->num_stages - 1; n >= ((int)(!this->left_endpoint)); --n ) {

        if ( this->c_stages[n] == nullptr ) {  continue;  }

        TransportOperator::Pmv( dt * Weight(0,n), 1.0, *this->c_stages[n], *this->phi );
        TransportOperator::Lmv( -dt * Weight(0,n), ( n == this->num_stages - 1 ? 0.0 : 1.0 ), *this->sigma_t,
                                *this->c_stages[n], c_result );
    }

    TransportOperator::Smv( 1.0, 1.0, *this->sigma_s, *this->phi, c_result );

    // Include collided initial condition if necessary (collided source is always zero).
    if ( c_initial_condition != nullptr )
        OrdinateFlux::AXPY( 1.0, 1.0, *c_initial_condition, c_result );

    Global::TMR_collided.Stop();
}


//------------------------------------------------------------------------------------------------------------
//! \brief  Computes residuals for the non-hybrid operator.
//!
//! For each substep value computed on the previous correction level, this routine evaluates the residual of
//! each substep. These residuals (which are stored in IDCIntegrator::residuals) are used to compute the
//! time-integral portion of the IDC augmented source term used to compute the next correction iteration.
//!
//! This routine is used by the non-hybrid and hybrid-II integrators.
//!
//! \see    IDCIntegrator::ComputeSource()
//! \see    IDCIntegrator::u_ComputeSource_HybridI()
//! \see    IDCIntegrator::c_ComputeSource_HybridI()
//------------------------------------------------------------------------------------------------------------
void IDCIntegrator::ComputeResiduals ( void ) {

    // Compute for each quadrature point.
    for ( int64_t n = ((int)(!this->left_endpoint)); n < this->num_stages; ++n ) {

        TransportOperator::Lmv( 1.0, 0.0, *this->sigma_t, *this->stages[n], *this->residuals[n] );
        TransportOperator::Pmv( 1.0, 0.0, *this->stages[n], *this->phi );
        this->phi->ZeroDensity( OpDomain::Boundary );
        TransportOperator::Smv( -1.0, 1.0, *this->sigma_s, *this->phi, *this->residuals[n] );
    }
}


//------------------------------------------------------------------------------------------------------------
//! \brief  Computes fictitious sources for non-hybrid corrections.
//!
//! For the IDC quadrature node with index \pp{n}, this routine computes the source term for solving the
//! non-hybrid error/correction equation. The result is stored in IDCIntegrator::error.
//!
//! This routine is used by the non-hybrid integrators.
//!
//! \param[in]      source      Source term of the transport system. <br>
//!                             The source term is assumed to be constant across the timestep interval.
//! \param[in]      dt          IDC timestep size.
//! \param[in]      n           Index denoting the time interval to integrate across.
//!
//! \see    IDCIntegrator::ComputeResiduals()
//------------------------------------------------------------------------------------------------------------
void IDCIntegrator::ComputeSource (

    const OrdinateFlux & source,
    const double dt,
    const int64_t n
) {

    if ( n <= 0 ) {  return;  }

    const double h = this->nodes[n] - this->nodes[n-1];
    OrdinateFlux * idc_source = nullptr;

    switch ( this->idc_type ) {

        case IDCType::Error:
        {
            idc_source = this->error;

            // - \frac{1}{h_n \Delta t} \psi^{n,[p-1]}
            OrdinateFlux::AXPY( -1.0/(h * dt), 0.0, *this->stages[n], *idc_source );

            // Use quadrature to integrate functional portion of residual.
            for ( int64_t l = ((int)(!this->left_endpoint)); l < this->num_stages; ++l ) {

                const double weight = -Weight(n,l) / h;
                OrdinateFlux::AXPY( weight, 1.0, *this->residuals[l], *idc_source );
            }
        } break;

        case IDCType::Update:
        {
            idc_source = this->stages[n];
            idc_source->ZeroDensity();

            // Use quadrature to integrate functional portion of residual.
            for ( int64_t l = ((int)(!this->left_endpoint)); l < this->num_stages; ++l ) {

                const double weight = ( n == l ? 1.0 : 0.0 ) - Weight(n,l) / h;
                OrdinateFlux::AXPY( weight, 1.0, *this->residuals[l], *idc_source );
            }
        } break;

        default:
        {   std::string error_message = "Invalid IDCType '" + IDCType_to_String.at( this->idc_type )
                                        + "' in '" + std::string(__func__) + "'.\n";

            PRINT_ERROR( error_message.c_str() )
            throw std::invalid_argument( error_message );
        }
    }

    // \frac{1}{h_n \Delta t} \psi^{n-1,[p]}
    OrdinateFlux::AXPY( 1.0/(h * dt), 1.0, *this->stages[n-1], *idc_source );

    if ( this->idc_type == IDCType::Update )
        idc_source->ZeroDensity( OpDomain::Boundary );

    // Include external source (assumed to be constant across timestep interval).
    OrdinateFlux::AXPY( 1.0, 1.0, source, *idc_source );

    // Boundary condition for the error is always zero.
    if ( this->idc_type == IDCType::Error )
        idc_source->ZeroDensity( OpDomain::Boundary );
}


//------------------------------------------------------------------------------------------------------------
//! \brief  Computes residuals for the hybrid-I operator.
//!
//! For each substep value computed on the previous correction level, this routine evaluates the residual of
//! the hybrid-I operator at each substep. These residuals (which are stored in IDCIntegrator::u_residuals
//! and IDCIntegrator::c_residuals) are used to compute the time-integral portion of the IDC augmented source
//! terms used to compute the next correction iteration.
//!
//! This routine is used by the hybrid-I integrators.
//!
//! \see    IDCIntegrator::u_ComputeSource_HybridI()
//! \see    IDCIntegrator::c_ComputeSource_HybridI()
//------------------------------------------------------------------------------------------------------------
void IDCIntegrator::ComputeResiduals_HybridI ( void ) {

    // --- Uncollided. -------------------------------------------------------------------------------- //

    Global::TMR_uncollided.Start();

    for ( int64_t n = ((int)(!this->left_endpoint)); n < this->num_stages; ++n )
        TransportOperator::Lmv( 1.0, 0.0, *this->sigma_t, *this->u_stages[n], *this->u_residuals[n] );

    Global::TMR_uncollided.Stop();

    // --- Collided. ---------------------------------------------------------------------------------- //

    Global::TMR_collided.Start();

    for ( int64_t n = ((int)(!this->left_endpoint)); n < this->num_stages; ++n ) {

        this->c_residuals[n]->ZeroDensity();

        TransportOperator::Pmv( 1.0, 0.0, *this->u_stages[n], *this->phi );

        if ( this->c_stages[n] != nullptr ) {

            TransportOperator::Lmv( 1.0, 1.0, *this->sigma_t, *this->c_stages[n], *this->c_residuals[n] );
            TransportOperator::Pmv( 1.0, 1.0, *this->c_stages[n], *this->phi );
        }

        this->phi->ZeroDensity( OpDomain::Boundary );
        TransportOperator::Smv( -1.0, 1.0, *this->sigma_s, *this->phi, *this->c_residuals[n] );
    }

    Global::TMR_collided.Stop();
}


//------------------------------------------------------------------------------------------------------------
//! \brief  Computes fictitious sources for the uncollided hybrid-I correction equation.
//!
//! For the IDC quadrature node stored in index \pp{n}, this routine computes the source term for solving the
//! uncollided error equation of the split time-dependent system.
//!
//! This routine is used by the hybrid-I integrators.
//!
//! \param[in]      source      Source term of the transport system. <br>
//!                             The source term is assumed to be constant across the timestep interval.
//! \param[in]      dt          IDC timestep size.
//! \param[in]      n           Index denoting the time interval to integrate across.
//!
//! \see    IDCIntegrator::ComputeResiduals_HybridI()
//! \see    IDCIntegrator::c_ComputeSource_HybridI()
//------------------------------------------------------------------------------------------------------------
void IDCIntegrator::u_ComputeSource_HybridI (

    const OrdinateFlux & source,
    const double dt,
    const int64_t n
) {

    if ( n <= 0 ) {  return;  }

    const double h = this->nodes[n] - this->nodes[n-1];
    OrdinateFlux * idc_source = nullptr;

    switch ( this->idc_type ) {

        case IDCType::Error:
        {
            // Set location of OrdinateFlux object in which to construct source.
            idc_source = this->u_error;

            // - \frac{1}{h_n \Delta t} \psi_u^{n,[p-1]}
            OrdinateFlux::AXPY( -1.0/(h * dt), 0.0, *this->u_stages[n], *idc_source );

            // Use quadrature to integrate functional portion of residual.
            for ( int64_t l = ((int)(!this->left_endpoint)); l < this->num_stages; ++l ) {

                const double weight = -Weight(n,l) / h;
                OrdinateFlux::AXPY( weight, 1.0, *this->u_residuals[l], *idc_source );
            }
        } break;

        case IDCType::Update:
        {
            // Set location of OrdinateFlux object in which to construct source.
            idc_source = this->u_stages[n];
            idc_source->ZeroDensity();

            // Use quadrature to integrate functional portion of residual.
            for ( int64_t l = ((int)(!this->left_endpoint)); l < this->num_stages; ++l ) {

                const double weight = ( n == l ? 1.0 : 0.0 ) - Weight(n,l) / h;
                OrdinateFlux::AXPY( weight, 1.0, *this->u_residuals[l], *idc_source );
            }
        } break;

        default:
        {   std::string error_message = "Invalid IDCType '" + IDCType_to_String.at( this->idc_type )
                                        + "' in '" + std::string(__func__) + "'.\n";

            PRINT_ERROR( error_message.c_str() )
            throw std::invalid_argument( error_message );
        }
    }

    // \frac{1}{h_n \Delta t} \psi_u^{n-1,[p]}
    OrdinateFlux::AXPY( 1.0/(h * dt), 1.0, *this->u_stages[n-1], *idc_source );

    if (    this->hybrid_method == HybridMethod::HybridIb
        && this->c_stages[n-1] != nullptr
    ) {
        // \frac{1}{h_n \Delta t} R \psi_c^{n-1,[p]}
        Global::TMR_uncollided.Stop();
        this->relabel_operator->Relabel( 1.0/(h * dt), *this->c_stages[n-1], *idc_source );
        Global::TMR_uncollided.Start();
    }

    if ( this->idc_type == IDCType::Update )
        idc_source->ZeroDensity( OpDomain::Boundary );

    // Include external source (assumed to be constant across timestep interval).
    OrdinateFlux::AXPY( 1.0, 1.0, source, *idc_source );

    // Boundary condition for the error is always zero.
    if ( this->idc_type == IDCType::Error )
        idc_source->ZeroDensity( OpDomain::Boundary );
}


//------------------------------------------------------------------------------------------------------------
//! \brief  Computes fictitious sources for the collided hybrid-I correction equation.
//!
//! For the IDC quadrature node stored in index \pp{n}, this routine computes the source term for solving the
//! collided error equation of the split time-dependent system.
//!
//! This routine is used by the hybrid-I integrators.
//!
//! \param[in]      dt          IDC timestep size.
//! \param[in]      n           Index denoting the time interval to integrate across.
//!
//! \see    IDCIntegrator::ComputeResiduals_HybridI()
//! \see    IDCIntegrator::u_ComputeSource_HybridI()
//------------------------------------------------------------------------------------------------------------
void IDCIntegrator::c_ComputeSource_HybridI (

    const double dt,
    const int64_t n
) {

    if ( n <= 0 ) {  return;  }

    const double h = this->nodes[n] - this->nodes[n-1];
    OrdinateFlux * idc_source = nullptr;

    switch ( this->idc_type ) {

        case IDCType::Error:
        {
            // Set location of OrdinateFlux object in which to construct source.
            idc_source = this->c_error;

            // - \frac{1}{h_n \Delta t} \psi_c^{n,[p-1]}
            OrdinateFlux::AXPY( -1.0/(h * dt), 0.0, *this->c_stages[n], *idc_source );

            // SP e_u^{n,[p-1]}
            TransportOperator::Pmv( 1.0, 0.0, *this->u_error, *this->phi );
            this->phi->ZeroDensity( OpDomain::Boundary );
            TransportOperator::Smv( 1.0, 1.0, *this->sigma_s, *this->phi, *idc_source );

            // Use quadrature to integrate functional portion of residual.
            for ( int64_t l = ((int)(!this->left_endpoint)); l < this->num_stages; ++l ) {

                const double weight = -Weight(n,l) / h;
                OrdinateFlux::AXPY( weight, 1.0, *this->c_residuals[l], *idc_source );
            }
        } break;

        case IDCType::Update:
        {
            // Set location of OrdinateFlux object in which to construct source.
            idc_source = this->c_stages[n];

            // SP \psi_u^{n,[p]}
            TransportOperator::Pmv( 1.0, 0.0, *this->u_stages[n], *this->phi );
            this->phi->ZeroDensity( OpDomain::Boundary );
            TransportOperator::Smv( 1.0, 0.0, *this->sigma_s, *this->phi, *idc_source );

            // Use quadrature to integrate functional portion of residual.
            for ( int64_t l = ((int)(!this->left_endpoint)); l < this->num_stages; ++l ) {

                const double weight = ( n == l ? 1.0 : 0.0 ) - Weight(n,l) / h;
                OrdinateFlux::AXPY( weight, 1.0, *this->c_residuals[l], *idc_source );
            }
        } break;

        default:
        {   std::string error_message = "Invalid IDCType '" + IDCType_to_String.at( this->idc_type )
                                        + "' in '" + std::string(__func__) + "'.\n";

            PRINT_ERROR( error_message.c_str() )
            throw std::invalid_argument( error_message );
        }
    }

    if (    this->hybrid_method == HybridMethod::HybridIa
         && this->c_stages[n-1] != nullptr
    ) {
        // \frac{1}{h_n \Delta t} \psi_c^{n-1,[p]}
        OrdinateFlux::AXPY( 1.0/(h * dt), 1.0, *this->c_stages[n-1], *idc_source );
    }
}


//------------------------------------------------------------------------------------------------------------
//! \brief  Computes fictitious sources for the uncollided portion of hybrid-I corrections.
//!
//! For the IDC quadrature node with index \pp{n}, this routine computes the source term for solving the
//! uncollided error/correction equation for the hybrid-I system.
//!
//! This routine is used by the hybrid-II integrators.
//!
//! \param[in]      source      Source term of the transport system. <br>
//!                             The source term is assumed to be constant across the timestep interval.
//! \param[in]      dt          IDC timestep size.
//! \param[in]      n           Index denoting the time interval to integrate across.
//!
//! \see    IDCIntegrator::ComputeResiduals()
//! \see    IDCIntegrator::c_ComputeSource_HybridII()
//------------------------------------------------------------------------------------------------------------
void IDCIntegrator::u_ComputeSource_HybridII (

    const OrdinateFlux & source,
    const double dt,
    const int64_t n
) {

    if ( n <= 0 ) {  return;  }

    const double h = this->nodes[n] - this->nodes[n-1];
    OrdinateFlux * idc_source = nullptr;

    if (    this->hybrid_method == HybridMethod::HybridIIc
         || this->idc_type == IDCType::Error
    ) {
        // Set location of OrdinateFlux object in which to construct source.
        idc_source = this->u_error;

        // \frac{1}{h_n \Delta t} ( \psi^{n-1,[p-1]} + e_u^{n-1,[p-1]} )
        OrdinateFlux::AXPY( 1.0/(h * dt), 0.0, *this->u_stages[n-1], *idc_source );

        // - \frac{1}{h_n \Delta t} \psi^{n,[p-1]}
        OrdinateFlux::AXPY( -1.0/(h * dt), 1.0, *this->u_stages[n], *idc_source );

    } else if ( this->idc_type == IDCType::Update ) {

        // Set location of OrdinateFlux object in which to construct source.
        idc_source = this->u_stages[n];
        /* Do not zero! */

        // - P \psi^{n,[p-1]}
        //
        // NOTE: Value stored in phi is used by c_ComputeSource_HybridII for collided correction.
        //
        TransportOperator::Pmv( -1.0, 0.0, *this->u_stages[n], *this->phi );

        // - R \psi_c^{n,[p-1]}
        Global::TMR_uncollided.Stop();
        this->relabel_operator->Relabel( -1.0, *this->c_stages[n], *idc_source );
        Global::TMR_uncollided.Start();

        // L ( \psi^{n,[p-1]} - R \psi_c^{n,[p-1]} )
        TransportOperator::Lmv( 1.0, 0.0, *this->sigma_t, *idc_source, *idc_source );

        // \frac{1}{h_n \Delta t} \psi{_u, }^{n-1,[p]}
        OrdinateFlux::AXPY( 1.0/(h * dt), 1.0, *this->u_stages[n-1], *idc_source );

        switch ( this->hybrid_method ) {

            case HybridMethod::HybridIIa:
            {
                if ( n == 1 )
                    this->c_error->ZeroDensity();

                // \psi_c^{n-1,[p-1]} - \psi_c^{n,[p-1]}
                //
                // NOTE: Value stored in c_error is used by c_ComputeSource_HybridII for collided
                //       correction.
                //
                OrdinateFlux::AXPY( -1.0, 1.0, *this->c_stages[n], *this->c_error );

                // \frac{1}{h_n \Delta t} R ( \psi_c^{n-1,[p-1]} - \psi_c^{n,[p-1]} )
                Global::TMR_uncollided.Stop();
                this->relabel_operator->Relabel( 1.0/(h * dt), *this->c_error, *idc_source );
                Global::TMR_uncollided.Start();

            } break;

            case HybridMethod::HybridIIb:
            {
                // - \frac{1}{h_n \Delta t} R \psic^{n,[p-1]}
                Global::TMR_uncollided.Stop();
                this->relabel_operator->Relabel( -1.0/(h * dt), *this->c_stages[n], *idc_source );
                Global::TMR_uncollided.Start();

            } break;

            default:
            {   std::string error_message = "Invalid hybrid splitting '"
                                            + HybridMethod_to_String.at( this->hybrid_method )
                                            + "' in '" + std::string(__func__) + "'.\n";

                PRINT_ERROR( error_message.c_str() )
                throw std::invalid_argument( error_message );
            }
        }

    } else {
        std::string error_message = "Invalid IDCType '" + IDCType_to_String.at( this->idc_type )
                                    + "' in '" + std::string(__func__) + "'.\n";

        PRINT_ERROR( error_message.c_str() )
        throw std::invalid_argument( error_message );
    }

    // Use quadrature to integrate functional portion of residual.
    for ( int64_t l = ((int)(!this->left_endpoint)); l < this->num_stages; ++l ) {

        const double weight = -Weight(n,l) / h;
        OrdinateFlux::AXPY( weight, 1.0, *this->u_residuals[l], *idc_source );
    }

    if ( this->idc_type == IDCType::Update )
        idc_source->ZeroDensity( OpDomain::Boundary );

    // Include external source (assumed to be constant across timestep interval).
    OrdinateFlux::AXPY( 1.0, 1.0, source, *idc_source );

    // Boundary condition for the error is always zero (note that Hybrid-IIc always solves for the error).
    if (    this->hybrid_method == HybridMethod::HybridIIc
         || this->idc_type == IDCType::Error
    ) {
        idc_source->ZeroDensity( OpDomain::Boundary );
    }
}


//------------------------------------------------------------------------------------------------------------
//! \brief  Computes the IDC augmented source term used to compute an approximation of the error in the
//!         collided component of the split time-dependent system at the \pp{n}th quadrature node.
//!
//! For the IDC quadrature node stored in index \pp{n}, this routine computes the source term for solving the
//! collided error equation of the split time-dependent system.
//!
//! This routine is used by the hybrid-II integrators.
//!
//! \param[in]      dt          IDC timestep size.
//! \param[in]      n           Index denoting the time interval to integrate across.
//!
//! \see    IDCIntegrator::ComputeResiduals()
//! \see    IDCIntegrator::u_ComputeSource_HybridII()
//------------------------------------------------------------------------------------------------------------
void IDCIntegrator::c_ComputeSource_HybridII (

    const double dt,
    const int64_t n
) {

    if ( n <= 0 ) {  return;  }

    const double h = this->nodes[n] - this->nodes[n-1];

    if (    this->hybrid_method == HybridMethod::HybridIIc
         || this->idc_type == IDCType::Error
    ) {
        // SP e_u^{n,[p-1]}
        TransportOperator::Pmv( 1.0, 0.0, *this->u_error, *this->phi );
        this->phi->ZeroDensity( OpDomain::Boundary );

        double coeff = 0.0;

        if (    this->hybrid_method == HybridMethod::HybridIIa
             && n > 1
        ) {
            // \frac{1}{h_n \Delta t} e_c^{n-1,[p-1]}
            coeff = 1.0/(h * dt);
        }

        TransportOperator::Smv( 1.0, coeff, *this->sigma_s, *this->phi, *this->c_error );

    } else if ( this->idc_type == IDCType::Update ) {
        //
        // NOTE: Hybrid-IIa and hybrid-IIb IDCType::Update methods require value of -P \psi^{n,[p-1]} stored
        //       in phi by u_ComputeSource_HybridII.
        //
        switch ( this->hybrid_method ) {

            case HybridMethod::HybridIIa:
            {
                // \frac{1}{h_n \Delta t} ( \psi_c^{n,[p-1]} - \psi_c^{n-1,[p-1]} )
                //
                // NOTE: Uses value stored in c_error by u_ComputeSource_HybridII.
                //
                this->c_error->Scale( -1.0/(h * dt) );

                if ( n > 1 )
                    // \frac{1}{h_n \Delta t} \psi_c^{n-1,[p]}
                    OrdinateFlux::AXPY( 1.0/(h * dt), 1.0, *this->c_stages[n-1], *this->c_error );

            } break;

            case HybridMethod::HybridIIb:
            {
                // \frac{1}{h_n \Delta t} \psi_c^{n,[p-1]}
                OrdinateFlux::AXPY( 1.0/(h * dt), 0.0, *this->c_stages[n], *this->c_error );

            } break;

            default:
            {   std::string error_message = "Invalid hybrid splitting '"
                                            + HybridMethod_to_String.at( this->hybrid_method )
                                            + "' in '" + std::string(__func__) + "'.\n";

                PRINT_ERROR( error_message.c_str() )
                throw std::invalid_argument( error_message );
            }
        }

        // L \psi_c^{n,[p-1]}
        TransportOperator::Lmv( 1.0, 1.0, *this->sigma_t, *this->c_stages[n], *this->c_error );

        // P \psi_u^{n,[p]}
        TransportOperator::Pmv( 1.0, 1.0, *this->u_stages[n], *this->phi );

        // SP ( \psi_u^{n,[p]} - \psi^{n,[p-1]} )
        this->phi->ZeroDensity( OpDomain::Boundary );
        TransportOperator::Smv( 1.0, 1.0, *this->sigma_s, *this->phi, *this->c_error );

    } else {
        std::string error_message = "Invalid IDCType '" + IDCType_to_String.at( this->idc_type )
                                    + "' in '" + std::string(__func__) + "'.\n";

        PRINT_ERROR( error_message.c_str() )
        throw std::invalid_argument( error_message );
    }
}


//------------------------------------------------------------------------------------------------------------
//! \brief  Computes source terms for the hybrid-IIc Nyström reconstruction.
//!
//! For the IDC quadrature node with index \pp{n}, this routine computes the source term for solving the
//! Nyström reconstruction system for the hybrid-IIc method.
//!
//! This routine is used by the hybrid-IIc integrators.
//!
//! \param[in]      source      Source term of the transport system. <br>
//!                             The source term is assumed to be constant across the timestep interval.
//! \param[in]      dt          IDC timestep size.
//! \param[in]      n           Index denoting the time interval to integrate across.
//!
//! \see    IDCIntegrator::ComputeResiduals()
//! \see    IDCIntegrator::u_ComputeSource_HybridII()
//! \see    IDCIntegrator::c_ComputeSource_HybridII()
//------------------------------------------------------------------------------------------------------------
void IDCIntegrator::ComputeNystromSource (

    const OrdinateFlux & source,
    const double dt,
    const int64_t n
) {

    if ( n <= 0 ) {  return;  }

    const double h = this->nodes[n] - this->nodes[n-1];
    OrdinateFlux * nystrom_source = nullptr;

    // P e_u^{n,[p-1]} + P e_c^{n,[p-1]}
    TransportOperator::Pmv( 1.0, 0.0, *this->u_error, *this->phi );
    TransportOperator::Pmv( 1.0, 1.0, *this->c_error, *this->phi );

    switch ( this->idc_type ) {

        case IDCType::Error:
        {
            // Set location of OrdinateFlux object in which to construct source.
            nystrom_source = this->u_error;

            // \frac{1}{h_n \Delta t} psi^{n,[p-1]}
            OrdinateFlux::AXPY( -1.0/(h * dt), 0.0, *this->u_stages[n], *nystrom_source );
        } break;

        case IDCType::Update:
        {
            // Set location of OrdinateFlux object in which to construct source.
            nystrom_source = this->u_stages[n];

            // L psi^{n,[p-1]}
            TransportOperator::Lmv( 1.0, 0.0, *this->sigma_t, *this->u_stages[n], *nystrom_source );
        } break;

        default:
        {   std::string error_message = "Invalid IDCType '" + IDCType_to_String.at( this->idc_type )
                                        + "' in '" + std::string(__func__) + "'.\n";

            PRINT_ERROR( error_message.c_str() )
            throw std::invalid_argument( error_message );
        }
    }

    // \frac{1}{h_n \Delta t} psi^{n-1,[p]}
    OrdinateFlux::AXPY( 1.0/(h * dt), 1.0, *this->u_stages[n-1], *nystrom_source );

    // SP e_u^{n,[p-1]} + SP e_c^{n,[p-1]}
    TransportOperator::Smv( 1.0, 1.0, *this->sigma_s, *this->phi, *nystrom_source );

    // Use quadrature to integrate functional portion of residual.
    for ( int64_t l = ((int)(!this->left_endpoint)); l < this->num_stages; ++l ) {

        const double weight = -Weight(n,l) / h;
        OrdinateFlux::AXPY( weight, 1.0, *this->u_residuals[l], *nystrom_source );
    }

    if ( this->idc_type == IDCType::Update )
        nystrom_source->ZeroDensity( OpDomain::Boundary );

    // Include external source (assumed to be constant across timestep interval).
    OrdinateFlux::AXPY( 1.0, 1.0, source, *nystrom_source );

    // Boundary condition for the error is always zero.
    if ( this->idc_type == IDCType::Error )
        nystrom_source->ZeroDensity( OpDomain::Boundary );
}


//------------------------------------------------------------------------------------------------------------
//! \brief  Sets up IDC quadrature nodes and weights.
//!
//! Reads node type and number of nodes from the global input file Global::kvr, allocates the required memory
//! at IDCIntegrator::nodes and IDCIntegrator::weights, sets IDCIntegrator::num_stages, and computes the node
//! values and quadrature rules required for the IDC integration scheme.
//------------------------------------------------------------------------------------------------------------
void IDCIntegrator::SetNodes ( void ) {

    int64_t num_nodes;

    Global::input_list.GetValue( "idc_nodes_type", this->nodes_type, Quadrule::String_to_NodesType );
    Global::input_list.GetValue( "idc_num_nodes", num_nodes );

    switch ( this->nodes_type ) {

        case NodesType::GaussLegendre:

            this->left_endpoint  = false;
            this->right_endpoint = false;
            break;

        case NodesType::GaussRadau:

            this->left_endpoint  = false;
            this->right_endpoint = true;
            break;

        case NodesType::GaussLobatto:

            this->left_endpoint  = true;
            this->right_endpoint = true;
            break;

        case NodesType::Chebyshev:

            this->left_endpoint  = false;
            this->right_endpoint = false;
            break;

        case NodesType::ChebyshevRadau:

            this->left_endpoint  = false;
            this->right_endpoint = true;
            break;

        default:
        {   std::string error_message = "Invalid NodesType '" + NodesType_to_String.at( this->nodes_type )
                                        + "' for initializing IDCIntegrator in '" + std::string(__func__)
                                        + "'.\n";

            PRINT_ERROR( error_message.c_str() )
            throw std::invalid_argument( error_message );
        }
    }

    // Determine number of internal stages.
    this->num_stages = num_nodes + ((int)(!this->left_endpoint));

    // Compute quadrature nodes and weights.
    this->nodes = new double[ this->num_stages ];
    this->weights = new double[ this->num_stages * this->num_stages ];

    std::memset( this->nodes, 0, this->num_stages * sizeof(double) );
    std::memset( this->weights, 0, this->num_stages * this->num_stages * sizeof(double) );

    ComputeQuadratureNodes( this->num_stages - ((int)(!this->left_endpoint)),
                            this->nodes + ((int)(!this->left_endpoint)),
                            this->nodes_type, 0.0, 1.0 );

    ComputeQuadratureWeights( this->num_stages - ((int)(!this->left_endpoint)),
                              this->nodes + ((int)(!this->left_endpoint)),
                              this->weights + ((int)(!this->left_endpoint)),
                              0.0, 1.0 );

    for ( int64_t i = 1; i < this->num_stages; ++i ) {

        ComputeQuadratureWeights( this->num_stages - ((int)(!this->left_endpoint)),
                                  this->nodes + ((int)(!this->left_endpoint)),
                                  this->weights + i*this->num_stages + ((int)(!this->left_endpoint)),
                                  this->nodes[i-1], this->nodes[i] );
    }
}

