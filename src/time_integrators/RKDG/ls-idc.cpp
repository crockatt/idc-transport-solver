//------------------------------------------------------------------------------------------------------------
//! \file   time_integrators/RKDG/ls-idc.cpp
//! \brief  Implementation of LS-IDC integrators using implicit Euler steps.
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
# include "time_integrators/RKDG/ls-idc.hpp"
# include "utils/CLog.hpp"
# include "utils/global.hpp"
# include "utils/Quadrule/Quadrule.hpp"


using namespace RKDG;
using namespace Quadrule;


//============================================================================================================
//=== CONSTRUCTORS AND DESTRUCTOR ============================================================================
//============================================================================================================


//------------------------------------------------------------------------------------------------------------
//! \brief  Creates an LSIDCIntegrator object using specified parameters and returns a pointer to the object.
//!
//! This is a wrapper function to ensure that LSIDCIntegrator::Allocate is called immediately after the object
//! is constructed.
//!
//! \param[in]  input_list      List of parameters to use for initialization.
//------------------------------------------------------------------------------------------------------------
LSIDCIntegrator * LSIDCIntegrator::Create (

    const ParameterList & input_list
) {

    LSIDCIntegrator * result = new LSIDCIntegrator( input_list );
    result->Allocate();

    return result;
}


//------------------------------------------------------------------------------------------------------------
//! \brief  Constructs an LSIDCIntegrator object using the provided values.
//!
//! \param[in]  input_list      List of parameters to use for initialization.
//------------------------------------------------------------------------------------------------------------
LSIDCIntegrator::LSIDCIntegrator (

    const ParameterList & input_list

) : IDCIntegrator( input_list ) {}


//------------------------------------------------------------------------------------------------------------
//! \brief  Allocates memory for internal stage values used by the LSIDCIntegrator object.
//------------------------------------------------------------------------------------------------------------
void LSIDCIntegrator::Allocate( void ) {

    // Hybrid/non-hybrid setup.
    switch ( this->hybrid_method ) {

        case HybridMethod::None:
        {
            this->phi = new ScalarFlux( *this, DG_degree );

            if ( idc_type == IDCType::Error )
                this->error = new OrdinateFlux( *this, this->ordinate_set, DG_degree );

            this->stages = new OrdinateFlux*[ this->num_stages ];

            this->stages[0] = nullptr;

            for ( int64_t i = 1; i < this->num_stages; ++i )
                this->stages[i] = new OrdinateFlux( *this, this->ordinate_set, DG_degree );

        } break;

        case HybridMethod::HybridIa:
        case HybridMethod::HybridIb:
        {
            this->phi = new ScalarFlux( *this, DG_degree );

            if ( idc_type == IDCType::Error ) {

                this->u_error = new OrdinateFlux( *this, this->u_ordinate_set, DG_degree );
                this->c_error = new OrdinateFlux( *this, this->c_ordinate_set, DG_degree );
            }

            this->u_stages = new OrdinateFlux*[ this->num_stages ];
            this->c_stages = new OrdinateFlux*[ this->num_stages ];

            this->u_stages[0] = nullptr;
            this->c_stages[0] = nullptr;

            for ( int64_t i = 1; i < this->num_stages; ++i ) {

                this->u_stages[i] = new OrdinateFlux( *this, this->u_ordinate_set, DG_degree );
                this->c_stages[i] = new OrdinateFlux( *this, this->c_ordinate_set, DG_degree );
            }
        } break;

        case HybridMethod::HybridIIb:
        case HybridMethod::HybridIIc:
        {
            this->phi = new ScalarFlux[2];

            this->phi[0].Reconfigure( *this, DG_degree );
            this->phi[1].Reconfigure( *this, DG_degree );

            if (    this->hybrid_method == HybridMethod::HybridIIc
                 || this->idc_type == IDCType::Error
            )
                this->u_error = new OrdinateFlux( *this, this->u_ordinate_set, DG_degree );

            this->c_error = new OrdinateFlux( *this, this->c_ordinate_set, DG_degree );

            this->u_stages = new OrdinateFlux*[ this->num_stages ];
            this->u_stages[0] = nullptr;

            this->stages = this->u_stages;  // Pointer aliasing for ComputeCollocation_Nonhybrid.

            if (    this->hybrid_method != HybridMethod::HybridIIc
                 && this->idc_type == IDCType::Update
            ) {
                this->c_stages = new OrdinateFlux*[ this->num_stages ];
                this->c_stages[0] = nullptr;
            }

            for ( int64_t i = 1; i < this->num_stages; ++i ) {

                this->u_stages[i] = new OrdinateFlux( *this, this->u_ordinate_set, DG_degree );

                if (    this->hybrid_method != HybridMethod::HybridIIc
                     && this->idc_type == IDCType::Update
                )
                    this->c_stages[i] = new OrdinateFlux( *this, this->c_ordinate_set, DG_degree );
            }
        } break;

        case HybridMethod::HybridIIa:
        {
            std::string error_message = "Hybrid-IIa LS-IDC integrators not implemented.\n";

            PRINT_ERROR( error_message.c_str() )
            throw std::invalid_argument( error_message );
        } break;

        default:
        {   std::string error_message = "Invalid hybrid splitting '"
                                        + HybridMethod_to_String.at( this->hybrid_method )
                                        + "' for initializing LS-IDC integrator in '" + std::string(__func__)
                                        + "'.\n";

            PRINT_ERROR( error_message.c_str() )
            throw std::invalid_argument( error_message );
        }
    }
}


//============================================================================================================
//=== PROTECTED HELPER ROUTINES ==============================================================================
//============================================================================================================


//------------------------------------------------------------------------------------------------------------
//! \brief  For LS-IDC integrators, this routine performs no action whatsoever.
//!
//! This routine is used by the non-hybrid and hybrid-II integrators.
//!
//! \see    LSIDCIntegrator::ComputeSource()
//! \see    LSIDCIntegrator::u_ComputeSource_HybridI()
//! \see    LSIDCIntegrator::c_ComputeSource_HybridI()
//------------------------------------------------------------------------------------------------------------
void LSIDCIntegrator::ComputeResiduals ( void ) {}


//------------------------------------------------------------------------------------------------------------
//! \brief  Computes fictitious sources for non-hybrid corrections.
//!
//! This routine is used by the non-hybrid integrators.
//!
//! \param[in]      source      Source term of the transport system. <br>
//!                             The source term is assumed to be constant across the timestep interval.
//! \param[in]      dt          IDC timestep size.
//! \param[in]      n           Index denoting the time interval to integrate across.
//!
//! \see    LSIDCIntegrator::ComputeResiduals()
//------------------------------------------------------------------------------------------------------------
void LSIDCIntegrator::ComputeSource (

    const OrdinateFlux & source,
    const double dt,
    const int64_t n
) {

    if ( n <= 0 ) {  return;  }

    const double h = this->nodes[n] - this->nodes[n-1];
    OrdinateFlux * idc_source = nullptr;
    this->phi->ZeroDensity();

    // \sum_{\ell=1}^{n-1} \frac{ \gamma_{n,\ell} }{ h_n } P \psi^{\ell,[p]}
    // + \sum_{\ell=n+1}^N \frac{ \gamma_{n,\ell} }{ h_n } P \psi^{\ell,[p-1]}
    for ( int64_t l = ((int)(!this->left_endpoint)); l < this->num_stages; ++l ) {

        if ( l == n ) {  continue;  }

        TransportOperator::Pmv( Weight(n,l) / h, 1.0, *this->stages[l], *this->phi );
    }

    switch ( this->idc_type ) {

        case IDCType::Error:
        {
            // Set location of OrdinateFlux object in which to construct source.
            idc_source = this->error;

            // - \frac{1}{h_n \Delta t} \psi^{n,[p-1]}
            OrdinateFlux::AXPY( -1.0/(h * dt), 0.0, *this->stages[n], *idc_source );

            if ( Weight(n,n) != 0.0 ) {

                // \frac{ \gamma_{n,n} }{ h_n } P \psi^{n,[p-1]}
                TransportOperator::Pmv( Weight(n,n) / h, 1.0, *this->stages[n], *this->phi );

                // - \frac{ \gamma_{n,n} }{ h_n } L \psi^{n,[p-1]}
                TransportOperator::Lmv( -Weight(n,n) / h, 1.0,
                                        *this->sigma_t, *this->stages[n], *idc_source );
            }
        } break;

        case IDCType::Update:
        {
            // Set location of OrdinateFlux object in which to construct source.
            idc_source = this->stages[n];

            // ( \frac{ \gamma_{n,n} }{ h_n } - 1.0 ) P \psi^{n,[p-1]}
            TransportOperator::Pmv( Weight(n,n) / h - 1.0, 1.0, *this->stages[n], *this->phi );

            // ( 1.0 - \frac{ \gamma_{n,n} }{ h_n } ) L \psi^{n,[p-1]}
            TransportOperator::Lmv( 1.0 - Weight(n,n) / h, 0.0,
                                    *this->sigma_t, *this->stages[n], *idc_source );
        } break;

        default:
        {   std::string error_message = "Invalid IDCType '"
                                        + IDCIntegrator::IDCType_to_String.at( this->idc_type )
                                        + "' in '" + std::string(__func__) + "'.\n";

            PRINT_ERROR( error_message.c_str() )
            throw std::invalid_argument( error_message );
        }
    }

    // \frac{1}{h_n \Delta t} \psi^{n-1,[p]}
    OrdinateFlux::AXPY( 1.0/(h * dt), 1.0, *this->stages[n-1], *idc_source );

    // S \sum_{\ell=1}^N ( \delta_{n,\ell} - \frac{ \gamma_{n,\ell} }{ h_n } ) P \psi^{\ell,[p??]}
    TransportOperator::Smv( 1.0, 1.0, *this->sigma_s, *this->phi, *idc_source );

    // - \sum_{\ell=1}^{n-1} \frac{ \gamma_{n,\ell} }{ h_n } L \psi^{\ell,[p]}
    // - \sum_{\ell=n+1}^N \frac{ \gamma_{n,\ell} }{ h_n } L \psi^{\ell,[p-1]}
    for ( int64_t l = ((int)(!this->left_endpoint)); l < this->num_stages; ++l ) {

        if ( l == n ) {  continue;  }

        TransportOperator::Lmv( -Weight(n,l) / h, 1.0, *this->sigma_t, *this->stages[l], *idc_source );
    }

    if ( this->idc_type == IDCType::Update )
        idc_source->ZeroDensity( OpDomain::Boundary );

    // Include external source (assumed to be constant across timestep interval).
    OrdinateFlux::AXPY( 1.0, 1.0, source, *idc_source );

    // Boundary condition for error is always zero.
    if ( this->idc_type == IDCType::Error )
        idc_source->ZeroDensity( OpDomain::Boundary );
}


//------------------------------------------------------------------------------------------------------------
//! \brief  For the LS-IDC integrators, this routine performs no action whatsoever.
//!
//! This routine is used by the hybrid-I integrators.
//!
//! \see    LSIDCIntegrator::u_ComputeSource_HybridI()
//! \see    LSIDCIntegrator::c_ComputeSource_HybridI()
//------------------------------------------------------------------------------------------------------------
void LSIDCIntegrator::ComputeResiduals_HybridI ( void ) {}


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
//! \see    LSIDCIntegrator::ComputeResiduals_HybridI()
//! \see    LSIDCIntegrator::c_ComputeSource_HybridI()
//------------------------------------------------------------------------------------------------------------
void LSIDCIntegrator::u_ComputeSource_HybridI (

    const OrdinateFlux & source,
    const double dt,
    const int64_t n
) {

    if ( n <= 0 ) {  return;  }

    const double h = this->nodes[n] - this->nodes[n-1];
    OrdinateFlux * idc_source = nullptr;
    this->phi->ZeroDensity();

    switch ( this->idc_type ) {

        case IDCType::Error:
        {
            // Set location of OrdinateFlux object in which to construct source.
            idc_source = this->u_error;

            // - \frac{1}{h_n \Delta t} \psi_u^{n,[p-1]}
            OrdinateFlux::AXPY( -1.0/(h * dt), 0.0, *this->u_stages[n], *idc_source );

            if ( Weight(n,n) != 0 ) {

                // \frac{ \gamma_{n,n} }{ h_n } P \psiu^{n,[p-1]}
                //
                // NOTE: Value stored in phi is used by c_ComputeSource_HybridI for collided correction.
                //
                Global::TMR_uncollided.Stop();
                Global::TMR_collided.Start();

                TransportOperator::Pmv( Weight(n,n) / h, 1.0, *this->u_stages[n], *this->phi );

                Global::TMR_collided.Stop();
                Global::TMR_uncollided.Start();

                // - \frac{ \gamma_{n,n} }{ h_n } L \psi_u^{n,[p-1]}
                TransportOperator::Lmv( -Weight(n,n) / h, 1.0,
                                        *this->sigma_t, *this->u_stages[n], *idc_source );
            }
        } break;

        case IDCType::Update:
        {
            // Set location of OrdinateFlux object in which to construct source.
            idc_source = this->u_stages[n];

            // ( \frac{ \gamma_{n,n} }{ h_n } - 1.0 ) P \psi_u^{n,[p-1]}
            //
            // NOTE: Value stored in phi is used by c_ComputeSource_HybridI for collided correction.
            //
            Global::TMR_uncollided.Stop();
            Global::TMR_collided.Start();

            TransportOperator::Pmv( Weight(n,n) / h - 1.0, 1.0, *this->u_stages[n], *this->phi );

            Global::TMR_collided.Stop();
            Global::TMR_uncollided.Start();

            // ( 1.0 - \frac{ \gamma_{n,n} }{ h_n } ) L \psi_u^{n,[p-1]}
            TransportOperator::Lmv( 1.0 - Weight(n,n) / h, 0.0,
                                    *this->sigma_t, *this->u_stages[n], *idc_source );
        } break;

        default:
        {   std::string error_message = "Invalid IDCType '"
                                        + IDCIntegrator::IDCType_to_String.at( this->idc_type )
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

    // - \sum_{\ell=1}^{n-1} \frac{ \gamma_{n,\ell} }{ h_n } L \psi_u^{\ell,[p]}
    // - \sum_{\ell=n+1}^N \frac{ \gamma_{n,\ell} }{ h_n } L \psi_u^{\ell,[p-1]}
    for ( int64_t l = ((int)(!this->left_endpoint)); l < this->num_stages; ++l ) {

        if ( l == n ) {  continue;  }

        TransportOperator::Lmv( -Weight(n,l) / h, 1.0, *this->sigma_t, *this->u_stages[l], *idc_source );
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
//! \see    LSIDCIntegrator::ComputeResiduals_HybridI()
//! \see    LSIDCIntegrator::u_ComputeSource_HybridI()
//------------------------------------------------------------------------------------------------------------
void LSIDCIntegrator::c_ComputeSource_HybridI (

    const double dt,
    const int64_t n
) {

    if ( n <= 0 ) {  return;  }

    const double h = this->nodes[n] - this->nodes[n-1];
    OrdinateFlux * idc_source = nullptr;
    /* Do not zero this->phi here! Contains values from u_ComputeSource_HybridI. */

    // \sum_{\ell=1}^{n-1} \frac{ \gamma_{n,\ell} }{ h_n } ( P \psi_c^{\ell,[p]} + P \psi_u^{\ell,[p]} )
    // + \sum_{\ell=n+1}^N \frac{ \gamma_{n,\ell} }{ h_n } ( P \psi_c^{\ell,[p-1]} + P \psi_u^{\ell,[p-1]} )
    for ( int64_t l = ((int)(!this->left_endpoint)); l < this->num_stages; ++l ) {

        if ( l == n ) {  continue;  }

        TransportOperator::Pmv( Weight(n,l) / h, 1.0, *this->c_stages[l], *this->phi );
        TransportOperator::Pmv( Weight(n,l) / h, 1.0, *this->u_stages[l], *this->phi );
    }

    switch ( this->idc_type ) {

        case IDCType::Error:
        {
            // Set location of OrdinateFlux object in which to construct source.
            idc_source = this->c_error;

            // - \frac{1}{h_n \Delta t} \psi_c^{n,[p-1]}
            OrdinateFlux::AXPY( -1.0/(h * dt), 0.0, *this->c_stages[n], *idc_source );

            if ( Weight(n,n) != 0.0 ) {

                // \frac{ \gamma_{n,n} }{ h_n } P \psi_c^{n,[p-1]}
                TransportOperator::Pmv( Weight(n,n) / h, 1.0, *this->c_stages[n], *this->phi );

                // - \frac{ \gamma_{n,n} }{ h_n } L \psi_c^{n,[p-1]}
                TransportOperator::Lmv( -Weight(n,n) / h, 1.0,
                                        *this->sigma_t, *this->c_stages[n], *idc_source );
            }

            // P e_u^{n,[p-1]}
            TransportOperator::Pmv( 1.0, 1.0, *this->u_error, *this->phi );
        } break;

        case IDCType::Update:
        {
            // Set location of OrdinateFlux object in which to construct source.
            idc_source = this->c_stages[n];

            // ( \frac{ \gamma_{n,n} }{ h_n } - 1.0 ) P \psi_c^{n,[p-1]}
            TransportOperator::Pmv( Weight(n,n) / h - 1.0, 1.0, *this->c_stages[n], *this->phi );

            // ( 1.0 - \frac{ \gamma_{n,n} }{ h_n } ) L \psi_c^{n,[p-1]}
            TransportOperator::Lmv( 1.0 - Weight(n,n) / h, 0.0,
                                    *this->sigma_t, *this->c_stages[n], *idc_source );

            // P \psi_u^{n,[p]}
            TransportOperator::Pmv( 1.0, 1.0, *this->u_stages[n], *this->phi );
        } break;

        default:
        {   std::string error_message = "Invalid IDCType '"
                                        + IDCIntegrator::IDCType_to_String.at( this->idc_type )
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

    // S (...)
    this->phi->ZeroDensity( OpDomain::Boundary );
    TransportOperator::Smv( 1.0, 1.0, *this->sigma_s, *this->phi, *idc_source );

    // - \sum_{\ell=1}^{n-1} \frac{ \gamma_{n,\ell} }{ h_n } L \psi_c^{\ell,[p]}
    // - \sum_{\ell=n+1}^N \frac{ \gamma_{n,\ell} }{ h_n } L \psi_c^{\ell,[p-1]}
    for ( int64_t l = ((int)(!this->left_endpoint)); l < this->num_stages; ++l ) {

        if ( l == n ) {  continue;  }

        TransportOperator::Lmv( -Weight(n,l) / h, 1.0, *this->sigma_t, *this->c_stages[l], *idc_source );
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
//! \see    LSIDCIntegrator::ComputeResiduals()
//! \see    LSIDCIntegrator::c_ComputeSource_HybridII()
//------------------------------------------------------------------------------------------------------------
void LSIDCIntegrator::u_ComputeSource_HybridII (

    const OrdinateFlux & source,
    const double dt,
    const int64_t n
) {

    if ( n <= 0 ) {  return;  }

    const double h = this->nodes[n] - this->nodes[n-1];
    OrdinateFlux * idc_source = nullptr;
    this->phi[1].ZeroDensity();

    // \sum_{\ell=1}^{n-1} \frac{ \gamma_{n,\ell} }{ h_n } P \psi^{\ell,[p]}
    // + \sum_{\ell=n+1}^N \frac{ \gamma_{n,\ell} }{ h_n } P \psi^{\ell,[p-1]}
    for ( int64_t l = ((int)(!this->left_endpoint)); l < this->num_stages; ++l )
        TransportOperator::Pmv( Weight(n,l) / h, 1.0, *this->u_stages[l], this->phi[1] );

    if (    this->hybrid_method == HybridMethod::HybridIIc
         || this->idc_type == IDCType::Error
    ) {
        // Set location of OrdinateFlux object in which to construct source.
        idc_source = this->u_error;

        // - \frac{ \gamma_{n,n} }{ h_n } L \psi^{n,[p-1]}
        TransportOperator::Lmv( -Weight(n,n) / h, 0.0, *this->sigma_t, *this->u_stages[n], *idc_source );

        // - \frac{1}{h_n \Delta t} \psi^{n,[p-1]}
        OrdinateFlux::AXPY( -1.0/(h * dt), 1.0, *this->u_stages[n], *idc_source );

    } else if ( this->idc_type == IDCType::Update ) {

        // Set location of OrdinateFlux object in which to construct source.
        idc_source = this->u_stages[n];
        /* Do not zero! */

        // - P \psi^{n,[p-1]}
        //
        // NOTE: Value stored in phi[0] is used by IDCIntegrator::c_ComputeSource_HybridII for collided
        //       correction.
        //
        TransportOperator::Pmv( -1.0, 0.0, *this->u_stages[n], this->phi[0] );

        // - \frac{ h_n }{ h_n - \gamma_{n,n} } R \psi_c^{n,[p-1]}
        Global::TMR_uncollided.Stop();
        this->relabel_operator->Relabel( -h / (h - Weight(n,n)), *this->c_stages[n], *idc_source );
        Global::TMR_uncollided.Start();

        // ( 1.0 - \frac{ \gamma_{n,n} }{ h_n } ) L ( \psi^{n,[p-1]} - \frac{ h_n }{ h_n - \gamma_{n,n} } R \psi_c^{n,[p-1]} )
        TransportOperator::Lmv( 1.0 - Weight(n,n) / h, 0.0, *this->sigma_t, *idc_source, *idc_source );

        // - \frac{1}{h_n \Delta t} R \psic^{n,[p-1]}
        Global::TMR_uncollided.Stop();
        this->relabel_operator->Relabel( -1.0/(h * dt), *this->c_stages[n], *idc_source );
        Global::TMR_uncollided.Start();

    } else {
        std::string error_message = "Invalid IDCType '"
                                    + IDCIntegrator::IDCType_to_String.at( this->idc_type )
                                    + "' in '" + std::string(__func__) + "'.\n";

        PRINT_ERROR( error_message.c_str() )
        throw std::invalid_argument( error_message );
    }

    // \frac{1}{h_n \Delta t} \psi^{n-1,[p]}
    OrdinateFlux::AXPY( 1.0/(h * dt), 1.0, *this->u_stages[n-1], *idc_source );

    // S \sum_{\ell=1}^N ( \delta_{n,\ell} - \frac{ \gamma_{n,\ell} }{ h_n } ) P \psi^{\ell,[p??]}
    TransportOperator::Smv( 1.0, 1.0, *this->sigma_s, this->phi[1], *idc_source );

    // - \sum_{\ell=1}^{n-1} \frac{ \gamma_{n,\ell} }{ h_n } L \psi^{\ell,[p]}
    // - \sum_{\ell=n+1}^N \frac{ \gamma_{n,\ell} }{ h_n } L \psi^{\ell,[p-1]}
    for ( int64_t l = ((int)(!this->left_endpoint)); l < this->num_stages; ++l ) {

        if ( l == n ) {  continue;  }

        TransportOperator::Lmv( -Weight(n,l) / h, 1.0, *this->sigma_t, *this->u_stages[l], *idc_source );
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


# if 0
//------------------------------------------------------------------------------------------------------------
// \brief  Computes the IDC augmented source term used to compute an approximation of the error in the
//         collided component of the split time-dependent system at the \pp{n}th quadrature node.
//
// For the IDC quadrature node stored in index \pp{n}, this routine computes the source term for solving the
// collided error equation of the split time-dependent system.
//
// This routine is used by the hybrid-II integrators.
//
// \param[in]      dt          IDC timestep size.
// \param[in]      n           Index denoting the time interval to integrate across.
//
// \see    LSIDCIntegrator::ComputeResiduals()
// \see    LSIDCIntegrator::u_ComputeSource_HybridII()
//------------------------------------------------------------------------------------------------------------
// void LSIDCIntegrator::c_ComputeSource_HybridII (
//
//     const double dt,
//     const int64_t n
// ) {
//
//     if ( n <= 0 ) {  return;  }
//
//     const double h = this->nodes[n] - this->nodes[n-1];
// }
# endif // if 0


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
//! \see    LSIDCIntegrator::u_ComputeSource_HybridII()
//! \see    IDCIntegrator::c_ComputeSource_HybridII()
//------------------------------------------------------------------------------------------------------------
void LSIDCIntegrator::ComputeNystromSource (

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

    // \sum_{\ell=1}^{n-1} \frac{ \gamma_{n,\ell} }{ h_n } P \psi^{\ell,[p]}
    // + \sum_{\ell=n+1}^N \frac{ \gamma_{n,\ell} }{ h_n } P \psi^{\ell,[p-1]}
    for ( int64_t l = ((int)(!this->left_endpoint)); l < this->num_stages; ++l )
        TransportOperator::Pmv( Weight(n,l) / h, 1.0, *this->u_stages[l], *this->phi );

    switch ( this->idc_type ) {

        case IDCType::Error:
        {
            // Set location of OrdinateFlux object in which to construct source.
            nystrom_source = this->u_error;

            // \frac{1}{h_n \Delta t} psi^{n,[p-1]}
            OrdinateFlux::AXPY( -1.0/(h * dt), 0.0, *this->u_stages[n], *nystrom_source );

            // \frac{ \gamma_{n,n} }{ h_n } L psi^{n,[p-1]}
            TransportOperator::Lmv( -Weight(n,n) / h, 1.0, *this->sigma_t, *this->u_stages[n],
                                    *nystrom_source );
        } break;

        case IDCType::Update:
        {
            // Set location of OrdinateFlux object in which to construct source.
            nystrom_source = this->u_stages[n];

            // ( 1.0 - \frac{ \gamma_{n,n} }{ h_n } ) L psi^{n,[p-1]}
            TransportOperator::Lmv( 1.0 - Weight(n,n) / h, 0.0, *this->sigma_t, *this->u_stages[n],
                                    *nystrom_source );
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

    // - \sum_{\ell=1}^{n-1} \frac{ \gamma_{n,\ell} }{ h_n } L \psi^{\ell,[p]}
    // - \sum_{\ell=n+1}^N \frac{ \gamma_{n,\ell} }{ h_n } L \psi^{\ell,[p-1]}
    for ( int64_t l = ((int)(!this->left_endpoint)); l < this->num_stages; ++l ) {

        if ( l == n ) {  continue;  }

        TransportOperator::Lmv( -Weight(n,l) / h, 1.0, *this->sigma_t, *this->u_stages[l], *nystrom_source );
    }

    if ( this->idc_type == IDCType::Update )
        nystrom_source->ZeroDensity( OpDomain::Boundary );

    // Include external source (assumed to be constant across timestep interval).
    OrdinateFlux::AXPY( 1.0, 1.0, source, *nystrom_source );

    // Boundary condition for the error is always zero.
    if ( this->idc_type == IDCType::Error )
        nystrom_source->ZeroDensity( OpDomain::Boundary );
}
