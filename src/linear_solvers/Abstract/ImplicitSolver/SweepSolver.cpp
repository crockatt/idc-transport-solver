//------------------------------------------------------------------------------------------------------------
//! \file   linear_solvers/Abstract/ImplicitSolver/SweepSolver.cpp
//! \brief  Contains implementations and instantiations of methods from SweepSolver class template.
//!
//! \author Michael M. Crockatt
//! \date   March 2018
//------------------------------------------------------------------------------------------------------------

# include <typeinfo>

# include "linear_solvers/Abstract/ImplicitSolver/SweepSolver.hpp"


//============================================================================================================
//=== MEMBER DEFINITIONS FOR SweepSolver CLASS ===============================================================
//============================================================================================================


//--------------------------------------------------------------------------------------------------------
//! \brief  Returns the string descriptor of the class type.
//--------------------------------------------------------------------------------------------------------
template< class AngularFlux >
const std::string SweepSolver<AngularFlux>::Descriptor( void ) {

    return "sweep";
}


//--------------------------------------------------------------------------------------------------------
//! \brief  Returns the string descriptor for an objects type.
//--------------------------------------------------------------------------------------------------------
template< class AngularFlux >
const std::string SweepSolver<AngularFlux>::GetDescriptor( void ) const {

    return Descriptor();
}


//------------------------------------------------------------------------------------------------------------
//! \brief  Initializes a SweepSolver object.
//!
//! \param[in]  domain_decomposition    Contains parameters of spatial discretization.
//! \param[in]  ordinate_set            Contains parameters of angular discretization.
//! \param[in]  input_list              Contains additional parameters.
//! \param[in]  dt                      Initial timestep size.
//------------------------------------------------------------------------------------------------------------
template< class AngularFlux >
SweepSolver<AngularFlux>::SweepSolver (

    const DomainDecomposition & domain_decomposition,
    const Quadrule::OrdinateSet & ordinate_set,
    const ParameterList & input_list,
    const double dt                 // = std::numeric_limits<double>::infinity()
) :
    ImplicitSolver<AngularFlux>( domain_decomposition, ordinate_set, input_list, dt )
{
    PRINT_STATUS( "Executing SweepSolver<%s>::%s.\n", typeid(AngularFlux).name(), __func__ )

    this->sweep_operator = std::shared_ptr<SweepOperator<AngularFlux>>(
            new SweepOperator<AngularFlux>( domain_decomposition, ordinate_set, input_list, dt )
        );
}


//------------------------------------------------------------------------------------------------------------
//! \brief  Virtual destructor.
//------------------------------------------------------------------------------------------------------------
template< class AngularFlux >
SweepSolver<AngularFlux>::~SweepSolver ( void ) {

    PRINT_STATUS( "Executing SweepSolver<%s>::%s.\n", typeid(AngularFlux).name(), __func__ )

    /* empty */
}


//------------------------------------------------------------------------------------------------------------
//! \brief  Prints a summary of the object's configuration to the logging interface.
//------------------------------------------------------------------------------------------------------------
template< class AngularFlux >
void SweepSolver<AngularFlux>::Print (

    const std::string & prefix // = "  "

) const {

    PRINT_STATUS( "Executing SweepSolver<%s>::%s.\n", typeid(AngularFlux).name(), __func__ )

    this->ImplicitSolver<AngularFlux>::Print( prefix );

    PRINT_LOG( "%s%-*s  %s\n", prefix.c_str(), Global::col_width, "Solve Type:", this->Descriptor().c_str() )

    this->sweep_operator->Print( prefix );
}


//------------------------------------------------------------------------------------------------------------
//! \brief  Sets parameters given a ParameterList.
//!
//! \param[in]  input_list  Contains parameters to set in the object.
//------------------------------------------------------------------------------------------------------------
template< class AngularFlux >
void SweepSolver<AngularFlux>::SetParameters (

    const ParameterList & input_list
) {
    PRINT_STATUS( "Executing SweepSolver<%s>::%s.\n", typeid(AngularFlux).name(), __func__ )

    ImplicitSolver<AngularFlux>::SetParameters( input_list );

    this->sweep_operator->SetParameters( input_list );
}


//------------------------------------------------------------------------------------------------------------
//! \brief  Solves a system of equations using the specified parameters.
//!
//! \param[in,out]  source      Initially contains the coefficients of the source term \f$ Q \f$ for the
//!                             transport system. Upon return, contains the coefficients \f$ \Psi \f$ that
//!                             satisfy the transport equation.
//! \param[in]      sigma_t     Total cross section \f$ \sigma_{\mathrm{t}} \f$.
//! \param[in]      dt          (optional) <br>
//!                             Timestep size. <br>
//!                             For solving steady-state problems on should use <tt>\pp{dt} = inf</tt>
//!                             (default value).
//! \param[in]      initial     (optional) <br>
//!                             Pointer to object containing the initial condition for the timestep.
//!                             This pointer is not dereferenced if it is null (default value).
//------------------------------------------------------------------------------------------------------------
template< class AngularFlux >
void SweepSolver<AngularFlux>::Solve (

    AngularFlux & source,
    const CrossSection & sigma_t,
    const CrossSection &,
    const double dt,                            // = std::numeric_limits<double>::infinity(),
    const AngularFluxInit * const initial    // = nullptr
) {

    PRINT_STATUS( "Executing SweepSolver<%s>::%s.\n", typeid(AngularFlux).name(), __func__ )

    this->sweep_operator->Apply( source, sigma_t, dt, initial );
}


//============================================================================================================
//=== EXPLICIT INSTANTIATIONS ================================================================================
//============================================================================================================

# include "linear_solvers/Abstract/ImplicitSolver/ImplicitSolverDerivedFactory.hpp"
# include "objects/RKDG/OrdinateFlux.hpp"
# include "objects/STDG/OrdinateFlux.hpp"


//!
//! \brief  Instantiation of SweepSolver<RKDG::OrdinateFlux> class.
//!
template class SweepSolver<RKDG::OrdinateFlux>;

//!
//! \brief  Instantiation of factory for SweepSolver<RKDG::OrdinateFlux>.
//!
template class
ImplicitSolver<RKDG::OrdinateFlux>::ImplicitSolverDerivedFactory<
        SweepSolver<RKDG::OrdinateFlux>
    >;


//!
//! \brief  Instantiation of SweepSolver<STDG::OrdinateFlux> class.
//!
template class SweepSolver<STDG::OrdinateFlux>;

//!
//! \brief  Instantiation of factory for SweepSolver<STDG::OrdinateFlux>.
//!
template class
ImplicitSolver<STDG::OrdinateFlux>::ImplicitSolverDerivedFactory<
        SweepSolver<STDG::OrdinateFlux>
    >;
