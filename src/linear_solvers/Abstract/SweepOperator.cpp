//------------------------------------------------------------------------------------------------------------
//! \file   linear_solvers/Abstract/SweepOperator.cpp
//! \brief  Contains implementations and instantiations of methods from SweepOperator class template.
//!
//! \author Michael M. Crockatt
//! \date   January 2018
//------------------------------------------------------------------------------------------------------------

# include <typeinfo>

# include "linear_solvers/Abstract/SweepOperator.hpp"
# include "linear_solvers/Abstract/SweepPattern/SweepPattern.hpp"
# include "linear_solvers/Abstract/SweepCellSolve/SweepCellSolveManager.hpp"


//------------------------------------------------------------------------------------------------------------
//! \brief  Initializes a SweepOperator object.
//!
//! \param[in]  domain_decomposition    Contains parameters of spatial discretization.
//! \param[in]  ordinate_set            Contains parameters of angular discretization.
//! \param[in]  input_list              Contains additional parameters.
//! \param[in]  dt                      Initial timestep size.
//------------------------------------------------------------------------------------------------------------
template< class OrdinateFlux >
SweepOperator<OrdinateFlux>::SweepOperator (

    const DomainDecomposition & domain_decomposition,
    const Quadrule::OrdinateSet & ordinate_set,
    const ParameterList & input_list,
    const double dt                 // = std::numeric_limits<double>::infinity()
) :
    DomainDecomposition( domain_decomposition ),
    Quadrule::OrdinateSet( ordinate_set )
{

    PRINT_STATUS( "Executing SweepOperator<%s>::%s.\n", typeid(OrdinateFlux).name(), __func__ )

    this->sweep_pattern =
        std::shared_ptr<SweepPattern>(
            SweepPattern::SweepPatternFactory::CreatePattern( *this, input_list )
        );

    this->sweep_pattern->SetParameters( input_list );

    this->sweep_cell_solve =
        std::shared_ptr<SweepCellSolveManager>(
            new SweepCellSolveManager( *this, input_list, dt )
        );
}


//------------------------------------------------------------------------------------------------------------
//! \brief  Virtual destructor.
//------------------------------------------------------------------------------------------------------------
template< class OrdinateFlux >
SweepOperator<OrdinateFlux>::~SweepOperator ( void ) {

    PRINT_STATUS( "Executing SweepOperator<%s>::%s.\n", typeid(OrdinateFlux).name(), __func__ )

    /* empty */
}


//------------------------------------------------------------------------------------------------------------
//! \brief  Prints a summary of the object's configuration to the logging interface.
//!
//! \param[in]  prefix  String prepended to each line of output, intended to be used for controlling
//!                     indentation of output.
//------------------------------------------------------------------------------------------------------------
template< class OrdinateFlux >
void SweepOperator<OrdinateFlux>::Print (

    const std::string & prefix // = "  "

) const {

    PRINT_STATUS( "Executing SweepOperator<%s>::%s.\n", typeid(OrdinateFlux).name(), __func__ )

    this->sweep_pattern->Print( prefix );
    this->sweep_cell_solve->Print( prefix );
}


//------------------------------------------------------------------------------------------------------------
//! \brief  Sets parameters given in a ParameterList.
//!
//! \param[in]  input_list  List of parameters to set.
//------------------------------------------------------------------------------------------------------------
template< class OrdinateFlux >
void SweepOperator<OrdinateFlux>::SetParameters (

    const ParameterList & input_list
) {

    PRINT_STATUS( "Executing SweepOperator<%s>::%s.\n", typeid(OrdinateFlux).name(), __func__ )

    this->sweep_pattern->SetParameters( input_list );
//     this->sweep_cell_solve->SetParameters( input_list );
}


//------------------------------------------------------------------------------------------------------------
//! \brief  Applies the sweep operator (\f$ \mathcal{L}^{-1} \f$).
//!
//! \param[in,out]  source              Initially contains the coefficients of the source term \f$ Q \f$ for
//!                                     the transport system. Upon return, contains the coefficients
//!                                     \f$ \Psi \f$ that satisfy the transport equation.
//! \param[in]      sigma               Total cross section \f$ \sigma_{\mathrm{t}} \f$.
//! \param[in]      dt                  (optional) <br>
//!                                     Timestep size. <br>
//!                                     For solving steady-state problems on should use
//!                                     <tt>\pp{dt} = inf</tt> (default value).
//! \param[in]      initial_condition   (optional) <br>
//!                                     Pointer to object containing the initial condition for the timestep.
//!                                     This pointer is not dereferenced if it is null (default value).
//------------------------------------------------------------------------------------------------------------
template< class OrdinateFlux >
void SweepOperator<OrdinateFlux>::Apply (

    OrdinateFlux & source,
    const RKDG::CrossSection & sigma,
    const double dt,                                    // = std::numeric_limits<double>::infinity()
    const RKDG::OrdinateFlux * const initial_condition  // = nullptr

) const {

    PRINT_STATUS( "Executing SweepOperator<%s>::%s.\n", typeid(OrdinateFlux).name(), __func__ )

    Global::TMR_AF_Lsv.Start();

    this->sweep_cell_solve->SetDt( dt );
    this->sweep_pattern->Sweep( source, source, sigma, initial_condition );

    Global::TMR_AF_Lsv.Stop();
}


//============================================================================================================
//=== EXPLICIT INSTANTIATIONS ================================================================================
//============================================================================================================

# include "objects/RKDG/OrdinateFlux.hpp"
# include "objects/STDG/OrdinateFlux.hpp"


//!
//! \brief  Instantiation of SweepOperator<RKDG::OrdinateFlux> class.
//!
template class SweepOperator<RKDG::OrdinateFlux>;


//!
//! \brief  Instantiation of SweepOperator<STDG::OrdinateFlux> class.
//!
template class SweepOperator<STDG::OrdinateFlux>;
