//------------------------------------------------------------------------------------------------------------
//! \file   time_integrators/RKDG/idc.hpp
//! \brief  Header file for implementation of IDC integrators using implicit Euler steps.
//!
//! \author Michael M. Crockatt
//! \date   February 2018
//------------------------------------------------------------------------------------------------------------

# ifndef __RKDG__IDC_INTEGRATOR_HPP__
# define __RKDG__IDC_INTEGRATOR_HPP__


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
# include "utils/Quadrule/Quadrule.hpp"


//------------------------------------------------------------------------------------------------------------
//! \brief  Class for IDC time integration schemes.
//!
//!
//------------------------------------------------------------------------------------------------------------
class IDCIntegrator : public Abstract::TimeIntegrator {

public:

    //--------------------------------------------------------------------------------------------------------
    //! \brief  Specifies the type of IDC method.
    //--------------------------------------------------------------------------------------------------------
    enum class IDCType {

        //!
        //! String descriptor: \e "none"
        //!
        //! Null value. Not a valid type of IDC method.
        //!
        None,

        //!
        //! String descriptor: \e "error"
        //!
        //! IDC method in which the error equation is solved directly, and then used to update the previous
        //! provisional solutions.
        //!
        Error,

        //!
        //! String descriptor: \e "update"
        //!
        //! IDC method in which a change of variable is applied to the error equation to yield an equation
        //! from which the corrected provisional solutions can be obtained directly.
        //!
        Update
    };


    //--------------------------------------------------------------------------------------------------------
    //! \brief  Specifies the type of hybrid splitting (if any) used by the IDC integrator.
    //--------------------------------------------------------------------------------------------------------
    enum class HybridMethod {

        //!
        //! String descriptor: \e "none"
        //!
        //! No hybrid splitting used. The transport system is solved using a standard discrete ordinates
        //! approximation.
        //!
        //! \see    IDCIntegrator::Step_Nonhybrid()
        //!
        None,

        //!
        //! String descriptor: \e "hybrid-ia"
        //!
        //! Applies a hybrid method derived by splitting the transport equation at the continuum level.
        //! Relabeling of the collided flux is (optionally) performed at the end of each timestep.
        //!
        //! This is the only hybrid method for IDCIntegrator objects that allows relabeling to be omitted.
        //!
        //! \see    IDCIntegrator::Step_HybridI()
        //!
        HybridIa,

        //!
        //! String descriptor: \e "hybrid-ib"
        //!
        //! Applies a hybrid method derived by splitting the transport equation at the continuum level.
        //! Relabeling of the collided flux is performed after each implicit Euler substep.
        //!
        //! \see    IDCIntegrator::Step_HybridI()
        //!
        HybridIb,

        //!
        //! String descriptor: \e "hybrid-iia"
        //!
        //! Applies a hybrid method derived by splitting the error equation used to perform correction
        //! iterations.
        //! Relabeling is applied between each correction iteration.
        //!
        //! \see    IDCIntegrator::Step_HybridII()
        //!
        HybridIIa,

        //!
        //! String descriptor: \e "hybrid-iib"
        //!
        //! Applies a hybrid method derived by splitting the error equation used to perform correction
        //! iterations.
        //! Relabeling is applied between each substep of each correction iteration.
        //!
        //! \see    IDCIntegrator::Step_HybridII()
        //!
        HybridIIb,

        //!
        //! String descriptor: \e "hybrid-iic"
        //!
        //! Applies a hybrid method derived by splitting the error equation used to perform correction
        //! iterations.
        //! Relabeling is applied between each substep of each correction iteration using a Nyström
        //! interpolation procedure.
        //!
        //! \see    IDCIntegrator::Step_HybridII()
        //!
        HybridIIc
    };


    //--------------------------------------------------------------------------------------------------------
    //! \brief  Used to convert IDCKIntegrator::IDCType values to their corresponding string descriptors.
    //--------------------------------------------------------------------------------------------------------
    static std::map< IDCIntegrator::IDCType, std::string > IDCType_to_String;


    //--------------------------------------------------------------------------------------------------------
    //! \brief  Used to convert string descriptors to IDCIntegrator::IDCType values.
    //--------------------------------------------------------------------------------------------------------
    static std::map< std::string, IDCIntegrator::IDCType > String_to_IDCType;


    //--------------------------------------------------------------------------------------------------------
    //! \brief  Used to convert IDCIntegrator::HybridMethod values to their corresponding string descriptors.
    //--------------------------------------------------------------------------------------------------------
    static std::map< IDCIntegrator::HybridMethod, std::string > HybridMethod_to_String;


    //--------------------------------------------------------------------------------------------------------
    //! \brief  Used to convert string descriptors to IDCIntegrator::HybridMethod values.
    //--------------------------------------------------------------------------------------------------------
    static std::map< std::string, IDCIntegrator::HybridMethod > String_to_HybridMethod;


    //========================================================================================================
    //=== CONSTRUCTORS AND DESTRUCTOR ========================================================================
    //========================================================================================================


    //--------------------------------------------------------------------------------------------------------
    //! \brief  Delete default constructor.
    //--------------------------------------------------------------------------------------------------------
    IDCIntegrator( void ) = delete;


    //--------------------------------------------------------------------------------------------------------
    //! \brief  Delete copy constructor.
    //--------------------------------------------------------------------------------------------------------
    IDCIntegrator( const IDCIntegrator & ) = delete;


    //--------------------------------------------------------------------------------------------------------
    //! \brief  Creates an IDCIntegrator object using the specified parameters and returns a pointer to the
    //!         object.
    //--------------------------------------------------------------------------------------------------------
    static IDCIntegrator * Create( const ParameterList & );


    //--------------------------------------------------------------------------------------------------------
    //! \brief  Destructor for class IDCIntegrator.
    //--------------------------------------------------------------------------------------------------------
    ~IDCIntegrator( void ) override;


    //========================================================================================================
    //=== OPERATOR OVERLOADS =================================================================================
    //========================================================================================================


    //--------------------------------------------------------------------------------------------------------
    //! \brief  Delete copy assignment operator.
    //--------------------------------------------------------------------------------------------------------
    IDCIntegrator & operator=( const IDCIntegrator & ) = delete;


    //========================================================================================================
    //=== INTERFACE ROUTINES =================================================================================
    //========================================================================================================


    //--------------------------------------------------------------------------------------------------------
    //! \brief  Interface routine for computing a single timestep with step size \pp{step_size}.
    //--------------------------------------------------------------------------------------------------------
    IDCIntegrator & Step (

        RKDG::OrdinateFlux & initial_condition,
        const RKDG::OrdinateFlux & source,
        const RKDG::CrossSection & sigma_t_in,
        const RKDG::CrossSection & sigma_s_in,
        const double step_size

    ) override;


    //--------------------------------------------------------------------------------------------------------
    //! \brief  Interface routine for computing a single timestep with step size \pp{step_size}.
    //--------------------------------------------------------------------------------------------------------
    IDCIntegrator & Step (

        RKDG::OrdinateFlux & u_initial_condition,
        RKDG::OrdinateFlux * c_initial_condition,
        const RKDG::OrdinateFlux & source,
        const RKDG::CrossSection & sigma_t_in,
        const RKDG::CrossSection & sigma_s_in,
        const double step_size

    ) override;


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
    TimeIntegratorType GetIntegratorType( void ) const override {  return TimeIntegratorType::IDC;  }


    //--------------------------------------------------------------------------------------------------------
    //! \brief  Prints the table of quadrature weights used by the IDC integrator to the logging interface.
    //--------------------------------------------------------------------------------------------------------
    void PrintQuadrature( void ) const;


protected:

    //========================================================================================================
    //=== PROTECTED MEMBER VARIABLES =========================================================================
    //========================================================================================================

    IDCType idc_type;                   //!< Specifies the type of IDC method used.
    Quadrule::NodesType nodes_type;     //!< Specifies the type of quadrature nodes used by the IDC method.
    double * nodes;                     //!< Pointer to array of IDC quadrature nodes.

    //!
    //! \brief  Pointer to array of IDC quadrature weights.
    //!
    //! Quadrature weights are stored in an integration matrix in row-major format, and can be accessed using
    //! IDCIntegrator::Weight. The quadrature weights for integrating across the entire timestep interval are
    //! stored in the first row of the matrix (with index zero).
    //!
    //! \see    IDCIntegrator::Weight()
    //!
    double * weights;

    bool left_endpoint;                 //!< Specifies whether or not the left endpoint of the timestep interval is a quadrature node.
    bool right_endpoint;                //!< Specifies whether or not the right endpoint of the timestep interval is a quadrature node.

    //!
    //! \brief  Specifies whether or not the collocation formula should be used to construct output values.
    //!
    //! Defaults to true.
    //!
    bool use_collocation;

    int64_t num_stages;                 //!< Number of stages in each correction level of the IDC scheme.
    int64_t num_corrections;            //!< Number of IDC correction iterations (in addition to implicit Euler prediction).

    int64_t DG_degree;

    RKDG::ScalarFlux * phi;             //!< Temporary storage for scalar fluxes.

    const RKDG::CrossSection * sigma_s;
    const RKDG::CrossSection * sigma_t;

    // Nonhybrid fluxes and solvers.
    Quadrule::OrdinateSet ordinate_set; //!< Nonhybrid ordinate set.

    RKDG::OrdinateFlux ** stages;       //!< Angular flux stage values for current correction level.
    RKDG::OrdinateFlux ** residuals;    //!< Residuals for stages values at previous correction level.
    RKDG::OrdinateFlux * error;         //!< Used as temporary storage for computing the current error.

    //! Solver for nonhybrid systems.
    std::shared_ptr<ImplicitSolver<RKDG::OrdinateFlux>> solver;

    // Variables for hybrid methods.
    HybridMethod hybrid_method;         //!< Specifies the type of hybrid splitting (if any) to use.

    Quadrule::OrdinateSet u_ordinate_set;  //!< Uncollided OrdinateSet.

    RKDG::OrdinateFlux ** u_stages;     //!< Uncollided flux stage values for current correction level.
    RKDG::OrdinateFlux ** u_residuals;  //!< Residuals for uncollided stage values at previous correction level.
    RKDG::OrdinateFlux * u_error;       //!< Used as temporary storage for computing the current error in the uncollided flux.

    //! Solver for uncollided systems.
    std::shared_ptr<ImplicitSolver<RKDG::OrdinateFlux>> u_solver;

    Quadrule::OrdinateSet c_ordinate_set;  //!< Collided OrdinateSet.

    RKDG::OrdinateFlux ** c_stages;     //!< Collided flux stage values for current correction level.
    RKDG::OrdinateFlux ** c_residuals;  //!< Residuals for collided stage values at previous correction level.
    RKDG::OrdinateFlux * c_error;       //!< Used as temporary storage for computing the current error in the collided flux.

    //! Solver for collided systems.
    std::shared_ptr<ImplicitSolver<RKDG::OrdinateFlux>> c_solver;

    bool relabel = true;                //!< Specifies whether relabeling should be performed or not. Defaults to true.
    RelabelOperator * relabel_operator; //!< RelabelOperator for mapping collided fluxes to uncollided.


    //========================================================================================================
    //=== PROTECTED CONSTRUCTOR ==============================================================================
    //========================================================================================================


    //--------------------------------------------------------------------------------------------------------
    //! \brief  Constructs an IDCIntegrator object using the provided values.
    //--------------------------------------------------------------------------------------------------------
    IDCIntegrator( const ParameterList & );


    //========================================================================================================
    //=== PROTECTED HELPER ROUTINES ==========================================================================
    //========================================================================================================


    //--------------------------------------------------------------------------------------------------------
    //! \brief  Indexing function for matrix of IDC integration weights.
    //!
    //! \param[in]      i       Subinterval to integrate over.
    //! \param[in]      j       Index of node corresponding to weight.
    //!
    //! \return     Returns a reference to the weight for the \pp{j}th quadrature node used for integrating
    //!             across the \pp{i}th subinterval.
    //--------------------------------------------------------------------------------------------------------
    inline double & Weight( const int64_t i, const int64_t j ) const {

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

    # endif // if defined (STRICT_CHECK)

        return this->weights[ j + this->num_stages*i ];
    }


    //--------------------------------------------------------------------------------------------------------
    //! \brief  Allocates memory for the internal stage values.
    //--------------------------------------------------------------------------------------------------------
    virtual void Allocate( void );


    //--------------------------------------------------------------------------------------------------------
    //! \brief  Sets up IDC quadrature nodes and weights.
    //--------------------------------------------------------------------------------------------------------
    void SetNodes( void );


    //--------------------------------------------------------------------------------------------------------
    //! \brief  Computes the final output value of the IDC method by applying the collocation formula to the
    //!         substep values of the non-hybrid or reconstructed hybrid approximation.
    //--------------------------------------------------------------------------------------------------------
    void ComputeCollocation_Nonhybrid( const RKDG::OrdinateFlux & initial_condition,
                                       const RKDG::OrdinateFlux & source, const double dt );


    //--------------------------------------------------------------------------------------------------------
    //! \brief  Computes the final output value of the IDC method by applying the collocation formula to the
    //!         substep values of the non-hybrid or reconstructed hybrid approximation.
    //--------------------------------------------------------------------------------------------------------
    void ComputeCollocation_Hybrid( const RKDG::OrdinateFlux & u_initial_condition,
                                    const RKDG::OrdinateFlux * const c_initial_condition,
                                    const RKDG::OrdinateFlux & source, const double dt );


    //--------------------------------------------------------------------------------------------------------
    //! \brief  Computes residuals for the non-hybrid operator.
    //--------------------------------------------------------------------------------------------------------
    virtual void ComputeResiduals( void );


    //--------------------------------------------------------------------------------------------------------
    //! \brief  Computes fictitious sources for non-hybrid corrections.
    //--------------------------------------------------------------------------------------------------------
    virtual void ComputeSource( const RKDG::OrdinateFlux & source, const double dt, const int64_t n );


    //--------------------------------------------------------------------------------------------------------
    //! \brief  Computes residuals for the hybrid-I operator.
    //--------------------------------------------------------------------------------------------------------
    virtual void ComputeResiduals_HybridI( void );


    //--------------------------------------------------------------------------------------------------------
    //! \brief  Computes fictitious sources for the uncollided portion of hybrid-I corrections.
    //--------------------------------------------------------------------------------------------------------
    virtual void u_ComputeSource_HybridI( const RKDG::OrdinateFlux & source,
                                          const double dt, const int64_t n );


    //--------------------------------------------------------------------------------------------------------
    //! \brief  Computes fictitious sources for the collided portion of hybrid-I corrections.
    //--------------------------------------------------------------------------------------------------------
    virtual void c_ComputeSource_HybridI( const double dt, const int64_t n );


    //--------------------------------------------------------------------------------------------------------
    //! \brief  Computes fictitious sources for the uncollided portion of hybrid-II corrections.
    //--------------------------------------------------------------------------------------------------------
    virtual void u_ComputeSource_HybridII( const RKDG::OrdinateFlux & source,
                                           const double dt, const int64_t n );


    //--------------------------------------------------------------------------------------------------------
    //! \brief  Computes fictitious sources for the collided portion of hybrid-II corrections.
    //--------------------------------------------------------------------------------------------------------
    virtual void c_ComputeSource_HybridII( const double dt, const int64_t n );


    //--------------------------------------------------------------------------------------------------------
    //! \brief  Computes source terms for the hybrid-IIc Nyström reconstruction.
    //--------------------------------------------------------------------------------------------------------
    virtual void ComputeNystromSource( const RKDG::OrdinateFlux & source, const double dt, const int64_t n );


    //========================================================================================================
    //=== PROTECTED TIME INTEGRATION ROUTINES ================================================================
    //========================================================================================================


    //--------------------------------------------------------------------------------------------------------
    //! \brief  Computes a single timestep of the non-hybrid IDC method.
    //--------------------------------------------------------------------------------------------------------
    void Step_Nonhybrid (

        RKDG::OrdinateFlux & initial_condition,
        const RKDG::OrdinateFlux & source,
        const double dt
    );


    //--------------------------------------------------------------------------------------------------------
    //! \brief  Computes a single timestep of the hybrid-I IDC method.
    //--------------------------------------------------------------------------------------------------------
    void Step_HybridI (

        RKDG::OrdinateFlux & u_initial_condition,
        RKDG::OrdinateFlux * const c_initial_condition,
        const RKDG::OrdinateFlux & source,
        const double dt
    );


    //--------------------------------------------------------------------------------------------------------
    //! \brief  Computes a single timestep of the hybrid-IIa or hybrid-IIb IDC method.
    //--------------------------------------------------------------------------------------------------------
    void Step_HybridII (

        RKDG::OrdinateFlux & initial_condition,
        const RKDG::OrdinateFlux & source,
        const double dt
    );
};


# endif // ifndef __RKDG__IDC_INTEGRATOR_HPP__
