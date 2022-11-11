//------------------------------------------------------------------------------------------------------------
//! \file   utils/Problems/Problem.hpp
//! \brief  Header for Problem class.
//!
//! \author Michael Crockatt
//! \date   June 2018
//------------------------------------------------------------------------------------------------------------

# ifndef __PROBLEM_HPP__
# define __PROBLEM_HPP__


# include <memory>
# include <string>

# include "objects/RKDG/DensityFunction.hpp"
# include "utils/Mollifier.hpp"
# include "utils/ParameterList.hpp"


// Forward declaration for typedef.
class Problem;


//!
//! \brief  Typedef for evaluation member functions.
//!
typedef double (Problem::*ProblemEvalFcn)(
    const double
# if SPACE_DIMS >= 2
  , const double
# endif
# if SPACE_DIMS == 3
  , const double
# endif
) const;


//------------------------------------------------------------------------------------------------------------
//! \brief  Interface for accessing parameters for test problems.
//!
//! The algorithm that is used to numerically compute mollified problem parameters is controlled by the
//! <i>feature size</i> of each parameter. The feature size for each problem is controlled by the following
//! virtual interface routines:
//!     - \c Problem::InitialConditionFtSize
//!     - \c Problem::SourceFtSize
//!     - \c Problem::TotalCrossFtSize
//!     - \c Problem::ScatterCrossFtSize
//!
//! The feature size for each parameter should be specified as zero (default) unless the function that
//! defines the parameter satisfies all of the following conditions:
//!     1. The parameter is a piecewise-constant function of position.
//!     2. All discontinuities in the function lie on a Cartesian grid aligned with the coordinate axes.
//!
//! If the above conditions are met, then the feature size of the parameter should be defined to be the
//! minimum distance in the \f$ \infty \f$-norm between any two points that lie on different edges.
//------------------------------------------------------------------------------------------------------------
class Problem {

public:

    //========================================================================================================
    //=== CONSTRUCTORS AND DESTRUCTOR ========================================================================
    //========================================================================================================

    //--------------------------------------------------------------------------------------------------------
    //! \brief  Delete default constructor.
    //--------------------------------------------------------------------------------------------------------
    Problem( void ) = delete;

    //--------------------------------------------------------------------------------------------------------
    //! \brief  Delete copy constructor.
    //--------------------------------------------------------------------------------------------------------
    Problem( const Problem & ) = delete;

    //--------------------------------------------------------------------------------------------------------
    //! \brief  Construct object from given parameters.
    //--------------------------------------------------------------------------------------------------------
    Problem( const ParameterList & );

    //--------------------------------------------------------------------------------------------------------
    //! \brief  Destructor for Problem class.
    //--------------------------------------------------------------------------------------------------------
    virtual ~Problem( void );


    //========================================================================================================
    //=== OPERATOR OVERLOADS =================================================================================
    //========================================================================================================

    //--------------------------------------------------------------------------------------------------------
    //! \brief  Delete copy assignment operator.
    //--------------------------------------------------------------------------------------------------------
    Problem & operator=( const Problem & ) = delete;


    //========================================================================================================
    //=== PUBLIC INTERFACE ROUTINES ==========================================================================
    //========================================================================================================

    //--------------------------------------------------------------------------------------------------------
    //! \brief  Returns the key used to search input lists for determining the derived type to construct.
    //--------------------------------------------------------------------------------------------------------
    static std::string GetInputKey( void );

    //--------------------------------------------------------------------------------------------------------
    //! \brief  Returns the strict descriptor for a problem.
    //--------------------------------------------------------------------------------------------------------
    virtual const std::string GetDescriptor( void ) const = 0;

    //--------------------------------------------------------------------------------------------------------
    //! \brief  Prints a summary of the problem configuration to the logging interface.
    //--------------------------------------------------------------------------------------------------------
    virtual void Print( const std::string & = "  " ) const;

    //--------------------------------------------------------------------------------------------------------
    //! \brief  Overrides values in given ParameterList with requirements imposed by test problem.
    //--------------------------------------------------------------------------------------------------------
    virtual void OverrideOptions( ParameterList & ) const {/* empty */}

    //--------------------------------------------------------------------------------------------------------
    //! \brief  Sets inflow boundary conditions for the test problem.
    //--------------------------------------------------------------------------------------------------------
    virtual void SetInflowBoundaries( RKDG::DensityFunction & ) const;

    //--------------------------------------------------------------------------------------------------------
    //! \brief  Returns the value of the initial condition at the specified position.
    //--------------------------------------------------------------------------------------------------------
    virtual double EvalInitialCondition(
            const double x
    # if SPACE_DIMS >= 2 || defined (DOXYCOMPILE)
        ,   const double y
    # endif
    # if SPACE_DIMS == 3
        ,   const double z
    # endif
    ) const;

    //--------------------------------------------------------------------------------------------------------
    //! \brief  Returns the value of the source term at the specified position.
    //--------------------------------------------------------------------------------------------------------
    virtual double EvalSource(
            const double x
    # if SPACE_DIMS >= 2 || defined (DOXYCOMPILE)
        ,   const double y
    # endif
    # if SPACE_DIMS == 3
        ,   const double z
    # endif
    ) const;

    //--------------------------------------------------------------------------------------------------------
    //! \brief  Returns the value of the total cross section at the specified position.
    //--------------------------------------------------------------------------------------------------------
    virtual double EvalTotalCross(
            const double x
    # if SPACE_DIMS >= 2 || defined (DOXYCOMPILE)
        ,   const double y
    # endif
    # if SPACE_DIMS == 3
        ,   const double z
    # endif
    ) const;

    //--------------------------------------------------------------------------------------------------------
    //! \brief  Returns the value of the scattering cross section at the specified position.
    //--------------------------------------------------------------------------------------------------------
    virtual double EvalScatterCross(
            const double x
    # if SPACE_DIMS >= 2 || defined (DOXYCOMPILE)
        ,   const double y
    # endif
    # if SPACE_DIMS == 3
        ,   const double z
    # endif
    ) const;


protected:

    //========================================================================================================
    //=== PROTECTED MEMBER VARIABLES =========================================================================
    //========================================================================================================

    double moll_radius;                 //!< Radius of mollifier.
    std::shared_ptr<Mollifier> moll;    //!< Mollifier object.

    int64_t GK_order;       //!< Order of Gauss-Kronrod quadrature used for computing mollifier integrals.
    double * GK_nodes;      //!< Stores Gauss-Kronrod quadrature nodes.
    double * G_weights;     //!< Stores Gaussian weights for Gauss-Kronrod quadrature.
    double * K_weights;     //!< Stores Kronrod weights for Gauss-Kronrod quadrature.


    //========================================================================================================
    //=== PROTECTED INTERFACE ROUTINES =======================================================================
    //========================================================================================================

    //--------------------------------------------------------------------------------------------------------
    //! \brief  Returns the value of the initial condition without mollification at the specified position.
    //--------------------------------------------------------------------------------------------------------
    virtual double InitialCondition(
        const double x
    # if SPACE_DIMS >= 2 || defined (DOXYCOMPILE)
      , const double y
    # endif
    # if SPACE_DIMS == 3
      , const double z
    # endif
    ) const;

    //--------------------------------------------------------------------------------------------------------
    //! \brief  Returns the value of the source term without mollification at the specified position.
    //--------------------------------------------------------------------------------------------------------
    virtual double Source(
        const double x
    # if SPACE_DIMS >= 2 || defined (DOXYCOMPILE)
      , const double y
    # endif
    # if SPACE_DIMS == 3
      , const double z
    # endif
    ) const;

    //--------------------------------------------------------------------------------------------------------
    //! \brief  Returns the value of the total cross section without mollification at the specified position.
    //--------------------------------------------------------------------------------------------------------
    virtual double TotalCross(
        const double x
    # if SPACE_DIMS >= 2 || defined (DOXYCOMPILE)
      , const double y
    # endif
    # if SPACE_DIMS == 3
      , const double z
    # endif
    ) const;

    //--------------------------------------------------------------------------------------------------------
    //! \brief  Returns the value of the scattering cross section without mollification at the specified
    //!         position.
    //--------------------------------------------------------------------------------------------------------
    virtual double ScatterCross(
        const double x
    # if SPACE_DIMS >= 2 || defined (DOXYCOMPILE)
      , const double y
    # endif
    # if SPACE_DIMS == 3
      , const double z
    # endif
    ) const;

    //--------------------------------------------------------------------------------------------------------
    //! \brief  Returns the feature size of the initial condition.
    //--------------------------------------------------------------------------------------------------------
    virtual double InitialConditionFtSize( void ) const {  return 0.0;  }

    //--------------------------------------------------------------------------------------------------------
    //! \brief  Returns the feature size of the source term.
    //--------------------------------------------------------------------------------------------------------
    virtual double SourceFtSize( void ) const {  return 0.0;  }

    //--------------------------------------------------------------------------------------------------------
    //! \brief  Returns the feature size of the total cross section.
    //--------------------------------------------------------------------------------------------------------
    virtual double TotalCrossFtSize( void ) const {  return 0.0;  }

    //--------------------------------------------------------------------------------------------------------
    //! \brief  Returns the feature size of the scattering cross section.
    //--------------------------------------------------------------------------------------------------------
    virtual double ScatterCrossFtSize( void ) const {  return 0.0;  }


    //========================================================================================================
    //=== PROTECTED INTEGRATION ROUTINES =====================================================================
    //========================================================================================================

    //--------------------------------------------------------------------------------------------------------
    //! \brief  Main mollifier routine. Evaluates the mollification of \pp eval_fcn at the specified point.
    //--------------------------------------------------------------------------------------------------------
    double EvalMoll(
        const ProblemEvalFcn eval_fcn,
        const double x_i,
    # if SPACE_DIMS >= 2
        const double y_j,
    # endif
    # if SPACE_DIMS == 3
        const double z_k,
    # endif
        const double ft_size
    ) const;

    //--------------------------------------------------------------------------------------------------------
    //! \brief  Used in the computation of numerical approximations of convolution integrals used to define
    //!         mollified problem parameters.
    //--------------------------------------------------------------------------------------------------------
    double IntegrateMollDiscontinuous(
        const ProblemEvalFcn eval_fcn,
        const double x_i,
    # if SPACE_DIMS >= 2
        const double y_j,
    # endif
    # if SPACE_DIMS == 3
        const double z_k,
    # endif
        const double ax,
        const double bx,
    # if SPACE_DIMS >= 2
        const double ay,
        const double by,
    # endif
    # if SPACE_DIMS == 3
        const double az,
        const double bz,
    # endif
        const double tol
    ) const;

    //--------------------------------------------------------------------------------------------------------
    //! \brief  Used in the computation of numerical approximations of convolution integrals used to define
    //!         mollified problem parameters.
    //--------------------------------------------------------------------------------------------------------
    double IntegrateMollSmooth(
        const ProblemEvalFcn eval_fcn,
        const double x_i,
    # if SPACE_DIMS >= 2
        const double y_j,
    # endif
    # if SPACE_DIMS == 3
        const double z_k,
    # endif
        const double ax,
        const double bx,
    # if SPACE_DIMS >= 2
        const double ay,
        const double by,
    # endif
    # if SPACE_DIMS == 3
        const double az,
        const double bz,
    # endif
        const double tol
    ) const;

};


# endif // ifndef __PROBLEM_HPP__
