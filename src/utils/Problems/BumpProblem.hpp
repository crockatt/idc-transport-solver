//------------------------------------------------------------------------------------------------------------
//! \file   utils/Problems/BumpProblem.hpp
//! \brief  Header for BumpProblem class.
//!
//! \author Michael M. Crockatt
//! \date   June 2018
//------------------------------------------------------------------------------------------------------------

# ifndef __BUMP_PROBLEM_HPP__
# define __BUMP_PROBLEM_HPP__


# include <memory>

# include "utils/Mollifier.hpp"
# include "utils/Problems/Problem.hpp"


//------------------------------------------------------------------------------------------------------------
//! \brief  Simple \f$ C^{\infty} \f$ initial condition relaxing in a uniform medium.
//!
//!
//------------------------------------------------------------------------------------------------------------
class BumpProblem :
    public Problem
{

public:

    //========================================================================================================
    //=== CONSTRUCTORS AND DESTRUCTOR ========================================================================
    //========================================================================================================

    //--------------------------------------------------------------------------------------------------------
    //! \brief  Delete default constructor.
    //--------------------------------------------------------------------------------------------------------
    BumpProblem( void ) = delete;

    //--------------------------------------------------------------------------------------------------------
    //! \brief  Delete copy constructor.
    //--------------------------------------------------------------------------------------------------------
    BumpProblem( const BumpProblem & ) = delete;

    //--------------------------------------------------------------------------------------------------------
    //! \brief  Construct object from given parameters.
    //--------------------------------------------------------------------------------------------------------
    BumpProblem( const ParameterList & );

    //--------------------------------------------------------------------------------------------------------
    //! \brief  Destructor for BumpProblem class.
    //--------------------------------------------------------------------------------------------------------
    ~BumpProblem( void ) override;


    //========================================================================================================
    //=== OPERATOR OVERLOADS =================================================================================
    //========================================================================================================

    //--------------------------------------------------------------------------------------------------------
    //! \brief  Delete copy assignment operator.
    //--------------------------------------------------------------------------------------------------------
    BumpProblem & operator=( const BumpProblem & ) = delete;


    //========================================================================================================
    //=== PUBLIC INTERFACE ROUTINES ==========================================================================
    //========================================================================================================

    //--------------------------------------------------------------------------------------------------------
    //! \brief  Returns the string descriptor of the class type.
    //--------------------------------------------------------------------------------------------------------
    static const std::string Descriptor( void );

    //--------------------------------------------------------------------------------------------------------
    //! \brief  Returns the string descriptor for an object's type.
    //--------------------------------------------------------------------------------------------------------
    const std::string GetDescriptor( void ) const override;

    //--------------------------------------------------------------------------------------------------------
    //! \brief  Prints a summary of the problem configuration to the logging interface.
    //--------------------------------------------------------------------------------------------------------
    void Print( const std::string & = "  " ) const override;


protected:

    //========================================================================================================
    //=== PROTECTED MEMBER VARIABLES =========================================================================
    //========================================================================================================

    double sigma_val;       //!< Cross section value.
    double bump_radius;     //!< Radius of initial condition.

    std::shared_ptr<Mollifier> bump_moll;   //!< Object used to evaluate initial condition.


    //========================================================================================================
    //=== PROTECTED INTERFACE ROUTINES =======================================================================
    //========================================================================================================

    //--------------------------------------------------------------------------------------------------------
    //! \brief  Returns the value of the initial condition without mollification at the specified position.
    //--------------------------------------------------------------------------------------------------------
    virtual double InitialCondition(
            const double x
    # if SPACE_DIMS >= 2 || defined (DOXYCOMPILE)
        ,   const double y
    # endif
    # if SPACE_DIMS == 3
        ,   const double z
    # endif
    ) const override;

    //--------------------------------------------------------------------------------------------------------
    //! \brief  Returns the value of the total cross section without mollification at the specified position.
    //--------------------------------------------------------------------------------------------------------
    virtual double TotalCross(
            const double x
    # if SPACE_DIMS >= 2 || defined (DOXYCOMPILE)
        ,   const double y
    # endif
    # if SPACE_DIMS == 3
        ,   const double z
    # endif
    ) const override;

    //--------------------------------------------------------------------------------------------------------
    //! \brief  Returns the value of the scattering cross section without mollification at the specified
    //!         position.
    //--------------------------------------------------------------------------------------------------------
    virtual double ScatterCross(
            const double x
    # if SPACE_DIMS >= 2 || defined (DOXYCOMPILE)
        ,   const double y
    # endif
    # if SPACE_DIMS == 3
        ,   const double z
    # endif
    ) const override;

};


# endif // ifndef __BUMP_PROBLEM_HPP__
