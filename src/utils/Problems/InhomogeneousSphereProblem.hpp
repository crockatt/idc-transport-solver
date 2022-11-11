//------------------------------------------------------------------------------------------------------------
//! \file   utils/Problems/InhomogeneousSphereProblem.hpp
//! \brief  Header for InhomogeneousSphereProblem class.
//!
//! \author Michael M. Crockatt
//! \date   June 2018
//------------------------------------------------------------------------------------------------------------

# ifndef __INHOMOGENEOUS_SPHERE_PROBLEM_HPP__
# define __INHOMOGENEOUS_SPHERE_PROBLEM_HPP__


# include <memory>

# include "utils/Mollifier.hpp"
# include "utils/Problems/Problem.hpp"


//------------------------------------------------------------------------------------------------------------
//! \brief  Radiating sphere with non-uniform density profile.
//!
//! \todo   Description of inhomogeneous sphere problem.
//------------------------------------------------------------------------------------------------------------
class InhomogeneousSphereProblem :
    public Problem
{

public:

    //========================================================================================================
    //=== CONSTRUCTORS AND DESTRUCTOR ========================================================================
    //========================================================================================================

    //--------------------------------------------------------------------------------------------------------
    //! \brief  Delete default constructor.
    //--------------------------------------------------------------------------------------------------------
    InhomogeneousSphereProblem( void ) = delete;

    //--------------------------------------------------------------------------------------------------------
    //! \brief  Delete copy constructor.
    //--------------------------------------------------------------------------------------------------------
    InhomogeneousSphereProblem( const InhomogeneousSphereProblem & ) = delete;

    //--------------------------------------------------------------------------------------------------------
    //! \brief  Construct object from given parameters.
    //--------------------------------------------------------------------------------------------------------
    InhomogeneousSphereProblem( const ParameterList & );

    //--------------------------------------------------------------------------------------------------------
    //! \brief  Destructor for InhomogeneousSphereProblem class.
    //--------------------------------------------------------------------------------------------------------
    ~InhomogeneousSphereProblem( void ) override;


    //========================================================================================================
    //=== OPERATOR OVERLOADS =================================================================================
    //========================================================================================================

    //--------------------------------------------------------------------------------------------------------
    //! \brief  Delete copy assignment operator.
    //--------------------------------------------------------------------------------------------------------
    InhomogeneousSphereProblem & operator=( const InhomogeneousSphereProblem & ) = delete;


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

    //!
    //! \brief  Specifies the radius at which the first transition in the sphere occurs.
    //!
    double inner_radius;

    //!
    //! \brief  Specifies the radius at which the second transition in the sphere occurs.
    //!
    double outer_radius;

    //!
    //! \brief  Specifies the length of the first transition in the sphere.
    //!
    double inner_trans_radius;

    //!
    //! \brief  Specifies the length of the second transition in the sphere.
    //!
    double outer_trans_radius;

    //!
    //! \brief  Specifies the cross section of the sphere in the innermost region.
    //!
    double inner_sigma_val;

    //!
    //! \brief  Specifies the cross section of the sphere in the middle region.
    //!
    double interm_sigma_val;

    //!
    //! \brief  Specifies the cross section of the sphere in the outermost region.
    //!
    double outer_sigma_val;


    //========================================================================================================
    //=== PROTECTED INTERFACE ROUTINES =======================================================================
    //========================================================================================================

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
    ) const override;

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
    ) const override;

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
    ) const override;

};


# endif // ifndef __INHOMOGENEOUS_SPHERE_PROBLEM_HPP__
