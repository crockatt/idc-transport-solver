//------------------------------------------------------------------------------------------------------------
//! \file   utils/Problems/Hohlraum2Problem.hpp
//! \brief  Header for Hohlraum2Problem class.
//!
//! \author Michael M. Crockatt
//! \date   June 2018
//------------------------------------------------------------------------------------------------------------

# ifndef __HOHLRAUM2_PROBLEM_HPP__
# define __HOHLRAUM2_PROBLEM_HPP__

# if SPACE_DIMS == 2 || defined (DOXYCOMPILE)


# include "utils/Problems/HohlraumProblem.hpp"


//------------------------------------------------------------------------------------------------------------
//! \brief  Hohlraum problem.
//!
//! \todo   Description of hohlraum test problem.
//------------------------------------------------------------------------------------------------------------
class Hohlraum2Problem :
    public HohlraumProblem
{

public:

    //========================================================================================================
    //=== CONSTRUCTORS AND DESTRUCTOR ========================================================================
    //========================================================================================================

    //--------------------------------------------------------------------------------------------------------
    //! \brief  Delete default constructor.
    //--------------------------------------------------------------------------------------------------------
    Hohlraum2Problem( void ) = delete;

    //--------------------------------------------------------------------------------------------------------
    //! \brief  Delete copy constructor.
    //--------------------------------------------------------------------------------------------------------
    Hohlraum2Problem( const Hohlraum2Problem & ) = delete;

    //--------------------------------------------------------------------------------------------------------
    //! \brief  Construct object from given parameters.
    //--------------------------------------------------------------------------------------------------------
    Hohlraum2Problem( const ParameterList & );

    //--------------------------------------------------------------------------------------------------------
    //! \brief  Destructor for Hohlraum2Problem class.
    //--------------------------------------------------------------------------------------------------------
    ~Hohlraum2Problem( void ) override;


    //========================================================================================================
    //=== OPERATOR OVERLOADS =================================================================================
    //========================================================================================================

    //--------------------------------------------------------------------------------------------------------
    //! \brief  Delete copy assignment operator.
    //--------------------------------------------------------------------------------------------------------
    Hohlraum2Problem & operator=( const Hohlraum2Problem & ) = delete;


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
    //=== PROTECTED INTERFACE ROUTINES =======================================================================
    //========================================================================================================

    //--------------------------------------------------------------------------------------------------------
    //! \brief  Returns the value of the scattering cross section without mollification at the specified
    //!         position.
    //--------------------------------------------------------------------------------------------------------
    virtual double ScatterCross(
        const double x,
        const double y
    ) const override;

};


# endif // if SPACE_DIMS == 2

# endif // ifndef __HOHLRAUM2_PROBLEM_HPP__