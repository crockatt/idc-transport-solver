//------------------------------------------------------------------------------------------------------------
//! \file   linear_solvers/Abstract/SweepOperator.hpp
//! \brief  Header file containing declaration of SweepOperator class template.
//!
//! \author Michael M. Crockatt
//! \date   January 2018
//------------------------------------------------------------------------------------------------------------

# ifndef __ABSTRACT__SWEEP_OPERATOR_HPP__
# define __ABSTRACT__SWEEP_OPERATOR_HPP__


# include <limits>
# include <memory>

# include "objects/DomainDecomposition.hpp"
# include "objects/RKDG/CrossSection.hpp"
# include "utils/ParameterList.hpp"
# include "utils/Quadrule/OrdinateSet.hpp"
# include "utils/SIMD.hpp"


//------------------------------------------------------------------------------------------------------------
//! \brief  Class template that wraps and manages components used to perform discrete ordinates transport
//!         sweeps.
//!
//------------------------------------------------------------------------------------------------------------
template< class OrdinateFlux >
class SweepOperator :
    public DomainDecomposition,
    public Quadrule::OrdinateSet
{

//
// Make sub-classes public only for deckmaker executable.
//
# if defined (DECKMAKER)
public:
# endif

    //
    // Declarations of nested SweepPattern classes.
    //
    class SweepPattern;
    class SweepPatternUniCart;
    class SweepPatternUniCartSequential;

# if SPACE_DIMS == 2

    class SweepPatternUniCartWavefront;

# if defined (ENABLE_MPI)

    class SweepPatternUniCartWavefrontMPI;
    class SweepPatternUniCartKBA;
    class SweepPatternUniCartKBA2;

# endif // if defined (ENABLE_MPI)
# endif // if SPACE_DIMS == 2

    //
    // Declarations of nested SweepCellSolve classes.
    //
    class SweepCellSolveManager;


    //========================================================================================================
    //=== CONSTRUCTORS AND DESTRUCTOR ========================================================================
    //========================================================================================================

public:

    //--------------------------------------------------------------------------------------------------------
    //! \brief  Delete default constructor.
    //--------------------------------------------------------------------------------------------------------
    SweepOperator( void ) = delete;

    //--------------------------------------------------------------------------------------------------------
    //! \brief  Delete copy constructor.
    //--------------------------------------------------------------------------------------------------------
    SweepOperator( const SweepOperator<OrdinateFlux> & ) = delete;

    //--------------------------------------------------------------------------------------------------------
    //! \brief  Initializes a SweepOperator object.
    //--------------------------------------------------------------------------------------------------------
    SweepOperator (
        const DomainDecomposition & domain_decomposition,
        const Quadrule::OrdinateSet & ordinate_set,
        const ParameterList & input_list,
        const double dt = std::numeric_limits<double>::infinity()
    );

    //--------------------------------------------------------------------------------------------------------
    //! \brief  Virtual destructor.
    //--------------------------------------------------------------------------------------------------------
    virtual ~SweepOperator( void );


    //========================================================================================================
    //=== OPERATOR OVERLOADS =================================================================================
    //========================================================================================================

    //--------------------------------------------------------------------------------------------------------
    //! \brief  Delete copy assignment operator.
    //--------------------------------------------------------------------------------------------------------
    SweepOperator<OrdinateFlux> & operator=( const SweepOperator<OrdinateFlux> & ) = delete;


    //========================================================================================================
    //=== PUBLIC INTERFACE ROUTINES ==========================================================================
    //========================================================================================================

    //--------------------------------------------------------------------------------------------------------
    //! \brief  Prints a summary of the object's configuration to the logging interface.
    //--------------------------------------------------------------------------------------------------------
    void Print( const std::string & = "  " ) const override;

    //--------------------------------------------------------------------------------------------------------
    //! \brief  Sets parameters given in a ParameterList.
    //--------------------------------------------------------------------------------------------------------
    virtual void SetParameters( const ParameterList & );

    //--------------------------------------------------------------------------------------------------------
    //! \brief  Applies the sweep operator.
    //--------------------------------------------------------------------------------------------------------
    virtual void Apply (
        OrdinateFlux & source,
        const RKDG::CrossSection & sigma,
        const double dt = std::numeric_limits<double>::infinity(),
        const RKDG::OrdinateFlux * const initial_condition = nullptr
    ) const;


protected:

    //========================================================================================================
    //=== INTERNAL VARIABLES =================================================================================
    //========================================================================================================

    //!
    //! \brief  Pointer to pattern object for sweeping mesh.
    //!
    std::shared_ptr<SweepPattern> sweep_pattern;

    //!
    //! \brief  Pointer to sweep solve object for local solves.
    //!
    std::shared_ptr<SweepCellSolveManager> sweep_cell_solve;

};


# endif // ifndef __ABSTRACT__SWEEP_OPERATOR_HPP__
