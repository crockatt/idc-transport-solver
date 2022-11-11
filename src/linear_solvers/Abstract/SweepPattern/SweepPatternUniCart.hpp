//------------------------------------------------------------------------------------------------------------
//! \file   linear_solvers/Abstract/SweepPattern/SweepPatternUniCart.hpp
//! \brief  Header file containing declaration of SweepPatternUniCart abstract class.
//!
//! \author Michael M. Crockatt
//! \date   January 2018
//------------------------------------------------------------------------------------------------------------

# ifndef __ABSTRACT__SWEEP_PATTERN_UNI_CART_HPP__
# define __ABSTRACT__SWEEP_PATTERN_UNI_CART_HPP__


# include "linear_solvers/Abstract/SweepPattern/SweepPattern.hpp"


//------------------------------------------------------------------------------------------------------------
//! \brief  Abstract class from which sweep algorithms for uniform Cartesian meshes are derived.
//!
//! \see    SweepPatternUniCartSequential
//! \see    SweepPatternUniCartWavefront
//! \see    SweepPatternUniCartKBA
//! \see    SweepPatternUniCartKBA2
//------------------------------------------------------------------------------------------------------------
template< class OrdinateFlux >
class SweepOperator<OrdinateFlux>::SweepPatternUniCart :
    public SweepOperator<OrdinateFlux>::SweepPattern
{

    //========================================================================================================
    //=== CONSTRUCTORS AND DESTRUCTOR ========================================================================
    //========================================================================================================

protected:

    //--------------------------------------------------------------------------------------------------------
    //! \brief  Delete default constructor.
    //--------------------------------------------------------------------------------------------------------
    SweepPatternUniCart( void ) = delete;

    //--------------------------------------------------------------------------------------------------------
    //! \brief  Delete copy constructor.
    //--------------------------------------------------------------------------------------------------------
    SweepPatternUniCart( const SweepPatternUniCart & ) = delete;

    //--------------------------------------------------------------------------------------------------------
    //! \brief  Initializes a SweepPatternUniCart object.
    //--------------------------------------------------------------------------------------------------------
    SweepPatternUniCart (
        const SweepOperator<OrdinateFlux> & enclosing,
        const ParameterList & input_list
    );

public:

    //--------------------------------------------------------------------------------------------------------
    //! \brief  Virtual destructor.
    //--------------------------------------------------------------------------------------------------------
    virtual ~SweepPatternUniCart( void );


    //========================================================================================================
    //=== OPERATOR OVERLOADS =================================================================================
    //========================================================================================================

    //--------------------------------------------------------------------------------------------------------
    //! \brief  Delete copy assignment operator.
    //--------------------------------------------------------------------------------------------------------
    SweepPatternUniCart & operator=( const SweepPatternUniCart & ) = delete;


    //========================================================================================================
    //=== PUBLIC INTERFACE ROUTINES ==========================================================================
    //========================================================================================================

    //--------------------------------------------------------------------------------------------------------
    //! \brief  Prints a summary of the object's configuration to the logging interface.
    //--------------------------------------------------------------------------------------------------------
    virtual void Print( const std::string & = "  " ) const override;

    //--------------------------------------------------------------------------------------------------------
    //! \brief  Sets parameters given in a ParameterList.
    //--------------------------------------------------------------------------------------------------------
    virtual void SetParameters( const ParameterList & ) override;


protected:

    //========================================================================================================
    //=== INTERNAL HELPER ROUTINES ===========================================================================
    //========================================================================================================

    //--------------------------------------------------------------------------------------------------------
    //! \brief  Computes the number of spatial cells in the specified diagonal of the mesh.
    //--------------------------------------------------------------------------------------------------------
    virtual int64_t ComputeDiagLength( const int64_t diag ) const;


    //========================================================================================================
    //=== INTERNAL VARIABLES =================================================================================
    //========================================================================================================

# if SPACE_DIMS >= 2 || defined (DOXYCOMPILE)

    //!
    //! \brief  Specifies whether quadrants/octants are swept simultaneously or sequentially.
    //!
    //! If true, quadrants/octants are swept sequentially in the order determined by
    //! SweepPatternUniCart::sweep_order.
    //! If false, all quadrants/octants are swept simultaneously.
    //!
    bool sweep_octants_in_sequence;

    //! Order in which to sweep the angular quadrants/octants.
    int64_t sweep_order
        # if SPACE_DIMS == 2 || defined (DOXYCOMPILE)
            [4];
        # elif SPACE_DIMS == 3
            [8];
        # endif

# endif // if SPACE_DIMS >= 2

};


# endif // ifndef __ABSTRACT__SWEEP_PATTERN_UNI_CART_HPP__
