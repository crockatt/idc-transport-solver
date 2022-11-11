//------------------------------------------------------------------------------------------------------------
//! \file   linear_solvers/Abstract/SweepPattern/SweepPatternUniCart.cpp
//! \brief  Contains implementations and instantiations of methods from SweepPatternUniCart class template.
//!
//! \author Michael M. Crockatt
//! \date   January 2018
//------------------------------------------------------------------------------------------------------------

# include <typeinfo>

# include "linear_solvers/Abstract/SweepPattern/SweepPatternUniCart.hpp"
# include "utils/DirectedGraph.hpp"


//------------------------------------------------------------------------------------------------------------
//! \brief  Initializes a SweepPatternUniCart object.
//!
//! \param[in]  enclosing   Reference to instance of enclosing SweepOperator class.
//! \param[in]  input_list  Contains parameters for initialization of object.
//------------------------------------------------------------------------------------------------------------
template< class OrdinateFlux >
SweepOperator<OrdinateFlux>::SweepPatternUniCart::SweepPatternUniCart (

    const SweepOperator<OrdinateFlux> & enclosing,
    const ParameterList & input_list
) :
    SweepPattern( enclosing, input_list )

# if SPACE_DIMS >= 2

    //
    // Default sweep order. Alternate order is computed only when one or more reflecting boundaries are used.
    //
,   sweep_order
        # if SPACE_DIMS == 2
            { 3, 2, 1, 0 }
        # elif SPACE_DIMS == 3
            { 0, 1, 6, 7, 4, 5, 2, 3 }
        # endif

# endif // if SPACE_DIMS >= 2
{
    PRINT_STATUS( "Executing SweepOperator<%s>::SweepPatternUniCart::%s.\n",
                  typeid(OrdinateFlux).name(), __func__ )

# if SPACE_DIMS >= 2

    // Determine whether quadrants/octants should be swept simultaneously or sequentially.
    if (
        // X dimension.
            (    BitmaskHasAll( Global::reflecting_boundaries, BoundaryEdge::X_Min )
            # if defined (ENABLE_MPI)
              && this->sw_op.MPI_block_coords(0) == 0
            # endif
            )
         || (    BitmaskHasAll( Global::reflecting_boundaries, BoundaryEdge::X_Max )
            # if defined (ENABLE_MPI)
              && this->sw_op.MPI_block_coords(0) == this->sw_op.MPI_num_blocks(0) - 1
            # endif
            )
    # if SPACE_DIMS >= 2
        // Y dimension.
         || (    BitmaskHasAll( Global::reflecting_boundaries, BoundaryEdge::Y_Min )
            # if defined (ENABLE_MPI)
              && this->sw_op.MPI_block_coords(1) == 0
            # endif
            )
         || (    BitmaskHasAll( Global::reflecting_boundaries, BoundaryEdge::Y_Max )
            # if defined (ENABLE_MPI)
              && this->sw_op.MPI_block_coords(1) == this->sw_op.MPI_num_blocks(1) - 1
            # endif
            )
    # endif // if SPACE_DIMS >= 2
    # if SPACE_DIMS == 3
        // Z dimension
         || (    BitmaskHasAll( Global::reflecting_boundaries, BoundaryEdge::Z_Min )
            # if defined (ENABLE_MPI)
              && this->sw_op.MPI_block_coords(2) == 0
            # endif
            )
         || (    BitmaskHasAll( Global::reflecting_boundaries, BoundaryEdge::Z_Max )
            # if defined (ENABLE_MPI)
              && this->sw_op.MPI_block_coords(2) == this->sw_op.MPI_block_coords(2) - 1
            # endif
            )
    # endif // if SPACE_DIMS == 3
    )    {  this->sweep_octants_in_sequence = true;   }
    else {  this->sweep_octants_in_sequence = false;  }

# if 0
    // Compute alternate sweep ordering if reflecting boundary conditions are used.
    if ( BitmaskHasAny( Global::reflecting_boundaries, BoundaryEdge::All ) ) {

        PRINT_STATUS( "Determining sweep ordering for reflecting boundary conditions.\n" )

        DirectedGraph<int64_t> quadrant_graph;

    # if SPACE_DIMS == 2

        if ( BitmaskHasAll( Global::reflecting_boundaries, BoundaryEdge::X_Min ) ) {

            quadrant_graph.AddEdge(2,0);    // II  ↦ I
            quadrant_graph.AddEdge(3,1);    // III ↦ IV
        }

        if ( BitmaskHasAll( Global::reflecting_boundaries, BoundaryEdge::X_Max ) ) {

            quadrant_graph.AddEdge(0,2);    // I  ↦ II
            quadrant_graph.AddEdge(1,3);    // IV ↦ III
        }

        if ( BitmaskHasAll( Global::reflecting_boundaries, BoundaryEdge::Y_Min ) ) {

            quadrant_graph.AddEdge(3,2);    // III ↦ II
            quadrant_graph.AddEdge(1,0);    // IV  ↦ I
        }

        if ( BitmaskHasAll( Global::reflecting_boundaries, BoundaryEdge::Y_Max ) ) {

            quadrant_graph.AddEdge(2,3);    // II ↦ III
            quadrant_graph.AddEdge(0,1);    // I  ↦ IV
        }

    # elif SPACE_DIMS == 3
        # warning "Implementation of SweepPatternUniCart::SweepPatternUniCart incomplete."
    # endif

        const auto sorted_list = quadrant_graph.TopologicalSort();

        int64_t quad_index = 0;

        for ( auto & quad : sorted_list ) {

            sweep_order[quad_index] = quad;
            ++quad_index;
        }
    }
# endif // if 0
# endif // if SPACE_DIMS >= 2
}


//------------------------------------------------------------------------------------------------------------
//! \brief  Virtual destructor.
//------------------------------------------------------------------------------------------------------------
template< class OrdinateFlux >
SweepOperator<OrdinateFlux>::SweepPatternUniCart::~SweepPatternUniCart ( void ) {

    PRINT_STATUS( "Executing SweepOperator<%s>::SweepPatternUniCart::%s.\n",
                  typeid(OrdinateFlux).name(), __func__ )

    /* empty */
}


//------------------------------------------------------------------------------------------------------------
//! \brief  Prints a summary of the object's configuration to the logging interface.
//!
//! \param[in]  prefix  String prepended to each line of output, intended to be used for controlling
//!                     indentation of output.
//------------------------------------------------------------------------------------------------------------
template< class OrdinateFlux >
void SweepOperator<OrdinateFlux>::SweepPatternUniCart::Print (

    const std::string & prefix // = "  "

) const {

    PRINT_STATUS( "Executing SweepOperator<%s>::SweepPatternUniCart::%s.\n",
                  typeid(OrdinateFlux).name(), __func__ )

# if SPACE_DIMS >= 2

    // Output sweep method for octants/quadrants.
    const std::string label = std::string("Sweep ") +
        # if SPACE_DIMS == 2
            "Quadrants"
        # elif SPACE_DIMS == 3
            "Octants"
        # endif
            + std::string(":");

    PRINT_LOG( "%s%-*s  %s\n", prefix.c_str(), Global::col_width, label.c_str(),
               ( this->sweep_octants_in_sequence ? "sequentially" : "simultaneously" ) )

    // Output sequence if sweeping octants/quadrants sequentially.
    if ( this->sweep_octants_in_sequence ) {

    # if SPACE_DIMS == 2

        PRINT_LOG( "%s%-*s  %" PRId64 " %" PRId64 " %" PRId64 " %" PRId64 "\n", prefix.c_str(),
                   Global::col_width, "Sweep Order:",
                   sweep_order[0], sweep_order[1], sweep_order[2], sweep_order[3] )

    # elif SPACE_DIMS == 3

        PRINT_LOG( "%s%-20s  %" PRId64 " %" PRId64 " %" PRId64 " %" PRId64 " %" PRId64 " %" PRId64 " %" PRId64
                   " %" PRId64 "\n", prefix.c_str(),
                   Global::col_width, "Sweep Order:",
                   sweep_order[0], sweep_order[1], sweep_order[2], sweep_order[3],
                   sweep_order[4], sweep_order[5], sweep_order[6], sweep_order[7] )

    # endif // if SPACE_DIMS == ?
    }

# endif // if SPACE_DIMS >= 2

    this->SweepPattern::Print( prefix );
}


//------------------------------------------------------------------------------------------------------------
//! \brief  Sets parameters given in a ParameterList.
//!
//! \param[in]  input_list  Contains list of parameters to set.
//------------------------------------------------------------------------------------------------------------
template< class OrdinateFlux >
void SweepOperator<OrdinateFlux>::SweepPatternUniCart::SetParameters (

    const ParameterList & input_list
) {

    PRINT_STATUS( "Executing SweepOperator<%s>::SweepPatternUniCart::%s.\n",
                  typeid(OrdinateFlux).name(), __func__ )

    this->SweepPattern::SetParameters( input_list );
}


//------------------------------------------------------------------------------------------------------------
//! \brief  Computes the number of spatial cells in the specified diagonal of the mesh.
//!
//! \param[in]  diag    Diagonal of the spatial mesh to compute length of.
//!
//! \return     Returns the number of spatial cells in the specified diagonal.
//------------------------------------------------------------------------------------------------------------
template< class OrdinateFlux >
int64_t SweepOperator<OrdinateFlux>::SweepPatternUniCart::ComputeDiagLength (

    const int64_t diag

) const {

    PRINT_STATUS( "Executing SweepOperator<%s>::SweepPatternUniCart::%s.\n",
                  typeid(OrdinateFlux).name(), __func__ )

    const int64_t num_diags = this->sw_op.nx(0) + this->sw_op.nx(1) - 1;

    return std::min( diag,
                     std::min( num_diags - diag + 1,
                               std::min( this->sw_op.nx(0), this->sw_op.nx(1) )
                             )
                   );
}


//============================================================================================================
//=== EXPLICIT INSTANTIATIONS ================================================================================
//============================================================================================================

# include "objects/RKDG/OrdinateFlux.hpp"
# include "objects/STDG/OrdinateFlux.hpp"


//!
//! \brief  Instantiation of SweepOperator<RKDG::OrdinateFlux>::SweepPatternUniCart class.
//!
template class SweepOperator<RKDG::OrdinateFlux>::SweepPatternUniCart;


//!
//! \brief  Instantiation of SweepOperator<STDG::OrdinateFlux>::SweepPatternUniCart class.
//!
template class SweepOperator<STDG::OrdinateFlux>::SweepPatternUniCart;
