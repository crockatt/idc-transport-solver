//------------------------------------------------------------------------------------------------------------
//! \file   objects/RKDG/OrdinateFlux.cpp
//! \brief  Implementation of RKDG::OrdinateFlux class.
//!
//! \author Michael M. Crockatt
//! \date   January 2018
//------------------------------------------------------------------------------------------------------------


# include <algorithm>
# include <cinttypes>
# include <cmath>
# include <cstring>
# include <limits>
# include <stdexcept>
# include <string>

# if defined (ENABLE_SIMD_BLOCKING)

    # ifndef LINF_SIMD_LEN
        # define LINF_SIMD_LEN 16
    # endif

    # include <vector>

# else // if defined (ENABLE_SIMD_BLOCKING)

    # define LINF_SIMD_LEN 1

# endif // if defined (ENABLE_SIMD_BLOCKING)

# if defined (ENABLE_MPI)
    # include <mpi.h>
# endif

# include "objects/RKDG/OrdinateFlux.hpp"
# include "objects/RKDG/ScalarFlux.hpp"
# include "operators/RelabelOperator.hpp"
# include "operators/TransportOperator.hpp"
# include "utils/CLog.hpp"
# include "utils/Quadrule/Quadrule.hpp"


using namespace RKDG;
using namespace Quadrule;


//============================================================================================================
//=== CONSTRUCTORS, DESTRUCTORS, AND ASSOCIATED HELPER ROUTINES ==============================================
//============================================================================================================


//------------------------------------------------------------------------------------------------------------
//! \brief  Constructs an empty RKDG::OrdinateFlux object with zero-initialized parameters.
//------------------------------------------------------------------------------------------------------------
OrdinateFlux::OrdinateFlux ( void ) :

    DomainDecomposition(),
    Abstract::DensityFunction(),
    Abstract::OrdinateFlux()
{}


//------------------------------------------------------------------------------------------------------------
//! \brief  Constructs an RKDG::OrdinateFlux object with the given set of parameters.
//!
//! \param[in]  spatial_params      Contains the parameters of the spatial discretization.
//! \param[in]  DG_degree_in        Maximum degree of DG basis polynomials.
//! \param[in]  ang_order           Order of the discrete ordinate set of type \pp{ordinate_type} to use.
//! \param[in]  symmetric_reduce    (optional) <br/>
//!                                 Specifies whether or not the number of ordinates should be reduced by
//!                                 symmetry. Only used for 2D problems. Defaults to false.
//! \param[in]  ordinate_type       (optional) <br/>
//!                                 OrdinateType specifying the type of discrete ordinate quadrature to use.
//!                                 Defaults to Gauss-Legendre quadrature in 1D and Chebyshev-Legendre
//!                                 quadrature in multi-D.
//!
//! \see    RKDG::OrdinateFlux::Allocate()
//------------------------------------------------------------------------------------------------------------
OrdinateFlux::OrdinateFlux (

    const DomainDecomposition & spatial_params,
    const int64_t DG_degree_in,
    const int64_t ang_order,
    const bool symmetric_reduce,        // = false
    const OrdinateType ordinate_type    // Defaulted.
) :
    DomainDecomposition( spatial_params ),
    Abstract::DensityFunction( spatial_params ),
    Abstract::OrdinateFlux( ang_order, symmetric_reduce, ordinate_type ),
    DG_degree{ DG_degree_in }
{
    Allocate();
}


//------------------------------------------------------------------------------------------------------------
//! \brief  Constructs an RKDG::OrdinateFlux object with the given set of parameters.
//!
//! Delegates to RKDG::OrdinateFlux::OrdinateFlux( const DomainDecomposition &, const int64_t, const int64_t, const OrdinateType ).
//!
//! \param[in]      spatial_params  Contains the parameters of the spatial discretization.
//! \param[in]      DG_degree_in    Maximum degree of DG basis polynomials.
//! \param[in]      ordinate_set    OrdinateSet containing the discrete ordinate quadrature to use.
//!
//! \see    RKDG::OrdinateFlux::Allocate()
//------------------------------------------------------------------------------------------------------------
OrdinateFlux::OrdinateFlux (

    const DomainDecomposition & spatial_params,
    const OrdinateSet & ordinate_set,
    const int64_t DG_degree_in

) : /* OrdinateFlux( spatial_params, DG_degree_in, ordinate_set.GetAngOrder(), ordinate_set.GetOrdinateSymmetry(),
                  ordinate_set.GetOrdinateType() )
{} */
    DomainDecomposition( spatial_params ),
    Abstract::DensityFunction( spatial_params ),
    Abstract::OrdinateFlux( ordinate_set.GetAngOrder(), ordinate_set.GetOrdinateSymmetry(),
                            ordinate_set.GetOrdinateType() ),
    DG_degree{ DG_degree_in }
{
    Allocate();
}


//------------------------------------------------------------------------------------------------------------
//! \brief  Sets appropriate values for object parameters
//------------------------------------------------------------------------------------------------------------
void OrdinateFlux::SetDensityDimensions ( void ) {

    // Set value for DOF per cell.
    this->DOF_per_cell = this->nq() * IntPow( this->DG_degree + 1, SPACE_DIMS );

    // Call up class hierarchy.
    this->Abstract::DensityFunction::SetDensityDimensions();
}


//------------------------------------------------------------------------------------------------------------
//! \brief  Reconfigures an RKDG::OrdinateFlux object with the given set of parameters.
//!
//! \param[in]  ang_order           Order of the discrete ordinate set of type \pp{ordinate_type} to use.
//! \param[in]  symmetric_reduce    (optional) <br/>
//!                                 Specifies whether or not the number of ordinates should be reduced by
//!                                 symmetry. Only used for 2D problems. Defaults to false.
//! \param[in]  ordinate_type       (optional) <br/>
//!                                 OrdinateType specifying the type of discrete ordinate quadrature to use.
//!                                 Defaults to Gauss-Legendre quadrature in 1D and Chebyshev-Legendre
//!                                 quadrature in multi-D.
//!
//! \see    RKDG::OrdinateFlux::Allocate()
//! \see    RKDG::OrdinateFlux::Deallocate()
//------------------------------------------------------------------------------------------------------------
OrdinateFlux & OrdinateFlux::Reconfigure (

    const int64_t ang_order,
    const bool symmetric_reduce,        // = false
    const OrdinateType ordinate_type    // Defaulted.
) {

    Deallocate();

    this->OrdinateSet::Reconfigure( ang_order, symmetric_reduce, ordinate_type );

    Allocate();

    return *this;
}


//------------------------------------------------------------------------------------------------------------
//! \brief  Reconfigures an RKDG::OrdinateFlux object with the parameters of a given OrdinateSet object.
//!
//! \param[in]      ordinate_set    OrdinateSet object containing parameters to use.
//!
//! \see    RKDG::OrdinateFlux::Allocate()
//! \see    RKDG::OrdinateFlux::Deallocate()
//------------------------------------------------------------------------------------------------------------
OrdinateFlux & OrdinateFlux::Reconfigure (

    const OrdinateSet & ordinate_set
) {
    return Reconfigure( ordinate_set.GetAngOrder(), ordinate_set.GetOrdinateSymmetry(),
                        ordinate_set.GetOrdinateType() );
}


//------------------------------------------------------------------------------------------------------------
//! \brief  Reconfigures an RKDG::OrdinateFlux object with the given set of parameters.
//!
//! \param[in]  spatial_params      Contains the parameters of the spatial discretization.
//! \param[in]  DG_degree_in        Maximum degree of DG basis polynomials.
//! \param[in]  ang_order           Order of the discrete ordinate set of type \pp{ordinate_type} to use.
//! \param[in]  symmetric_reduce    (optional) <br/>
//!                                 Specifies whether or not the number of ordinates should be reduced by
//!                                 symmetry. Only used for 2D problems. Defaults to false.
//! \param[in]  ordinate_type       (optional) <br/>
//!                                 OrdinateType specifying the type of discrete ordinate quadrature to use.
//!                                 Defaults to Gauss-Legendre quadrature in 1D and Chebyshev-Legendre
//!                                 quadrature in multi-D.
//!
//! \see    RKDG::OrdinateFlux::Allocate()
//! \see    RKDG::OrdinateFlux::Deallocate()
//------------------------------------------------------------------------------------------------------------
OrdinateFlux & OrdinateFlux::Reconfigure (

    const DomainDecomposition & spatial_params,
    const int64_t DG_degree_in,
    const int64_t ang_order,
    const bool symmetric_reduce,        // = false
    const OrdinateType ordinate_type    // Defaulted.
) {

    Deallocate();

    DomainDecomposition::operator=( spatial_params );
    this->DG_degree = DG_degree_in;
    this->OrdinateSet::Reconfigure( ang_order, symmetric_reduce, ordinate_type );

    Allocate();

    return *this;
}


//------------------------------------------------------------------------------------------------------------
//! \brief  Reconfigures an RKDG::OrdinateFlux object with the given set of parameters.
//!
//! \param[in]      spatial_params  Contains the parameters of the spatial discretization.
//! \param[in]      DG_degree_in    Maximum degree of DG basis polynomials.
//! \param[in]      ordinate_set    OrdinateSet object containing parameters to use.
//!
//! \see    RKDG::OrdinateFlux::Allocate()
//! \see    RKDG::OrdinateFlux::Deallocate()
//------------------------------------------------------------------------------------------------------------
OrdinateFlux & OrdinateFlux::Reconfigure (

    const DomainDecomposition & spatial_params,
    const OrdinateSet & ordinate_set,
    const int64_t DG_degree_in
) {
    return Reconfigure( spatial_params, DG_degree_in, ordinate_set.GetAngOrder(),
                        ordinate_set.GetOrdinateSymmetry(), ordinate_set.GetOrdinateType() );
}


//============================================================================================================
//=== ONE SPATIAL DIMENSION ==================================================================================
//============================================================================================================

# if SPACE_DIMS == 1


//------------------------------------------------------------------------------------------------------------
//! \brief  Routine to obtain the value of the density function represented in an #OrdinateFlux object at a
//!         specified point and ordinate.
//!
//! \param[in]      x           Point to obtain value at.
//! \param[in]      q           Index of ordinate to obtain value at.
//!
//! \return     Returns the value of the density at the given point.
//------------------------------------------------------------------------------------------------------------
double OrdinateFlux::EvaluateAtPoint (

    const double x,
    const int64_t q

) const {

# if defined (STRICT_CHECK)

    if (    x < this->global_ax(0)
         || x > this->global_bx(0)
    ) {
        std::ostringstream error_message;
        error_message << std::scientific
                      << "Evaluation point "
                      << x
                      << " is outside of OrdinateFlux's domain of [ "
                      << this->global_ax(0) << ", " << this->global_bx(0)
                      << " ].\n";

        PRINT_ERROR( error_message.str().c_str() )
        throw std::range_error( error_message.str() );
    }

# endif // if defined (STRICT_CHECK)

    double ret_val = 0.0;
    int block_coords [SPACE_DIMS];
    int local_coords [SPACE_DIMS];

    PointToCell( {x}, block_coords, local_coords );

# if defined (ENABLE_MPI)

    // If point is inside MPI block of current process, compute the value at that point.
    if ( block_coords[0] == this->MPI_block_coords(0) )
# endif
    {
        // The point 'x' lies inside the spatial cell with index 'i'.
        const int64_t i = local_coords[0];

        // Compute \bar{x} to evaluate Legendre polynomials with.
        const double x_bar = 2.0 * (x - this->ax(0)) / this->dx(0) - (2*i + 1);

        // Evaluate density function at given point.
        for ( int64_t d = 0; d <= this->DG_degree; ++d )
            ret_val += (*this)(q,i+1,d) * Legendre( d, x_bar );

    # if defined (ENABLE_MPI)

        if ( Global::MPI_rank != 0 ) {

            MPI_Send( &ret_val, 1, MPI_DOUBLE, 0, 0, Global::MPI_cart_comm );
            ret_val = std::nan("0");
        }
    # endif // if defined (ENABLE_MPI)
    }

# if defined (ENABLE_MPI)

    else {

        if ( Global::MPI_rank == 0 ) {

            int recv_from_rank;
            MPI_Cart_rank( Global::MPI_cart_comm, block_coords, &recv_from_rank );

            MPI_Recv( &ret_val, 1, MPI_DOUBLE, recv_from_rank, 0, Global::MPI_cart_comm, MPI_STATUS_IGNORE );

        } else {  ret_val = std::nan("0");  }
    }

    MPI_Barrier( Global::MPI_cart_comm );

# endif // if defined (ENABLE_MPI)

    return ret_val;
}


//------------------------------------------------------------------------------------------------------------
//! \brief  Routine to obtain the value of the density function represented in an #OrdinateFlux object at a
//!         specified point and ordinate. The given point must lie within the local domain decomposition of
//!         the current MPI rank.
//!
//! \param[in]      x           Point to obtain value at.
//! \param[in]      q           Index of ordinate to obtain value at.
//!
//! \return     Returns the value of the density at the given point.
//------------------------------------------------------------------------------------------------------------
double OrdinateFlux::EvaluateAtLocalPoint (

    const double x,
    const int64_t q

) const {

# if defined (STRICT_CHECK)

    if (    x < this->global_ax(0)
         || x > this->global_bx(0)
    ) {
        std::ostringstream error_message;
        error_message << std::scientific
                      << "Evaluation point "
                      << x
                      << " is outside of OrdinateFlux's domain of [ "
                      << this->global_ax(0) << ", " << this->global_bx(0)
                      << " ].\n";

        PRINT_ERROR( error_message.str().c_str() )
        throw std::range_error( error_message.str() );
    }

# endif // if defined (STRICT_CHECK)

    double ret_val = 0.0;
    int block_coords [SPACE_DIMS];
    int local_coords [SPACE_DIMS];

    PointToCell( {x}, block_coords, local_coords );

# if defined (ENABLE_MPI)

    // Point needs to be inside the current MPI rank.
    if ( block_coords[0] != this->MPI_block_coords(0) ) {

        std::ostringstream error_message;
        error_message << std::scientific
                      << "Evaluation point " << x << " is inside MPI block ("
                      << block_coords[0]
                      << ") but the current MPI rank has block coordinates ("
                      << this->MPI_block_coords(0)
                      << ").\n";

        PRINT_ERROR( error_message.str().c_str() )
        throw std::range_error( error_message.str() );
    }
# endif  // if defined (ENABLE_MPI)

    {
        // The point 'x' lies inside the spatial cell with index 'i'.
        const int64_t i = local_coords[0];

        // Compute \bar{x} to evaluate Legendre polynomials with.
        const double x_bar = 2.0 * (x - this->ax(0)) / this->dx(0) - (2*i + 1);

        // Evaluate density function at given point.
        for ( int64_t d = 0; d <= this->DG_degree; ++d )
            ret_val += (*this)(q,i+1,d) * Legendre( d, x_bar );
    }

    return ret_val;
}


//------------------------------------------------------------------------------------------------------------
//! \brief  Computes the difference between two RKDG::OrdinateFlux objects.
//!
//! This operation is lossless in the sense that the result is upscaled to the coarsest finite element space
//! which is guaranteed to represent the difference between the two density functions exactly.
//!
//! \attention  If the angular resolutions of the two fluxes differ, then the solution with the coarser
//!             resolution is interpolated to the higher resolution before computing the difference.
//!
//! \attention  \pp{result} is either <tt>(\pp{data1} - \pp{data2})</tt> or <tt>(\pp{data2} - \pp{data1})</tt>
//!             depending on the structure of each solution; i.e., this routine is only useful when absolute
//!             differences are computed using one of the norm routines.
//!
//! \attention  The output object \pp{result} will be reconfigured as necessary.
//!
//! \param[in]      data1       First RKDG::OrdinateFlux object in difference.
//! \param[in]      data2       Second RKDG::OrdianteFlux object in difference.
//! \param[out]     result      Upon return, contains the difference between \pp{first} and \pp{second}.
//------------------------------------------------------------------------------------------------------------
void OrdinateFlux::ComputeDifference (

    const OrdinateFlux & data1,
    const OrdinateFlux & data2,
    OrdinateFlux & result
) {

    const double TOL = 1e-14;

    if (    std::abs( data1.ax(0) - data2.ax(0) ) > TOL
         || std::abs( data1.bx(0) - data2.bx(0) ) > TOL
    ) {
        std::string error_message
            = "Domains of RKDG::OrdainteFlux objects in '"
              + std::string(__func__)
              + "' are not compatible for computing difference.\n";

        PRINT_ERROR( error_message.c_str() )
        throw std::invalid_argument( error_message );
    }

    const OrdinateFlux * first = nullptr, * second = nullptr;
    const OrdinateFlux * coarse = nullptr, * fine = nullptr;

    double * GL_nodes = nullptr, * GL_weights = nullptr;
    double * fine_project = nullptr, * coarse_project = nullptr;

    int64_t GL_order = 0, fine_num_new_cells = 0, coarse_num_new_cells = 0;


    // --- STEP 1: angular projection. ---------------------------------------------------------------- //

    // Project data1 to higher quadrature if data2 has more angles.
    if ( data1.nq() < data2.nq() ) {

        second = &data2;
        first = new OrdinateFlux( data1, data2, data1.DG_degree );

        RelabelOperator relabel_operator( data1, *first );

        relabel_operator.Relabel( 1.0, data1, const_cast<OrdinateFlux &>( *first ) );

    // Project data2 to higher quadrature if data1 has more angles.
    } else if ( data1.nq() > data2.nq() ) {

        first = &data1;
        second = new OrdinateFlux( data2, data1, data2.DG_degree );

        RelabelOperator relabel_operator( data2, *second );

        relabel_operator.Relabel( 1.0, data2, const_cast<OrdinateFlux &>( *second ) );

    // Otherwise they have the same number of angles, so no projection is necessary.
    } else {

        first = &data1;
        second = &data2;
    }


    // --- STEP 2A: If solutions use the same spatial mess, just take the difference of their       --- //
    // ---          coefficients.                                                                   --- //

    if ( first->nx(0) == second->nx(0) ) {

        if ( first->DG_degree < second->DG_degree ) {

            coarse = first;
            fine = second;

        } else {

            fine = first;
            coarse = second;
        }

        result.Reconfigure( *fine, *fine, fine->DG_degree );

        #pragma omp parallel for
        for ( int64_t i = 1; i <= result.nx(0); ++i ) {
        for ( int64_t q = 0; q <  result.nq();  ++q ) {

            for ( int64_t d = 0; d <= coarse->DG_degree; ++d )
                result(q,i,d) = (*fine)(q,i,d) - (*coarse)(q,i,d);

            for ( int64_t d = coarse->DG_degree + 1; d <= fine->DG_degree; ++d )
                result(q,i,d) = (*fine)(q,i,d);
        }}

        goto cleanup;


    // --- STEP 2B: If spatial meshes are different, the compute spatial projection onto suitable   --- //
    // ---          mesh.                                                                           --- //

    } else if ( second->nx(0) > first->nx(0) ) {

        coarse = first;
        fine = second;

    } else {

        fine = first;
        coarse = second;
    }

    // Configure result object.
    {
        SpatialDomain result_domain( { lcm( coarse->nx(0), fine->nx(0) ) }, first->ax(), first->bx() );

        DomainDecomposition result_decomposition( result_domain );

        result.Reconfigure( result_decomposition, *fine, std::max( coarse->DG_degree, fine->DG_degree ) );
    }

    // Retrieve Gauss-Legendre quadrature of optimal efficiency.
    GL_order = result.DG_degree + 1;

    GL_nodes   = new double[ GL_order ];
    GL_weights = new double[ GL_order ];

    ComputeQuadrature( GL_order, GL_nodes, GL_weights, NodesType::GaussLegendre );


    // Macro for indexing into matrices of spatial projection operators.
    # define IPROJ(in,out,i,d,a) ( d + ((in)->DG_degree + 1)*( a + ((out).DG_degree + 1)*(i) ) )


    // Compute spatial projection operator for fine solution.
    fine_num_new_cells = result.nx(0) / fine->nx(0);

    fine_project = new double[ (fine->DG_degree + 1) * fine_num_new_cells * (result.DG_degree + 1) ];

    std::memset( fine_project, 0,
                 (fine->DG_degree + 1) * fine_num_new_cells * (result.DG_degree + 1) * sizeof(double) );

    for ( int64_t i = 0; i <  fine_num_new_cells; ++i ) { // Index of new subcell.
    for ( int64_t a = 0; a <= result.DG_degree ;  ++a ) { // Degree of new basis function.
    for ( int64_t d = 0; d <= fine->DG_degree;    ++d ) { // Degree of original basis function.
    for ( int64_t l = 0; l <  GL_order;           ++l ) { // Index of Gaussian quadrature node.

        fine_project[ IPROJ(fine,result,i,d,a) ]
            += GL_weights[l] * (2*a + 1) / 2.0 * Legendre( a, GL_nodes[l] )
               * Legendre( d, ( GL_nodes[l] + 2*i + 1 ) / fine_num_new_cells - 1.0 );
    }}}}

    // Compute spatial projection operator for coarse solution.
    coarse_num_new_cells = result.nx(0) / coarse->nx(0);

    coarse_project = new double[ (coarse->DG_degree + 1) * coarse_num_new_cells * (result.DG_degree + 1) ];

    std::memset( coarse_project, 0,
                 (coarse->DG_degree + 1) * coarse_num_new_cells * (result.DG_degree + 1) * sizeof(double) );

    for ( int64_t i = 0; i <  coarse_num_new_cells; ++i ) { // Index of new subcell.
    for ( int64_t a = 0; a <= result.DG_degree;     ++a ) { // Degree of new basis function.
    for ( int64_t d = 0; d <= coarse->DG_degree;    ++d ) { // Degree of original basis function.
    for ( int64_t l = 0; l <  GL_order;             ++l ) { // Index of Gaussian quadrature node.

        coarse_project[ IPROJ(coarse,result,i,d,a) ]
            += GL_weights[l] * (2*a + 1) / 2.0 * Legendre( a, GL_nodes[l] )
               * Legendre( d, ( GL_nodes[l] + 2*i + 1 ) / coarse_num_new_cells - 1.0 );
    }}}}


    // --- Do projections and compute difference. ----------------------------------------------------- //

    {
    const int64_t (& nx) [SPACE_DIMS] = result.nx();
    const int64_t nq = result.nq();

    # pragma omp parallel for collapse(2)
    for ( int64_t i = 0; i <  nx[0];            ++i ) {
    for ( int64_t q = 0; q <  nq;               ++q ) {
    for ( int64_t a = 0; a <= result.DG_degree; ++a ) {

        // Add fine projection.
        for ( int64_t d = 0; d <= fine->DG_degree; ++d )
            result(q,i+1,a) += fine_project[ IPROJ(fine,result, i % fine_num_new_cells, d,a) ]
                               * (*fine)( q, (i / fine_num_new_cells) + 1, d);

        // Subtract coarse projection.
        for ( int64_t d = 0; d <= coarse->DG_degree; ++d )
            result(q,i+1,a) -= coarse_project[ IPROJ(coarse,result, i % coarse_num_new_cells, d,a) ]
                               * (*coarse)( q, (i / coarse_num_new_cells) + 1, d);
    }}}
    }


    // --- Cleanup and return. ------------------------------------------------------------------------ //

cleanup:

    # undef IPROJ

    delete [] GL_nodes;
    delete [] GL_weights;
    delete [] coarse_project;
    delete [] fine_project;
}


//------------------------------------------------------------------------------------------------------------
//! \brief  Computes the angular moments of the ordinate flux. The moments are returned as
//!         RKDG::DensityFunction elements in an std::vector container, in order of moment degree.
//!
//! \return     Returns an std::vector container with OrdinateFlux::nq() DensityFunction objects, where each
//!             DensityFunction contains the spatial distribution of the given moment. The moments are sorted
//!             in degree order, so that, e.g., array access of the form \c [n] yields the angular moment of
//!             degree \f$ n \f$.
//------------------------------------------------------------------------------------------------------------
std::vector< RKDG::DensityFunction > OrdinateFlux::ComputeAngularMoments ( void ) const {

    std::vector< RKDG::DensityFunction > moments( this->nq() );

    for ( int64_t m = 0; m < this->nq(); ++m ) {

        moments[m].Reconfigure( *this, this->DG_degree );

        const double legendre_norm = std::sqrt( 2.0 / (2*m + 1) );

        const int64_t (& nx) [SPACE_DIMS] = this->nx();

        # pragma omp parallel for collapse(2)
        for ( int64_t i = 1; i <= nx[0];           ++i ) {
        for ( int64_t d = 0; d <= this->DG_degree; ++d ) {
        for ( int64_t q = 0; q <  this->nq();      ++q ) {

            (moments[m])(i,d) += this->w(q) * (*this)(q,i,d) * Legendre( m, this->xi(q) ) / legendre_norm;
        }}}
    }

    return moments;
}


# endif // if SPACE_DIMS == 1


//============================================================================================================
//=== TWO SPATIAL DIMENSIONS =================================================================================
//============================================================================================================

# if SPACE_DIMS == 2 || defined (DOXYCOMPILE)


//------------------------------------------------------------------------------------------------------------
//! \brief  Routine to obtain the value of the density function represented in an #OrdinateFlux object at a
//!         specified point and ordinate.
//!
//! \param[in]      x           \f$ x \f$ coordinate of point to obtain value at.
//! \param[in]      y           \f$ y \f$ coordinate of point to obtain value at.
//! \param[in]      q           Index of ordinate to obtain value at.
//!
//! \return     Returns the value of the density at the given point.
//------------------------------------------------------------------------------------------------------------
double OrdinateFlux::EvaluateAtPoint (

    const double x,
    const double y,
    const int64_t q

) const {

# if defined (STRICT_CHECK)

    if (    x < this->global_ax(0)
         || x > this->global_bx(0)
         || y < this->global_ax(1)
         || y > this->global_bx(1)
    ) {
        std::ostringstream error_message;
        error_message << std::scientific
                      << "Evaluation point ("
                      << x << "," << y
                      << ") is outside of RKDG::OrdinateFlux's domain of [ "
                      << this->global_ax(0) << ", " << this->global_bx(0)
                      << " ] × [ "
                      << this->global_ax(1) << ", " << this->global_bx(1)
                      << " ].\n";

        PRINT_ERROR( error_message.str().c_str() )
        throw std::range_error( error_message.str() );
    }

# endif // if defined (STRICT_CHECK)

    double ret_val = 0.0;
    int block_coords [SPACE_DIMS];
    int local_coords [SPACE_DIMS];

    PointToCell( {x,y}, block_coords, local_coords );

# if defined (ENABLE_MPI)

    // If point is inside MPI block of current process, compute the value at that point.
    if (    block_coords[0] == this->MPI_block_coords(0)
         && block_coords[1] == this->MPI_block_coords(1)
    )
# endif // if defined (ENABLE_MPI)
    {
        // The point (x,y) lies inside the spatial cell with indices 'i,j' of the MPI block.
        const int64_t i = local_coords[0];
        const int64_t j = local_coords[1];

        // Compute \bar{x} and \bar{y} to evaluate Legendre polynomials with.
        const double xBar = 2.0 * (x - this->ax(0)) / this->dx(0) - (2*i + 1);
        const double yBar = 2.0 * (y - this->ax(1)) / this->dx(1) - (2*j + 1);

        // Evaluate density function at given point.
        for ( int64_t d = 0; d <= this->DG_degree; ++d ) {
        for ( int64_t e = 0; e <= this->DG_degree; ++e ) {

            ret_val += (*this)(q,i+1,j+1,d,e) * Legendre( d, xBar ) * Legendre( e, yBar );
        }}

    # if defined (ENABLE_MPI)

        if ( Global::MPI_rank != 0 ) {

            MPI_Send( &ret_val, 1, MPI_DOUBLE, 0, 0, Global::MPI_cart_comm );
            ret_val = std::nan("0");
        }
    # endif // if defined (ENABLE_MPI)
    }

# if defined (ENABLE_MPI)

    else {

        if ( Global::MPI_rank == 0 ) {

            int recv_from_rank;
            MPI_Cart_rank( Global::MPI_cart_comm, block_coords, &recv_from_rank );

            MPI_Recv( &ret_val, 1, MPI_DOUBLE, recv_from_rank, 0, Global::MPI_cart_comm, MPI_STATUS_IGNORE );

        } else {  ret_val = std::nan("0");  }
    }

    MPI_Barrier( MPI_COMM_WORLD );

# endif // if defined (ENABLE_MPI)

    return ret_val;
}


//------------------------------------------------------------------------------------------------------------
//! \brief  Routine to obtain the value of the density function represented in an #OrdinateFlux object at a
//!         specified point and ordinate. The given point must lie within the local domain decomposition of
//!         the current MPI rank.
//!
//! \param[in]      x           \f$ x \f$ coordinate of point to obtain value at.
//! \param[in]      y           \f$ y \f$ coordinate of point to obtain value at.
//! \param[in]      q           Index of ordinate to obtain value at.
//!
//! \return     Returns the value of the density at the given point.
//------------------------------------------------------------------------------------------------------------
double OrdinateFlux::EvaluateAtLocalPoint (

    const double x,
    const double y,
    const int64_t q

) const {

# if defined (STRICT_CHECK)

    if (    x < this->global_ax(0)
         || x > this->global_bx(0)
         || y < this->global_ax(1)
         || y > this->global_bx(1)
    ) {
        std::ostringstream error_message;
        error_message << std::scientific
                      << "Evaluation point ("
                      << x << "," << y
                      << ") is outside of RKDG::OrdinateFlux's domain of [ "
                      << this->global_ax(0) << ", " << this->global_bx(0)
                      << " ] × [ "
                      << this->global_ax(1) << ", " << this->global_bx(1)
                      << " ].\n";

        PRINT_ERROR( error_message.str().c_str() )
        throw std::range_error( error_message.str() );
    }

# endif // if defined (STRICT_CHECK)

    double ret_val = 0.0;
    int block_coords [SPACE_DIMS];
    int local_coords [SPACE_DIMS];

    PointToCell( {x,y}, block_coords, local_coords );

# if defined (ENABLE_MPI)

    // Point needs to be inside the current MPI rank.
    if (    block_coords[0] == this->MPI_block_coords(0)
         && block_coords[1] == this->MPI_block_coords(1)
    ) {
        std::ostringstream error_message;
        error_message << std::scientific
                      << "Evaluation point " << x << " is inside MPI block ("
                      << block_coords[0] << "," << block_coords[1]
                      << ") but the current MPI rank has block coordinates ("
                      << this->MPI_block_coords(0) << "," << this->MPI_block_coords(1)
                      << ").\n";

        PRINT_ERROR( error_message.str().c_str() )
        throw std::range_error( error_message.str() );
    }
# endif // if defined (ENABLE_MPI)

    {
        // The point (x,y) lies inside the spatial cell with indices 'i,j' of the MPI block.
        const int64_t i = local_coords[0];
        const int64_t j = local_coords[1];

        // Compute \bar{x} and \bar{y} to evaluate Legendre polynomials with.
        const double xBar = 2.0 * (x - this->ax(0)) / this->dx(0) - (2*i + 1);
        const double yBar = 2.0 * (y - this->ax(1)) / this->dx(1) - (2*j + 1);

        // Evaluate density function at given point.
        for ( int64_t d = 0; d <= this->DG_degree; ++d ) {
        for ( int64_t e = 0; e <= this->DG_degree; ++e ) {

            ret_val += (*this)(q,i+1,j+1,d,e) * Legendre( d, xBar ) * Legendre( e, yBar );
        }}
    }

    return ret_val;
}


# endif // if SPACE_DIMS == 2


//============================================================================================================
//=== THREE SPATIAL DIMENSIONS ===============================================================================
//============================================================================================================

# if SPACE_DIMS == 3


//------------------------------------------------------------------------------------------------------------
//! \brief  Routine to obtain the value of the density function represented in an #OrdinateFlux object at a
//!         specified point and ordinate.
//!
//! \param[in]      x           \f$ x \f$ coordinate of point to obtain value at.
//! \param[in]      y           \f$ y \f$ coordinate of point to obtain value at.
//! \param[in]      z           \f$ z \f$ coordinate of point to obtain value at.
//! \param[in]      q           Index of ordinate to obtain value at.
//!
//! \return     Returns the value of the density at the given point.
//------------------------------------------------------------------------------------------------------------
double OrdinateFlux::EvaluateAtPoint (

    const double x,
    const double y,
    const double z,
    const int64_t q

) const {

# if defined (STRICT_CHECK)

    if (    x < this->global_ax(0)
         || x > this->global_bx(0)
         || y < this->global_ax(1)
         || y > this->global_bx(1)
         || z < this->global_ax(2)
         || z > this->global_bx(2)
    ) {
        std::ostringstream error_message;
        error_message << std::scientific
                      << "Evaluation point ("
                      << x << "," << y << "," << z
                      << ") is outside of RKDG::OrdinateFlux's domain of [ "
                      << this->global_ax(0) << ", " << this->global_bx(0)
                      << " ] × [ "
                      << this->global_ax(1) << ", " << this->global_bx(1)
                      << " ] × [ "
                      << this->global_ax(2) << ", " << this->global_bx(2)
                      << " ].\n";

        PRINT_ERROR( error_message.str().c_str() )
        throw std::range_error( error_message.str() );
    }

# endif // if defined (STRICT_CHECK)

    double ret_val = 0.0;
    int block_coords [SPACE_DIMS];
    int local_coords [SPACE_DIMS];

    PointToCell( {x,y,z}, block_coords, local_coords );

# if defined (ENABLE_MPI)

    // If point is inside MPI block of current process, compute the value at that point.
    if (    block_coords[0] == this->MPI_block_coords(0)
         && block_coords[1] == this->MPI_block_coords(1)
         && block_coords[2] == this->MPI_block_coords(2) )
# endif
    {
        // The point (x,y,z) lies inside the spatial cell with indices 'i,j,k' of the MPI block.
        const int64_t i = local_coords[0];
        const int64_t j = local_coords[1];
        const int64_t k = local_coords[2];

        // Compute \bar{x}, \bar{y}, and \bar{z} to evaluate Legendre polynomials with.
        const double xBar = 2.0 * (x - this->ax(0)) / this->dx(0) - (2*i + 1);
        const double yBar = 2.0 * (y - this->ax(1)) / this->dx(1) - (2*j + 1);
        const double zBar = 2.0 * (z - this->ax(2)) / this->dx(2) - (2*k + 1);

        // Evaluate density function at given point.
        for ( int64_t d = 0; d <= this->DG_degree; ++d ) {
        for ( int64_t e = 0; e <= this->DG_degree; ++e ) {
        for ( int64_t f = 0; f <= this->DG_degree; ++f ) {

            ret_val += (*this)(q,i+1,j+1,k+1,d,e,f)
                      * Legendre( d, xBar ) * Legendre( e, yBar ) * Legendre( f, zBar );
        }}}

    # if defined (ENABLE_MPI)

        if ( Global::MPI_rank != 0 ) {

            MPI_Send( &ret_val, 1, MPI_DOUBLE, 0, 0, Global::MPI_cart_comm );
            ret_val = std::nan("0");
        }
    # endif // if defined (ENABLE_MPI)
    }

# if defined (ENABLE_MPI)

    else {

        if ( Global::MPI_rank == 0 ) {

            int recv_from_rank;
            MPI_Cart_rank( Global::MPI_cart_comm, block_coords, &recv_from_rank );

            MPI_Recv( &ret_val, 1, MPI_DOUBLE, recv_from_rank, 0, Global::MPI_cart_comm, MPI_STATUS_IGNORE );

        } else {  ret_val = std::nan("0");  }
    }

    MPI_Barrier( MPI_COMM_WORLD );

# endif // if defined (ENABLE_MPI)

    return ret_val;
}


# endif // if SPACE_DIMS == 3


//============================================================================================================
//=== STATIC MEMBER FUNCTIONS ================================================================================
//============================================================================================================


//------------------------------------------------------------------------------------------------------------
//! \brief  Creates and returns an OrdinateFlux object.
//!
//! \param[in]  domain_decomposition    Contains parameters of spatial discretization.
//! \param[in]  ordinate_set            Contains parameters of angular discretization.
//! \param[in]  input_list              List of additional parameters.
//!
//! \return     Returns the newly created object.
//------------------------------------------------------------------------------------------------------------
OrdinateFlux * OrdinateFlux::Create (

    const DomainDecomposition & domain_decomposition,
    const OrdinateSet & ordinate_set,
    const ParameterList & input_list
) {

    int64_t DG_degree = GetDGDegreeX( input_list );

    return new OrdinateFlux( domain_decomposition, ordinate_set, DG_degree );
}


//------------------------------------------------------------------------------------------------------------
//! \brief  Two RKDG::OrdinateFlux objects "are matching" if all of their parameters are equal.
//!
//! \param[in]      first       The first of the two RKDG::OrdinateFlux objects to compare.
//! \param[in]      second      The second of the two RKDG::OrdinateFlux objects to compare.
//!
//! \return     Returns true if the RKDG::OrdinateFlux objects "are matching" and false otherwise.
//------------------------------------------------------------------------------------------------------------
bool OrdinateFlux::AreMatching (

    const OrdinateFlux & first,
    const OrdinateFlux & second
) {

    return      DensityFunction::AreMatching( first, second )
            &&  OrdinateSet::AreMatching( first, second )
            &&  first.DG_degree == second.DG_degree;
}


//------------------------------------------------------------------------------------------------------------
//! \brief  Swaps the contents of two RKDG::OrdinateFlux objects.
//!
//! Throws an error if the two objects are not matching.
//!
//! \todo   Implement RKDG::OrdinateFlux::swap for case when objects !AreMatching. (Print warning when this
//!         happens.)
//------------------------------------------------------------------------------------------------------------
void OrdinateFlux::swap (

    OrdinateFlux & thing1,
    OrdinateFlux & thing2
) {

    if ( AreMatching( thing1, thing2 ) ) {

        std::swap( thing1.density,                  thing2.density                  );
        std::swap( thing1.halo_cells_dirty,         thing2.halo_cells_dirty         );
        std::swap( thing1.upwind_halo_cells_dirty,  thing2.upwind_halo_cells_dirty  );

    } else {

        std::string error_message = "Parameter mismatch between RKDG::OrdinateFlux objects in '"
                                    + std::string(__func__) + "'.\n";

        PRINT_ERROR( error_message.c_str() )
        throw std::invalid_argument( error_message );
    }
}


# if SPACE_DIMS == 1 || defined (DOXYCOMPILE)

//------------------------------------------------------------------------------------------------------------
//! \brief  Computes the distances (with respect to the spatial measure) between the angular moments of two
//!         OrdinateFlux objects.
//!
//! The number of moments computed is determined by the angular resolution of the OrdinateFlux objects.
//!
//! \param[in]  data1   First OrdinateFlux object in difference.
//! \param[in]  data2   Second OrdinateFlux object in difference.
//! \param[in]  norm    (optional) <br/>
//!                     DensityFunction norm that induces the metric with which the distances between moments
//!                     are to be computed. Default is
//!
//! \return Returns an std::vector object containing the distances between each of the angular moments.
//------------------------------------------------------------------------------------------------------------
std::vector<double> OrdinateFlux::ComputeMomentDistances (

    const OrdinateFlux & data1,
    const OrdinateFlux & data2,
    double (RKDG::DensityFunction::*norm)() const   // = &DensityFunction::L2Norm
) {

    // Determine which OrdinateFlux object has the higher angular resolution.
    const OrdinateFlux & coarse = ( data1.nq() > data2.nq() ? data2 : data1 );
    const OrdinateFlux & fine   = ( data1.nq() > data2.nq() ? data1 : data2 );

    // Compute angular moments.
    std::vector< RKDG::DensityFunction > coarse_moments = coarse.ComputeAngularMoments();
    std::vector< RKDG::DensityFunction > fine_moments = fine.ComputeAngularMoments();

    // Compute difference between angular moments.
    std::vector< RKDG::DensityFunction > diff_moments( fine.nq() );
    std::vector<double> distances( fine.nq() );

    RKDG::DensityFunction zero( fine_moments.at(0), fine_moments.at(0).DG_degree );

    for ( uint64_t n = 0; n < coarse_moments.size(); ++n )
        RKDG::DensityFunction::ComputeDifference( fine_moments.at(n), coarse_moments.at(n),
                                                  diff_moments.at(n) );

    for ( uint64_t n = coarse_moments.size(); n < fine_moments.size(); ++n )
        RKDG::DensityFunction::ComputeDifference( zero, fine_moments.at(n), diff_moments.at(n) );

    // Compute norms of differences to obtain distances.
    for ( uint64_t n = 0; n < diff_moments.size(); ++n )
        distances.at(n) = (diff_moments.at(n).*norm)();

    return distances;
}

# endif // if SPACE_DIMS == 1


//============================================================================================================
//=== ADDITIONAL MEMBER FUNCTIONS ============================================================================
//============================================================================================================


//------------------------------------------------------------------------------------------------------------
//! \brief  Prints a summary of the object configuration to the logging interface.
//!
//! \param[in]  prefix  String prepended to each line of output, intended to be used for controlling
//!                     indentation of output.
//------------------------------------------------------------------------------------------------------------
void OrdinateFlux::Print (

    const std::string & prefix // = "  "

) const {

    PRINT_STATUS( "Executing RKDG::OrdinateFlux::%s.\n", __func__ )

    PRINT_LOG( "\n" )

    PRINT_LOG( "%s%-*s % " PRId64 "\n", prefix.c_str(),
               Global::col_width, "DG degree:", this->DG_degree )

    this->Abstract::OrdinateFlux::Print( prefix );

    PRINT_LOG( "\n" )

    PRINT_LOG( "%s%-*s % .2e\n", prefix.c_str(),
               Global::col_width, "L1 norm:", this->L1Norm() )

    PRINT_LOG( "%s%-*s % .2e\n", prefix.c_str(),
               Global::col_width, "L2 norm:", this->L2Norm() )

    PRINT_LOG( "%s%-*s % .2e\n", prefix.c_str(),
               Global::col_width, "Linf norm:", this->LinfNorm() )

    PRINT_LOG( "\n" )

    PRINT_LOG( "%s%-*s % .2e\n", prefix.c_str(),
               Global::col_width, "Total mass:", this->TotalMass() )
}


//------------------------------------------------------------------------------------------------------------
//! \brief  Updates ghost cells on reflecting boundaries.
//!
//! \todo   3D Implementation.
//------------------------------------------------------------------------------------------------------------
OrdinateFlux & OrdinateFlux::ReflectBoundaries ( void ) {

    PRINT_STATUS( "Updating ghost cells with reflecting boundary conditions.\n" )

    # pragma omp master
    Global::TMR_AF_boundary.Start();

    /*  Get thread count and thread IDs before entering parallel region below.
     *
     *  If this function is called from within a parallel region (as is the case for the mesh-sweeping
     *  routines) then each thread executes the function call in its own thread-local scope. Because of this,
     *  the loops below cannot be divided amongst the threads using standard worksharing constructs such as
     *  "omp for". Hence this implementation manually divides the loop iteration amongst the threads.
     *
     *  The manual division is performed using the thread count and thread IDs. There are two cases:
     *
     *      (i) When this function is called from a sequential portion of the code, the parallel region
     *          below is used to create a team of threads to execute the reflect operation. In this case
     *          the thread count and thread IDs must be determined from within the parallel region.
     *
     *     (ii) When this function is called from within a parallel region (the "outer" parallel region)
     *          the parallel region below (the "inner" parallel region) does not spawn additional threads.
     *          However, function calls to the OpenMP runtime library are always bound to the innermost
     *          parallel region, regardless of whether additional threads are spawned there. In order to
     *          obtain the correct thread counts and thread IDs, OpenMP library calls must be made before
     *          entering the innermost parallel region.
     */
    const bool in_parallel =
        # if defined (_OPENMP)
            omp_in_parallel();
        # else
            false;
        # endif

    const int64_t num_threads_outer =
        # if defined (_OPENMP)
            omp_get_num_threads();
        # else
            1;
        # endif

    const int64_t tid_outer =
        # if defined (_OPENMP)
            omp_get_thread_num();
        # else
            0;
        # endif

    // Spawn parallel threads here only if this function was called from a sequential section of the code.
    # pragma omp parallel if( !in_parallel )
    {
        /*  Get thread count and thread IDs for innermost parallel region. */
        const int64_t num_threads_inner =
            # if defined (_OPENMP)
                omp_get_num_threads();
            # else
                1;
            # endif

        const int64_t tid_inner =
            # if defined (_OPENMP)
                omp_get_thread_num();
            # else
                0;
            # endif

        /*  Determine "real" thread count and thread IDs.
         *
         *  If this function was called from within a parallel region, then the thread count and thread IDs
         *  that we want are num_threads_outer and tid_outer, respectively. If this function was called from
         *  a sequential portion of the code, then the thread count and thread IDs that we want are
         *  num_threads_inner and tid_inner, respectively.
         */
        const int64_t num_threads = ( in_parallel ? num_threads_outer : num_threads_inner );
        const int64_t tid         = ( in_parallel ? tid_outer         : tid_inner         );

    # if SPACE_DIMS == 1
        ((void)(num_threads));  // To avoid compiler warning in 1D code.
    # endif


        // --- X_Min ---------------------------------------------------------------------------------- //

        if (    BitmaskHasAll( Global::reflecting_boundaries, BoundaryEdge::X_Min )
        # if defined (ENABLE_MPI)
             && this->MPI_block_coords(0) == 0
        # endif
        ) {

        # if SPACE_DIMS == 1

            if ( tid == 0 ) {

                for ( int64_t q = 0; q <  this->nq();      ++q ) {
                for ( int64_t d = 0; d <= this->DG_degree; ++d ) {

                    (*this)(q,0,d) = neg1pow(d) * (*this)( this->nq() - q, 1,d);
                }}
            }

        # elif SPACE_DIMS == 2

            const int64_t Q = this->nx(1) / num_threads;    // Quotient and remainder, respectively.
            const int64_t R = this->nx(1) % num_threads;    //

            const int64_t j_start = 1 + tid * Q + std::min( tid, R );
            const int64_t j_stop  = j_start + Q - ( tid < R ? 0 : 1 );

            for ( int64_t quad : { 0, 1 } ) {
            for ( int64_t q = this->Quadrants(quad); q < this->Quadrants(quad + 1); ++q ) {

                for ( int64_t j = j_start; j <= j_stop; ++j ) {
                for ( int64_t d = 0; d <= this->DG_degree; ++d ) {
                for ( int64_t e = 0; e <= this->DG_degree; ++e ) {

                    (*this)(q, 0, j,d,e) = neg1pow(d) * (*this)( ReflectMap(0,q), 1, j,d,e);
                }}}
            }}

        # elif SPACE_DIMS == 3
            # warning "RKDG::OrdinateFlux::ReflectBoundaries not implemented in 3D."
        # endif // if SPACE_DIMS == ?
        }

        // --- X_Max ---------------------------------------------------------------------------------- //

        if (    BitmaskHasAll( Global::reflecting_boundaries, BoundaryEdge::X_Max )
        # if defined (ENABLE_MPI)
             && this->MPI_block_coords(0) == this->MPI_num_blocks(0) - 1
        # endif
        ) {

        # if SPACE_DIMS == 1

            if ( tid == 0 ) {

                for ( int64_t q = 0; q <  this->nq();      ++q ) {
                for ( int64_t d = 0; d <= this->DG_degree; ++d ) {

                    (*this)(q, this->nx(0) + 1, d) = neg1pow(d) * (*this)( this->nq() - q, this->nx(0), d);
                }}
            }

        # elif SPACE_DIMS == 2

            const int64_t Q = this->nx(1) / num_threads;    // Quotient and remainder, respectively.
            const int64_t R = this->nx(1) % num_threads;    //

            const int64_t j_start = 1 + tid * Q + std::min( tid, R );
            const int64_t j_stop  = j_start + Q - ( tid < R ? 0 : 1 );

            for ( int64_t quad : { 2, 3 } ) {
            for ( int64_t q = this->Quadrants(quad); q < this->Quadrants(quad + 1); ++q ) {

                for ( int64_t j = j_start; j <= j_stop; ++j ) {
                for ( int64_t d = 0; d <= this->DG_degree; ++d ) {
                for ( int64_t e = 0; e <= this->DG_degree; ++e ) {

                    (*this)(q, this->nx(0) + 1, j,d,e)
                        = neg1pow(d) * (*this)( ReflectMap(0,q), this->nx(0), j,d,e);
                }}}
            }}

        # elif SPACE_DIMS == 3
            # warning "RKDG::OrdinateFlux::ReflectBoundaries not implemented in 3D."
        # endif // if SPACE_DIMS == ?
        }

    # if SPACE_DIMS >= 2

        // --- Y_Min ---------------------------------------------------------------------------------- //

        if (    BitmaskHasAll( Global::reflecting_boundaries, BoundaryEdge::Y_Min )
        # if defined (ENABLE_MPI)
             && this->MPI_block_coords(1) == 0
        # endif
        ) {

        # if SPACE_DIMS == 2

            const int64_t Q = this->nx(0) / num_threads;    // Quotient and remainder, respectively.
            const int64_t R = this->nx(0) % num_threads;    //

            const int64_t i_start = 1 + tid * Q + std::min( tid, R );
            const int64_t i_stop  = i_start + Q - ( tid < R ? 0 : 1 );

            for ( int64_t quad : { 0, 2 } ) {
            for ( int64_t q = this->Quadrants(quad); q < this->Quadrants(quad + 1); ++q ) {

                for ( int64_t i = i_start; i <= i_stop; ++i ) {
                for ( int64_t d = 0; d <= this->DG_degree; ++d ) {
                for ( int64_t e = 0; e <= this->DG_degree; ++e ) {

                    (*this)(q,i, 0, d,e) = neg1pow(e) * (*this)( ReflectMap(1,q), i, 1, d,e);
                }}}
            }}

        # elif SPACE_DIMS == 3
            # warning "RKDG::OrdinateFlux::ReflectBoundaries not implemented in 3D."
        # endif // if SPACE_DIMS == ?
        }

        // --- Y_Max ---------------------------------------------------------------------------------- //

        if (    BitmaskHasAll( Global::reflecting_boundaries, BoundaryEdge::Y_Max )
        # if defined (ENABLE_MPI)
             && this->MPI_block_coords(1) == this->MPI_num_blocks(1) - 1
        # endif
        ) {

        # if SPACE_DIMS == 2

            const int64_t Q = this->nx(0) / num_threads;    // Quotient and remainder, respectively.
            const int64_t R = this->nx(0) % num_threads;    //

            const int64_t i_start = 1 + tid * Q + std::min( tid, R );
            const int64_t i_stop  = i_start + Q - ( tid < R ? 0 : 1 );

            for ( int64_t quad : { 0, 1, 2, 3 } ) {
            for ( int64_t q = this->Quadrants(quad); q < this->Quadrants(quad + 1); ++q ) {

                for ( int64_t i = i_start; i <= i_stop; ++i ) {
                for ( int64_t d = 0; d <= this->DG_degree; ++d ) {
                for ( int64_t e = 0; e <= this->DG_degree; ++e ) {

                    (*this)(q,i, this->nx(1) + 1, d,e)
                        = neg1pow(e) * (*this)( ReflectMap(1,q), i, this->nx(1), d,e);
                }}}
            }}

        # elif SPACE_DIMS == 3
            # warning "RKDG::OrdinateFlux::ReflectBoundaries not implemented in 3D."
        # endif // if SPACE_DIMS == ?
        }

    # endif // if SPACE_DIMS >= 2
    } // end parallel region.

    # pragma omp master
    Global::TMR_AF_boundary.Stop();

    return *this;
}


//------------------------------------------------------------------------------------------------------------
//! \brief  Returns the stride in the global array between subsequent spatial cells, moving in the specified
//!         dimension.
//!
//! \param[in]  dim     Dimension to compute stride for.
//!
//! \return     Returns the stride in the global array for the specified dimension.
//------------------------------------------------------------------------------------------------------------
int64_t OrdinateFlux::CellStride (

    const int64_t dim

) const {

    switch ( dim ) {

        case 0:
            return
            # if SPACE_DIMS == 1
                Index(0,1,0) - Index(0,0,0);
            # elif SPACE_DIMS == 2
                Index(0,1,0,0,0) - Index(0,0,0,0,0);
            # elif SPACE_DIMS == 3
                Index(0,1,0,0,0,0,0) - Index(0,0,0,0,0,0,0);
            # endif

    # if SPACE_DIMS >= 2

        case 1:
            return
            # if SPACE_DIMS == 2
                Index(0,0,1,0,0) - Index(0,0,0,0,0);
            # elif SPACE_DIMS == 3
                Index(0,0,1,0,0,0,0) - Index(0,0,0,0,0,0,0);
            # endif

    # endif // if SPACE_DIMS >= 2

    # if SPACE_DIMS == 3

        case 2:
            return Index(0,0,0,1,0,0,0) - Index(0,0,0,0,0,0,0);

    # endif // if SPACE_DIMS == 3

        default:
        {   std::string error_message =   "Spatial dimension "
                                        + std::to_string(dim)
                                        + " not in valid range [0,"
                                        + std::to_string( SPACE_DIMS )
                                        + ") for "
                                        + std::string(__func__)
                                        + ".\n";

            PRINT_ERROR( error_message.c_str() )
            throw std::out_of_range( error_message );
        }
    }
}


//------------------------------------------------------------------------------------------------------------
//! \brief  Returns a pointer to the location within DensityFunction::density where the coefficients for the
//!         given spatial cell begin.
//------------------------------------------------------------------------------------------------------------
double * OrdinateFlux::PointerAtCell (

    const int64_t i
# if SPACE_DIMS >= 2
,   const int64_t j
# endif
# if SPACE_DIMS == 3
,   const int64_t k
# endif

) const {

    return this->density + Index
                            # if SPACE_DIMS == 1
                                (0,i,0)
                            # elif SPACE_DIMS == 2
                                (0,i,j,0,0)
                            # elif SPACE_DIMS == 3
                                (0,i,j,k,0,0,0)
                            # endif
                            ;
}


//------------------------------------------------------------------------------------------------------------
//! \brief  Returns a pointer to the location within DensityFunction::density where the coefficients for the
//!         given spatial cell begin.
//------------------------------------------------------------------------------------------------------------
double * OrdinateFlux::PointerAtOrdinateCell (

    const int64_t q
,   const int64_t i
# if SPACE_DIMS >= 2
,   const int64_t j
# endif
# if SPACE_DIMS == 3
,   const int64_t k
# endif

) const {

    return this->density + Index
                            # if SPACE_DIMS == 1
                                (q,i,0)
                            # elif SPACE_DIMS == 2
                                (q,i,j,0,0)
                            # elif SPACE_DIMS == 3
                                (q,i,j,k,0,0,0)
                            # endif
                            ;
}


# if SPACE_DIMS == 2

//------------------------------------------------------------------------------------------------------------
//! \brief  Returns a pointer to the location within DensityFunction::density where the coefficients for the
//!         given spatial cell begin.
//------------------------------------------------------------------------------------------------------------
double * OrdinateFlux::PointerAtQuadrantCell (

    const int64_t quad,
    const int64_t i,
    const int64_t j

) const {

    return this->density + Index( Quadrants(quad), i,j,0,0);
}

# elif SPACE_DIMS == 3

//------------------------------------------------------------------------------------------------------------
//! \brief  Returns a pointer to the location within DensityFunction::density where the coefficients for the
//!         given spatial cell begin.
//------------------------------------------------------------------------------------------------------------
double * OrdinateFlux::PointerAtOctantCell (

    const int64_t oct,
    const int64_t i,
    const int64_t j

) const {

    return this->density + Index( Octants(oct), i,j,k,0,0,0);
}

# endif // if SPACE_DIMS == ?


//------------------------------------------------------------------------------------------------------------
//! \brief  Computes the total mass of the density function stored in an RKDG::OrdinateFlux object.
//!
//! The total mass is given by the integral of the density function over the spatial variables. Note that
//! due to the use of Legendre polynomial bases for the DG spatial approximation, the total mass is obtained
//! by integrating over the angular variable with the discrete ordinates quadrature, adding up the degree zero
//! coefficients of the spatial discretization, and multiplying by the mesh size.
//!
//! \return     Returns the total mass in the given object.
//------------------------------------------------------------------------------------------------------------
double OrdinateFlux::TotalMass ( void ) const {

    PRINT_STATUS( "Computing total mass of RKDG::OrdinateFlux object.\n" )

    double mass = 0.0;

# if SPACE_DIMS == 1

    for ( int64_t i = 1; i <= this->nx(0); ++i ) {
    for ( int64_t q = 0; q <  this->nq();  ++q ) {

        mass += this->w(q) * (*this)(q,i,0);
    }}

# elif SPACE_DIMS == 2

    for ( int64_t i = 1; i <= this->nx(0); ++i ) {
    for ( int64_t j = 1; j <= this->nx(1); ++j ) {
    for ( int64_t q = 0; q <  this->nq();  ++q ) {

        mass += this->w(q) * (*this)(q,i,j,0,0);
    }}}

# elif SPACE_DIMS == 3

    for ( int64_t i = 1; i <= this->nx(0); ++i ) {
    for ( int64_t j = 1; j <= this->nx(1); ++j ) {
    for ( int64_t k = 1; k <= this->nx(2); ++k ) {
    for ( int64_t q = 0; q <  this->nq();  ++q ) {

        mass += this->w(q) * (*this)(q,i,j,k,0,0,0);
    }}}}

# endif // if SPACE_DIMS == ?

# if defined (ENABLE_MPI)
    MPI_Allreduce( MPI_IN_PLACE, &mass, 1, MPI_DOUBLE, MPI_SUM, Global::MPI_cart_comm );
# endif

    for ( int64_t i = 0; i < SPACE_DIMS; ++i ) {  mass *= this->dx(i);  }

    return mass;
}


//------------------------------------------------------------------------------------------------------------
//! \brief  Computes the \f$ L^1( d\vec{x} \, d\vec{\Omega} ) \f$ norm of the density function stored in an
//!         RKDG::OrdinateFlux object.
//!
//! \return     Returns the \f$ L^1( d\vec{x} \, d\vec{\Omega} ) \f$ norm of the density function.
//------------------------------------------------------------------------------------------------------------
double OrdinateFlux::L1Norm ( void ) const {

    PRINT_STATUS( "Computing L1 norm of RKDG::OrdinateFlux object.\n" )

    // Retrieve Gauss-Legendre quadrature.
    const int64_t GL_order = this->DG_degree + 1;

    double * const GL_nodes   = new double[ GL_order ];
    double * const GL_weights = new double[ GL_order ];

    ComputeQuadrature( GL_order, GL_nodes, GL_weights, NodesType::GaussLegendre );


    // Compute norm using quadrature.
    double norm = 0.0;

    const int64_t (& nx) [SPACE_DIMS] = this->nx();

# if SPACE_DIMS == 1

    # pragma omp parallel for schedule(static) reduction(+:norm)
    for ( int64_t i = 1; i <= nx[0];    ++i ) {
    for ( int64_t l = 0; l <  GL_order; ++l ) {

        double temp = 0.0;

        // Value of quadrature point scaled to current cell.
        const double x_l = this->ax(0) + this->dx(0) * (i - 0.5) + this->dx(0) * GL_nodes[l] / 2.0;

        for ( int64_t q = 0; q < this->nq(); ++q )
            temp += this->w(q) * std::abs( EvaluateAtLocalPoint( x_l, q ) );

        norm += GL_weights[l] * this->dx(0) * temp / 2.0;
    }}

# elif SPACE_DIMS == 2

    # pragma omp parallel for collapse(2) schedule(static) reduction(+:norm)
    for ( int64_t i = 1; i <= nx[0];    ++i ) {
    for ( int64_t j = 1; j <= nx[1];    ++j ) {
    for ( int64_t l = 0; l <  GL_order; ++l ) {
    for ( int64_t m = 0; m <  GL_order; ++m ) {

        double temp = 0.0;

        // Value of quadrature point scaled to current cell.
        const double x_l = this->ax(0) + this->dx(0) * (i - 0.5) + this->dx(0) * GL_nodes[l] / 2.0;
        const double y_m = this->ax(1) + this->dx(1) * (j - 0.5) + this->dx(1) * GL_nodes[m] / 2.0;

        for ( int64_t q = 0; q < this->nq(); ++q )
            temp += this->w(q) * std::abs( EvaluateAtLocalPoint( x_l, y_m, q ) );

        norm += GL_weights[l] * GL_weights[m] * this->dx(0) * this->dx(1) * temp / 4.0;
    }}}}

# elif SPACE_DIMS == 3

    # pragma omp parallel for collapse(3) schedule(static) reduction(+:norm)
    for ( int64_t i = 1; i <= nx[0];    ++i ) {
    for ( int64_t j = 1; j <= nx[1];    ++j ) {
    for ( int64_t k = 1; k <= nx[2];    ++k ) {
    for ( int64_t l = 0; l <  GL_order; ++l ) {
    for ( int64_t m = 0; m <  GL_order; ++m ) {
    for ( int64_t n = 0; n <  GL_order; ++n ) {

        double temp = 0.0;

        // Value of quadrature point scaled to current cell.
        const double x_l = this->ax(0) + this->dx(0) * (i - 0.5) + this->dx(0) * GL_nodes[l] / 2.0;
        const double y_m = this->ax(1) + this->dx(1) * (j - 0.5) + this->dx(1) * GL_nodes[m] / 2.0;
        const double z_n = this->ax(2) + this->dx(2) * (k - 0.5) + this->dx(2) * GL_nodes[n] / 2.0;

        for ( int64_t q = 0; q < this->nq(); ++q )
            temp += this->w(q) * std::abs( EvaluateAtLocalPoint( x_l, y_m, z_n, q ) );

        norm += GL_weights[l] * GL_weights[m] * GL_weights[n]
                * this->dx(0) * this->dx(1) * this->dx(2) * temp / 8.0;
    }}}}}}

# endif // if SPACE_DIMS == ?

# if defined (ENABLE_MPI)
    MPI_Allreduce( MPI_IN_PLACE, &norm, 1, MPI_DOUBLE, MPI_SUM, Global::MPI_cart_comm );
# endif

    // Cleanup Gauss-Legendre quadrature.
    delete [] GL_nodes;
    delete [] GL_weights;

    return norm;
}


//------------------------------------------------------------------------------------------------------------
//! \brief  Computes the \f$ L^2( d\vec{x} \, d\vec{\Omega} ) \f$ norm of the density function stored in an
//!         RKDG::OrdinateFlux object.
//!
//! \todo   Use Parseval's identity.
//!
//! \return     Returns the \f$ L^2( d\vec{x} \, d\vec{\Omega} ) \f$ norm of the density function.
//------------------------------------------------------------------------------------------------------------
double OrdinateFlux::L2Norm ( void ) const {

    PRINT_STATUS( "Computing L2 norm of RKDG::OrdinateFlux object.\n" )

    // Retrieve Gauss-Legendre quadrature.
    const int64_t GL_order = this->DG_degree + 1;

    double * const GL_nodes   = new double[ GL_order ];
    double * const GL_weights = new double[ GL_order ];

    ComputeQuadrature( GL_order, GL_nodes, GL_weights, NodesType::GaussLegendre );


    // Compute norm using quadrature.
    double norm = 0.0;

    const int64_t (& nx) [SPACE_DIMS] = this->nx();

# if SPACE_DIMS == 1

    # pragma omp parallel for schedule(static) reduction(+:norm)
    for ( int64_t i = 1; i <= nx[0];    ++i ) {
    for ( int64_t l = 0; l <  GL_order; ++l ) {

        double temp = 0.0;

        // Value of quadrature point scaled to current cell.
        const double x_l = this->ax(0) + this->dx(0) * (i - 0.5) + this->dx(0) * GL_nodes[l] / 2.0;

        for ( int64_t q = 0; q < this->nq(); ++q ) {

            const double pt_val = EvaluateAtLocalPoint( x_l, q );
            temp += this->w(q) * pt_val * pt_val;
        }

        norm += GL_weights[l] * this->dx(0) * temp / 2.0;
    }}

# elif SPACE_DIMS == 2

    # pragma omp parallel for collapse(2) schedule(static) reduction(+:norm)
    for ( int64_t i = 1; i <= nx[0];    ++i ) {
    for ( int64_t j = 1; j <= nx[1];    ++j ) {
    for ( int64_t l = 0; l <  GL_order; ++l ) {
    for ( int64_t m = 0; m <  GL_order; ++m ) {

        double temp = 0.0;

        // Value of quadrature point scaled to current cell.
        const double x_l = this->ax(0) + this->dx(0) * (i - 0.5) + this->dx(0) * GL_nodes[l] / 2.0;
        const double y_m = this->ax(1) + this->dx(1) * (j - 0.5) + this->dx(1) * GL_nodes[m] / 2.0;

        for ( int64_t q = 0; q < this->nq(); ++q ) {

            const double pt_val = EvaluateAtLocalPoint( x_l, y_m, q );
            temp += this->w(q) * pt_val * pt_val;
        }

        norm += GL_weights[l] * GL_weights[m] * this->dx(0) * this->dx(1) * temp / 4.0;
    }}}}

# elif SPACE_DIMS == 3

    # pragma omp parallel for collapse(3) schedule(static) reduction(+:norm)
    for ( int64_t i = 1; i <= nx[0];    ++i ) {
    for ( int64_t j = 1; j <= nx[1];    ++j ) {
    for ( int64_t k = 1; k <= nx[2];    ++k ) {
    for ( int64_t l = 0; l <  GL_order; ++l ) {
    for ( int64_t m = 0; m <  GL_order; ++m ) {
    for ( int64_t n = 0; n <  GL_order; ++n ) {

        double temp = 0.0;

        // Value of quadrature point scaled to current cell.
        const double x_l = this->ax(0) + this->dx(0) * (i - 0.5) + this->dx(0) * GL_nodes[l] / 2.0;
        const double y_m = this->ax(1) + this->dx(1) * (j - 0.5) + this->dx(1) * GL_nodes[m] / 2.0;
        const double z_n = this->ax(2) + this->dx(2) * (k - 0.5) + this->dx(2) * GL_nodes[n] / 2.0;

        for ( int64_t q = 0; q < this->nq; ++q ) {

            const double pt_val = EvaluateAtLocalPoint( x_l, y_m, z_n, q );
            temp += this->w[q] * pt_val * pt_val;
        }

        norm += GL_weights[l] * GL_weights[m] * GL_weights[n]
                * this->dx(0) * this->dx(1) * this->dx(2) * temp / 8.0;
    }}}}}}

# endif // if SPACE_DIMS == ?

# if defined (ENABLE_MPI)
    MPI_Allreduce( MPI_IN_PLACE, &norm, 1, MPI_DOUBLE, MPI_SUM, Global::MPI_cart_comm );
# endif

    // Cleanup Gauss-Legendre quadrature.
    delete [] GL_nodes;
    delete [] GL_weights;

    return std::sqrt( norm );
}


//------------------------------------------------------------------------------------------------------------
//! \brief  Computes the \f$ L^{\infty}( d\vec{x} \, d\vec{\Omega} ) \f$ norm of the density function stored
//!         in an RKDG::OrdinateFlux object.
//!
//! \param[in]  tol     (optional)
//!                     Relative tolerance with with to compute \f$ L^{\infty} \f$ norm.
//!                     The default value of 1e-3 computes three digits of accuracy.
//!
//! \return     Returns the \f$ L^{\infty}( d\vec{x} \, d\vec{\Omega} ) \f$ norm of the density function.
//------------------------------------------------------------------------------------------------------------
double OrdinateFlux::LinfNorm (

    const double tol    // = 1e-3

) const {

    PRINT_STATUS( "Computing L-inf norm of RKDG::OrdinateFlux object.\n" )

    // Stores approximation of norm.
    double norm = -1.0;

    const int64_t (& nx) [SPACE_DIMS] = this->nx();
    const int64_t nq = this->nq();

# if SPACE_DIMS == 1

    # if _OPENMP >= 201107
        # pragma omp parallel for collapse(2) schedule(dynamic) reduction(max:norm)
    # endif
    for ( int64_t i = 1; i <= nx[0]; ++i ) {
    for ( int64_t q = 0; q <  nq;    ++q ) {

        int64_t part_size = 1;                                  // Number of points per dimension to sample.
        double c_norm     = std::numeric_limits<double>::max(); // Coarse approximation of norm.
        double f_norm     = std::numeric_limits<double>::max(); // Fine approximation of norm.
        double norm_diff  = std::numeric_limits<double>::max(); // Difference between coarse and fine approximations.

        do {
            c_norm = f_norm;
            f_norm = 0.0;

            // Use twice as many sample points per dimension as the previous iteration.
            part_size *= 2;
            const double pt_width = 2.0 / (part_size - 1.0);

            // Find maximum magnitude of point values in cell.
            for ( int64_t l = 0; l < part_size; ++l ) {

                // Point to evaluate Legendre polynomials at.
                const double x_L = -1.0 + pt_width * l;

                double temp = 0.0;

                // Evaluate point by going through Legendre recurrence relation only once.
                double pn  = 1.0;
                double pn1 = 0.0;
                double pn2;

                for ( int64_t d = 0; d <= this->DG_degree; ++d ) {

                    temp += (*this)(q,i,d) * pn;

                    pn2 = pn1;
                    pn1 = pn;
                    pn  = ( (2*d + 1) * x_L * pn1 - d * pn2 ) / (d + 1);
                }

                f_norm = std::max( f_norm, std::abs(temp) );
            }

            norm_diff = std::abs( f_norm - c_norm );

        //
        // Continue looping, doubling the mesh at each iteration, until first three significant digits of
        // successive approximations to norm of current cell agree.
        //
        } while ( (norm_diff / f_norm) > tol );

        norm = std::max( norm, f_norm );
    }}

# elif SPACE_DIMS == 2

    const size_t bytes_allocd = (this->DG_degree + 1)*(this->DG_degree + 1) * sizeof(double) * LINF_SIMD_LEN;

    # if _OPENMP >= 201107
        # pragma omp parallel reduction(max:norm)
    # endif
    {
        // Allocate temporary memory for traversing recurrence relations.
    # if defined (ENABLE_HWLOC) && defined (USE_HWLOC_ALLOC)

        hwloc_cpuset_t thread_mask = hwloc_bitmap_alloc();
        hwloc_get_cpubind( Global::machine_topology, thread_mask, HWLOC_CPUBIND_THREAD );

        double * const recurr_coeff = (double *)
            hwloc_alloc_membind( Global::machine_topology, bytes_allocd, thread_mask, HWLOC_MEMBIND_BIND,
                                 HWLOC_MEMBIND_THREAD );

        hwloc_bitmap_free( thread_mask );

    # else // if defined (ENABLE_HWLOC) && defined (USE_HWLOC_ALLOC)

        double * const recurr_coeff = (double *)
            # if defined (USE_ALIGNED_ALLOC)
                aligned_alloc( ALIGNED_ALLOC_ALIGNMENT, bytes_allocd );
            # else
                std::malloc( bytes_allocd );
            # endif

    # endif // if defined (ENABLE_HWLOC) && defined (USE_HWLOC_ALLOC)

        # if defined (ENABLE_SIMD_BLOCKING)

            # define IRECURR(d,e,l) recurr_coeff[ (l) + LINF_SIMD_LEN*( (e) + (this->DG_degree + 1)*(d) ) ]

            std::vector<std::tuple<double,double>> sample_points;

        # else // if defined (ENABLE_SIMD_BLOCKING)

            # define IRECURR(d,e) recurr_coeff[ (e) + (this->DG_degree + 1)*(d) ]

        # endif // if defined (ENABLE_SIMD_BLOCKING)

        # if _OPENMP >= 201107
            # pragma omp for collapse(3) schedule(dynamic)
        # endif
        for ( int64_t i = 1; i <= nx[0]; ++i ) {
        for ( int64_t j = 1; j <= nx[1]; ++j ) {
        for ( int64_t q = 0; q <  nq;    ++q ) {

            int64_t part_size = 1;                                  // Number of points per dimension to sample.
            double c_norm     = std::numeric_limits<double>::max(); // Coarse approximation of norm.
            double f_norm     = std::numeric_limits<double>::max(); // Fine approximation of norm.
            double norm_diff  = std::numeric_limits<double>::max(); // Difference between coarse and fine approximations.

            do {
                c_norm = f_norm;
                f_norm = 0.0;

                // Use twice as many sample points per dimension as the previous iteration.
                part_size *= 2;
                const double pt_width = 2.0 / (part_size - 1.0);

            # if defined (ENABLE_SIMD_BLOCKING)

                sample_points.clear();
                sample_points.reserve( part_size * part_size );

                for ( int64_t l = 0; l < part_size; ++l ) {
                for ( int64_t m = 0; m < part_size; ++m ) {

                    const double x = -1.0 + pt_width * l;
                    const double y = -1.0 + pt_width * m;

                    sample_points.push_back( std::make_tuple(x,y) );
                }}

                for ( int64_t blk = 0; blk < part_size * part_size; blk += LINF_SIMD_LEN ) {

                    const int64_t blk_len = std::min( ((int64_t)(LINF_SIMD_LEN)), part_size * part_size - blk );

                    alignas( LINF_SIMD_LEN * sizeof(double) ) double x [LINF_SIMD_LEN] = {};
                    alignas( LINF_SIMD_LEN * sizeof(double) ) double y [LINF_SIMD_LEN] = {};

                    for ( int64_t l = 0; l < blk_len; ++l ) {

                        x[l] = std::get<0>( sample_points[blk + l] );
                        y[l] = std::get<1>( sample_points[blk + l] );
                    }

                    //
                    // Evaluate Legendre polynomials using dynamic programming approach to tensor product of
                    // recurrence relations.
                    //

                    // Handle piecewise constant approximations as a special case.
                    //
                    // Note: Do not need to consider vector of points since they all yield the same value.
                    //
                    if ( this->DG_degree == 0 ) {

                        f_norm = std::max( f_norm, std::abs( (*this)(q,i,j,0,0) ) );
                        continue;
                    }

                    // Phase 1: base cases.
                    # pragma omp simd aligned( recurr_coeff, x, y : LINF_SIMD_LEN * sizeof(double) )
                    for ( int64_t l = 0; l < LINF_SIMD_LEN; ++l ) {

                        IRECURR(0,0,l) = 1.0;
                        IRECURR(1,0,l) = x[l];
                        IRECURR(0,1,l) = y[l];
                        IRECURR(1,1,l) = x[l] * y[l];
                    }

                    // Phase 2: edge cases.
                    for ( int64_t d = 2; d <= this->DG_degree; ++d ) {
                    # pragma omp simd aligned( recurr_coeff, x, y : LINF_SIMD_LEN * sizeof(double) )
                    for ( int64_t l = 0; l < LINF_SIMD_LEN; ++l ) {

                        IRECURR(d,0,l) = ( (2*d - 1) * x[l] * IRECURR(d-1,0,l) - (d-1) * IRECURR(d-2,0,l) ) / d;
                        IRECURR(d,1,l) = y[l] * IRECURR(d,0,l);
                    }}

                    for ( int64_t e = 2; e <= this->DG_degree; ++e ) {
                    # pragma omp simd aligned( recurr_coeff, x, y : LINF_SIMD_LEN * sizeof(double) )
                    for ( int64_t l = 0; l < LINF_SIMD_LEN; ++l ) {

                        IRECURR(0,e,l) = ( (2*e - 1) * y[l] * IRECURR(0,e-1,l) - (e-1) * IRECURR(0,e-2,l) ) / e;
                        IRECURR(1,e,l) = x[l] * IRECURR(0,e,l);
                    }}

                    // Phase 3: general case.
                    for ( int64_t d = 2; d <= this->DG_degree; ++d ) {
                    for ( int64_t e = 2; e <= this->DG_degree; ++e ) {
                    # pragma omp simd aligned( recurr_coeff, x, y : LINF_SIMD_LEN * sizeof(double) )
                    for ( int64_t l = 0; l < LINF_SIMD_LEN; ++l ) {

                        IRECURR(d,e,l) = (   (2*d - 1)*(2*e - 1) * x[l] * y[l] * IRECURR(d-1,e-1,l)
                                           - (d-1)    *(2*e - 1)        * y[l] * IRECURR(d-2,e-1,l)
                                           - (2*d - 1)*(e-1)     * x[l]        * IRECURR(d-1,e-2,l)
                                           + (d-1)    *(e-1)                   * IRECURR(d-2,e-2,l) ) / (d*e);
                    }}}

                    // Phase 4: evaluate expansion.
                    # pragma omp simd aligned( x : LINF_SIMD_LEN * sizeof(double) )
                    for ( int64_t l = 0; l < LINF_SIMD_LEN; ++l )
                        x[l] = 0.0;

                    for ( int64_t d = 0; d <= this->DG_degree; ++d ) {
                    for ( int64_t e = 0; e <= this->DG_degree; ++e ) {
                    # pragma omp simd aligned( recurr_coeff : LINF_SIMD_LEN * sizeof(double) )
                    for ( int64_t l = 0; l < LINF_SIMD_LEN; ++l ) {

                        x[l] += (*this)(q,i,j,d,e) * IRECURR(d,e,l);
                    }}}

                    for ( int64_t l = 0; l < blk_len; ++l )
                        f_norm = std::max( f_norm, std::abs( x[l] ) );
                }

            # else // if defined (ENABLE_SIMD_BLOCKING)

                for ( int64_t l = 0; l < part_size; ++l ) {
                for ( int64_t m = 0; m < part_size; ++m ) {

                    // Point to evaluate Legendre polynomials at.
                    const double x_L = -1.0 + pt_width * l;
                    const double y_L = -1.0 + pt_width * m;

                    double temp = 0.0;

                    // Evaluate point using Legendre function.
                # if 0
                    for ( int64_t d = 0; d <= this->DG_degree; ++d ) {
                    for ( int64_t e = 0; e <= this->DG_degree; ++e ) {

                        temp += (*this)(i,j,d,e) * Legendre( d, x_L ) * Legendre( e, y_L );
                    }}
                # endif

                    //
                    // Evaluate Legendre polynomials using dynamic programming approach to tensor product
                    // of recurrence relations.
                    //

                    // Handle piecewise constant approximations as a special case.
                    if ( this->DG_degree == 0 ) {

                        f_norm = std::max( f_norm, std::abs( (*this)(q,i,j,0,0) ) );
                        continue;
                    }

                    // Phase 1: base cases.
                    IRECURR(0,0) = 1.0;
                    IRECURR(1,0) = x_L;
                    IRECURR(0,1) = y_L;
                    IRECURR(1,1) = x_L * y_L;

                    // Phase 2: edge cases.
                    for ( int64_t d = 2; d <= this->DG_degree; ++d ) {

                        IRECURR(d,0) = ( (2*d - 1) * x_L * IRECURR(d-1,0) - (d-1) * IRECURR(d-2,0) ) / d;
                        IRECURR(d,1) = y_L * IRECURR(d,0);
                    }

                    for ( int64_t e = 2; e <= this->DG_degree; ++e ) {

                        IRECURR(0,e) = ( (2*e - 1) * y_L * IRECURR(0,e-1) - (e-1) * IRECURR(0,e-2) ) / e;
                        IRECURR(1,e) = x_L * IRECURR(0,e);
                    }

                    // Phase 3: general case.
                    for ( int64_t d = 2; d <= this->DG_degree; ++d ) {
                    for ( int64_t e = 2; e <= this->DG_degree; ++e ) {

                        IRECURR(d,e) = (   (2*d - 1)*(2*e - 1) * x_L * y_L * IRECURR(d-1,e-1)
                                         - (d-1)    *(2*e - 1)       * y_L * IRECURR(d-2,e-1)
                                         - (2*d - 1)*(e-1)     * x_L       * IRECURR(d-1,e-2)
                                         + (d-1)    *(e-1)                 * IRECURR(d-2,e-2) ) / (d*e);
                    }}

                    // Phase 4: evaluate expansion.
                    for ( int64_t d = 0; d <= this->DG_degree; ++d ) {
                    for ( int64_t e = 0; e <= this->DG_degree; ++e ) {

                        temp += (*this)(q,i,j,d,e) * IRECURR(d,e);
                    }}

                    f_norm = std::max( f_norm, std::abs(temp) );
                }}

        # endif // if defined (ENABLE_SIMD_BLOCKING)

                norm_diff = std::abs( f_norm - c_norm );

            //
            // Continue looping, doubling the sample mesh at each iteration, until first three significant digits
            // of successive approximations to norm of current cell agree.
            //
            } while ( (norm_diff / f_norm) > tol );

            norm = std::max( norm, f_norm );
        }}}

    # undef IRECURR

    # if defined (ENABLE_HWLOC) && defined (USE_HWLOC_ALLOC)
        hwloc_free( Global::machine_topology, recurr_coeff, bytes_allocd );
    # else
        std::free( recurr_coeff );
    # endif
    }

# elif SPACE_DIMS == 3
    # warning "3D version of RKDG::DensityFunction LinfNorm not implemented."
# endif // if SPACE_DIMS == ?

# if defined (ENABLE_MPI)
    MPI_Allreduce( MPI_IN_PLACE, &norm, 1, MPI_DOUBLE, MPI_MAX, Global::MPI_cart_comm );
# endif

    return norm;
}


//------------------------------------------------------------------------------------------------------------
//! \brief  Routine to output a grid of sample points from the OrdinateFlux density function for plotting.
//!
//! In one spatial dimensions, this routine outputs points from the phase-space distribution (i.e., dependent
//! on space and angle). In higher dimensions, this routine outputs point values of the scalar flux
//! distribution.
//!
//! \attention  Output is (destructively!) written to the file with name given by \pp{filename}.
//!
//! \param[in]      filename    Name of file in which to write output.
//------------------------------------------------------------------------------------------------------------
void OrdinateFlux::OutputPlot (

    const std::string filename

) const {

# if SPACE_DIMS == 1

    const char mode [] = "w";
    std::FILE * outfp = nullptr;

    // Open output file.
# if defined (ENABLE_MPI)
    if ( Global::MPI_rank == 0 )
# endif
    {
        outfp = std::fopen( filename.c_str(), mode );

        if ( outfp == nullptr ) {

            std::string error_message = "Failed to open file '" + filename + "' in '" + std::string(__func__)
                                        + "'.\n";

            PRINT_ERROR( error_message.c_str() )
            throw std::runtime_error( error_message );
        }
    }

    // Number of sample points in spatial dimension.
    const int64_t nx = 4096;

    // Distance between sample points.
    const double dx = (this->global_bx(0) - this->global_ax(0)) / nx;

    for ( int64_t q = 0; q < this->nq(); ++q ) {

        for ( int64_t i = 0; i < nx; ++i ) {

            // Compute value of spatial point to sample density at.
            const double x_i = this->global_ax(0) + dx * (i + 0.5);

            // Evaluate density function at sample spatial point.
            double point_val = this->EvaluateAtPoint( x_i, q );

        # if defined (ENABLE_MPI)
            if ( Global::MPI_rank == 0 )
        # endif
                std::fprintf( outfp, "% .4e % .4e % .4e\n", x_i, this->xi(q), point_val );
        }

        std::fprintf( outfp, "\n" );
    }

# if defined (ENABLE_MPI)
    MPI_Barrier( Global::MPI_cart_comm );

    if ( Global::MPI_rank == 0 )
# endif
    std::fclose( outfp );

# else // if SPACE_DIMS == 1

    ScalarFlux * temp = new ScalarFlux( *this, this->DG_degree );

    TransportOperator::Pmv( 1.0, 0.0, *this, *temp );
    temp->OutputPlot( filename );

    delete temp;

# endif // if SPACE_DIMS == ?
}


//------------------------------------------------------------------------------------------------------------
//! \brief  Routine to output a grid of points sampling the difference between two OrdinateFlux objects for
//!         plotting.
//!
//! \attention  Output is (destructively!) written to the file with name given by \pp{filename}.
//!
//! \todo   Check if OrdinateFlux objects are matching.
//!
//! \param[in]      data1       First RKDG::OrdinateFlux object in difference.
//! \param[in]      data2       Second RKDG::OrdinateFlux object in difference.
//! \param[in]      filename    Name of file in which to write output.
//------------------------------------------------------------------------------------------------------------
void OrdinateFlux::OutputDiffPlot (

    const OrdinateFlux & data1,
    const OrdinateFlux & data2,
    const std::string filename

) const {

# if SPACE_DIMS == 1

    const char mode [] = "w";
    FILE * outfp = nullptr;

    // Open output file.
# if defined (ENABLE_MPI)
    if ( Global::MPI_rank == 0 )
# endif
    {
        outfp = std::fopen( filename.c_str(), mode );

        if ( outfp == nullptr ) {

            std::string error_message = "Failed to open file '" + filename + "' in '" + std::string(__func__)
                                        + "'.\n";

            PRINT_ERROR( error_message.c_str() )
            throw std::runtime_error( error_message );
        }
    }

    // Number of sample points in spatial dimension.
    const int64_t nx = 4096;

    // Distance between sample points.
    const double dx = (data1.global_bx(0) - data1.global_ax(0)) / nx;

    for ( int64_t q = 0; q < data1.nq(); ++q ) {

        for ( int64_t i = 0; i < nx; ++i ) {

            // Compute value of spatial point to sample density at.
            const double x_i = data1.global_ax(0) + dx * (i + 0.5);

            // Evaluate density function at sample spatial point.
            double diff_val = data1.EvaluateAtPoint( x_i, q ) - data2.EvaluateAtPoint( x_i, q );

        # if defined (ENABLE_MPI)
            if ( Global::MPI_rank == 0 )
        # endif
                std::fprintf( outfp, "% .4e % .4e % .4e\n", x_i, data1.xi(q), diff_val );
        }

        std::fprintf( outfp, "\n" );
    }

# if defined (ENABLE_MPI)
    if ( Global::MPI_rank == 0 )
# endif
        std::fclose( outfp );

# else // if SPACE_DIMS == 1

    ScalarFlux * scalar1 = new ScalarFlux( data1, data1.DG_degree );
    ScalarFlux * scalar2 = new ScalarFlux( data2, data2.DG_degree );

    TransportOperator::Pmv( 1.0, 0.0, data1, *scalar1 );
    TransportOperator::Pmv( 1.0, 0.0, data2, *scalar2 );
    RKDG::DensityFunction::OutputDiffPlot( *scalar1, *scalar2, filename );

    delete scalar1;
    delete scalar2;

# endif // if SPACE_DIMS == ?
}


//------------------------------------------------------------------------------------------------------------
//! \brief  Reconstructs an RKDG::OrdinateFlux object using data from a binary file on disk.
//!
//! Data is read from the file pointed to by \pp{fp} starting at the initial value of \pp{fp}. Upon return,
//! \pp{fp} has been incremented to the first position past the data for the object read.
//!
//! \param[in,out]  fp                  File pointer to read binary data from.
//!
//! \see    RKDG::OrdinateFlux::WriteToDisk()
//------------------------------------------------------------------------------------------------------------
void OrdinateFlux::ReadFromDisk (

# if defined (ENABLE_MPI)
    MPI_File & fp
# else
    std::FILE * & fp
# endif
) {

# if defined (ENABLE_MPI)

//!
//! \brief  Macro to check return values from MPI commands.
//!
# define CHECK_READ \
    if ( ret != MPI_SUCCESS ) {  goto read_error;  }


    OrdinateFlux condensed {};
    char err_str [MPI_MAX_ERROR_STRING];
    int err_len, ret = MPI_SUCCESS;

    int64_t sizeof_density_read, dimof_density_read, ang_order;
    std::string ordinate_string;
    bool symmetric_reduce = false;

    Zero();

    // Read global spatial domain parameters from disk, but do not apply domain decomposition.
    condensed.DomainDecomposition::ReadFromDisk( fp, false );

    ReadStringFromFile( ordinate_string, fp );
    ret = MPI_File_read( fp, &ang_order, 1, MPI_INT64_T, MPI_STATUS_IGNORE );   CHECK_READ

    if ( ordinate_string.find( "#" ) != std::string::npos ) {

        if ( ordinate_string.back() == 'R' )
            symmetric_reduce = true;

        ordinate_string.erase( ordinate_string.size() -2, 2 );
    }

    condensed.OrdinateSet::Reconfigure( ang_order, symmetric_reduce,
                                        String_to_OrdinateType.at( ordinate_string ) );

    ret = MPI_File_read( fp, &condensed.DG_degree, 1, MPI_INT64_T, MPI_STATUS_IGNORE ); CHECK_READ
    ret = MPI_File_read( fp, &sizeof_density_read, 1, MPI_INT64_T, MPI_STATUS_IGNORE ); CHECK_READ
    ret = MPI_File_read( fp, &dimof_density_read,  1, MPI_INT64_T, MPI_STATUS_IGNORE ); CHECK_READ

    // Perform sanity checks.
    condensed.SetDensityDimensions();

    if ( condensed.sizeof_density != sizeof_density_read ) {  goto read_error;  }
    if ( condensed.dimof_density  != dimof_density_read  ) {  goto read_error;  }

    // Construct decomposed object.
    DomainDecomposition::operator=( DomainDecomposition( condensed.global_domain ) );

    this->OrdinateSet::Reconfigure( condensed.GetAngOrder(),
                                    condensed.GetOrdinateSymmetry(),
                                    condensed.GetOrdinateType() );

    this->DG_degree = condensed.DG_degree;

    // density offset contains the position in the file where this->density begins.
    MPI_Offset density_offset;
    MPI_File_get_position( fp, &density_offset );

    Allocate();

    // Read density coefficients.
# if SPACE_DIMS == 1

    {
        int64_t ii = 1;

        for ( int64_t blk_i = 0; blk_i < this->MPI_block_coords(0); ++blk_i )
            ii += this->Subdomain(blk_i).nx(0);

        MPI_Offset offset = density_offset + sizeof(double) * condensed.Index(0,ii,0);

        ret = MPI_File_seek( fp, offset, MPI_SEEK_SET );
        CHECK_READ

        ret = MPI_File_read( fp,
                            this->PointerAt(0,1,0),
                            (this->DG_degree + 1) * this->nq() * this->nx(0),
                            MPI_DOUBLE, MPI_STATUS_IGNORE );
        CHECK_READ
    }

# elif SPACE_DIMS == 2

    {
        int64_t jj = 1;

        for ( int64_t blk_j = 0; blk_j < this->MPI_block_coords(1); ++blk_j )
            jj += this->Subdomain(0,blk_j).nx(1);

        for ( int64_t i = 1; i <= this->nx(0); ++i ) {

            int64_t ii = i;

            for ( int64_t blk_i = 0; blk_i < this->MPI_block_coords(0); ++blk_i )
                ii += this->Subdomain(blk_i,0).nx(0);

            MPI_Offset offset = density_offset + sizeof(double) * condensed.Index(0,ii,jj,0,0);

            ret = MPI_File_seek( fp, offset, MPI_SEEK_SET );
            CHECK_READ

            ret = MPI_File_read( fp,
                                this->PointerAt(0,i,1,0,0),
                                (this->DG_degree + 1)*(this->DG_degree + 1) * this->nq() * this->nx(1),
                                MPI_DOUBLE, MPI_STATUS_IGNORE );
            CHECK_READ
        }
    }

# elif SPACE_DIMS == 3

    # warning "RKDG::OrdinateFlux::WriteToDisk not implemented with MPI for 3 spatial dimensions."

# if 0
    for ( int64_t i = 1; i <= this->nx(0); ++i ) {
    for ( int64_t j = 1; j <= this->nx(1); ++j ) {
    for ( int64_t k = 1; k <= this->nx(2); ++k ) {

        int64_t ii = this->MPI_block_coords(0) * this->nx(0) + i;
        int64_t jj = this->MPI_block_coords(1) * this->nx(1) + j;
        int64_t kk = this->MPI_block_coords(2) * this->nx(2) + k;

        MPI_Offset offset = density_offset + sizeof(double) * condensed.Index(0,ii,jj,kk,0,0,0);

        MPI_File_seek( fp, offset, MPI_SEEK_SET );
        CHECK_READ

        ret = MPI_File_read( fp,
                             this->PointerAt(0,i,j,k,0,0,0),
                             (this->DG_degree + 1)*(this->DG_degree + 1)*(this->DG_degree + 1) * this->nq(),
                             MPI_DOUBLE, MPI_STATUS_IGNORE );
        CHECK_READ
    }}}
# endif // if 0

# endif // if SPACE_DIMS == ?

    // Advance file pointer past end of current object (for proper return value).
    {
        MPI_Offset offset = density_offset + condensed.sizeof_density;
        ret = MPI_File_seek( fp, offset, MPI_SEEK_SET );
        CHECK_READ
    }

    this->SynchronizeHalos();

    return;


# undef CHECK_READ

read_error:

    MPI_Error_string( ret, err_str, &err_len );
    std::string error_string = std::string( err_str );
    std::replace( error_string.begin(), error_string.end(), '\n', ' ' );

    Zero();

    std::string error_message =   "Failed to read RKDG::OrdinateFlux from file pointer with error '"
                                + error_string
                                + "'.\n";

    PRINT_ERROR( error_message.c_str() )
    throw std::runtime_error( error_message );

# else // if defined (ENABLE_MPI)

//!
//! \brief  Macro to check return values from \c std::fread for errors.
//!
# define CHECK_READ \
    success &= ((bool)ret); \
    if ( !success ) {  goto read_error;  }


    Zero();

    bool success = true;
    size_t ret;

    std::string ordinate_string;
    bool symmetric_reduce = false;
    int64_t ang_order;
    int64_t sizeof_density_read, dimof_density_read;

    // First read DomainDecomposition parameters.
    DomainDecomposition::ReadFromDisk( fp );

    // Next read angular quadrature information.
    ReadStringFromFile( ordinate_string, fp );
    ret = std::fread( &ang_order, sizeof(int64_t), 1, fp ); CHECK_READ

    if ( ordinate_string.find( "#" ) != std::string::npos ) {

        if ( ordinate_string.back() == 'R' )
            symmetric_reduce = true;

        ordinate_string.erase( ordinate_string.size() -2, 2 );
    }

    OrdinateSet::Reconfigure( ang_order, symmetric_reduce, String_to_OrdinateType.at( ordinate_string ) );

    // Next read density function parameters.
    ret = std::fread( &this->DG_degree,     sizeof(int64_t), 1, fp );  CHECK_READ
    ret = std::fread( &sizeof_density_read, sizeof(int64_t), 1, fp );  CHECK_READ
    ret = std::fread( &dimof_density_read,  sizeof(int64_t), 1, fp );  CHECK_READ

    //
    // Allocate memory. Note that this computes this->sizeof_density and this->dimof_density from the
    // DomainDecomposition parameters, permitting the sanity checks performed below.
    //
    Allocate();

    // Perform sanity checks.
    if ( this->sizeof_density != sizeof_density_read ) {  goto read_error;  }
    if ( this->dimof_density  != dimof_density_read  ) {  goto read_error;  }

    // Read density coefficients.
    ret = std::fread( this->density, sizeof(double), this->dimof_density, fp ); CHECK_READ

    return;


# undef CHECK_READ

read_error:

    Zero();

    std::string error_message = "Failed to read RKDG::OrdinateFlux from file pointer.\n";
    PRINT_ERROR( error_message.c_str() )
    throw std::runtime_error( error_message );

# endif // if defined (ENABLE_MPI)
}


//------------------------------------------------------------------------------------------------------------
//! \brief  Reconstructs an RKDG::OrdinateFlux object using data from a binary file on disk.
//!
//! \param[in]  filename            Name of the file from which to read binary data from.
//!
//! \see    RKDG::OrdinateFlux::WriteToDisk()
//------------------------------------------------------------------------------------------------------------
void OrdinateFlux::ReadFromDisk (

    const std::string filename
) {

    Global::TMR_io.Start();

# if defined (ENABLE_MPI)

//!
//! \brief  Macro to check return values from MPI commands.
//!
# define CHECK_READ \
    if ( ret != MPI_SUCCESS ) {  goto read_error;  }


    MPI_File fp;
    char err_str [MPI_MAX_ERROR_STRING];
    int err_len;

    int ret = MPI_File_open( Global::MPI_cart_comm, filename.c_str(), MPI_MODE_RDONLY, MPI_INFO_NULL, &fp );
    CHECK_READ

    OrdinateFlux::ReadFromDisk( fp );

    MPI_File_close( &fp );
    Global::TMR_io.Stop();
    return;


# undef CHECK_READ

read_error:

    MPI_Error_string( ret, err_str, &err_len );
    std::string error_string = std::string( err_str );
    std::replace( error_string.begin(), error_string.end(), '\n', ' ' );

    Global::TMR_io.Stop();

    std::string error_message =   "Failed to read RKDG::OrdinateFlux from file '"
                                + filename
                                + "' with error '"
                                + error_string
                                + "'.\n";

    PRINT_ERROR( error_message.c_str() )
    throw std::runtime_error( error_message );

# else // if defined (ENABLE_MPI)

    std::FILE * fp = std::fopen( filename.c_str(), "rb" );

    if ( fp == nullptr ) {

        PRINT_ERROR( "Failed to open file '%s' in '%s'.\n", filename.c_str(), __func__ )
        goto read_error;
    }

    OrdinateFlux::ReadFromDisk( fp );

    std::fclose( fp );
    Global::TMR_io.Stop();
    return;


read_error:

    Global::TMR_io.Stop();

    std::string error_message = "Failed to read RKDG::OrdinateFlux from file '" + filename + "'.\n";

    PRINT_ERROR( error_message.c_str() )
    throw std::runtime_error( error_message );

# endif // if defined (ENABLE_MPI)
}


//------------------------------------------------------------------------------------------------------------
//! \brief  Writes the contents of an RKDG::OrdinateFlux object to a binary file on disk.
//!
//! Data is written to the file pointer \pp{fp} starting at its initial value. Upon return, \pp{fp} has been
//! incremented to the first position past the data for the object written.
//!
//! \param[in,out]  fp      File pointer to write binary data to.
//!
//! \todo   Use more efficient method for MPI 3D code.
//!
//! \see    RKDG::OrdinateFlux::ReadFromDisk()
//------------------------------------------------------------------------------------------------------------
void OrdinateFlux::WriteToDisk (

# if defined (ENABLE_MPI)
    MPI_File & fp
# else
    std::FILE * & fp
# endif

) const {

# if defined (ENABLE_MPI)

    OrdinateFlux condensed {};
    std::string ordinate_string;

    // Write condensed DomainDecomposition and OrdinateSet parameters to disk.
    if ( Global::MPI_rank == 0 ) {

        this->DomainDecomposition::WriteToDisk( fp );

        ordinate_string = OrdinateType_to_String.at( this->GetOrdinateType() );

        if ( this->GetOrdinateSymmetry() )
            ordinate_string += "#R";

        WriteStringToFile( ordinate_string, fp );
        MPI_File_write( fp, &this->GetAngOrder(), 1, MPI_INT64_T, MPI_STATUS_IGNORE );
    }

    //
    // Copy parameters into a condensed OrdinateFlux object. Memory is not allocated to the pointer
    // OrdinateFlux::density in the condensed object. The condensed object is only used to leverage its
    // indexing function.
    //
    condensed.global_domain = this->global_domain;
    condensed.DomainDecomposition::SetTopology( false /* Do not apply domain decomposition. */ );
    condensed.DomainDecomposition::ComputeSubdomains();

    condensed.DG_degree = this->DG_degree;
    condensed.OrdinateSet::Reconfigure( this->GetAngOrder(), this->GetOrdinateSymmetry(),
                                        this->GetOrdinateType() );
    condensed.SetDensityDimensions();

    // Write additional parameters.
    if ( Global::MPI_rank == 0 ) {

        int64_t sizeof_density_write = condensed.sizeof_density;

        MPI_File_write( fp, &condensed.DG_degree,     1, MPI_INT64_T, MPI_STATUS_IGNORE );
        MPI_File_write( fp, &sizeof_density_write,    1, MPI_INT64_T, MPI_STATUS_IGNORE );
        MPI_File_write( fp, &condensed.dimof_density, 1, MPI_INT64_T, MPI_STATUS_IGNORE );
    }

    // density_offset contains the offset position in the file where this->density begins.
    MPI_Offset density_offset;
    MPI_File_get_position( fp, &density_offset );
    MPI_Bcast( &density_offset, 1, MPI_OFFSET, 0, Global::MPI_cart_comm );

# if SPACE_DIMS == 1

    {
        int64_t ii = 1;

        for ( int64_t blk_i = 0; blk_i < this->MPI_block_coords(0); ++blk_i )
            ii += this->Subdomain(blk_i).nx(0);

        MPI_Offset offset = density_offset + sizeof(double) * condensed.Index(0,ii,0);

        MPI_File_seek( fp, offset, MPI_SEEK_SET );

        MPI_File_write( fp,
                        this->PointerAt(0,1,0),
                        (this->DG_degree + 1) * this->nq() * this->nx(0),
                        MPI_DOUBLE, MPI_STATUS_IGNORE );
    }

# elif SPACE_DIMS == 2

    {
        int64_t jj = 1;

        for ( int64_t blk_j = 0; blk_j < this->MPI_block_coords(1); ++blk_j )
            jj += this->Subdomain(0,blk_j).nx(1);

        for ( int64_t i = 1; i <= this->nx(0); ++i ) {

            int64_t ii = i;

            for ( int64_t blk_i = 0; blk_i < this->MPI_block_coords(0); ++blk_i )
                ii += this->Subdomain(blk_i,0).nx(0);

            MPI_Offset offset = density_offset + sizeof(double) * condensed.Index(0,ii,jj,0,0);

            MPI_File_seek( fp, offset, MPI_SEEK_SET );

            MPI_File_write( fp,
                            this->PointerAt(0,i,1,0,0),
                            (this->DG_degree + 1)*(this->DG_degree + 1) * this->nq() * this->nx(1),
                            MPI_DOUBLE, MPI_STATUS_IGNORE );
        }
    }

# elif SPACE_DIMS == 3

    # warning "RKDG::OrdinateFlux::ReadFromDisk not implemented with MPI for 3 spatial dimensions."

# if 0
    for ( int64_t i = 1; i <= this->nx(0); ++i ) {
    for ( int64_t j = 1; j <= this->nx(1); ++j ) {
    for ( int64_t k = 1; k <= this->nx(2); ++k ) {

        int64_t ii = this->MPI_block_coords(0) * this->nx(0) + i;
        int64_t jj = this->MPI_block_coords(1) * this->nx(1) + j;
        int64_t kk = this->MPI_block_coords(2) * this->nx(2) + k;

        MPI_Offset offset = density_offset + sizeof(double) * condensed.Index(0,ii,jj,kk,0,0,0);

        MPI_File_seek( fp, offset, MPI_SEEK_SET );

        MPI_File_write( fp,
                        this->PointerAt(0,i,j,k,0,0,0),
                        (this->DG_degree + 1)*(this->DG_degree + 1)*(this->DG_degree + 1) * this->nq(),
                        MPI_DOUBLE, MPI_STATUS_IGNORE );
    }}}
# endif // if 0

# endif // if SPACE_DIMS == ?

    // Advance file pointer past end of current object (for proper return value).
    MPI_Offset offset = density_offset + condensed.sizeof_density;
    MPI_File_seek( fp, offset, MPI_SEEK_SET );

# else // if defined (ENABLE_MPI)

    // First write DomainDecomposition parameters.
    DomainDecomposition::WriteToDisk( fp );

    // Next write angular quadrature information.
    std::string ordinate_string = OrdinateType_to_String.at( GetOrdinateType() );

    if ( GetOrdinateSymmetry() )
        ordinate_string += "#R";

    WriteStringToFile( ordinate_string, fp );
    std::fwrite( &GetAngOrder(), sizeof(int64_t), 1, fp );

    // Next write density function parameters.
    int64_t sizeof_density_write = this->sizeof_density;

    std::fwrite( &this->DG_degree,      sizeof(int64_t), 1, fp );
    std::fwrite( &sizeof_density_write, sizeof(int64_t), 1, fp );
    std::fwrite( &this->dimof_density,  sizeof(int64_t), 1, fp );

    std::fwrite( this->density, sizeof(double), this->dimof_density, fp );

# endif // if defined (ENABLE_MPI)
}


//------------------------------------------------------------------------------------------------------------
//! \brief  Writes the contents of an RKDG::OrdinateFlux object to a binary file on disk.
//!
//! \param[in]      filename    Name of the file to write data to.
//!
//! \see    RKDG::OrdinateFlux::ReadFromDisk()
//------------------------------------------------------------------------------------------------------------
void OrdinateFlux::WriteToDisk (

    const std::string filename

) const {

    Global::TMR_io.Start();

# if defined (ENABLE_MPI)

    MPI_File fp;

    // Open file.
    MPI_File_open( Global::MPI_cart_comm,
                   filename.c_str(),
                   MPI_MODE_CREATE | MPI_MODE_RDWR,
                   MPI_INFO_NULL,
                   &fp );

    OrdinateFlux::WriteToDisk( fp );

    MPI_File_close( &fp );

# else // if defined (ENABLE_MPI)

    std::FILE * fp = std::fopen( filename.c_str(), "wb" );

    if ( fp == nullptr ) {

        std::string error_message =   "Failed to open file '"
                                    + filename
                                    + "' in '"
                                    + std::string(__func__)
                                    + "'.\n";

        PRINT_ERROR( error_message.c_str() )
        throw std::runtime_error( error_message );
    }

    OrdinateFlux::WriteToDisk( fp );

    std::fclose( fp );

# endif // if defined (ENABLE_MPI)

    Global::TMR_io.Stop();
}
