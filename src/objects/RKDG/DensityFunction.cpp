//------------------------------------------------------------------------------------------------------------
//! \file   objects/RKDG/DensityFunction.cpp
//! \brief  Implementation of RKDG::DensityFunction class.
//!
//! \author Michael M. Crockatt
//! \date   December 2017
//------------------------------------------------------------------------------------------------------------


# include <algorithm>
# include <cinttypes>
# include <cmath>
# include <cstdint>
# include <cstdlib>
# include <cstring>
# include <limits>
# include <sstream>
# include <string>

# include "objects/RKDG/DensityFunction.hpp"
# include "utils/CLog.hpp"
# include "utils/global.hpp"
# include "utils/Quadrule/Quadrule.hpp"

# if defined (ENABLE_SIMD_BLOCKING)

    # ifndef LINF_SIMD_LEN
        # define LINF_SIMD_LEN 16
    # endif

    # include <vector>

# else // if defined (ENABLE_SIMD_BLOCKING)

    # define LINF_SIMD_LEN 1

# endif // if defined (ENABLE_SIMD_BLOCKING)

//
// PETSc redefines these as macros. This avoids compiler warnings caused by their macro expansion.
//
# if defined (MPI_Allreduce)
    # undef MPI_Allreduce
# endif
# if defined (MPI_Recv)
    # undef MPI_Recv
# endif
# if defined (MPI_Send)
    # undef MPI_Send
# endif


using namespace RKDG;
using namespace Quadrule;


//============================================================================================================
//=== CONSTRUCTORS, DESTRUCTORS, AND ASSOCIATED HELPER ROUTINES ==============================================
//============================================================================================================


//------------------------------------------------------------------------------------------------------------
//! \brief  Constructs an empty RKDG::DensityFunction object with zero-initialized parameters.
//------------------------------------------------------------------------------------------------------------
DensityFunction::DensityFunction ( void ) :

    DomainDecomposition(),
    Abstract::DensityFunction(),
    DG_degree {}
{}


//------------------------------------------------------------------------------------------------------------
//! \brief  Constructs an RKDG::DensityFunction object with given parameters.
//!
//! The elements of internal data vectors are initialized to zero.
//!
//! \param[in]      spatial_params  Contains the parameters of the spatial discretization.
//! \param[in]      DG_degree_in    Maximum degree of DG basis polynomials.
//!
//! \see    RKDG::DensityFunction::Allocate()
//------------------------------------------------------------------------------------------------------------
DensityFunction::DensityFunction (

    const DomainDecomposition & spatial_params,
    const int64_t DG_degree_in
) :
    DomainDecomposition( spatial_params ),
    Abstract::DensityFunction( spatial_params ),
    DG_degree( DG_degree_in )
{
    Allocate();
}


//------------------------------------------------------------------------------------------------------------
//! \brief  Sets appropriate values for object parameters.
//------------------------------------------------------------------------------------------------------------
void DensityFunction::SetDensityDimensions ( void ) {

    // Set value for DOF per cell.
    this->DOF_per_cell = IntPow( this->DG_degree + 1, SPACE_DIMS );

    // Call up class hierarchy.
    this->Abstract::DensityFunction::SetDensityDimensions();
}


//------------------------------------------------------------------------------------------------------------
//! \brief  Reconfigures an RKDG::DensityFunction object with a different set of parameters.
//!
//! \param[in]      DG_degree_in    Maximum degree of DG basis polynomials.
//!
//! \see    RKDG::DensityFunction::Allocate()
//! \see    RKDG::DensityFunction::Deallocate()
//------------------------------------------------------------------------------------------------------------
DensityFunction & DensityFunction::Reconfigure (

    const int64_t DG_degree_in
) {

    Deallocate();

    this->DG_degree = DG_degree_in;

    Allocate();

    return *this;
}


//------------------------------------------------------------------------------------------------------------
//! \brief  Reconfigures an RKDG::DensityFunction object with a different set of parameters.
//!
//! \param[in]      spatial_params  Contains the parameters of the spatial discretization.
//! \param[in]      DG_degree_in    Maximum degree of DG basis polynomials.
//!
//! \see    RKDG::DensityFunction::Allocate()
//! \see    RKDG::DensityFunction::Deallocate()
//------------------------------------------------------------------------------------------------------------
DensityFunction & DensityFunction::Reconfigure (

    const DomainDecomposition & spatial_params,
    const int64_t DG_degree_in
) {

    Deallocate();

    DomainDecomposition::operator=( spatial_params );
    this->DG_degree = DG_degree_in;

    Allocate();

    return *this;
}


//============================================================================================================
//=== ONE SPATIAL DIMENSION ==================================================================================
//============================================================================================================

# if SPACE_DIMS == 1


//------------------------------------------------------------------------------------------------------------
//! \brief  Evaluates the density function stored in an RKDG::DensityFunction object at the specified spatial
//!         point.
//!
//! \param[in]  x   Point to obtain value at.
//!
//! \return     Returns the value of the density at the given point.
//------------------------------------------------------------------------------------------------------------
double RKDG::DensityFunction::EvaluateAtPoint (

    const double x

) const {

    if (    x < this->global_ax(0)
         || x > this->global_bx(0)
    ) {
        std::ostringstream error_message;
        error_message << std::scientific
                      << "Evaluation point "
                      << x
                      << " is outside of RKDG::DensityFunction's domain of [ "
                      << this->global_ax(0) << ", " << this->global_bx(0)
                      << " ].\n";

        PRINT_ERROR( error_message.str().c_str() )
        throw std::range_error( error_message.str() );
    }

    double ret_val = 0.0;
    int block_coords [SPACE_DIMS];
    int local_coords [SPACE_DIMS];

    PointToCell( {x}, block_coords, local_coords );

# if defined (ENABLE_MPI)

    // If point is inside MPI block of current process, compute the value at that point.
    if ( block_coords[0] == this->MPI_block_coords(0) )
# endif
    {
        // The point 'x' lies inside the spatial cell with index 'i' of the MPI block.
        const int64_t i = local_coords[0];

        // Compute \bar{x} to evaluate Legendre polynomials with.
        const double x_bar = 2.0 * (x - this->ax(0)) / this->dx(0) - (2*i + 1);

        // Evaluate density function at given point.
        for ( int64_t d = 0; d <= this->DG_degree; ++d )
            ret_val += (*this)(i+1,d) * Legendre( d, x_bar );

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
//! \brief  Computes the difference between two RKDG::DensityFunction objects.
//!
//! This operation is lossless in the sense that the result is upscaled to the coarsest finite element space
//! which is guaranteed to represent the difference between the two density functions exactly.
//!
//! \attention  The output object \pp{result} will be reconfigured as necessary.
//!
//! \param[in]      first       First RKDG::DensityFunction object in difference.
//! \param[in]      second      Second RKDG::DensityFunction object in difference.
//! \param[out]     result      Upon return, contains the difference \pp{first} minus \pp{second}.
//------------------------------------------------------------------------------------------------------------
void DensityFunction::ComputeDifference (

    const DensityFunction & first,
    const DensityFunction & second,
    DensityFunction & result
) {

    const double TOL = 1e-14;

    if (    std::abs( first.ax(0) - second.ax(0) ) > TOL
         || std::abs( first.bx(0) - second.bx(0) ) > TOL
    ) {
        std::string error_message
            = "Domains of RKDG::DensityFunction objects in '"
              + std::string(__func__)
              + "' are not compatible for computing difference.\n";

        PRINT_ERROR( error_message.c_str() )
        throw std::invalid_argument( error_message );
    }


    // --- If solutions use the same spatial mesh, just take the difference of their coefficients. ---- //

    if ( first.nx(0) == second.nx(0) ) {

        result.Reconfigure( first, std::max( first.DG_degree, second.DG_degree ) );

        const int64_t (& nx) [SPACE_DIMS] = result.nx();

        # pragma omp parallel for collapse(2)
        for ( int64_t i = 1; i <= nx[0];           ++i ) {
        for ( int64_t d = 0; d <= first.DG_degree; ++d ) {

            result(i,d) += first(i,d);
        }}

        # pragma omp parallel for collapse(2)
        for ( int64_t i = 1; i <= nx[0];            ++i ) {
        for ( int64_t d = 0; d <= second.DG_degree; ++d ) {

            result(i,d) -= second(i,d);
        }}

        return;
    }


    // --- Otherwise, solutions require projection onto a mesh defined by least-common-multiple. ------ //

    // Configure result.
    {
        SpatialDomain result_domain( { lcm( first.nx(0), second.nx(0) ) }, first.ax(), first.bx() );

        DomainDecomposition result_decomposition( result_domain );

        result.Reconfigure( result_decomposition, std::max( first.DG_degree, second.DG_degree ) );
    }

    // Retrieve Gauss-Legendre quadrature of optimal efficiency.
    const int64_t GL_order = result.DG_degree + 1;

    double * const GL_nodes   = new double[ GL_order ];
    double * const GL_weights = new double[ GL_order ];

    ComputeQuadrature( GL_order, GL_nodes, GL_weights, NodesType::GaussLegendre );


    // Macro for indexing into matrices of spatial projection operators.
    # define IPROJ(in,out,i,d,a) ( d + ((in).DG_degree + 1)*( a + ((out).DG_degree + 1)*(i) ) )


    // Compute spatial projection operator for first solution.
    const int64_t first_num_new_cells = result.nx(0) / first.nx(0);

    const size_t dimof_first_project = (first.DG_degree + 1) * first_num_new_cells * (result.DG_degree + 1);

    double * const first_project = new double[ dimof_first_project ];
    std::memset( first_project, 0, dimof_first_project * sizeof(double) );

    for ( int64_t i = 0; i <  first_num_new_cells; ++i ) { // Index of new subcell.
    for ( int64_t a = 0; a <= result.DG_degree;    ++a ) { // Degree of new basis function.
    for ( int64_t d = 0; d <= first.DG_degree;     ++d ) { // Degree of original basis function.
    for ( int64_t l = 0; l <  GL_order;            ++l ) { // Index of Gaussian quadrature node.

        first_project[ IPROJ(first,result,i,d,a) ]
            += GL_weights[l] * (2*a + 1) / 2.0 * Legendre( a, GL_nodes[l] )
               * Legendre( d, ( GL_nodes[l] + 2*i + 1 ) / first_num_new_cells - 1.0 );
    }}}}

    // Compute spatial projection operator for second solution.
    const int64_t second_num_new_cells = result.nx(0) / second.nx(0);

    const size_t dimof_second_project = (second.DG_degree + 1) * second_num_new_cells * (result.DG_degree + 1);

    double * const second_project = new double[ dimof_second_project ];
    std::memset( second_project, 0, dimof_second_project * sizeof(double) );

    for ( int64_t i = 0; i <  second_num_new_cells; ++i ) { // Index of new subcell.
    for ( int64_t a = 0; a <= result.DG_degree;     ++a ) { // Degree of new basis function.
    for ( int64_t d = 0; d <= second.DG_degree;     ++d ) { // Degree of original basis function.
    for ( int64_t l = 0; l <  GL_order;             ++l ) { // Index of Gaussian quadrature node.

        second_project[ IPROJ(second,result,i,d,a) ]
            += GL_weights[l] * (2*a + 1) / 2.0 * Legendre( a, GL_nodes[l] )
               * Legendre( d, ( GL_nodes[l] + 2*i + 1 ) / second_num_new_cells - 1.0 );
    }}}}


    // Do projections and compute difference.
    {
    const int64_t (& nx) [SPACE_DIMS] = result.nx();

    # pragma omp parallel for collapse(2)
    for ( int64_t i = 0; i <  nx[0];            ++i ) {
    for ( int64_t a = 0; a <= result.DG_degree; ++a ) {

        // Add projection of first solution.
        for ( int64_t d = 0; d <= first.DG_degree; ++d )
            result(i+1,a) += first_project[ IPROJ(first,result, i % first_num_new_cells, d,a) ]
                             * first( (i / first_num_new_cells) + 1, d);

        // Subtract projection of second solution.
        for ( int64_t d = 0; d <= second.DG_degree; ++d )
            result(i+1,a) -= second_project[ IPROJ(second,result, i % second_num_new_cells, d,a) ]
                             * second( (i / second_num_new_cells) + 1, d);
    }}
    }

    // Clean up.
    # undef IPROJ

    delete [] GL_nodes;
    delete [] GL_weights;
    delete [] first_project;
    delete [] second_project;
}


# endif // if SPACE_DIMS == 1


//============================================================================================================
//=== TWO SPATIAL DIMENSIONS =================================================================================
//============================================================================================================

# if SPACE_DIMS == 2 || defined (DOXYCOMPILE)


//------------------------------------------------------------------------------------------------------------
//! \brief  Evaluates the density function stored in an RKDG::DensityFunction object at the specified spatial
//!         point.
//!
//! \param[in]      x           First coordinate of point to obtain value at.
//! \param[in]      y           Second coordinate of point to obtain value at.
//!
//! \return     Returns the value of the density at the given point.
//------------------------------------------------------------------------------------------------------------
double RKDG::DensityFunction::EvaluateAtPoint (

    const double x,
    const double y

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
                      << ") is outside of RKDG::DensityFunction's domain of [ "
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
        const double x_bar = 2.0 * (x - this->ax(0)) / this->dx(0) - (2*i + 1);
        const double y_bar = 2.0 * (y - this->ax(1)) / this->dx(1) - (2*j + 1);

        // Evaluate density function at given point.
        for ( int64_t d = 0; d <= this->DG_degree; ++d ) {
        for ( int64_t e = 0; e <= this->DG_degree; ++e ) {

            ret_val += (*this)(i+1,j+1,d,e) * Legendre( d, x_bar ) * Legendre( e, y_bar );
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

    MPI_Barrier( Global::MPI_cart_comm );

# endif // if defined (ENABLE_MPI)

    return ret_val;
}


//------------------------------------------------------------------------------------------------------------
//! \brief  Computes the difference between two RKDG::DensityFunction objects.
//!
//! This operation is lossless in the sense that the result is upscaled to the coarsest finite element space
//! which is guaranteed to represent the difference between the two density functions exactly.
//!
//! \attention  \pp{result} is either <tt>(\pp{data1} - \pp{data2})</tt> or <tt>(\pp{data2} - \pp{data1})</tt>
//!             depending on the structure of each solution; i.e., this routine is only useful when absolute
//!             differences are computed using one of the norm routines.
//!
//! \attention  The output object \pp{result} will be reconfigured as necessary.
//!
//! \param[in]      first       First RKDG::DensityFunction object in difference.
//! \param[in]      second      Second RKDG::DensityFunction object in difference.
//! \param[out]     result      Upon return, contains the difference between \pp{first} and \pp{second}.
//------------------------------------------------------------------------------------------------------------
void DensityFunction::ComputeDifference (

    const DensityFunction & first,
    const DensityFunction & second,
    DensityFunction & result
) {

    const double TOL = 1e-14;

    if (    std::abs( first.ax(0) - second.ax(0) ) > TOL
         || std::abs( first.ax(1) - second.ax(1) ) > TOL
         || std::abs( first.bx(0) - second.bx(0) ) > TOL
         || std::abs( first.bx(1) - second.bx(1) ) > TOL
    ) {
        std::string error_message
            = "Domains of RKDG::DensityFunction objects in "
              + std::string(__func__)
              + "' are not compatible for computing difference.\n";

        PRINT_ERROR( error_message.c_str() )
        throw std::invalid_argument( error_message );
    }


    // --- If solutions use the same spatial mesh, just take the difference of their coefficients. ---- //

    if (    first.nx(0) == second.nx(0)
         && first.nx(1) == second.nx(1)
    ) {

        result.Reconfigure( first, std::max( first.DG_degree, second.DG_degree ) );

        const int64_t (& nx) [SPACE_DIMS] = result.nx();

        # pragma omp parallel for collapse(2) schedule(static)
        for ( int64_t i = 1; i <= nx[0];           ++i ) {
        for ( int64_t j = 1; j <= nx[1];           ++j ) {
        for ( int64_t d = 0; d <= first.DG_degree; ++d ) {
        for ( int64_t e = 0; e <= first.DG_degree; ++e ) {

            result(i,j,d,e) += first(i,j,d,e);
        }}}}

        # pragma omp parallel for collapse(2) schedule(static)
        for ( int64_t i = 1; i <= nx[0];            ++i ) {
        for ( int64_t j = 1; j <= nx[1];            ++j ) {
        for ( int64_t d = 0; d <= second.DG_degree; ++d ) {
        for ( int64_t e = 0; e <= second.DG_degree; ++e ) {

            result(i,j,d,e) -= second(i,j,d,e);
        }}}}

        return;
    }


    // --- Otherwise, solutions require projection onto a mesh defined by least-common-multiple. ------ //

    // Configure result.
    {
        SpatialDomain result_domain( { lcm( first.nx(0), second.nx(0) ), lcm( first.nx(1), second.nx(1) ) },
                                     first.ax(),
                                     first.bx() );

        DomainDecomposition result_decomposition( result_domain /*, true */ );

        result.Reconfigure( result_decomposition, std::max( first.DG_degree, second.DG_degree ) );
    }

    // Retrieve Gauss-Legendre quadrature of optimal efficiency.
    const int64_t GL_order = result.DG_degree + 1;

    double * const GL_nodes   = new double[ GL_order ];
    double * const GL_weights = new double[ GL_order ];

    ComputeQuadrature( GL_order, GL_nodes, GL_weights, NodesType::GaussLegendre );


    // Macro for indexing into matrices of spatial projection operators.
    # define IPROJ(in,out,i,j,d,e,a,b) (                        \
        (e) + ((in).DG_degree + 1)*(                            \
            (d) + ((in).DG_degree + 1)*(                        \
                (b) + ((out).DG_degree + 1)*(                   \
                    (a) + ((out).DG_degree + 1)*(               \
                        (j) + ((out).nx(1))/((in).nx(1))*(i)    \
                    )                                           \
                )                                               \
            )                                                   \
        )                                                       \
    )


    // Compute spatial projection operator for first solution.
    const int64_t first_num_new_cells [] = { result.nx(0) / first.nx(0), result.nx(1) / first.nx(1) };

    const size_t dimof_first_project = (first.DG_degree + 1)*(first.DG_degree + 1)
                                       * first_num_new_cells[0] * first_num_new_cells[1]
                                       * (result.DG_degree + 1)*(result.DG_degree + 1);

    double * const first_project = new double[ dimof_first_project ];
    std::memset( first_project, 0, dimof_first_project * sizeof(double) );

    for ( int64_t i = 0; i <  first_num_new_cells[0]; ++i ) { // Index of new subcell wrt x.
    for ( int64_t j = 0; j <  first_num_new_cells[1]; ++j ) { // Index of new subcell wrt y.
    for ( int64_t a = 0; a <= result.DG_degree;       ++a ) { // Degree of new basis function wrt x.
    for ( int64_t b = 0; b <= result.DG_degree;       ++b ) { // Degree of new basis function wrt y.
    for ( int64_t d = 0; d <= first.DG_degree;        ++d ) { // Degree of original basis function wrt x.
    for ( int64_t e = 0; e <= first.DG_degree;        ++e ) { // Degree of original basis function wrt y.
    for ( int64_t l = 0; l <  GL_order;               ++l ) { // Index of Gaussian quadrature node wrt x.
    for ( int64_t m = 0; m <  GL_order;               ++m ) { // Index of Gaussian quadrature node wrt y.

        first_project[ IPROJ(first,result,i,j,d,e,a,b) ]
            += GL_weights[l] * (2*a + 1) / 2.0 * Legendre( a, GL_nodes[l] )
               * Legendre( d, ( GL_nodes[l] + 2*i + 1 ) / first_num_new_cells[0] - 1.0 )

               * GL_weights[m] * (2*b + 1) / 2.0 * Legendre( b, GL_nodes[m] )
               * Legendre( e, ( GL_nodes[m] + 2*j + 1 ) / first_num_new_cells[1] - 1.0 );
    }}}}}}}}

    // Compute spatial projection operator for second solution.
    const int64_t second_num_new_cells [] = { result.nx(0) / second.nx(0), result.nx(1) / second.nx(1) };

    const size_t dimof_second_project = (second.DG_degree + 1)*(second.DG_degree + 1)
                                        * second_num_new_cells[0] * second_num_new_cells[1]
                                        * (result.DG_degree + 1)*(result.DG_degree + 1);

    double * const second_project = new double[ dimof_second_project ];
    std::memset( second_project, 0, dimof_second_project * sizeof(double) );

    for ( int64_t i = 0; i <  second_num_new_cells[0]; ++i ) { // Index of new subcell wrt x.
    for ( int64_t j = 0; j <  second_num_new_cells[1]; ++j ) { // Index of new subcell wrt y.
    for ( int64_t a = 0; a <= result.DG_degree;        ++a ) { // Degree of new basis function wrt x.
    for ( int64_t b = 0; b <= result.DG_degree;        ++b ) { // Degree of new basis function wrt y.
    for ( int64_t d = 0; d <= second.DG_degree;        ++d ) { // Degree of original basis function wrt x.
    for ( int64_t e = 0; e <= second.DG_degree;        ++e ) { // Degree of original basis function wrt y.
    for ( int64_t l = 0; l <  GL_order;                ++l ) { // Index of Gaussian quadrature node wrt x.
    for ( int64_t m = 0; m <  GL_order;                ++m ) { // Index of Gaussian quadrature node wrt y.

        second_project[ IPROJ(second,result,i,j,d,e,a,b) ]
            += GL_weights[l] * (2*a + 1) / 2.0 * Legendre( a, GL_nodes[l] )
               * Legendre( d, ( GL_nodes[l] + 2*i + 1 ) / second_num_new_cells[0] - 1.0 )

               * GL_weights[m] * (2*b + 1) / 2.0 * Legendre( b, GL_nodes[m] )
               * Legendre( e, ( GL_nodes[m] + 2*j + 1 ) / second_num_new_cells[1] - 1.0 );
    }}}}}}}}


    // Do projections and compute difference.
    {
    const int64_t (& nx) [SPACE_DIMS] = result.nx();

    # pragma omp parallel for collapse(2)
    for ( int64_t i = 0; i <  nx[0];            ++i ) {
    for ( int64_t j = 0; j <  nx[1];            ++j ) {
    for ( int64_t a = 0; a <= result.DG_degree; ++a ) {
    for ( int64_t b = 0; b <= result.DG_degree; ++b ) {

        // Add projection of first solution.
        for ( int64_t d = 0; d <= first.DG_degree; ++d ) {
        for ( int64_t e = 0; e <= first.DG_degree; ++e ) {

            result(i+1,j+1,a,b) += first_project[ IPROJ(first,result, i % first_num_new_cells[0],
                                                                      j % first_num_new_cells[1], d,e,a,b) ]
                                   * first( (i / first_num_new_cells[0]) + 1,
                                            (j / first_num_new_cells[1]) + 1, d,e);
        }}

        // Subtract projection of second solution.
        for ( int64_t d = 0; d <= second.DG_degree; ++d ) {
        for ( int64_t e = 0; e <= second.DG_degree; ++e ) {

            result(i+1,j+1,a,b) -= second_project[ IPROJ(second,result, i % second_num_new_cells[0],
                                                                        i % second_num_new_cells[1], d,e,a,b) ]
                                   * second( (i / second_num_new_cells[0]) + 1,
                                             (j / second_num_new_cells[1]) + 1, d,e);
        }}
    }}}}
    }

    // Clean up.
    # undef IPROJ

    delete [] GL_nodes;
    delete [] GL_weights;
    delete [] second_project;
    delete [] first_project;
}


# endif // if SPACE_DIMS == 2


//============================================================================================================
//=== THREE SPATIAL DIMENSIONS ===============================================================================
//============================================================================================================

# if SPACE_DIMS == 3


//------------------------------------------------------------------------------------------------------------
//! \brief  Evaluates the density function stored in an RKDG::DensityFunction object at the specified spatial
//!         point.
//!
//! \param[in]      x           First coordinate of point to obtain value at.
//! \param[in]      y           Second coordinate of point to obtain value at.
//! \param[in]      z           Third coordinate of point to obtain value at.
//!
//! \return     Returns the value of the density at the given point.
//------------------------------------------------------------------------------------------------------------
double RKDG::DensityFunction::EvaluateAtPoint (

    const double x,
    const double y,
    const double z

) const {

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
                      << ") is outside of RKDG::DensityFunction's domain of [ "
                      << this->global_ax(0) << ", " << this->global_bx(0)
                      << " ] × [ "
                      << this->global_ax(1) << ", " << this->global_bx(1)
                      << " ] × [ "
                      << this->global_ax(2) << ", " << this->global_bx(2)
                      << " ].\n";

        PRINT_ERROR( error_message.str().c_str() )
        throw std::range_error( error_message.str() );
    }

    double ret_val = 0.0;
    int block_coords [SPACE_DIMS];
    int local_coords [SPACE_DIMS];

    PointToCell( {x,y,z}, block_coords, local_coords );

# if defined (ENABLE_MPI)

    // If point is inside MPI block of current process, compute the value at that point.
    if (    block_coords[0] == this->MPI_block_coords[0]
         && block_coords[1] == this->MPI_block_coords[1]
         && block_coords[2] == this->MPI_block_coords[2] )
# endif
    {
        // The point (x,y,z) lies inside the spatial cell with indices 'i,j,k' of the MPI block.
        const int64_t i = local_coords[0];
        const int64_t j = local_coords[1];
        const int64_t k = local_coords[2];

        // Compute \bar{x}, \bar{y}, and \bar{z} to evaluate Legendre polynomials with.
        const double x_bar = 2.0 * (x - this->ax(0)) / this->dx(0) - (2*i + 1);
        const double y_bar = 2.0 * (y - this->ax(1)) / this->dx(1) - (2*j + 1);
        const double z_bar = 2.0 * (z - this->ax(2)) / this->dx(2) - (2*k + 1);

        // Evaluate density function at given point.
        for ( int64_t d = 0; d <= this->DG_degree; ++d ) {
        for ( int64_t e = 0; e <= this->DG_degree; ++e ) {
        for ( int64_t f = 0; f <= this->DG_degree; ++f ) {

            ret_val += (*this)(i+1,j+1,k+1,d,e,f)
                      * Legendre( d, x_bar ) * Legendre( e, y_bar ) * Legendre( f, z_bar );
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
//! \brief  Creates and returns a DensityFunction object.
//!
//! \param[in]  domain_decomposition    Contains parameters of spatial discretization.
//! \param[in]  input_list              List of additional parameters.
//!
//! \return     Returns the newly created object.
//------------------------------------------------------------------------------------------------------------
DensityFunction * DensityFunction::Create (

    const DomainDecomposition & domain_decomposition,
    const ParameterList & input_list
) {

    int64_t DG_degree = GetDGDegreeX( input_list );

    return new DensityFunction( domain_decomposition, DG_degree );
}


//------------------------------------------------------------------------------------------------------------
//! \brief  Two RKDG::DensityFunction objects "are matching" if all of their parameters are equal.
//!
//! \param[in]      first       The first of the two RKDG::DensityFunction objects to compare.
//! \param[in]      second      The second of the two RKDG::DensityFunction objects to compare.
//!
//! \return     Returns true if the RKDG::DensityFunction objects "are matching" and false otherwise.
//------------------------------------------------------------------------------------------------------------
bool DensityFunction::AreMatching (

    const DensityFunction & first,
    const DensityFunction & second
) {

    return      DomainDecomposition::AreMatching( first, second )
            &&  first.DG_degree == second.DG_degree;
}


//============================================================================================================
//=== ADDITIONAL MEMBER FUNCTIONS ============================================================================
//============================================================================================================


//------------------------------------------------------------------------------------------------------------
//! \brief  Prints a summary of the object configuration to the logging interface.
//!
//! \param[in]  prefix  String prepended to each line of output, intended to be used for controlling
//!                     indentation of output.
//------------------------------------------------------------------------------------------------------------
void DensityFunction::Print (

    const std::string & prefix // = "  "

) const {

    PRINT_STATUS( "Executing RKDG::DensityFunction::%s.\n", __func__ )

    PRINT_LOG( "\n" )

    PRINT_LOG( "%s%-*s % " PRId64 "\n", prefix.c_str(),
               Global::col_width, "DG degree:", this->DG_degree )

    this->Abstract::DensityFunction::Print( prefix );

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
//! \brief  Deallocates memory at all internal pointers and zeros all parameters of an RKDG::DensityFunction
//!         object.
//!
//! \see    Abstract::DensityFunction::Zero()
//------------------------------------------------------------------------------------------------------------
DensityFunction & DensityFunction::Zero ( void ) {

    this->DG_degree = 0;
    Abstract::DensityFunction::Zero();

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
int64_t DensityFunction::CellStride (

    const int64_t dim

) const {

    switch ( dim ) {

        case 0:
            return
            # if SPACE_DIMS == 1
                Index(1,0) - Index(0,0);
            # elif SPACE_DIMS == 2
                Index(1,0,0,0) - Index(0,0,0,0);
            # elif SPACE_DIMS == 3
                Index(1,0,0,0,0,0) - Index(0,0,0,0,0,0);
            # endif

    # if SPACE_DIMS >= 2

        case 1:
            return
            # if SPACE_DIMS == 2
                Index(0,1,0,0) - Index(0,0,0,0);
            # elif SPACE_DIMS == 3
                Index(0,1,0,0,0,0) - Index(0,0,0,0,0,0);
            # endif

    # endif // if SPACE_DIMS >= 2

    # if SPACE_DIMS == 3

        case 2:
            return Index(0,0,1,0,0,0) - Index(0,0,0,0,0,0);

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
double * DensityFunction::PointerAtCell (

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
                                (i,0)
                            # elif SPACE_DIMS == 2
                                (i,j,0,0)
                            # elif SPACE_DIMS == 3
                                (i,j,k,0,0,0)
                            # endif
                            ;
}


//------------------------------------------------------------------------------------------------------------
//! \brief  Computes the total mass of the density function stored in an RKDG::DensityFunction object.
//!
//! The total mass is given by the integral of the density function over the spatial variables. Note that
//! due to the use of Legendre polynomial bases for the DG spatial approximation, the total mass is obtained
//! by adding up the degree zero coefficients and multiplying by the mesh size.
//!
//! \return     Returns the total mass in the given object.
//------------------------------------------------------------------------------------------------------------
double DensityFunction::TotalMass ( void ) const {

    PRINT_STATUS( "Computing total mass of RKDG::DensityFunction object.\n" )

    double mass = 0.0;

# if SPACE_DIMS == 1

    for ( int64_t i = 1; i <= this->nx(0); ++i )
        mass += (*this)(i,0);

# elif SPACE_DIMS == 2

    for ( int64_t i = 1; i <= this->nx(0); ++i ) {
    for ( int64_t j = 1; j <= this->nx(1); ++j ) {

        mass += (*this)(i,j,0,0);
    }}

# elif SPACE_DIMS == 3

    for ( int64_t i = 1; i <= this->nx(0); ++i ) {
    for ( int64_t j = 1; j <= this->nx(1); ++j ) {
    for ( int64_t k = 1; k <= this->nx(2); ++k ) {

        mass += (*this)(i,j,k,0,0,0);
    }}}

# endif // if SPACE_DIMS == ?

# if defined (ENABLE_MPI)
    MPI_Allreduce( MPI_IN_PLACE, &mass, 1, MPI_DOUBLE, MPI_SUM, Global::MPI_cart_comm );
# endif

    for ( int64_t i = 0; i < SPACE_DIMS; ++i ) {  mass *= this->dx(i);  }

    return mass;
}


//------------------------------------------------------------------------------------------------------------
//! \brief  Computes the \f$ L^1( d\vec{x} ) \f$ norm of the density function stored in an
//!         RKDG::DensityFunction object.
//!
//! \todo   Check this. DensityFunction::Print() produces nan when using 2x2 grid of MPI tasks.
//!
//! \return     Returns the \f$ L^1( d\vec{x} ) \f$ norm of the density function.
//------------------------------------------------------------------------------------------------------------
double DensityFunction::L1Norm ( void ) const {

    PRINT_STATUS( "Computing L1 norm of RKDG::DensityFunction object.\n" )

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

        // Value of x point for evaluating spatial quadrature rule.
        const double x_l = this->ax(0) + this->dx(0) * (i - 0.5) + this->dx(0) * GL_nodes[l] / 2.0;

        norm += GL_weights[l] * this->dx(0) * std::abs( EvaluateAtPoint( x_l ) ) / 2.0;
    }}

# elif SPACE_DIMS == 2

    # pragma omp parallel for collapse(2) schedule(static) reduction(+:norm)
    for ( int64_t i = 1; i <= nx[0];    ++i ) {
    for ( int64_t j = 1; j <= nx[1];    ++j ) {
    for ( int64_t l = 0; l <  GL_order; ++l ) {
    for ( int64_t m = 0; m <  GL_order; ++m ) {

        // Value of quadrature point scaled to current cell.
        const double x_l = this->ax(0) + this->dx(0) * (i - 0.5) + this->dx(0) * GL_nodes[l] / 2.0;
        const double y_m = this->ax(1) + this->dx(1) * (j - 0.5) + this->dx(1) * GL_nodes[m] / 2.0;

        norm += GL_weights[l] * GL_weights[m]
                * this->dx(0) * this->dx(1)
                * std::abs( EvaluateAtPoint( x_l, y_m ) )
                / 4.0;
    }}}}

# elif SPACE_DIMS == 3

    # pragma omp parallel for collapse(3) schedule(static) reduction(+:norm)
    for ( int64_t i = 1; i <= nx[0];    ++i ) {
    for ( int64_t j = 1; j <= nx[1];    ++j ) {
    for ( int64_t k = 1; k <= nx[2];    ++k ) {
    for ( int64_t l = 0; l <  GL_order; ++l ) {
    for ( int64_t m = 0; m <  GL_order; ++m ) {
    for ( int64_t n = 0; n <  GL_order; ++n ) {

        // Value of quadrature point scaled to current cell.
        const double x_l = this->ax(0) + this->dx(0) * (i - 0.5) + this->dx(0) * GL_nodes[l] / 2.0;
        const double y_m = this->ax(1) + this->dx(1) * (j - 0.5) + this->dx(1) * GL_nodes[m] / 2.0;
        const double z_n = this->ax(2) + this->dx(2) * (k - 0.5) + this->dx(2) * GL_nodes[n] / 2.0;

        norm += GL_weights[l] * GL_weights[m] * GL_weights[n]
                * this->dx(0) * this->dx(1) * this->dx(2)
                * std::abs( EvaluateAtPoint( x_l, y_m, z_n ) )
                / 8.0;
    }}}}}}

# endif // if SPACE_DIMS == ?

# if defined (ENABLE_MPI)
    MPI_Allreduce( MPI_IN_PLACE, &norm, 1, MPI_DOUBLE, MPI_SUM, Global::MPI_cart_comm );
# endif

    // Cleanup quadrature.
    delete [] GL_nodes;
    delete [] GL_weights;

    return norm;
}


//------------------------------------------------------------------------------------------------------------
//! \brief  Computes the \f$ L^2( d\vec{x} ) \f$ norm of the density function stored in an
//!         RKDG::DensityFunction object.
//!
//! The norm is computed using
//! <a href="https://en.wikipedia.org/wiki/Parseval%27s_identity">Parseval's identity</a> applied to the
//! Legendre polynomial basis.
//!
//! \return     Returns the \f$ L^2( d\vec{x} ) \f$ norm of the density function
//------------------------------------------------------------------------------------------------------------
double DensityFunction::L2Norm ( void ) const {

    PRINT_STATUS( "Computing L2 norm of RKDG::DensityFunction object.\n" )

    double norm = 0.0;

    const int64_t (& nx) [SPACE_DIMS] = this->nx();

# if SPACE_DIMS == 1

    # pragma omp parallel for schedule(static) reduction(+:norm)
    for ( int64_t i = 1; i <= nx[0];           ++i ) {
    for ( int64_t d = 0; d <= this->DG_degree; ++d ) {

        const double temp = (*this)(i,d);
        norm += temp*temp * this->dx(0) / (2*d + 1);
    }}

# elif SPACE_DIMS == 2

    # pragma omp parallel for collapse(2) schedule(static) reduction(+:norm)
    for ( int64_t i = 1; i <= nx[0];           ++i ) {
    for ( int64_t j = 1; j <= nx[1];           ++j ) {
    for ( int64_t d = 0; d <= this->DG_degree; ++d ) {
    for ( int64_t e = 0; e <= this->DG_degree; ++e ) {

        const double temp = (*this)(i,j,d,e);
        norm += temp*temp * this->dx(0) * this->dx(1) / ( (2*d + 1)*(2*e + 1) );
    }}}}

# elif SPACE_DIMS == 3

    # pragma omp parallel for collapse(3) schedule(static) reduction(+:norm)
    for ( int64_t i = 1; i <= nx[0];           ++i ) {
    for ( int64_t j = 1; j <= nx[1];           ++j ) {
    for ( int64_t k = 1; k <= nx[2];           ++k ) {
    for ( int64_t d = 0; d <= this->DG_degree; ++d ) {
    for ( int64_t e = 0; e <= this->DG_degree; ++e ) {
    for ( int64_t f = 0; f <= this->DG_degree; ++f ) {

        const double temp = (*this)(i,j,k,d,e,f);
        norm += temp*temp * this->dx(0) * this->dx(1) * this->dx(2) / ( (2*d + 1)*(2*e + 1)*(2*f + 1) );
    }}}}}}

# endif // if SPACE_DIMS == ?

# if defined (ENABLE_MPI)
    MPI_Allreduce( MPI_IN_PLACE, &norm, 1, MPI_DOUBLE, MPI_SUM, Global::MPI_cart_comm );
# endif

    return std::sqrt( norm );
}


//------------------------------------------------------------------------------------------------------------
//! \brief  Computes the \f$ L^{\infty}( d\vec{x} ) \f$ norm of the density function stored in an
//!         RKDG::DensityFunction object.
//!
//! \param[in]  tol     (optional)
//!                     Relative tolerance with with to compute \f$ L^{\infty} \f$ norm.
//!                     The default value of 1e-3 computes three digits of accuracy.
//!
//! \return     Returns the \f$ L^{\infty}( d\vec{x} ) \f$ norm of the density function.
//------------------------------------------------------------------------------------------------------------
double DensityFunction::LinfNorm (

    const double tol    // = 1e-3

) const {

    PRINT_STATUS( "Computing L-inf norm of RKDG::DensityFunction object.\n" )

    // Stores approximation of norm.
    double norm = -1.0;

# if SPACE_DIMS == 1

    # if _OPENMP >= 201107
        # pragma omp parallel for schedule(dynamic) reduction(max:norm)
    # endif
    for ( int64_t i = 1; i <= this->nx(0); ++i ) {

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

                    temp += (*this)(i,d) * pn;

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
    }

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

        const int64_t (& nx) [SPACE_DIMS] = this->nx();

        # if _OPENMP >= 201107
            # pragma omp for collapse(2) schedule(dynamic)
        # endif
        for ( int64_t i = 1; i <= nx[0]; ++i ) {
        for ( int64_t j = 1; j <= nx[1]; ++j ) {

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

                        f_norm = std::max( f_norm, std::abs( (*this)(i,j,0,0) ) );
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

                        x[l] += (*this)(i,j,d,e) * IRECURR(d,e,l);
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

                        f_norm = std::max( f_norm, std::abs( (*this)(i,j,0,0) ) );
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

                        temp += (*this)(i,j,d,e) * IRECURR(d,e);
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
        }}

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
//! \brief  Outputs a grid of sample points from the density function stored in an RKDG::DensityFunction
//!         object for plotting.
//!
//! \attention  Output is (destructively!) written to the file with name given by \pp{filename}.
//!
//! \param[in]      filename    Name of file in which to write output.
//------------------------------------------------------------------------------------------------------------
void DensityFunction::OutputPlot (

    const std::string filename

) const {

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

# if SPACE_DIMS == 1

    // Number of sample points.
    const int64_t nx = 4096;

    // Distance between sample points.
    const double dx = (this->global_bx(0) - this->global_ax(0)) / nx;

    for ( int64_t i = 0; i < nx; ++i ) {

        // Compute value of spatial point to sample density at.
        const double x_i = this->global_ax(0) + dx * (i + 0.5);

        // Evaluate density function at sample point.
        double point_val = this->EvaluateAtPoint( x_i );

    # if defined (ENABLE_MPI)
        if ( Global::MPI_rank == 0 )
    # endif
            std::fprintf( outfp, "%.4e %.4e\n", x_i, point_val );
    }

# elif SPACE_DIMS == 2

    // Number of sample points per dimension.
    const int64_t nx = 1024;
    const int64_t ny = 1024;

    // Distance between sample points.
    const double dx = (this->global_bx(0) - this->global_ax(0)) / nx;
    const double dy = (this->global_bx(1) - this->global_ax(1)) / ny;

    for ( int64_t i = 0; i < nx; ++i ) {
    for ( int64_t j = 0; j < ny; ++j ) {

        // Compute value of spatial point to sample density at.
        const double x_i = this->global_ax(0) + dx * (i + 0.5);
        const double y_j = this->global_ax(1) + dy * (j + 0.5);

        // Evaluate density function at sample point.
        double point_val = this->EvaluateAtPoint( x_i, y_j );

    # if defined (ENABLE_MPI)
        if ( Global::MPI_rank == 0 )
    # endif
            std::fprintf( outfp, "%.4e %.4e %.4e\n", x_i, y_j, point_val );
    }}

# elif SPACE_DIMS == 3

    if ( Global::MPI_rank == 0 ) {  std::fprintf( outfp, "X Y Z scalar_density\n" );  }

    // Number of sample points per dimension.
    const int64_t nx = 256;
    const int64_t ny = 256;
    const int64_t nz = 256;

    // Distance between sample points.
    const double dx = (this->global_bx(0) - this->global_ax(0)) / nx;
    const double dy = (this->global_bx(1) - this->global_ax(1)) / ny;
    const double dz = (this->global_bx(2) - this->global_ax(2)) / nz;

    for ( int64_t i = 0; i < nx; ++i ) {
    for ( int64_t j = 0; j < ny; ++j ) {
    for ( int64_t k = 0; k < nz; ++k ) {

        const double x_i = this->global_ax(0) + dx * (i + 0.5);
        const double y_j = this->global_ax(1) + dy * (j + 0.5);
        const double z_k = this->global_ax(2) + dz * (k + 0.5);

        // Evaluate density function at sample point.
        double point_val = EvaluateAtPoint( x_i, y_j, z_k );

        if ( Global::MPI_rank == 0 )
            std::fprintf( outfp, "%.4e %.4e %.4e %.4e\n", x_i, y_j, z_k, std::abs( point_val ) );
    }}}

# endif // if SPACE_DIMS == ?

# if defined (ENABLE_MPI)
    MPI_Barrier( Global::MPI_cart_comm );

    if ( Global::MPI_rank == 0 )
# endif
        std::fclose( outfp );
}


# ifndef ENABLE_MPI

//------------------------------------------------------------------------------------------------------------
//! \brief  Outputs cell averages of the density function to a binary file.
//!
//! \attention  Output is (destructively!) written to the file with name given by \pp{filename}.
//!
//! \param[in]      filename    Name of file in which to write output.
//------------------------------------------------------------------------------------------------------------
void DensityFunction::OutputCellAveragesToBin (

    const std::string filename

) const {

    const char mode [] = "w";
    std::FILE * outfp = nullptr;

    // Open output file.
    {
        outfp = std::fopen( filename.c_str(), mode );

        if ( outfp == nullptr ) {

            std::string error_message = "Failed to open file '" + filename + "' in '" + std::string(__func__)
                                        + "'.\n";

            PRINT_ERROR( error_message.c_str() )
            throw std::runtime_error( error_message );
        }
    }

# if SPACE_DIMS == 1

    for ( int64_t i = 1; i <= this->nx(0); ++i )
        std::fwrite( this->PointerAt(i,0), sizeof(double), 1, outfp );

# elif SPACE_DIMS == 2

    for ( int64_t i = 1; i <= this->nx(0); ++i ) {
    for ( int64_t j = 1; j <= this->nx(1); ++j ) {

        std::fwrite( this->PointerAt(i,j,0,0), sizeof(double), 1, outfp );
    }}

# elif SPACE_DIMS == 3

    for ( int64_t i = 1; i <= this->nx(0); ++i ) {
    for ( int64_t j = 1; j <= this->nx(1); ++j ) {
    for ( int64_t k = 1; k <= this->nx(2); ++k ) {

        std::fwrite( this->PointerAt(i,j,k,0,0,0), sizeof(double), 1, outfp );
    }}}

# endif // if SPACE_DIMS == ?

    std::fclose( outfp );
}

# endif // ifndef ENABLE_MPI


//------------------------------------------------------------------------------------------------------------
//! \brief  Outputs a grid of points sampling the difference between two density functions for plotting.
//!
//! \attention  Output is (destructively!) written to the file with name given by \pp{filename}.
//!
//! \todo   Update 1D code for MPI.
//!
//! \param[in]      first       First RKDG::DensityFunction object in difference.
//! \param[in]      second      Second RKDG::DensityFunction object in difference.
//! \param[in]      filename    Name of file in which to write output.
//------------------------------------------------------------------------------------------------------------
void DensityFunction::OutputDiffPlot (

    const DensityFunction & first,
    const DensityFunction & second,
    const std::string filename
) {

    const char mode[] = "w";
    FILE * outfp = nullptr;

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

# if SPACE_DIMS == 1

    // Number of sample points.
    const int64_t nx = 4096;

    // Distance between sample points.
    const double dx = (first.global_bx(0) - first.global_ax(0)) / nx;

    for ( int64_t i = 0; i < nx; ++i ) {

        // Compute value of spatial point to sample density at.
        const double x_i = first.global_ax(0) + dx * (i + 0.5);

        // Evaluate density function at sample spatial point.
        std::fprintf( outfp, "%.4e %.4e\n", x_i,
                      first.EvaluateAtPoint( x_i ) - second.EvaluateAtPoint( x_i ) );
    }

# elif SPACE_DIMS == 2

    // Number of sample points per dimension.
    const int64_t nx = 1024;
    const int64_t ny = 1024;

    // Distance between sample points.
    const double dx = (first.global_bx(0) - first.global_ax(0)) / nx;
    const double dy = (first.global_bx(1) - first.global_ax(1)) / ny;

    for ( int64_t i = 0; i < nx; ++i ) {
    for ( int64_t j = 0; j < ny; ++j ) {

        // Compute value of spatial point to sample density at.
        const double x_i = first.global_ax(0) + dx * (i + 0.5);
        const double y_j = first.global_ax(1) + dy * (j + 0.5);

        // Evaluate density function at sample point.
        double diff_val = first.EvaluateAtPoint( x_i, y_j ) - second.EvaluateAtPoint( x_i, y_j );

    # if defined (ENABLE_MPI)
        if ( Global::MPI_rank == 0 )
    # endif
            std::fprintf( outfp, "%.4e %.4e %.4e\n", x_i, y_j, diff_val );
    }}

# elif SPACE_DIMS == 3

    if ( Global::MPI_rank == 0 ) {  std::fprintf( outfp, "X Y Z scalar_density\n" );  }

    // Number of sample points per dimension.
    const int64_t nx = 256;
    const int64_t ny = 256;
    const int64_t nz = 256;

    // Distance between sample points.
    const double dx = (first.global_bx(0) - first.global_ax(0)) / nx;
    const double dy = (first.global_bx(1) - first.global_ax(1)) / ny;
    const double dz = (first.global_bx(2) - first.global_ax(2)) / nz;

    for ( int64_t i = 0; i < nx; ++i ) {
    for ( int64_t j = 0; j < ny; ++j ) {
    for ( int64_t k = 0; k < nz; ++k ) {

        const double x_i = first.global_ax(0) + dx * (i + 0.5);
        const double y_j = first.global_ax(1) + dy * (j + 0.5);
        const double z_k = first.global_ax(2) + dz * (k + 0.5);

        // Evaluate density function at sample point.
        double diff_val = first.EvaluateAtPoint( x_i, y_j, z_k ) - second.EvaluateAtPoint( x_i, y_j, z_k );

        if ( Global::MPI_rank == 0 )
            std::fprintf( outfp, "%.4e %.4e %.4e %.4e\n", x_i, y_j, z_k, std::abs( diff_val ) );
    }}}

# endif // if SPACE_DIMS == ?

# if defined (ENABLE_MPI)
    MPI_Barrier( MPI_COMM_WORLD );

    if ( Global::MPI_rank == 0 )
# endif
        std::fclose( outfp );
}


//------------------------------------------------------------------------------------------------------------
//! \brief  Reconstructs an RKDG::DensityFunction object using data from a binary file on disk.
//!
//! Data is read from the file pointed to by \pp{fp} starting at the initial value of \pp{fp}. Upon return,
//! \pp{fp} has been incremented to the first position past the data for the object read.
//!
//! \param[in,out]  fp                  File pointer to read binary data from.
//! \param[in]      version             Integer for file format to use. This option in here only to enable
//!                                     conversion of old output files to the current format. Generally the
//!                                     default value should be used.
//!
//! \see    RKDG::DensityFunction::WriteToDisk()
//------------------------------------------------------------------------------------------------------------
void DensityFunction::ReadFromDisk (

# if defined (ENABLE_MPI)
    MPI_File & fp,
# else
    std::FILE * & fp,
# endif
    const int64_t version // = 0
) {

# if defined (ENABLE_MPI)

//!
//! \brief  Macro to check return values from MPI commands.
//!
# define CHECK_READ \
    if ( ret != MPI_SUCCESS ) {  goto read_error;  }


    DensityFunction condensed {};
    char err_str [MPI_MAX_ERROR_STRING];
    int err_len, ret = MPI_SUCCESS;

    int64_t sizeof_density_read, dimof_density_read;

    Zero();

    // Read global spatial domain parameters from disk, but do not apply domain decomposition.
    condensed.DomainDecomposition::ReadFromDisk( fp, false, version );

    MPI_File_read( fp, &condensed.DG_degree, 1, MPI_INT64_T, MPI_STATUS_IGNORE );  CHECK_READ
    MPI_File_read( fp, &sizeof_density_read, 1, MPI_INT64_T, MPI_STATUS_IGNORE );  CHECK_READ
    MPI_File_read( fp, &dimof_density_read,  1, MPI_INT64_T, MPI_STATUS_IGNORE );  CHECK_READ

    // Perform sanity checks.
    condensed.SetDensityDimensions();

    if ( condensed.sizeof_density != sizeof_density_read ) {  goto read_error;  }
    if ( condensed.dimof_density  != dimof_density_read  ) {  goto read_error;  }

    // Construct decomposed object.
    DomainDecomposition::operator=( DomainDecomposition( condensed.global_domain ) );
    this->DG_degree = condensed.DG_degree;

    // density_offset contains the position in the file where this->density begins.
    MPI_Offset density_offset;
    MPI_File_get_position( fp, &density_offset );

    Allocate();

    // Read density coefficients.
# if SPACE_DIMS == 1

    {
        int64_t ii = 1;

        for ( int64_t blk_i = 0; blk_i < this->MPI_block_coords(0); ++blk_i )
            ii += this->Subdomain(blk_i).nx(0);

        MPI_Offset offset = density_offset + sizeof(double) * condensed.Index(ii,0);

        ret = MPI_File_seek( fp, offset, MPI_SEEK_SET );
        CHECK_READ

        ret = MPI_File_read( fp,
                            this->PointerAt(1,0),
                            (this->DG_degree + 1) * this->nx(0),
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

            MPI_Offset offset = density_offset + sizeof(double) * condensed.Index(ii,jj,0,0);

            ret = MPI_File_seek( fp, offset, MPI_SEEK_SET );
            CHECK_READ

            ret = MPI_File_read( fp,
                                this->PointerAt(i,1,0,0),
                                (this->DG_degree + 1)*(this->DG_degree + 1) * this->nx(1),
                                MPI_DOUBLE, MPI_STATUS_IGNORE );
            CHECK_READ
        }
    }

# elif SPACE_DIMS == 3

    # warning "RKDG::DensityFunction::WriteToDisk not implemented with MPI for 3 spatial dimensions."

# if 0
    for ( int64_t i = 1; i <= this->nx(0); ++i ) {
    for ( int64_t j = 1; j <= this->nx(1); ++j ) {
    for ( int64_t k = 1; k <= this->nx(2); ++k ) {

        int64_t ii = this->MPI_block_coords(0) * this->nx(0) + i;
        int64_t jj = this->MPI_block_coords(1) * this->nx(1) + j;
        int64_t kk = this->MPI_block_coords(2) * this->nx(2) + k;

        MPI_Offset offset = density_offset + sizeof(double) * condensed.Index(ii,jj,kk,0,0,0);

        MPI_File_seek( fp, offset, MPI_SEEK_SET );
        CHECK_READ

        ret = MPI_File_read( fp,
                             this->PointerAt(i,j,k,0,0,0),
                             (this->DG_degree + 1)*(this->DG_degree + 1)*(this->DG_degree + 1),
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

    if (    Global::MPI_num_ranks > 1
         || Global::periodic
    )
        this->halo_cells_dirty = true;
    else
        this->halo_cells_dirty = false;

//     this->SynchronizeHalos();

    return;


# undef CHECK_READ

read_error:

    MPI_Error_string( ret, err_str, &err_len );
    std::string error_string = std::string( err_str );
    std::replace( error_string.begin(), error_string.end(), '\n', ' ' );

    Zero();

    std::string error_message =   "Failed to read RKDG::DensityFunction from file pointer with error '"
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

    // First read DomainDecomposition parameters.
    DomainDecomposition::ReadFromDisk( fp, version );

    // Next read density function parameters.
    int64_t sizeof_density_read, dimof_density_read;

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

    if ( Global::periodic )
        this->halo_cells_dirty = true;
    else
        this->halo_cells_dirty = false;

//     this->SynchronizeHalos();

    return;


# undef CHECK_READ

read_error:

    Zero();

    std::string error_message = "Failed to read RKDG::DensityFunction from file pointer.\n";
    PRINT_ERROR( error_message.c_str() )
    throw std::runtime_error( error_message );

# endif // if defined (ENABLE_MPI)
}


//------------------------------------------------------------------------------------------------------------
//! \brief  Reconstructs an RKDG::DensityFunction object using data from a binary file on disk.
//!
//! \param[in]  filename            Name of the file from which to read binary data from.
//! \param[in]  version             Integer for file format to use. This option in here only to enable
//!                                 conversion of old output files to the current format. Generally the
//!                                 default value should be used.
//!
//! \see    RKDG::DensityFunction::WriteToDisk()
//------------------------------------------------------------------------------------------------------------
void DensityFunction::ReadFromDisk (

    const std::string filename,
    const int64_t version // = 0
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

    DensityFunction::ReadFromDisk( fp, version );

    MPI_File_close( &fp );
    Global::TMR_io.Stop();
    return;


# undef CHECK_READ

read_error:

    MPI_Error_string( ret, err_str, &err_len );
    std::string error_string = std::string( err_str );
    std::replace( error_string.begin(), error_string.end(), '\n', ' ' );

    Global::TMR_io.Stop();

    std::string error_message =   "Failed to read RKDG::DensityFunction from file '"
                                + filename
                                + "' with error '"
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


    std::FILE * fp = std::fopen( filename.c_str(), "rb" );

    if ( fp == nullptr ) {

        PRINT_ERROR( "Failed to open file '%s' in '%s'.\n", filename.c_str(), __func__ )
        goto read_error;
    }

    DensityFunction::ReadFromDisk( fp, version );

    std::fclose( fp );
    Global::TMR_io.Stop();
    return;


# undef CHECK_READ

read_error:

    Global::TMR_io.Stop();

    std::string error_message = "Failed to read RKDG::DensityFunction from file '" + filename + "'.\n";

    PRINT_ERROR( error_message.c_str() )
    throw std::runtime_error( error_message );

# endif // if defined (ENABLE_MPI)
}


//------------------------------------------------------------------------------------------------------------
//! \brief  Writes the contents of an RKDG::DensityFunction object to a binary file on disk.
//!
//! Data is written to the file pointer \pp{fp} starting at its initial value. Upon return, \pp{fp} has been
//! incremented to the first position past the data for the object written.
//!
//! \param[in,out]  fp      File pointer to write binary data to.
//!
//! \todo   Use more efficient method for MPI 3D code.
//!
//! \see    RKDG::DensityFunction::ReadFromDisk()
//------------------------------------------------------------------------------------------------------------
void DensityFunction::WriteToDisk (

# if defined (ENABLE_MPI)
    MPI_File & fp
# else
    std::FILE * & fp
# endif

) const {

# if defined (ENABLE_MPI)

    DensityFunction condensed {};

    // Write condensed DomainDecomposition parameters to disk.
    if ( Global::MPI_rank == 0 )
        DomainDecomposition::WriteToDisk( fp );

    //
    // Copy parameters into a condensed DensityFunction object. Memory is not allocated to the pointer
    // DensityFunction::density in the condensed object. The condensed object is only used to leverage its
    // indexing function.
    //
    condensed.global_domain = this->global_domain;
    condensed.DomainDecomposition::SetTopology( false /* Do not apply domain decomposition. */ );
    condensed.DomainDecomposition::ComputeSubdomains();

    condensed.DG_degree = this->DG_degree;
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

        MPI_Offset offset = density_offset + sizeof(double) * condensed.Index(ii,0);

        MPI_File_seek( fp, offset, MPI_SEEK_SET );

        MPI_File_write( fp,
                        this->PointerAt(1,0),
                        (this->DG_degree + 1) * this->nx(0),
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

            MPI_Offset offset = density_offset + sizeof(double) * condensed.Index(ii,jj,0,0);

            MPI_File_seek( fp, offset, MPI_SEEK_SET );

            MPI_File_write( fp,
                            this->PointerAt(i,1,0,0),
                            (this->DG_degree + 1)*(this->DG_degree + 1) * this->nx(1),
                            MPI_DOUBLE, MPI_STATUS_IGNORE );
        }
    }

# elif SPACE_DIMS == 3

    # warning "RKDG::DensityFunction::ReadFromDisk not implemented with MPI for 3 spatial dimensions."

# if 0
    for ( int64_t i = 1; i <= this->nx(0); ++i ) {
    for ( int64_t j = 1; j <= this->nx(1); ++j ) {
    for ( int64_t k = 1; k <= this->nx(2); ++k ) {

        int64_t ii = this->MPI_block_coords(0) * this->nx(0) + i;
        int64_t jj = this->MPI_block_coords(1) * this->nx(1) + j;
        int64_t kk = this->MPI_block_coords(2) * this->nx(2) + k;

        MPI_Offset offset = density_offset + sizeof(double) * condensed.Index(ii,jj,kk,0,0,0);

        MPI_File_seek( fp, offset, MPI_SEEK_SET );

        MPI_File_write( fp,
                        this->PointerAt(i,j,k,0,0,0),
                        (this->DG_degree + 1)*(this->DG_degree + 1)*(this->DG_degree + 1),
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

    // Next write density function parameters.
    int64_t sizeof_density_write = this->sizeof_density;

    std::fwrite( &this->DG_degree,      sizeof(int64_t), 1, fp );
    std::fwrite( &sizeof_density_write, sizeof(int64_t), 1, fp );
    std::fwrite( &this->dimof_density,  sizeof(int64_t), 1, fp );

    std::fwrite( this->density, sizeof(double), this->dimof_density, fp );

# endif // if defined (ENABLE_MPI)
}


//------------------------------------------------------------------------------------------------------------
//! \brief  Writes the contents of an RKDG::DensityFunction object to a binary file on disk.
//!
//! \param[in]      filename    Name of the file to write data to.
//!
//! \see    RKDG::DensityFunction::ReadFromDisk()
//------------------------------------------------------------------------------------------------------------
void DensityFunction::WriteToDisk (

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

    DensityFunction::WriteToDisk( fp );

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

    DensityFunction::WriteToDisk( fp );

    std::fclose( fp );

# endif // if defined (ENABLE_MPI)

    Global::TMR_io.Stop();
}
