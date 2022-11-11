//------------------------------------------------------------------------------------------------------------
//! \file   main_deckmaker.cpp
//! \brief  Creates input files for ranges of parameters.
//!
//! \author Michael M. Crockatt
//! \date   February 2018
//------------------------------------------------------------------------------------------------------------


# include <cinttypes>
# include <cmath>
# include <cstdio>
# include <cstdlib>
# include <cstring>
# include <iomanip>
# include <iostream>
# include <fstream>
# include <map>
# include <set>
# include <sstream>
# include <vector>

# define DECKMAKER

# include "time_integrators/Abstract/TimeIntegrator.hpp"
# include "time_integrators/RKDG/dirk.hpp"
# include "time_integrators/RKDG/idc.hpp"
# include "time_integrators/STDG/stdg.hpp"

# include "linear_solvers/Abstract/ImplicitSolver/AztecGMRESSolver.hpp"
# include "linear_solvers/Abstract/ImplicitSolver/AztecPeierlsSolver.hpp"
# include "linear_solvers/Abstract/ImplicitSolver/BelosPeierlsSolver.hpp"
# include "linear_solvers/Abstract/ImplicitSolver/PETScPeierlsSolver.hpp"
# include "linear_solvers/Abstract/ImplicitSolver/SISolver.hpp"
# include "linear_solvers/Abstract/ImplicitSolver/SweepSolver.hpp"

# include "linear_solvers/Abstract/SweepCellSolve/SweepGESolve.hpp"
# include "linear_solvers/Abstract/SweepCellSolve/SweepHessSolve.hpp"

# include "linear_solvers/Abstract/SweepPattern/SweepPatternUniCartKBA.hpp"
# include "linear_solvers/Abstract/SweepPattern/SweepPatternUniCartKBA2.hpp"
# include "linear_solvers/Abstract/SweepPattern/SweepPatternUniCartSequential.hpp"
# include "linear_solvers/Abstract/SweepPattern/SweepPatternUniCartWavefront.hpp"

# include "utils/CLog.hpp"
# include "utils/init.hpp"
# include "utils/ParameterList.hpp"
# include "utils/job_scripts.hpp"
# include "utils/Quadrule/OrdinateSet.hpp"


using namespace Abstract;
using namespace Quadrule;


# if SPACE_DIMS == 1
    const std::string spacedims_str = "1";
# elif SPACE_DIMS == 2
    const std::string spacedims_str = "2";
# elif SPACE_DIMS == 3
    const std::string spacedims_str = "3";
# else
    # error "Invalid number of spatial dimensions specified."
# endif


//============================================================================================================
//=== HELPER FUNCTIONS =======================================================================================
//============================================================================================================


//------------------------------------------------------------------------------------------------------------
//! \brief  Constructs a map from integers to hybrid method strings for an arbitrary integrator. Nonhybrid
//!         methods are excluded.
//!
//! \param[out] hybrid_methods          Map in which to store results.
//! \param[in]  HybridMethod_to_String  Map from *::HybridType to std::string values for a time integrator
//!                                     class.
//------------------------------------------------------------------------------------------------------------
template< typename HybridType >
void InsertHybridStrings (

    std::map< int64_t, std::string > & hybrid_methods,
    const std::map< HybridType, std::string > HybridMethod_to_String
) {

    for ( auto hybrid_method : HybridMethod_to_String )
        hybrid_methods.insert( std::make_pair( (int64_t) hybrid_method.first, hybrid_method.second ) );
}


//------------------------------------------------------------------------------------------------------------
//! \brief  Constructs a map from integers to hybrid method strings for an arbitrary integrator. Nonhybrid
//!         methods are excluded.
//!
//! \param[in]  time_integrator     Type of time integrator to generate list of hybrid methods for
//!
//! \return Returns a map from integers to string descriptors for hybrid methods. Integers are assigned in
//!         the order determined by the appropriate HybridMethod_to_String object, starting at 1.
//------------------------------------------------------------------------------------------------------------
std::map< int64_t, std::string >
ConstructHybridStringsMap (

    const TimeIntegratorType time_integrator
) {

    std::map< int64_t, std::string > hybrid_methods;

    switch ( time_integrator )
    {
        case TimeIntegratorType::STDG:
            InsertHybridStrings( hybrid_methods, STDGIntegrator::HybridMethod_to_String );
            break;

        case TimeIntegratorType::DIRK:
            InsertHybridStrings( hybrid_methods, DIRKIntegrator::HybridMethod_to_String );
            break;

        case TimeIntegratorType::IDC:
        case TimeIntegratorType::LSIDC:
            InsertHybridStrings( hybrid_methods, IDCIntegrator::HybridMethod_to_String );
            break;

        default:
        {   std::string error_message = "Invalid time integrator '"
                                      + TimeIntegratorType_to_String.at( time_integrator )
                                      + "' in '"
                                      + std::string(__func__)
                                      + "'.\n";

            PRINT_ERROR( error_message.c_str() )
            throw std::invalid_argument( error_message );
        }
    }

    return hybrid_methods;
}


//------------------------------------------------------------------------------------------------------------
//! \brief  Factory method for job script template strings.
//------------------------------------------------------------------------------------------------------------
std::string GetJobScriptTemplate (

    const std::string template_name
) {

    if ( template_name == "omp_parallel_restart" )
    {
        return omp_parallel_restart;
    }
    if ( template_name == "mpi_parallel_restart" )
    {
        return mpi_parallel_restart;
    }
    else
    {
        std::string error_message = "Unknown job script template '" + template_name + "' in '"
                                  + std::string(__func__) + "'.\n";

        PRINT_ERROR( error_message.c_str() )
        throw std::invalid_argument( error_message );
    }
}


//------------------------------------------------------------------------------------------------------------
//! \brief  Generate the string for setting \pp{OMP_NUM_THREADS}.
//!
//! \param[in]  num_threads     The total number of threads.
//!
//! \return     Returns the appropriate \pp{OMP_NUM_THREADS} string for the given number of threads.
//------------------------------------------------------------------------------------------------------------
std::string GenerateOMPNumThreads (

    const int64_t num_threads
) {

    if ( num_threads <= 0 )
    {
        std::string error_message = "Number of threads passed to " + std::string(__func__) 
                                  + " must be positive (" + std::to_string( num_threads ) + " was given).\n";

        PRINT_ERROR( error_message.c_str() )
        throw std::invalid_argument( error_message );
    }

    if ( num_threads == 1 )
        return "1,1";
    else
        return std::to_string( num_threads ) + "," + std::to_string( num_threads - 1 );
}


//------------------------------------------------------------------------------------------------------------
//! \brief  Generate the string for setting \pp{OMP_PLACES}
//!
//! \param[in]  num_threads     The total number of threads.
//!
//! \return     Returns the appropriate \pp{OMP_PLACES} string for the given number of threads.
//------------------------------------------------------------------------------------------------------------
std::string GenerateOMPPlaces (

    const int64_t num_threads
) {

    if ( num_threads <= 0 )
    {
        std::string error_message = "Number of threads passed to " + std::string(__func__) 
                                  + " must be positive (" + std::to_string( num_threads ) + " was given).\n";

        PRINT_ERROR( error_message.c_str() )
        throw std::invalid_argument( error_message );
    }

    // Start with 1 thread.
    std::string omp_places = "{1}";

    for ( int64_t k = 1; k < num_threads; ++k )
        omp_places += ",{" + std::to_string( k+1 ) + "}";

    return omp_places;
}


//------------------------------------------------------------------------------------------------------------
//! \brief  Replaces all occurences of \pp{placeholder} with \pp{replacement} in \pp{target}.
//!
//! \param[in,out]  target          String to perform replacements in.
//! \param[in]      placeholder     Placeholder string to replace.
//! \param[in]      replacement     String to replace placeholder with.
//------------------------------------------------------------------------------------------------------------
void ReplaceAllPlaceholder (

    std::string & target,
    const std::string & placeholder,
    const std::string & replacement
) {

    std::string::size_type substr_pos = 0;
    while ( ( substr_pos = target.find( placeholder, substr_pos ) ) != std::string::npos )
    {
        target.replace( substr_pos, placeholder.size(), replacement );
        substr_pos += placeholder.size();
    }
}


//------------------------------------------------------------------------------------------------------------
//! \brief  Writes a job script to a file with name \pp{filename} using the template \pp{job_template} and
//!         replacing each occurrence of #jobname_placeholder in \pp{job_template} with the given
//!         \pp{job_prefix}. Other substitutions (job account, thread/process counts, etc.) are also made into
//!         the generated job script.
//!
//! \param[in]  job_template    String containing the template for the job script.
//! \param[in]  job_prefix      Filename prefix for the job.
//! \param[in]  filename        Filename to output jobscript to.
//! \param[in]  procs_per_node  Number of MPI ranks or OpenMP threads per node, depending on the job script
//!                             configuration.
//! \param[in]  num_nodes       Number of nodes to use for the run.
//! \param[in]  job_account     Account string to substitute into the job script, if needed.
//------------------------------------------------------------------------------------------------------------
void WriteJobScript (

    const std::string & job_template,
    const std::string & job_prefix,
    const std::string & filename,
    const int64_t procs_per_node,
    const int64_t num_nodes = 1,
    const std::string job_account = "NOACCOUNT"
) {

    // Open output file for writing.
    std::ofstream out_stream( filename, std::ios::out );

    if ( !out_stream.is_open() )
    {
        std::string error_message = "Failed to open file '" + filename + "' in '" + std::string(__func__)
                                    + "'.\n";

        PRINT_ERROR( error_message.c_str() )
        throw std::runtime_error( error_message );
    }

    // Make copy of job template and perform replacements for jobname and number of spatial dimensions.
    std::string job_script = job_template;

    ReplaceAllPlaceholder( job_script, jobname_placeholder, job_prefix );
    ReplaceAllPlaceholder( job_script, spacedims_placeholder, spacedims_str );

    // Replace node count based on num_nodes.
    const std::string nodes_in = "nodes=1";
    const std::string nodes_out = "nodes=" + std::to_string( num_nodes );

    ReplaceAllPlaceholder( job_script, nodes_in, nodes_out );

    // Handle account used in jobscript, if specified.
    if ( job_account == "NOACCOUNT" )
        ReplaceAllPlaceholder( job_script, accountname_placeholder, job_account );
    else
        ReplaceAllPlaceholder( job_script, accountname_placeholder, "SBATCH --account=" + job_account );

    // Substitute for OMP_NUM_THREADS.
    ReplaceAllPlaceholder( job_script, ompnumthreads_placeholder, GenerateOMPNumThreads( procs_per_node ) );

    // Substitute for OMP_PLACES.
    ReplaceAllPlaceholder( job_script, ompplaces_placeholder, GenerateOMPPlaces( procs_per_node ) );

    // Substitute for --npernode.
    ReplaceAllPlaceholder( job_script, npernode_placeholder, std::to_string( procs_per_node ) );

    // Write resulting job script to file and close it.
    out_stream << job_script;

    out_stream.close();
}


//------------------------------------------------------------------------------------------------------------
//! \brief  Splits a list of comma-separated values into a vector of individual strings.
//!
//! \param[in]  csv_string  CSV string to split.
//!
//! \return Returns a vector containing the components of the split string.
//------------------------------------------------------------------------------------------------------------
std::vector<std::string>
SplitCSVString (

    const std::string & csv_string
) {

    std::vector<std::string> str_vec;

    size_t first = 0;
    size_t second = 0;

    while ( second < csv_string.size() )
    {
        while (     ( second < csv_string.size() )
                &&  ( csv_string.at(second) != ',' )
        ) {
            ++second;
        }

        str_vec.push_back( csv_string.substr( first, second - first ) );
        first = ++second;
    }

    return str_vec;
}


//------------------------------------------------------------------------------------------------------------
//! \brief  Formats the angle string portion of the file naming pattern. For the nonhybrid case.
//!
//! \param[in]  ang_order   Integer specifying the angular resolution.
//!
//! \return     Returns the formatted angle string.
//------------------------------------------------------------------------------------------------------------
std::string
FormatAngleString (

    const int64_t ang_order
) {

    std::ostringstream angle_str;

    angle_str << "q"
              << std::setw(3) << std::setfill('0') << std::right << ang_order;

    return angle_str.str();
}


//------------------------------------------------------------------------------------------------------------
//! \brief  Formats the angle string portion of the file naming pattern. For the hybrid case.
//!
//! \param[in]  u_ang_order   Integer specifying the uncollided angular resolution.
//! \param[in]  c_ang_order   Integer specifying the collided angular resolution.
//!
//! \return     Returns the formatted angle string.
//------------------------------------------------------------------------------------------------------------
std::string
FormatAngleString (

    const int64_t u_ang_order,
    const int64_t c_ang_order
) {

    std::ostringstream angle_str;

    angle_str << "q"
              << std::setw(2) << std::setfill('0') << std::right << u_ang_order
              << "-"
              << std::setw(2) << std::setfill('0') << std::right << c_ang_order;

    return angle_str.str();
}


//------------------------------------------------------------------------------------------------------------
//! \brief  Formats the mesh string portion of the file naming pattern.
//!
//! \param[in]  nx_int  Integer specifying the mesh resolution.
//!
//! \return     Returns the formatted mesh string.
//------------------------------------------------------------------------------------------------------------
std::string
FormatMeshString (

    const int64_t nx_int
) {

    std::ostringstream mesh_str;

    mesh_str << "nx";
    mesh_str << std::setw(6) << std::setfill('0') << std::right << nx_int;

    return mesh_str.str();
}


//------------------------------------------------------------------------------------------------------------
//! \brief  Formats the mesh string portion of the file naming pattern.
//!
//! \param[in]  nx_str  Integer string specifying the mesh resolution.
//!
//! \return     Returns the formatted mesh string.
//------------------------------------------------------------------------------------------------------------
std::string
FormatMeshString (

    const std::string & nx_str
) {

    // First convert to integer value.
    int64_t nx_int = 0;
    if ( std::stringstream( nx_str ) >> nx_int ) {/* empty */}
    else {
        std::string error_message = "Error converting value '" + nx_str + "' to type int.\n";
        PRINT_ERROR( error_message.c_str() )
        throw std::runtime_error( error_message );
    }

    return FormatMeshString( nx_int );
}


//------------------------------------------------------------------------------------------------------------
//! \brief  Formats the order string portion of the file naming pattern.
//!
//! \param[in]  order_int   Integer specifying the method order.
//!
//! \return     Returns the formatted order string.
//------------------------------------------------------------------------------------------------------------
std::string
FormatOrderString (

    const int64_t order_int
) {

    return "ord" + std::to_string( order_int );
}


//------------------------------------------------------------------------------------------------------------
//! \brief  Formats the order string portion of the file naming pattern.
//!
//! \param[in]  order_str   Integer string specifying the method order.
//!
//! \return     Returns the formatted order string.
//------------------------------------------------------------------------------------------------------------
std::string
FormatOrderString (

    const std::string & order_str
) {

    // First convert to integer value.
    int64_t order_int = 0;
    if ( std::stringstream( order_str ) >> order_int ) {/* empty */}
    else {
        std::string error_message = "Error converting value '" + order_str + "' to type int.\n";
        PRINT_ERROR( error_message.c_str() )
        throw std::runtime_error( error_message );
    }

    return FormatOrderString( order_int );
}


//============================================================================================================
//=== HELPERS FOR SPACE-TIME REFINEMENT ======================================================================
//============================================================================================================


//! \brief  Specifies which strategy to use for refinement in space and time.
enum class RefinementType {

    //! No refinement in space or time. Assumes a constant CFL and handles different mesh resolutions for
    //! methods of different orders of accuracy.
    None,

    //! Refinement in space and time with a constant CFL. Fixed CFL is assumed to be given, so we only need to
    //! handle generating the specification for all of the mesh options.
    SpaceTime

};

//! \brief  Used to convert RefinementType values to their corresponding string descriptors.
const std::map< RefinementType, std::string > RefinementType_to_String =
{
    { RefinementType::None,             "none"              },
    { RefinementType::SpaceTime,        "spacetime"         }
};

//! \brief  Used to convert string descriptors to RefinementType values.
const std::map< std::string, RefinementType > String_to_RefinementType =
{
    { "none",           RefinementType::None            },
    { "spacetime",      RefinementType::SpaceTime       }
};

// Key values for list parameters.
const std::string jobscript_template_key = "jobscript_template";
const std::string jobscript_nodes_key = "jobscript_nodes";
const std::string jobscript_ppn_key = "jobscript_ppn";
const std::string jobscript_account_key = "jobscript_account";

const std::string nx_list_key = "nx_vals";
const std::string kba_block_list_key = "kba_block_sizes";

// Key values for input parameters.
const std::string hybrid_str_key = "hybrid_string";
const std::string order_str_key = "order_string";
const std::string angle_str_key = "angle_string";
const std::string mesh_str_key = "mesh_string";

const std::string hybrid_method_key = "hybrid_method";


//------------------------------------------------------------------------------------------------------------
//! \brief  Parses order options for arbitrary methods, and applies some error checking on certain inputs.
//!
//! \param[in]      maker_list      ParameterList containing maker specifications.
//------------------------------------------------------------------------------------------------------------
void
ParseOrderOptions (

    ParameterList & maker_list,
    std::vector<std::string> & orders_vec,
    std::vector<std::string> & dg_degrees_vec,
    const bool remove_values = true
) {

    const std::string orders_list_key = "orders";
    const std::string dg_degrees_list_key = "dg_degrees";

    orders_vec = SplitCSVString( maker_list.GetValue<std::string>( orders_list_key ) );
    if ( remove_values )
        maker_list.RemoveValue( orders_list_key );

    dg_degrees_vec = SplitCSVString( maker_list.GetValue<std::string>( dg_degrees_list_key ) );
    if ( remove_values )
        maker_list.RemoveValue( dg_degrees_list_key );

    if ( dg_degrees_vec.size() != orders_vec.size() )
    {
        std::string error_message = "Number of items in '" + dg_degrees_list_key + "' list ("
                                  + std::to_string( dg_degrees_vec.size() )
                                  + ") does not match the number of items in '" + orders_list_key
                                  + "' list ("
                                  + std::to_string( orders_vec.size() )
                                  + ").\n";

        PRINT_ERROR( error_message.c_str() )
        throw std::invalid_argument( error_message );
    }
}


//------------------------------------------------------------------------------------------------------------
//! \brief  Parses order options for DIRK methods, and applies some error checking on certain inputs.
//!
//! \param[in]      maker_list      ParameterList containing maker specifications.
//------------------------------------------------------------------------------------------------------------
void
ParseDIRKOptions (

    ParameterList & maker_list,
    std::vector<std::string> & dirk_methods_vec,
    const bool remove_values = true
) {

    const std::string dirk_methods_key = "dirk_methods";

    dirk_methods_vec = SplitCSVString( maker_list.GetValue<std::string>( dirk_methods_key ) );
    if ( remove_values )
        maker_list.RemoveValue( dirk_methods_key );

    for ( const std::string & dirk_method_str : dirk_methods_vec )
    {
        try {
            const auto dirk_method = DIRKIntegrator::String_to_DIRKMethod.at( dirk_method_str );
            ((void)(dirk_method));
        }
        catch (...) {
            std::string error_message = "Failed to verify DIRK method string '"
                                      + dirk_method_str + "'.\n";

            PRINT_ERROR( error_message.c_str() )
            throw std::invalid_argument( error_message );
        }
    }
}


//------------------------------------------------------------------------------------------------------------
//! \brief  Parses order options for STDG methods.
//!
//! \param[in]      maker_list      ParameterList containing maker specifications.
//------------------------------------------------------------------------------------------------------------
void
ParseSTDGOptions (

    ParameterList & maker_list,
    std::vector<std::string> & stdg_degrees_vec,
    const bool remove_values = true
) {

    const std::string stdg_degrees_key = "stdg_degrees";

    stdg_degrees_vec = SplitCSVString( maker_list.GetValue<std::string>( stdg_degrees_key ) );
    if ( remove_values )
        maker_list.RemoveValue( stdg_degrees_key );
}


//------------------------------------------------------------------------------------------------------------
//! \brief  Generates a list of the necessary space-time discretizations.
//!
//! \param[in]      maker_list      ParameterList containing maker specifications.
//!
//! \return Returns a vector of ParameterList objects, each of which specifies a space-time mesh
//!         configuration.
//------------------------------------------------------------------------------------------------------------
std::vector<ParameterList>
MakeRefinementList_None (

    ParameterList & maker_list
) {

    std::vector<ParameterList> mesh_vec;

    // Determine (and verify) the type of time integrator.
    TimeIntegratorType time_integrator
        = maker_list.GetValue( TimeIntegrator::GetInputKey(), String_to_TimeIntegratorType );

    // Get common order options.
    std::vector<std::string> orders_vec;
    std::vector<std::string> dg_degrees_vec;
    ParseOrderOptions( maker_list, orders_vec, dg_degrees_vec );

    // Get list of spatial cell counts, and verify length of the list.
    const std::vector<std::string> nx_vec
        = SplitCSVString( maker_list.GetValue<std::string>( nx_list_key ) );
    maker_list.RemoveValue( nx_list_key );

    if ( nx_vec.size() != orders_vec.size() )
    {
        std::string error_message = "Number of items in '" + nx_list_key + "' list ("
                                  + std::to_string( nx_vec.size() )
                                  + ") does not match the expected number of items based on the length of the"
                                    " orders list ("
                                  + std::to_string( orders_vec.size() )
                                  + ").\n";

        PRINT_ERROR( error_message.c_str() )
        throw std::invalid_argument( error_message );
    }

    // Consider this value "optional" and don't cause any errors if it is missing.
    std::vector<std::string> kba_block_vec;
    try {
        kba_block_vec = SplitCSVString( maker_list.GetValue<std::string>( kba_block_list_key ) );
        maker_list.RemoveValue( kba_block_list_key );
    }
    catch (...) {/* empty */}

    // Generate parameter lists for space-time discretizations. Need to handle each time integrator type as a
    // special case.
    switch ( time_integrator )
    {
        case TimeIntegratorType::DIRK:
        {
            std::vector<std::string> dirk_methods_vec;
            ParseDIRKOptions( maker_list, dirk_methods_vec );

            for ( size_t k = 0; k < orders_vec.size(); ++k )
            {
                // Get string arguments for this case.
                const std::string & order_str = orders_vec.at(k);
                const std::string & dirk_method_str = dirk_methods_vec.at(k);
                const std::string & dg_degree_str = dg_degrees_vec.at(k);
                const std::string & nx_str = nx_vec.at(k);

                // Generate the ParameterList for the given discretization.
                ParameterList plist = {{
                    // These are temporary values used for naming input files.
                    { order_str_key, FormatOrderString( order_str ) },
                    { mesh_str_key, FormatMeshString( nx_str ) },

                    // These are the actual input values for the code.
                    { TimeIntegrator::GetInputKey(), TimeIntegratorType_to_String.at( time_integrator ) },
                    { DIRKIntegrator::GetInputKey(), dirk_method_str },
                    { "dg_degree", dg_degree_str },
                    { "nx", nx_str }
                # if SPACE_DIMS >= 2
                ,   { "ny", nx_str }
                # endif
                # if SPACE_DIMS == 3
                ,   { "nz", nx_str }
                # endif
                }};

                // Handle optional KBA parameter.
                try {
                    const std::string & kba_block_str = kba_block_vec.at(k);
                    plist.SetValue( "kba_block_size", kba_block_str );
                }
                catch (...) {/* empty */}

                mesh_vec.push_back( plist );
            }
            break;
        }

        case TimeIntegratorType::STDG:
        {
            std::vector<std::string> stdg_degrees_vec;
            ParseSTDGOptions( maker_list, stdg_degrees_vec );

            for ( size_t k = 0; k < orders_vec.size(); ++k )
            {
                // Get string arguments for this case.
                const std::string & order_str = orders_vec.at(k);
                const std::string & stdg_degree_str = stdg_degrees_vec.at(k);
                const std::string & dg_degree_str = dg_degrees_vec.at(k);
                const std::string & nx_str = nx_vec.at(k);

                // Generate the ParameterList for the given discretization.
                ParameterList plist = {{
                    // These are temporary values used for naming input files.
                    { order_str_key, FormatOrderString( order_str ) },
                    { mesh_str_key, FormatMeshString( nx_str ) },

                    // These are the actual input values for the code.
                    { TimeIntegrator::GetInputKey(), TimeIntegratorType_to_String.at( time_integrator ) },
                    { "dg_degree_t", stdg_degree_str },
                    { "dg_degree_x", dg_degree_str },
                    { "nx", nx_str }
                # if SPACE_DIMS >= 2
                ,   { "ny", nx_str }
                # endif
                # if SPACE_DIMS == 3
                ,   { "nz", nx_str }
                # endif
                }};

                // Handle optional KBA parameter.
                try {
                    const std::string & kba_block_str = kba_block_vec.at(k);
                    plist.SetValue( "kba_block_size", kba_block_str );
                }
                catch (...) {/* empty */}

                mesh_vec.push_back( plist );
            }
            break;
        }

        default:
        {   std::string error_message = "Invalid time integrator '"
                                      + TimeIntegratorType_to_String.at( time_integrator )
                                      + "' in '"
                                      + std::string(__func__)
                                      + "'.\n";

            PRINT_ERROR( error_message.c_str() )
            throw std::invalid_argument( error_message );
        }
    }

    return mesh_vec;
}


//------------------------------------------------------------------------------------------------------------
//! \brief  Generates input options for mesh resolutions used as reference solutions.
//!
//! Want to handle this as a special case, since we only want a minimum number of simulations at this
//! resolution, and only for the high-resolution nonhybrid angular discretizations.
//!
//! \param[in]      maker_list      ParameterList containing maker specifications.
//!
//! \return Returns a vector of ParameterList objects, each of which specifies a set of space/time/angle
//!         options for reference solutions.
//------------------------------------------------------------------------------------------------------------
std::vector<ParameterList>
MakeReferenceList_SpaceTime (

    ParameterList & maker_list
) {
    std::vector<ParameterList> mesh_vec;

    // Determine (and verify) the type of time integrator.
    TimeIntegratorType time_integrator
        = maker_list.GetValue( TimeIntegrator::GetInputKey(), String_to_TimeIntegratorType );

    // Get common order options (use false below so that speficiations are not removed from maker_list yet).
    std::vector<std::string> orders_vec;
    std::vector<std::string> dg_degrees_vec;
    ParseOrderOptions( maker_list, orders_vec, dg_degrees_vec, false );

    const std::string nx_str = maker_list.GetValue<std::string>( "ref_nx" );
    maker_list.RemoveValue( "ref_nx" );

    // Determine angular resolutions. Here just use the u_ resolution specifiers.
    const int64_t u_ang_order_min = maker_list.GetValue<int64_t>( "u_ang_order_min" );
    const int64_t u_ang_order_max = maker_list.GetValue<int64_t>( "u_ang_order_max" );

    // Generate lists.
    switch ( time_integrator )
    {
        case TimeIntegratorType::DIRK:
        {
            std::vector<std::string> dirk_methods_vec;
            ParseDIRKOptions( maker_list, dirk_methods_vec, false );

            for ( size_t order_idx = 0; order_idx < orders_vec.size(); ++order_idx ) {
            for ( int64_t u_ang_order = u_ang_order_min; u_ang_order <= u_ang_order_max; u_ang_order *= 2 ) {

                const std::string & order_str = orders_vec.at( order_idx );
                const std::string & dirk_method_str = dirk_methods_vec.at( order_idx );
                const std::string & dg_degree_str = dg_degrees_vec.at( order_idx );

                ParameterList plist = {{
                    // These are temporary values used for naming input files.
                    { order_str_key, FormatOrderString( order_str ) },
                    { mesh_str_key, FormatMeshString( nx_str ) },
                    { hybrid_str_key, "C" },
                    { angle_str_key, FormatAngleString( u_ang_order ) },

                    // These are the actual input values for the code.
                    { TimeIntegrator::GetInputKey(), TimeIntegratorType_to_String.at( time_integrator ) },
                    { DIRKIntegrator::GetInputKey(), dirk_method_str },
                    { "ang_order", std::to_string( u_ang_order ) },
                    { "hybrid_method", "none" },        // No error checking, but these should all be none.
                # if SPACE_DIMS == 2
                    { "ordinate_sym_reduce", "true" },  // Add this here because why not.
                # endif
                    { "dg_degree", dg_degree_str },
                    { "nx", nx_str }
                # if SPACE_DIMS >= 2
                ,   { "ny", nx_str }
                # endif
                # if SPACE_DIMS == 3
                ,   { "nz", nx_str }
                # endif
                }};

                mesh_vec.push_back( plist );
            }}
            break;
        }

        case TimeIntegratorType::STDG:
        {
            std::vector<std::string> stdg_degrees_vec;
            ParseSTDGOptions( maker_list, stdg_degrees_vec, false );

            for ( size_t order_idx = 0; order_idx < orders_vec.size(); ++order_idx ) {
            for ( int64_t u_ang_order = u_ang_order_min; u_ang_order <= u_ang_order_max; u_ang_order *= 2 ) {

                const std::string & order_str = orders_vec.at( order_idx );
                const std::string & stdg_degree_str = stdg_degrees_vec.at( order_idx );
                const std::string & dg_degree_str = dg_degrees_vec.at( order_idx );

                ParameterList plist = {{
                    // These are temporary values used for naming input files.
                    { order_str_key, FormatOrderString( order_str ) },
                    { mesh_str_key, FormatMeshString( nx_str ) },
                    { hybrid_str_key, "C" },
                    { angle_str_key, FormatAngleString( u_ang_order ) },

                    // These are the actual input values for the code.
                    { TimeIntegrator::GetInputKey(), TimeIntegratorType_to_String.at( time_integrator ) },
                    { "dg_degree_t", stdg_degree_str },
                    { "ang_order", std::to_string( u_ang_order ) },
                    { "hybrid_method", "none" },        // No error checking, but these should all be none.
                # if SPACE_DIMS == 2
                    { "ordinate_sym_reduce", "true" },  // Add this here because why not.
                # endif
                    { "dg_degree_x", dg_degree_str },
                    { "nx", nx_str }
                # if SPACE_DIMS >= 2
                ,   { "ny", nx_str }
                # endif
                # if SPACE_DIMS == 3
                ,   { "nz", nx_str }
                # endif
                }};

                mesh_vec.push_back( plist );
            }}
            break;
        }

        default:
        {   std::string error_message = "Invalid time integrator '"
                                      + TimeIntegratorType_to_String.at( time_integrator )
                                      + "' in '"
                                      + std::string(__func__)
                                      + "'.\n";

            PRINT_ERROR( error_message.c_str() )
            throw std::invalid_argument( error_message );
        }
    }

    return mesh_vec;
}


//------------------------------------------------------------------------------------------------------------
//! \brief  Generates a list of the necessary space-time discretizations.
//!
//! The use case for this refinement type at present is for one-dimensional convergence tests. It is generally
//! sufficient in this setting to run each test case on a single node/workstation, and options relating to
//! scaling on multiple cluster nodes (as would be required for multidimensional refinement tests) are
//! therefore neglected at present.
//!
//! \param[in]      maker_list      ParameterList containing maker specifications.
//!
//! \return Returns a vector of ParameterList objects, each of which specifies a space-time mesh
//!         configuration.
//------------------------------------------------------------------------------------------------------------
std::vector<ParameterList>
MakeRefinementList_SpaceTime (

    ParameterList & maker_list
) {

    std::vector<ParameterList> mesh_vec;

    // Determine (and verify) the type of time integrator.
    TimeIntegratorType time_integrator
        = maker_list.GetValue( TimeIntegrator::GetInputKey(), String_to_TimeIntegratorType );

    // Get common order options.
    std::vector<std::string> orders_vec;
    std::vector<std::string> dg_degrees_vec;
    ParseOrderOptions( maker_list, orders_vec, dg_degrees_vec );

    // Generate list of spatial cell counts, using a factor of two.
    std::vector<std::string> nx_vec;
    {
        const int64_t nx_min = maker_list.GetValue<int64_t>( "nx_min" );
        const int64_t nx_max = maker_list.GetValue<int64_t>( "nx_max" );

        maker_list.RemoveValue( "nx_min" );
        maker_list.RemoveValue( "nx_max" );

        for ( int64_t nx = nx_min; nx <= nx_max; nx *= 2 )
            nx_vec.push_back( std::to_string( nx ) );
    }

    // Generate parameter lists for space-time discretizations. Need to handle each time integrator type as a
    // special case.
    switch ( time_integrator )
    {
        case TimeIntegratorType::DIRK:
        {
            std::vector<std::string> dirk_methods_vec;
            ParseDIRKOptions( maker_list, dirk_methods_vec );

            for ( size_t order_idx = 0; order_idx < orders_vec.size(); ++order_idx )
            {
                // Get string arguments for this case.
                const std::string & order_str = orders_vec.at( order_idx );
                const std::string & dirk_method_str = dirk_methods_vec.at( order_idx );
                const std::string & dg_degree_str = dg_degrees_vec.at( order_idx );

                // Now we loop over the the number of spatial mesh cells and generate each set of inputs.
                for ( const std::string & nx_str : nx_vec )
                {
                    ParameterList plist = {{
                        // These are temporary values used for naming input files.
                        { order_str_key, FormatOrderString( order_str ) },
                        { mesh_str_key, FormatMeshString( nx_str ) },

                        // These are the actual input values for the code.
                        { TimeIntegrator::GetInputKey(), TimeIntegratorType_to_String.at( time_integrator ) },
                        { DIRKIntegrator::GetInputKey(), dirk_method_str },
                        { "dg_degree", dg_degree_str },
                        { "nx", nx_str }
                    # if SPACE_DIMS >= 2
                    ,   { "ny", nx_str }
                    # endif
                    # if SPACE_DIMS == 3
                    ,   { "nz", nx_str }
                    # endif
                    }};

                    mesh_vec.push_back( plist );
                }
            }
            break;
        }

        case TimeIntegratorType::STDG:
        {
            std::vector<std::string> stdg_degrees_vec;
            ParseSTDGOptions( maker_list, stdg_degrees_vec );

            for ( size_t order_idx = 0; order_idx < orders_vec.size(); ++order_idx )
            {
                // Get string arguments for this case.
                const std::string & order_str = orders_vec.at( order_idx );
                const std::string & stdg_degree_str = stdg_degrees_vec.at( order_idx );
                const std::string & dg_degree_str = dg_degrees_vec.at( order_idx );

                // Now we loop over the number of spatial mesh cells and generate each set of inputs.
                for ( const std::string & nx_str : nx_vec )
                {
                    ParameterList plist = {{
                        // These are temporary values used for naming input files.
                        { order_str_key, FormatOrderString( order_str ) },
                        { mesh_str_key, FormatMeshString( nx_str ) },

                        // These are the actual input values for the code.
                        { TimeIntegrator::GetInputKey(), TimeIntegratorType_to_String.at( time_integrator ) },
                        { "dg_degree_t", stdg_degree_str },
                        { "dg_degree_x", dg_degree_str },
                        { "nx", nx_str }
                    # if SPACE_IMDS >= 2
                    ,   { "ny", nx_str }
                    # endif
                    # if SPACE_DIMS == 3
                    ,   { "nz", nx_str }
                    # endif
                    }};

                    mesh_vec.push_back( plist );
                }
            }
            break;
        }

        default:
        {   std::string error_message = "Invalid time integrator '"
                                      + TimeIntegratorType_to_String.at( time_integrator )
                                      + "' in '"
                                      + std::string(__func__)
                                      + "'.\n";

            PRINT_ERROR( error_message.c_str() )
            throw std::invalid_argument( error_message );
        }
    }

    return mesh_vec;
}



//------------------------------------------------------------------------------------------------------------
//! \brief  Generates a list of the necessary angular/hybrid discretizations.
//!
//! \param[in]      maker_list      ParameterList containing maker specifications.
//!
//! \return Returns a vector of ParameterList objects, each of which specifies a configuration for the angular
//!         discretization
//------------------------------------------------------------------------------------------------------------
std::vector<ParameterList>
MakeAnglesList (

    ParameterList & maker_list
) {

    std::vector<ParameterList> angles_vec;

    // Determine hybrid methods.
    const TimeIntegratorType time_integrator
        = maker_list.GetValue( TimeIntegrator::GetInputKey(), String_to_TimeIntegratorType );

    const auto hybrid_methods = ConstructHybridStringsMap( time_integrator );

    // Ordinate configuration options.
    const bool c_include_nonhybrid = maker_list.GetValue<bool>( "c_include_nonhybrid" );
    const bool c_include_equal = maker_list.GetValue<bool>( "c_include_equal" );
    maker_list.RemoveValue( "c_include_nonhybrid" );
    maker_list.RemoveValue( "c_include_equal" );

    const int64_t u_ang_order_min = maker_list.GetValue<int64_t>( "u_ang_order_min" );
    const int64_t u_ang_order_max = maker_list.GetValue<int64_t>( "u_ang_order_max" );
    maker_list.RemoveValue( "u_ang_order_min" );
    maker_list.RemoveValue( "u_ang_order_max" );

    const int64_t c_ang_order_min = maker_list.GetValue<int64_t>( "c_ang_order_min" );
    const int64_t c_ang_order_max = maker_list.GetValue<int64_t>( "c_ang_order_max" );
    maker_list.RemoveValue( "c_ang_order_min" );
    maker_list.RemoveValue( "c_ang_order_max" );

    // Use this option to only output nonhybrid input files for u_* angular resolutions.
    bool include_hybrid = true;
    try {
        include_hybrid = maker_list.GetValue<bool>( "include_hybrid" );
        maker_list.RemoveValue( "include_hybrid" );
    }
    catch (...) {/* empty */}

    //
    // Handle "reference" solution as a special case. This angular resolution is excluded from the hybrid
    // cases, and we specifically override the ordinate_sym_reduce parameter to true (regardless of what is
    // used for the other input files).
    //
    try {
        const int64_t ref_angles = maker_list.GetValue<int64_t>( "ref_angles" );
        maker_list.RemoveValue( "ref_angles" );

        ParameterList plist = {{
            // These are temporary values used for naming input files.
            { hybrid_str_key, "C" },
            { angle_str_key, FormatAngleString( ref_angles ) },

            // These are the actual input values for the code.
            { "ang_order", std::to_string(ref_angles) },
            { "hybrid_method", hybrid_methods.at(0) },
            { "ordinate_sym_reduce", "true" }
        }};

        angles_vec.push_back( plist );
    }
    catch (...) {/* empty */}

    //
    // First do the nonhybrid methods.
    //

    std::set<int64_t> nonhybrid_ordinate_resolutions;

    // Include the uncollided ordinate resolutions.
    for ( int64_t u_ang_order = u_ang_order_min; u_ang_order <= u_ang_order_max; u_ang_order *= 2 )
        nonhybrid_ordinate_resolutions.insert( u_ang_order );

    // Include the collided ordinate resolutions if desired.
    if ( c_include_nonhybrid && include_hybrid )
        for ( int64_t c_ang_order = c_ang_order_min; c_ang_order <= c_ang_order_max; c_ang_order *= 2 )
            nonhybrid_ordinate_resolutions.insert( c_ang_order );

    // Now we loop over the set and avoid duplication.
    for ( const int64_t ang_order : nonhybrid_ordinate_resolutions )
    {
        const std::string angle_str = FormatAngleString( ang_order );

        ParameterList plist = {{
            // These are temporary values used for naming input files.
            { hybrid_str_key, "C" },
            { angle_str_key, angle_str },

            // These are the actual input values for the code.
            { "ang_order", std::to_string(ang_order) },
            { "hybrid_method", hybrid_methods.at(0) }
        }};

        angles_vec.push_back( plist );
    }

    //
    // Next we handle the hybrid cases. Here the collided angles are always no larger than the uncollided,
    // but we can make them equal if that option is set.
    //

    if ( include_hybrid )
    {
        for ( int64_t u_ang_order = u_ang_order_min; u_ang_order <= u_ang_order_max; u_ang_order *= 2 )
        {
            const int64_t c_iter_max = std::min( u_ang_order, c_ang_order_max );

            for ( int64_t c_ang_order = c_ang_order_min; c_ang_order <= c_iter_max; c_ang_order *= 2 )
            {
                if (    (u_ang_order == c_ang_order)
                     && (!c_include_equal)
                ) {  continue;  }

                const std::string & angle_str = FormatAngleString( u_ang_order, c_ang_order );

                for ( const auto & hybrid : hybrid_methods )
                {
                    // Skip nonhybrid methods.
                    if ( hybrid.first == 0 ) {  continue;  }

                    ParameterList plist = {{
                        // These are temporary values used for naming input files.
                        { hybrid_str_key, "H" + std::to_string( hybrid.first ) },
                        { angle_str_key, FormatAngleString( u_ang_order, c_ang_order ) },

                        // These are the actual input values for the code.
                        { "u_ang_order", std::to_string(u_ang_order) },
                        { "c_ang_order", std::to_string(c_ang_order) },
                        { "hybrid_method", hybrid.second }
                    }};

                    angles_vec.push_back( plist );
                }
            }
        }
    }

    return angles_vec;
}


//============================================================================================================
//! \brief  Main driver routine for deckmaker.
//============================================================================================================

int main (

    int argc,
    char ** argv
) {

    if ( argc < 2 )
    {
        const std::string error_message = "Deckmaker expects 1 argument: "
                                          "Filename for deck specifications file.";

        throw std::invalid_argument( error_message );
    }

    ParameterList maker_list;
    maker_list.ReadFromFile( argv[1] );

    std::cout << std::endl;
    std::cout << "Contents of deck specifications file:" << std::endl;
    maker_list.Print();
    std::cout << std::endl;

    // Determine jobscript options.
    const std::string jobscript_template
        = GetJobScriptTemplate( maker_list.GetValue<std::string>( jobscript_template_key ) );
    maker_list.RemoveValue( jobscript_template_key );

    const int64_t jobscript_nodes = maker_list.GetValue<int64_t>( jobscript_nodes_key );
    maker_list.RemoveValue( jobscript_nodes_key );

    const int64_t jobscript_ppn = maker_list.GetValue<int64_t>( jobscript_ppn_key );
    maker_list.RemoveValue( jobscript_ppn_key );

    std::string jobscript_account = "NOACCOUNT";
    try {
        jobscript_account = maker_list.GetValue<std::string>( jobscript_account_key );
        maker_list.RemoveValue( jobscript_account_key );
    }
    catch (...) {/* empty */}

    // Determine how to refine in space and time. Note that for some cases, the resolution for the reference
    // solution(s) requires special handling.
    const RefinementType refinement = maker_list.GetValue( "refinement", String_to_RefinementType );
    maker_list.RemoveValue( "refinement" );

    std::vector<ParameterList> mesh_angle_list_vec;
    std::vector<ParameterList> mesh_list_vec;
    switch ( refinement )
    {
        case RefinementType::None:
            mesh_list_vec = MakeRefinementList_None( maker_list );
            break;

        case RefinementType::SpaceTime:
        {
            try {
                mesh_angle_list_vec = MakeReferenceList_SpaceTime( maker_list );
            }
            catch (...)
            {
                PRINT_WARNING( "Failed to construct list of reference inputs.\n" )
            }

            mesh_list_vec = MakeRefinementList_SpaceTime( maker_list );
            break;
        }

        default:
        {   std::string error_message = "Invalid refinement type '"
                                        + RefinementType_to_String.at( refinement )
                                        + "' in deckmaker.\n";

            PRINT_ERROR( error_message.c_str() )
            throw std::invalid_argument( error_message );
        }
    }

    // Generate lists for angular discretization options.
    const std::vector<ParameterList> angle_list_vec = MakeAnglesList( maker_list );

    std::cout << "Contents of deck specifications after processing:" << std::endl;
    maker_list.Print();
    std::cout << std::endl;

    // Normal output is the outer product of the space-time and angular discretization options.
    for ( const ParameterList & mesh_list : mesh_list_vec ) {
    for ( const ParameterList & angle_list : angle_list_vec ) {

        ParameterList plist = mesh_list + angle_list;
        mesh_angle_list_vec.push_back( plist );
    }}

    // Now we output the outer product of the space-time and angular discretization options.
    std::cout << "Generating output files..." << std::endl;

    for ( const ParameterList & mesh_angle_list : mesh_angle_list_vec )
    {
        ParameterList plist = maker_list + mesh_angle_list;

        const std::string hybrid_string = plist.GetValue<std::string>( hybrid_str_key );
        const std::string order_string = plist.GetValue<std::string>( order_str_key );
        const std::string angle_string = plist.GetValue<std::string>( angle_str_key );
        const std::string mesh_string = plist.GetValue<std::string>( mesh_str_key );

        plist.RemoveValue( hybrid_str_key );
        plist.RemoveValue( order_str_key );
        plist.RemoveValue( angle_str_key );
        plist.RemoveValue( mesh_str_key );

        const std::string filename = hybrid_string  + "_"
                                   + order_string   + "_"
                                   + angle_string   + "_"
                                   + mesh_string;

        std::cout << "\t" << filename << std::endl;

        plist.WriteToFile( filename + ".deck" );
        WriteJobScript( jobscript_template, filename, filename + ".sub", jobscript_ppn, jobscript_nodes, jobscript_account );
    }
}
