//------------------------------------------------------------------------------------------------------------
//! \file   main_solver.cpp
//! \brief  Contains the main driver function to run a simulation of collisional or non-collisional
//!         transport problems.
//!
//! \authors    Michael M. Crockatt, Kris Garrett
//! \date       Febrary 2018
//------------------------------------------------------------------------------------------------------------


# include <cinttypes>
# include <cmath>
# include <cstdint>
# include <cstdlib>
# include <cstring>
# include <fstream>
# include <limits>
# include <memory>

# if defined (__linux__) || defined (DOXYCOMPILE)
    # include <sys/types.h>
    # include <unistd.h>
# endif

# if defined (ENABLE_MPI) || defined (DOXYCOMPILE)
    # include <mpi.h>
# endif

# if defined (ENABLE_HWLOC) || defined (DOXYCOMPILE)
    # include <hwloc.h>
# endif

# if defined (_OPENMP)
    # include <omp.h>
# endif

# include "sha.h"

# include "objects/RKDG/CrossSection.hpp"
# include "objects/RKDG/DensityFunction.hpp"
# include "objects/RKDG/OrdinateFlux.hpp"
# include "operators/TransportOperator.hpp"
# include "time_integrators/RKDG/dirk.hpp"
# include "time_integrators/RKDG/erk.hpp"
# include "time_integrators/RKDG/idc.hpp"
# include "time_integrators/RKDG/ls-idc.hpp"
# include "time_integrators/STDG/stdg.hpp"
# include "utils/CLog.hpp"
# include "utils/global.hpp"
# include "utils/init.hpp"
# include "utils/ParameterList.hpp"
# include "utils/Quadrule/OrdinateSet.hpp"
# include "utils/Quadrule/Quadrule.hpp"


using namespace RKDG;
using namespace Quadrule;


//------------------------------------------------------------------------------------------------------------
//! \brief  Main driver routine for transport solver.
//------------------------------------------------------------------------------------------------------------
int main (

    int argc,
    char ** argv
) {

# if defined (ENABLE_MPI)

    int required = MPI_THREAD_MULTIPLE;
    int provided = -1;

    MPI_Init_thread( &argc, &argv, required, &provided );

    MPI_Comm_size( MPI_COMM_WORLD, &Global::MPI_num_ranks );
    MPI_Comm_rank( MPI_COMM_WORLD, &Global::MPI_rank );

    if ( provided < required )
        PRINT_WARNING( "Thread support provided by MPI library (%d) is less than requested (%d).\n",
                       provided, required )

# endif // if defined (ENABLE_MPI)

# if defined (ENABLE_PETSC)
    PetscInitialize( &argc, &argv, NULL, NULL );
# endif

# if defined (__linux__)
    const std::string status_match = "VmHWM";
# endif

    const std::string jobid_env_variable = "SLURM_JOB_ID";

    char filename[128];
    char prefix_buff[128];
    char input_file[128];

    bool output_af, output_sf, output_plot, output_mem;

    double final_time;
    double initial_mass, final_mass;

# if defined (DO_SWEEP_SUBTIMING) || defined (ENABLE_HWLOC)

    int num_threads
        # if defined (_OPENMP)
            = omp_get_max_threads();
        # else
            = 1;
        # endif

# endif // if defined (DO_SWEEP_SUBTIMING) || defined (ENABLE_HWLOC)

    DomainDecomposition spatial_params;

    OrdinateFlux & psi    = *new OrdinateFlux();
    OrdinateFlux & source = *new OrdinateFlux();
    CrossSection sigma_t, sigma_s;

    Walltimer init_time;

# if defined (DO_SWEEP_SUBTIMING)

    Global::TMR_calcA       = new Walltimer[ num_threads ];
    Global::TMR_calcB       = new Walltimer[ num_threads ];
    Global::TMR_linearSolve = new Walltimer[ num_threads ];

# endif // if defined (DO_SWEEP_SUBTIMING)

    init_time.Start();


    // --- Determine naming of input and output files. ------------------------------------------------ //

    // If no input filename given, use default.
    if ( argc == 1 ) {  std::strcpy( input_file, "input.deck" );  }
    else             {  std::strcpy( input_file, argv[1] );       }

    // Get prefix from input file.
    char * last_period = std::strrchr( input_file, '.' );

    if ( last_period != NULL ) {

        size_t offset = last_period - input_file;

        std::strncpy( prefix_buff, input_file, offset );
        prefix_buff[offset] = '\0';

    } else {

        // If the input file does not contain a suffix, then use the entire filename as the output prefix.
        std::strcpy( prefix_buff, input_file );
    }

    Global::file_prefix = prefix_buff;


    // --- Read parameters for run from input file. --------------------------------------------------- //

    // Read input file.
    Global::input_list.ReadFromFile( input_file );

    // Get run parameters from input keys.
    Global::DG_degree = GetDGDegreeX( Global::input_list );

    Global::input_list.GetValue( "tfinal", final_time );

    try {  Global::input_list.GetValue( "output_af",    output_af    );  } catch (...) {}
    try {  Global::input_list.GetValue( "output_sd",    output_sf    );  } catch (...) {}
    try {  Global::input_list.GetValue( "output_plot",  output_plot  );  } catch (...) {}
    try {  Global::input_list.GetValue( "output_mem",   output_mem   );  } catch (...) {}


    // --- Begin text output. ------------------------------------------------------------------------- //

    bool output_log;

    Global::input_list.GetValue( "output_log", output_log );

    if (    output_log
    # if defined (ENABLE_MPI)
         && Global::MPI_rank == 0
    # endif
    ) {

        /*  Set filename for text output log file.
         *
         *  Include the job ID in the name if it is found in the execution environment.
         */
        try {

            // Get job ID from environment.
            std::string job_id_str;

            if ( const char * job_id_cstr = std::getenv( jobid_env_variable.c_str() ) ) {

                job_id_str = std::string( job_id_cstr );

            } else {

                throw std::runtime_error( "Variable " + jobid_env_variable
                                          + " not found in environment." );
            }

            // Find first period in job ID string.
            auto first_period = job_id_str.find_first_of( "." );

            // Remove characters after first period if found.
            if ( first_period != std::string::npos )
                job_id_str.erase( first_period, job_id_str.size() );

            std::sprintf( filename, "%s.%s.log", Global::file_prefix.c_str(), job_id_str.c_str() );

        } catch (...) {

            std::sprintf( filename, "%s.log", Global::file_prefix.c_str() );
        }

        Global::logfp = std::fopen( filename, "w" );

        if ( Global::logfp == nullptr )
            PRINT_ERROR( "Failed to open file '%s'.\n", filename )

        // Turn off buffering to log file.
        std::setbuf( stdout, NULL );
    }

    // Output code version and name of input file.
    PRINT_LOG( "\n" )
    PRINT_LOG( "%.*s\n", 75,
               "================================================================================" )

    {
        const std::string sha( reinterpret_cast<const char *>( src_sha ), src_sha_len );

        PRINT_LOG( "    %-13s %s\n", "Version:", sha.c_str() )
    }

    PRINT_LOG( "    %-13s %s\n", "Input file:", input_file )
    PRINT_LOG( "%.*s\n", 75,
               "================================================================================" )
    PRINT_LOG( "\n" )


    // Print heading for problem parameter descriptions.
    PRINT_LOG( "---------------------------------------------------------------------------\n" )
    PRINT_LOG( "--- PROBLEM PARAMETERS ----------------------------------------------------\n" )
    PRINT_LOG( "---------------------------------------------------------------------------\n" )

    // --- Setup hwloc configuration and print core bindings for the threads. ------------------------- //
# if defined (ENABLE_HWLOC)

    // Get machine topology configuration.
    hwloc_topology_init( &Global::machine_topology );
    hwloc_topology_load( Global::machine_topology );

    Global::active_core_mask = hwloc_bitmap_alloc();
    Global::thread_masks = (hwloc_bitmap_s **) std::malloc( num_threads * sizeof(hwloc_cpuset_t) );

    // Get per-thread core bindings.
    # pragma omp parallel
    {
        int tid
            # if defined (_OPENMP)
                = omp_get_thread_num();
            # else
                = 0;
            # endif

        Global::thread_masks[tid] = hwloc_bitmap_alloc();
        hwloc_get_cpubind( Global::machine_topology, Global::thread_masks[tid], HWLOC_CPUBIND_THREAD );

        # pragma omp critical
        hwloc_bitmap_or( Global::active_core_mask, Global::active_core_mask, Global::thread_masks[tid] );
    }

    // Print core bindings to logging interface.
    char * cpubinding;

    PRINT_LOG( "\n" )
    PRINT_LOG( "  OpenMP Thread binding:\n" )
    PRINT_LOG( "\n" )
    PRINT_LOG( "     %s | %s \n", "Thread ID", "Binding" )
    PRINT_LOG( "    -%s-|-%s-\n", "---------", "----------" )

    for ( int tid = 0; tid < num_threads; ++tid ) {

        hwloc_bitmap_asprintf( &cpubinding, Global::thread_masks[tid] );

        PRINT_LOG( "     %7d   | %s\n", tid, cpubinding )
        std::free( cpubinding );
    }

    hwloc_bitmap_asprintf( &cpubinding, Global::active_core_mask );

    PRINT_LOG( "\n" )
    PRINT_LOG( "  %-20s  %-20s\n", "Active cores:", cpubinding )
    std::free( cpubinding );

# endif // if defined (ENABLE_HWLOC)

    // --- Initialize problem. ------------------------------------------------------------------------ //

    Init( psi, source, sigma_t, sigma_s, spatial_params );

    PRINT_LOG( "  %-20s % .4e\n", "Final Time:", final_time )
    PRINT_LOG( "\n" )

    // --- Initialize solvers and contexts. ----------------------------------------------------------- //

    Abstract::TimeIntegrator * integrator = Abstract::TimeIntegrator::CreateIntegrator( Global::input_list );

    integrator->Print();

    // Create object for collided flux if necessary.
    OrdinateFlux * psic = nullptr;
    bool relabel;

    try {

        Global::input_list.GetValue<bool>( "relabel", relabel );

    } catch (...) {  relabel = true;  }

    if (    relabel == false
         && integrator->IsHybrid()
    ) {

        OrdinateSet::OrdinateType c_ordinate_type;
        bool symmetric_reduce = false;

        int64_t c_ang_order = Global::input_list.GetValue<int64_t>( "c_ang_order" );

        try {

            Global::input_list.GetValue( "c_ordinate_type", c_ordinate_type, OrdinateSet::String_to_OrdinateType );

        } catch (...) {

        # if SPACE_DIMS == 1
            c_ordinate_type = OrdinateSet::OrdinateType::GaussLegendre;
        # elif SPACE_DIMS >= 2
            c_ordinate_type = OrdinateSet::OrdinateType::ChebyshevLegendre;
        # endif
        }

    # if SPACE_DIMS == 2

        try         {  Global::input_list.GetValue( "ordinate_sym_reduce", symmetric_reduce );  }
        catch (...) {  symmetric_reduce = false;                                                }

    # endif

        psic = new OrdinateFlux( psi, psi.DG_degree, c_ang_order, symmetric_reduce, c_ordinate_type );
    }

    // --- Compute and print the total mass of the initial condition. --------------------------------- //

    initial_mass = psi.TotalMass();

    PRINT_LOG( "\n" )
    PRINT_LOG( "  %-20s % .4e\n", "Initial Mass:", initial_mass )

    PRINT_LOG( "---------------------------------------------------------------------------\n" )
    PRINT_LOG( "\n" )

    init_time.Stop();

    // Clear timers so that timing of calls during the initialization phase are ignored.
    Global::TMR_AF_Rmv.Clear();
    Global::TMR_Pmv.Clear();
    Global::TMR_AF_copy.Clear();
    Global::TMR_AF_zero.Clear();
    Global::TMR_SD_zero.Clear();
    Global::TMR_SynchronizeHalos.Clear();

    // --- Compute the solution. ---------------------------------------------------------------------- //

    integrator->Integrate( psi, psic, source, sigma_t, sigma_s, final_time );


    // --- Begin outputting results and stats about the solve. ---------------------------------------- //

    // Print heading for result summary.
# if LOGLEVEL == 0
    PRINT_LOG( "\n" )
# endif
    PRINT_LOG( "---------------------------------------------------------------------------\n" )
    PRINT_LOG( "--- SUMMARY ---------------------------------------------------------------\n" )
    PRINT_LOG( "---------------------------------------------------------------------------\n" )


    // Output timing results.
    PRINT_LOG( "  %-20s % .4e\n", "Solve Seconds:", integrator->GetSolveTime() )
    PRINT_LOG( "\n" )

    int total_hours   = integrator->GetSolveTime() / 3600;
    int total_minutes = ( integrator->GetSolveTime() - ( total_hours * 3600 ) ) / 60;
    double total_seconds = ( integrator->GetSolveTime() - ( total_hours * 3600 ) - ( total_minutes * 60 ) );

    PRINT_LOG( "  %-20s % d\n", "Hours:",   total_hours )
    PRINT_LOG( "  %-20s % d\n", "Minutes:", total_minutes )
    PRINT_LOG( "  %-20s % d\n", "Seconds:", ((int) total_seconds) )

    PRINT_LOG( "\n" )
    PRINT_LOG( "  Operation Summary:\n" )
    PRINT_LOG( "\n" )

    double unacc_seconds = integrator->GetSolveTime();

    total_seconds = Global::TMR_AF_Lsv.GetElapsedTime();
    unacc_seconds -= total_seconds;

    PRINT_LOG( "        %-20s % .4e    %05.2f%%\n", "AF Lsv:",
                 total_seconds, 100.0 * total_seconds / integrator->GetSolveTime() )

    total_seconds = Global::TMR_AF_Lmv.GetElapsedTime();
    unacc_seconds -= total_seconds;

    PRINT_LOG( "        %-20s % .4e    %05.2f%%\n", "AF Lmv:",
                 total_seconds, 100.0 * total_seconds / integrator->GetSolveTime() )

    total_seconds = Global::TMR_AF_Rmv.GetElapsedTime();
    unacc_seconds -= total_seconds;

    PRINT_LOG( "        %-20s % .4e    %05.2f%%\n", "AF Rmv:",
                 total_seconds, 100.0 * total_seconds / integrator->GetSolveTime() )

    total_seconds = Global::TMR_Pmv.GetElapsedTime();
    unacc_seconds -= total_seconds;

    PRINT_LOG( "        %-20s % .4e    %05.2f%%\n", "AF Pmv:",
                 total_seconds, 100.0 * total_seconds / integrator->GetSolveTime() )

    total_seconds = Global::TMR_AF_Smv.GetElapsedTime();
    unacc_seconds -= total_seconds;

    PRINT_LOG( "        %-20s % .4e    %05.2f%%\n", "AF Smv:",
               total_seconds, 100.0 * total_seconds / integrator->GetSolveTime() )

    total_seconds = Global::TMR_AF_scal.GetElapsedTime();
    unacc_seconds -= total_seconds;

    PRINT_LOG( "        %-20s % .4e    %05.2f%%\n", "AF scal:",
               total_seconds, 100.0 * total_seconds / integrator->GetSolveTime() )

    total_seconds = Global::TMR_AF_axpy.GetElapsedTime();
    unacc_seconds -= total_seconds;

    PRINT_LOG( "        %-20s % .4e    %05.2f%%\n", "AF axpy:",
               total_seconds, 100.0 * total_seconds / integrator->GetSolveTime() )

    total_seconds = Global::TMR_AF_copy.GetElapsedTime();
    unacc_seconds -= total_seconds;

    PRINT_LOG( "        %-20s % .4e    %05.2f%%\n", "AF copy:",
               total_seconds, 100.0 * total_seconds / integrator->GetSolveTime() )

    total_seconds = Global::TMR_AF_zero.GetElapsedTime();
    unacc_seconds -= total_seconds;

    PRINT_LOG( "        %-20s % .4e    %05.2f%%\n", "AF zero:",
               total_seconds, 100.0 * total_seconds / integrator->GetSolveTime() )

    total_seconds = Global::TMR_AF_boundary.GetElapsedTime();
    unacc_seconds -= total_seconds;

    PRINT_LOG( "        %-20s % .4e    %05.2f%%\n", "AF boundary:",
               total_seconds, 100.0 * total_seconds / integrator->GetSolveTime() )

    PRINT_LOG( "\n" )

    total_seconds = Global::TMR_SD_zero.GetElapsedTime();
    unacc_seconds -= total_seconds;

    PRINT_LOG( "        %-20s % .4e    %05.2f%%\n", "SD zero:",
               total_seconds, 100.0 * total_seconds / integrator->GetSolveTime() )

    total_seconds = Global::TMR_SD_pack.GetElapsedTime();
    unacc_seconds -= total_seconds;

    PRINT_LOG( "        %-20s % .4e    %05.2f%%\n", "SD pack/unpack:",
               total_seconds, 100.0 * total_seconds / integrator->GetSolveTime() )

    total_seconds = Global::TMR_SD_axpy.GetElapsedTime();
    unacc_seconds -= total_seconds;

    PRINT_LOG( "        %-20s % .4e    %05.2f%%\n", "SD axpy:",
               total_seconds, 100.0 * total_seconds / integrator->GetSolveTime() )

    PRINT_LOG( "\n" )

    total_seconds = Global::TMR_PETSc.GetElapsedTime();
    unacc_seconds -= total_seconds;

    PRINT_LOG( "        %-20s % .4e    %05.2f%%\n", "PETSc:",
               total_seconds, 100.0 * total_seconds / integrator->GetSolveTime() )

    total_seconds = Global::TMR_SynchronizeHalos.GetElapsedTime();
    unacc_seconds -= total_seconds;

    PRINT_LOG( "        %-20s % .4e    %05.2f%%\n", "Halo Sync:",
               total_seconds, 100.0 * total_seconds / integrator->GetSolveTime() )

# if SPACE_DIMS == 1

    total_seconds = Global::TMR_DSA_Assemble.GetElapsedTime();
    unacc_seconds -= total_seconds;   // Don't include assembly of diffusion matrix in total solve time.

    PRINT_LOG( "        %-20s % .4e    %05.2f%%\n", "DSA Assemble:",
               total_seconds, 100.0 * total_seconds / integrator->GetSolveTime() )

    total_seconds = Global::TMR_DSA_Solve.GetElapsedTime();
    unacc_seconds -= total_seconds;

    PRINT_LOG( "        %-20s % .4e    %05.2f%%\n", "DSA Solve:",
               total_seconds, 100.0 * total_seconds / integrator->GetSolveTime() )

# endif // if SPACE_DIMS == 1

    PRINT_LOG( "\n" )

    PRINT_LOG( "        %-20s % .4e    %05.2f%%\n", "Unaccounted:",
               unacc_seconds, 100.0 * unacc_seconds / integrator->GetSolveTime() )

    if ( integrator->IsHybrid() ) {

        PRINT_LOG( "\n" )
        PRINT_LOG( "  Hybrid Splitting Summary:\n" )
        PRINT_LOG( "\n" )

        unacc_seconds = integrator->GetSolveTime();

        total_seconds = Global::TMR_collided.GetElapsedTime();
        unacc_seconds -= total_seconds;

        PRINT_LOG( "        %-20s % .4e    %05.2f%%\n", "Collided:",
                   total_seconds, 100.0 * total_seconds / integrator->GetSolveTime() )

        total_seconds = Global::TMR_uncollided.GetElapsedTime();
        unacc_seconds -= total_seconds;

        PRINT_LOG( "        %-20s % .4e    %05.2f%%\n", "Uncollided:",
                   total_seconds, 100.0 * total_seconds / integrator->GetSolveTime() )

        total_seconds = Global::TMR_AF_Rmv.GetElapsedTime()
                        + Global::TMR_relabel.GetElapsedTime();
        unacc_seconds -= total_seconds;

        PRINT_LOG( "        %-20s % .4e    %05.2f%%\n", "Relabel:",
                   total_seconds, 100.0 * total_seconds / integrator->GetSolveTime() )

        PRINT_LOG( "\n" )

        PRINT_LOG( "        %-20s % .4e    %05.2f%%\n", "Unaccounted:",
                   unacc_seconds, 100.0 * unacc_seconds / integrator->GetSolveTime() )
    }

# if defined (DO_SWEEP_SUBTIMING)

    PRINT_LOG( "\n" )
    PRINT_LOG( "  Sweep Sub-timings:\n" )
    PRINT_LOG( "\n" )

    double calcA_seconds = 0.0;
    double calcB_seconds = 0.0;
    double linearSolve_seconds = 0.0;

    // Per-thread subtimings.
    for ( int i = 0; i < num_threads; ++i ) {

        total_seconds = Global::TMR_calcA[i].GetElapsedTime()
                        + Global::TMR_calcB[i].GetElapsedTime()
                        + Global::TMR_linearSolve[i].GetElapsedTime();

        calcA_seconds += Global::TMR_calcA[i].GetElapsedTime();
        calcB_seconds += Global::TMR_calcB[i].GetElapsedTime();
        linearSolve_seconds += Global::TMR_linearSolve[i].GetElapsedTime();

        PRINT_LOG( " %2d     %-20s % .4e    %05.2f%%\n", i, "Calc A:",
                   Global::TMR_calcA[i].GetElapsedTime(),
                   100.0 * Global::TMR_calcA[i].GetElapsedTime() / total_seconds )

        PRINT_LOG( " %2d     %-20s % .4e    %05.2f%%\n", i, "Calc B:",
                   Global::TMR_calcB[i].GetElapsedTime(),
                   100.0 * Global::TMR_calcB[i].GetElapsedTime() / total_seconds )

        PRINT_LOG( " %2d     %-20s % .4e    %05.2f%%\n", i, "Linear Solve:",
                   Global::TMR_linearSolve[i].GetElapsedTime(),
                   100.0 * Global::TMR_linearSolve[i].GetElapsedTime() / total_seconds )

        PRINT_LOG( "\n" )
    }

    // Aggregate timings.
    total_seconds = calcA_seconds + calcB_seconds + linearSolve_seconds;

    PRINT_LOG( "        %-20s % .4e    %05.2f%%\n", "Calc A:",
               calcA_seconds, 100.0 * calcA_seconds / total_seconds )

    PRINT_LOG( "        %-20s % .4e    %05.2f%%\n", "Calc B:",
               calcB_seconds, 100.0 * calcB_seconds / total_seconds )

    PRINT_LOG( "        %-20s % .4e    %05.2f%%\n", "Linear Solve:",
               linearSolve_seconds, 100.0 * linearSolve_seconds / total_seconds )

    PRINT_LOG( "\n" )

# endif // if defined (DO_SWEEP_SUBTIMING)

# if defined (ENABLE_MPI) && defined (ENABLE_KBA_WAIT_TIMING)

    // Wait time of computation threads.
    {
        double sendbuf [] = { Global::TMR_KBA_comp_wait.GetElapsedTime(),
                              100.0 * Global::TMR_KBA_comp_wait.GetElapsedTime() / Global::TMR_AF_Lsv.GetElapsedTime() };

        double * recvbuf = new double[ 2 * Global::MPI_num_ranks ];

        MPI_Gather( sendbuf, 2, MPI_DOUBLE, recvbuf, 2, MPI_DOUBLE, 0, Global::MPI_cart_comm );

        if ( Global::MPI_rank == 0 ) {

            PRINT_LOG( "\n" )

            for ( int rank = 0; rank < Global::MPI_num_ranks; ++rank ) {

                std::string rank_string = "(" + std::to_string(rank) + "/"
                                          + std::to_string(Global::MPI_num_ranks) + ")";

                PRINT_LOG( "%-7s %-20s % .4e    %05.2f%%\n",
                           rank_string.c_str(), "KBA Comp Wait:",
                           recvbuf[ 2*rank ], recvbuf[ 2*rank + 1 ] );
            }
        }

        delete [] recvbuf;
    }

    // Wait time of communication thread.
    {
        double sendbuf [] = { Global::TMR_KBA_comm_wait.GetElapsedTime(),
                              100.0 * Global::TMR_KBA_comm_wait.GetElapsedTime() / Global::TMR_AF_Lsv.GetElapsedTime() };

        double * recvbuf = new double[ 2 * Global::MPI_num_ranks ];

        MPI_Gather( sendbuf, 2, MPI_DOUBLE, recvbuf, 2, MPI_DOUBLE, 0, Global::MPI_cart_comm );

        if ( Global::MPI_rank == 0 ) {

            PRINT_LOG( "\n" )

            for ( int rank = 0; rank < Global::MPI_num_ranks; ++rank ) {

                std::string rank_string = "(" + std::to_string(rank) + "/"
                                          + std::to_string(Global::MPI_num_ranks) + ")";

                PRINT_LOG( "%-7s %-20s % .4e    %05.2f%%\n",
                           rank_string.c_str(), "KBA Comm Wait:",
                           recvbuf[ 2*rank ], recvbuf[ 2*rank + 1 ] );
            }
        }

        delete [] recvbuf;
    }

# endif // if defined (ENABLE_KBA_WAIT_TIMING)

    PRINT_LOG( "\n" )

    PRINT_LOG( "  %-20s % .4e\n", "Init Seconds:", init_time.GetElapsedTime() )

    // Output final mass and change from initial mass.
    final_mass = psi.TotalMass();

    if ( psic != nullptr )
        final_mass += psic->TotalMass();

    PRINT_LOG( "\n" )

    PRINT_LOG( "  %-20s % .4e\n", "Final Mass:", final_mass )
    PRINT_LOG( "  %-20s % .4e\n", "Mass Change:", initial_mass - final_mass )


    // --- Output data file for plotting the total density function (if desired). --------------------- //

    if ( output_plot ) {

        RKDG::DensityFunction temp( psi, psi.DG_degree );
        TransportOperator::Pmv( 1.0, 0.0, psi, temp );

        if ( psic != nullptr )
            TransportOperator::Pmv( 1.0, 1.0, *psic, temp );

        temp.OutputPlot(
            # if SPACE_DIMS == 1
                "idc.out"
            # elif SPACE_DIMS == 2
                "heatmap.out"
            # elif SPACE_DIMS == 3
                "volume.out"
            # endif
                );
    }

    // Output scalar and angular fluxes as desired.
    integrator->OutputSolution( psi, psic );


    // --- Print the peak memory usage as recorded by the OS. ----------------------------------------- //

# if defined (__linux__)

    if ( output_mem ) {

        double mem = std::nan("0");

        const std::string status_filename = "/proc/" + std::to_string( getpid() ) + "/status";

        std::ifstream status_file( status_filename, std::ios_base::in );

        // Print error if /proc/<pid>/status cannot be opened.
        if ( ! status_file.is_open() ) {

            PRINT_ERROR( "Failed to open file '%s'.\n", status_filename.c_str() )

        } else {

            std::string line;
            while ( std::getline( status_file, line ) ) {

                if ( line.find( status_match ) != std::string::npos ) {

                    std::size_t start = line.find_first_of( "0123456789" );
                    std::size_t stop = line.find_last_of( "0123456789" );

                    std::stringstream( line.substr( start, stop +1 ) ) >> mem;
                    break;
                }
            }

            status_file.close();
        }

        mem /= (1 << 20);   // Convert memory from KB to GB.

    # if defined (ENABLE_MPI)

        double * recvbuf = new double[ Global::MPI_num_ranks ];

        MPI_Gather( &mem, 1, MPI_DOUBLE, recvbuf, 1, MPI_DOUBLE, 0, Global::MPI_cart_comm );

        if ( Global::MPI_rank == 0 ) {

            double total_mem = 0.0;

            PRINT_LOG( "\n" )

            for ( int rank = 0; rank < Global::MPI_num_ranks; ++rank ) {

                total_mem += recvbuf[rank];

                std::string rank_string =   status_match
                                          + " ("
                                          + std::to_string(rank)
                                          + "/"
                                          + std::to_string(Global::MPI_num_ranks)
                                          + "):";

                PRINT_LOG( "  %-20s  %10.3f GB\n", rank_string.c_str(), recvbuf[rank] )
            }

            PRINT_LOG( "  %-20s  %10.3f GB\n", "Total:", total_mem )
        }

        delete [] recvbuf;

    # else // if defined (ENABLE_MPI)

        PRINT_LOG( "\n" )
        PRINT_LOG( "  %-20s  %10.3f GB\n", status_match.c_str(), mem )

    # endif // if defined (ENABLE_MPI)
    }

# endif // if defined (__linux__)

    PRINT_LOG( "\n" )
    PRINT_LOG( "  %-20s % .4e\n", "I/O Seconds:", Global::TMR_io.GetElapsedTime() )

    // --- Cleanup memory allocated in main and call cleanup routines for solvers. -------------------- //

    delete [] Global::TPI;
    delete integrator;
    delete psic;

# if defined (ENABLE_PETSC)
    PetscFinalize();
# endif

# if defined (ENABLE_MPI)
    MPI_Finalize();
# endif

# if defined (DO_SWEEP_SUBTIMING)

    delete [] Global::TMR_calcA;
    delete [] Global::TMR_calcB;
    delete [] Global::TMR_linearSolve;

# endif // if defined (DO_SWEEP_SUBTIMING)

# if defined (ENABLE_HWLOC)

    for ( int tid = 0; tid < num_threads; ++tid )
        hwloc_bitmap_free( Global::thread_masks[tid] );

    hwloc_bitmap_free( Global::active_core_mask );
    hwloc_topology_destroy( Global::machine_topology );

# endif // if defined (ENABLE_HWLOC)


    // --- Finish and close text output. -------------------------------------------------------------- //

    // Terminate results block.
    PRINT_LOG( "---------------------------------------------------------------------------\n" )

    // Close output file.
    if ( Global::logfp != nullptr ) {  std::fclose( Global::logfp );  }
}
