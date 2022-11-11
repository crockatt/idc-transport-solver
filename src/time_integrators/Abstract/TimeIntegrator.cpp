//------------------------------------------------------------------------------------------------------------
//! \file   time_integrators/Abstract/TimeIntegrator.cpp
//! \brief  Implementation of Abstract::TimeIntegrator class.
//!
//! \author Michael M. Crockatt
//! \date   February 2018
//------------------------------------------------------------------------------------------------------------

# include <algorithm>
# include <cmath>
# include <csignal>
# include <iostream>

# if defined (ENABLE_MPI)
    # include <mpi.h>
# endif

# include "objects/RKDG/ScalarFlux.hpp"
# include "operators/TransportOperator.hpp"
# include "time_integrators/Abstract/TimeIntegrator.hpp"
# include "time_integrators/RKDG/dirk.hpp"
# include "time_integrators/RKDG/erk.hpp"
# include "time_integrators/RKDG/idc.hpp"
# include "time_integrators/RKDG/ls-idc.hpp"
# include "time_integrators/STDG/stdg.hpp"
# include "utils/CLog.hpp"
# include "utils/global.hpp"


using namespace Abstract;
using namespace RKDG;


//============================================================================================================
//=== IMPLEMENTATIONS FOR SIGNAL-BASED CHECKPOINT CAPABILITIES ===============================================
//============================================================================================================


//!
//! \namespace  SignalCheckpointData
//!
//! \brief  Contains static variables with data to write during a checkpoint triggered by a received signal.
//!
namespace SignalCheckpointData {

    //!
    //! \brief  Pointer to integrator to checkpoint from.
    //!
    TimeIntegrator * integrator = nullptr;

    //!
    //! \brief  Pointer to uncollided flux.
    //!
    RKDG::OrdinateFlux * uncollided = nullptr;

    //!
    //! \brief  Pointer to collided flux.
    //!
    RKDG::OrdinateFlux * collided = nullptr;
}


//------------------------------------------------------------------------------------------------------------
//! \brief  Signal handler to write checkpoint file on receipt of signal.
//------------------------------------------------------------------------------------------------------------
extern "C"
void CheckpointSignalHandler (

    int sig_num
) {

    PRINT_STATUS( "Executing CheckpointSignalHandler.\n" )
    PRINT_WARNING( "***Received signal %d: Checkpointing program state.\n", sig_num )

    // Stop timers.
    Global::TMR_AF_Lsv.Stop();
    Global::TMR_AF_Lmv.Stop();
    Global::TMR_AF_Rmv.Stop();
    Global::TMR_Pmv.Stop();
    Global::TMR_AF_Smv.Stop();
    Global::TMR_AF_axpy.Stop();
    Global::TMR_AF_scal.Stop();
    Global::TMR_AF_copy.Stop();
    Global::TMR_AF_zero.Stop();
    Global::TMR_AF_boundary.Stop();
    Global::TMR_SD_zero.Stop();
    Global::TMR_SD_axpy.Stop();
    Global::TMR_SD_pack.Stop();
    Global::TMR_PETSc.Stop();
    Global::TMR_SynchronizeHalos.Stop();
    Global::TMR_DSA_Solve.Stop();
    Global::TMR_collided.Stop();
    Global::TMR_uncollided.Stop();
    Global::TMR_relabel.Stop();

    SignalCheckpointData::integrator->solve_time.Stop();

    // Checkpoint program state.
    SignalCheckpointData::integrator->WriteCheckpoint( *SignalCheckpointData::uncollided,
                                                        SignalCheckpointData::collided );

    std::exit( ExitCode::Checkpoint );
}


//============================================================================================================
//=== CONSTRUCTORS, DESTRUCTOR, AND RECONFIGURATION ROUTINES =================================================
//============================================================================================================


//------------------------------------------------------------------------------------------------------------
//! \brief  Constructs an Abstract::TimeIntegrator object with the given set of parameters.
//!
//! \param[in]  input_list      List of parameters to use for initialization.
//------------------------------------------------------------------------------------------------------------
TimeIntegrator::TimeIntegrator (

    const ParameterList & input_list
) :
    DomainDecomposition( input_list ),

    output_steps{false},
    output_sf{false},
    output_af{false},

    current_step {0},
    current_time {0.0},

    CFL { input_list.GetValue<double>( "cfl" ) },
    max_step_size { CFL * dx(0) },

    checkpoint_steps {0},

    checkpoint_wtime {0.0},
    last_wtime {0.0},
    est_wtime {0.0},
    A_wtime {0.3},
    first_step {true}
{
    // Determine output settings.
    try {  this->output_steps = input_list.GetValue<bool>( "output_steps" );  } catch (...) {/* empty */}
    try {  this->output_af    = input_list.GetValue<bool>( "output_af"    );  } catch (...) {/* empty */}
    try {  this->output_sf    = input_list.GetValue<bool>( "output_sd"    );  } catch (...) {/* empty */}

    // Determine checkpoint settings.
    try         {  this->checkpoint_steps = input_list.GetValue<int64_t>( "checkpoint_steps" );  }
    catch (...) {/* empty */}

    if ( this->checkpoint_steps < 0 ) {

        PRINT_WARNING( "Value of checkpoint_steps is negative: setting to zero.\n" )
        this->checkpoint_steps = 0;
    }

    try         {  this->checkpoint_wtime = input_list.GetValue<double>( "checkpoint_walltime" );  }
    catch (...) {/* empty */}

    if ( this->checkpoint_wtime < 0.0 ) {

        PRINT_WARNING( "Value of checkpoint_walltime is negative: setting to zero.\n" )
        this->checkpoint_wtime = 0.0;
    }
}


//------------------------------------------------------------------------------------------------------------
//! \brief  Creates an integrator of the specified type.
//!
//! \param[in]  input_list      List of parameters to use for initialization.
//!
//! \return     Returns a pointer to the newly created integrator object.
//------------------------------------------------------------------------------------------------------------
TimeIntegrator * TimeIntegrator::CreateIntegrator (

    const ParameterList & input_list
) {

    const TimeIntegratorType integrator_type
        = input_list.GetValue<TimeIntegratorType>( GetInputKey(), String_to_TimeIntegratorType );

    switch ( integrator_type ) {

        case TimeIntegratorType::IDC:

            return IDCIntegrator::Create( input_list );

        case TimeIntegratorType::LSIDC:

            return LSIDCIntegrator::Create( input_list );

        case TimeIntegratorType::DIRK:

            return new DIRKIntegrator( input_list );

        case TimeIntegratorType::STDG:

            return new STDGIntegrator( input_list );

        case TimeIntegratorType::ERK:

            return new ERKIntegrator( input_list );

        default:
        {   std::string error_message =   "Invalid time integrator '"
                                        + TimeIntegratorType_to_String.at( integrator_type )
                                        + "' in '"
                                        + std::string(__func__)
                                        + "'.\n";

            PRINT_ERROR( error_message.c_str() )
            throw std::invalid_argument( error_message );
        }
    }
}



//============================================================================================================
//=== INTERFACE ROUTINES FOR TIME MARCHING ===================================================================
//============================================================================================================


//------------------------------------------------------------------------------------------------------------
//! \brief  Integrates from time zero to time \pp{final_time} using a step size \pp{max_step_size}.
//!
//! \attention  If the step size is on the order of \f$ 1.0e-14 \f$ or smaller, this routine is not guaranteed
//!             to be numerically stable in the sense that the number of timesteps used may not be correct.
//!
//! \param[in,out]  initial_condition   Initially contains the initial condition.
//!                                     Upon return, contains the solution at the final time \pp{final_time}.
//! \param[in]      source              RKDG::OrdinateFlux object containing the source term of the transport
//!                                     system.
//! \param[in]      sigma_t             RKDG::CrossSection object containing the total cross section
//!                                     \f$ \sigma_{\mathrm{t}} \f$.
//! \param[in]      sigma_s             RKDG::CrossSection object containing the scattering cross section
//!                                     \f$ \sigma_{\mathrm{s}} \f$.
//! \param[in]      final_time          The final time to integrate the transport system to.
//------------------------------------------------------------------------------------------------------------
TimeIntegrator & TimeIntegrator::Integrate (

    RKDG::OrdinateFlux & initial_condition,
    const RKDG::OrdinateFlux & source,
    const RKDG::CrossSection & sigma_t,
    const RKDG::CrossSection & sigma_s,
    const double final_time
) {

    const double TOL = 1.0e-14;

    solve_time.Start();

    // Output initial condition if desired.
    if ( this->output_steps )
        OutputSolution( initial_condition, nullptr );

    // Compute number of steps.
    double current_step_size = max_step_size;
    const int64_t num_steps = final_time / max_step_size
        + ( std::fmod( final_time, max_step_size ) > TOL ? 1 : 0 );

    // Timestepping loop.
    for ( /* this->current_step set by constructor and ResumeCheckpoint */ ;
          this->current_step < num_steps;
          this->current_step++
    ) {

        // Compute step size.
        if ( (this->current_step + 1) * max_step_size > final_time + TOL )
            current_step_size = final_time - (this->current_step * max_step_size);

        if ( current_step_size != max_step_size ) {

            PRINT_WARNING( "Reducing timestep size for final timestep.\n" )
        }

        PRINT_LOG( "% 4d%%    t = % .4e    dt = % .4e\n",
                   100 * this->current_step / num_steps, this->current_time, current_step_size )

        // Perform step.
        Step( initial_condition, source, sigma_t, sigma_s, current_step_size );

        this->current_time += current_step_size;

        // Output solution if desired.
        if ( this->output_steps )
            OutputSolution( initial_condition, nullptr );
    }

    solve_time.Stop();

    // Enforce strong equality in floating point values at final time.
    if ( this->current_step == num_steps )
        this->current_time = final_time;

    PRINT_LOG( "% 4d%%    t = % .4e\n",
               (int) std::round(100 * this->current_time / final_time), this->current_time )
    PRINT_LOG( "\n" )

    return *this;
}


//------------------------------------------------------------------------------------------------------------
//! \brief  Integrates from time zero to time \pp{final_time} using a step size \pp{max_step_size}.
//!
//! \attention  If the step size is on the order of \f$ 1.0e-14 \f$ or smaller, this routine is not guaranteed
//!             to be numerically stable in the sense that the number of timesteps used may not be correct.
//!
//! \param[in,out]  u_initial_condition Initially contains the initial condition for the uncollided flux.
//!                                     Upon return, contains either the total flux (if relabeling is used)
//!                                     or contains only the uncollided flux (if relabeling is not used) at
//!                                     the final time \pp{final_time}.
//! \param[in,out]  c_initial_condition Initially contains the initial condition for the collided flux.
//!                                     Upon return, contains either a zero flux (if relabeling is used) or
//!                                     contains the collided flux (if relabeling is not used) at the final
//!                                     time \pp{final_time}. <br/>
//!                                     This pointer may be \c null, in which case relabeling is forced if a
//!                                     hybrid integrator is used.
//! \param[in]      source              RKDG::OrdinateFlux object containing the source term of the transport
//!                                     system.
//! \param[in]      sigma_t             RKDG::CrossSection object containing the total cross section
//!                                     \f$ \sigma_{\mathrm{t}} \f$.
//! \param[in]      sigma_s             RKDG::CrossSection object containing the scattering cross section
//!                                     \f$ \sigma_{\mathrm{s}} \f$.
//! \param[in]      final_time          The final time to integrate the transport system to.
//------------------------------------------------------------------------------------------------------------
TimeIntegrator & TimeIntegrator::Integrate (

    RKDG::OrdinateFlux & u_initial_condition,
    RKDG::OrdinateFlux * c_initial_condition,
    const RKDG::OrdinateFlux & source,
    const RKDG::CrossSection & sigma_t,
    const RKDG::CrossSection & sigma_s,
    const double final_time
) {

    const double TOL = 1.0e-14;

    if (    ( this->checkpoint_steps > 0 )
         || ( this->checkpoint_wtime > 0.0 )
    ) {
        try {
            ResumeCheckpoint( u_initial_condition, c_initial_condition );
        }
        catch (...)
        {
            PRINT_WARNING( "Failed to resume from checkpoint file: Starting from beginning.\n" )
            PRINT_LOG( "\n" )
        }

        /*  Only enable signal-handling capabilities for sequential (non-mpi) compilation.
         *
         *  The propagation of signals between MPI ranks on a system is implementation defined. In order to
         *  ensure signal handling would work as intended, it would be necessary for the processes to
         *  synchronize and notify each other of signals received. Since this has not been implemented,
         *  signal handling has been disabled when MPI is enabled.
         */
    # if ! defined (ENABLE_MPI)
        this->SetSignalHandler( u_initial_condition, c_initial_condition );
    # endif
    }

    solve_time.Start();

    // Output initial condition if desired.
    if ( this->output_steps )
        OutputSolution( u_initial_condition, c_initial_condition );

    // Compute number of steps.
    double current_step_size = max_step_size;
    const int64_t num_steps = final_time / max_step_size
        + ( std::fmod( final_time, max_step_size ) > TOL ? 1 : 0 );

    while ( this->current_step < num_steps ) {

        // Compute step size.
        if ( (this->current_step + 1) * max_step_size > final_time + TOL )
            current_step_size = final_time - (this->current_step * max_step_size);

        if ( current_step_size != max_step_size ) {

            PRINT_WARNING( "Reducing timestep size for final timestep.\n" )
        }

        PRINT_LOG( "% 4d%%    t = % .4e    dt = % .4e\n",
                   100 * this->current_step / num_steps, this->current_time, current_step_size )

        // Perform step.
        Step( u_initial_condition, c_initial_condition, source, sigma_t, sigma_s, current_step_size );

        // Increment counters.
        this->current_time += current_step_size;
        this->current_step++;

        // Output solution if desired.
        if ( this->output_steps )
            OutputSolution( u_initial_condition, c_initial_condition );

        // Write checkpoint and exit if necessary. Need to stop the time so that we get the correct estimate
        // from CheckWalltimeCheckpoint.
        solve_time.Stop();

        const bool checkpoint_steps_flag = this->checkpoint_steps
                                        && this->current_step < num_steps
                                        && this->current_step % this->checkpoint_steps == 0;

        if (    checkpoint_steps_flag
             || CheckWalltimeCheckpoint()
        ) {
            PRINT_LOG( "% 4d%%    t = % .4e    dt = % .4e\n",
                       100 * this->current_step / num_steps, this->current_time, current_step_size )
            PRINT_LOG( "\n" )

            WriteCheckpoint( u_initial_condition, c_initial_condition );

            std::exit( ExitCode::Checkpoint );
        }

        solve_time.Start();
    }

    solve_time.Stop();

    // Enforce strong equality in floating point values at final time.
    if ( this->current_step == num_steps )
        this->current_time = final_time;

    PRINT_LOG( "% 4d%%    t = % .4e\n",
               (int) std::round(100 * this->current_time / final_time), this->current_time )
    PRINT_LOG( "\n" )

    return *this;
}


//============================================================================================================
//=== STATIC MEMBER FUNCTIONS ================================================================================
//============================================================================================================


//------------------------------------------------------------------------------------------------------------
//! \brief  Splits keys into three categories: "collided", "uncollided", "none"
//!
//! - Keys prefixed with "c_" are "collided".
//! - Keys prefixed with "u_" are "uncollided".
//! - Keys prefixed with neither "c_" nor "u_" are "none".
//!
//! \return     Returns a map from the strings "collided", "uncollided", "none" to ParameterList objects
//!             containing each category.
//------------------------------------------------------------------------------------------------------------
static std::map< std::string, ParameterList > SplitList (

    const ParameterList & input_list
) {

    ParameterList collided, uncollided, none;

    for ( auto & item : input_list ) {

        if ( item.first.find( "c_" ) == 0 )
            collided.SetValue( item );
        else if ( item.first.find( "u_" ) == 0 )
            uncollided.SetValue( item );
        else
            none.SetValue( item );
    }

    return std::map< std::string, ParameterList >{
            { "collided",   collided },
            { "uncollided", uncollided },
            { "none",       none },
        };
}


//------------------------------------------------------------------------------------------------------------
//! \brief  Prepares a list for use when constructing objects for the collided component.
//!
//! Rules for constructing list:
//! - All pairs with keys that do not begin with "u_" or "c_" are copied exactly.
//! - All pairs with keys that begin with "u_" are removed.
//! - All pairs with keys that begin with "c_" are added back to the list with the prefix "c_" removed from
//!   the key.
//!
//------------------------------------------------------------------------------------------------------------
ParameterList TimeIntegrator::MakeCollidedList (

    const ParameterList & input_list
) {

    auto split_lists = SplitList( input_list );

    ParameterList output_list( split_lists.at( "none" ) );

    for ( auto & item : split_lists.at( "collided" ) ) {

        std::string new_key = item.first;
        new_key.erase(0,2);

        output_list.SetValue( new_key, item.second );
    }

    return output_list;
}


//------------------------------------------------------------------------------------------------------------
//! \brief  Prepares a list for use when constructing objects for the uncollided component.
//!
//! Rules for constructing list:
//! - All pairs with keys that do not begin with "u_" or "c_" are copied exactly.
//! - All pairs with keys that begin with "c_" are removed.
//! - All pairs with keys that begin with "u_" are added back to the list with the prefix "u_" removed from
//!   the key.
//! - The following values are explicitly overridden:
//!     - "solve_type" is set to "sweep".
//!
//------------------------------------------------------------------------------------------------------------
ParameterList TimeIntegrator::MakeUncollidedList (

    const ParameterList & input_list
) {

    auto split_lists = SplitList( input_list );

    ParameterList output_list( split_lists.at( "none" ) );

    for ( auto item : split_lists.at( "uncollided" ) ) {

        std::string new_key = item.first;
        new_key.erase(0,2);

        output_list.SetValue( new_key, item.second );
    }

    output_list.SetValue( "solve_type", "sweep" );

    return output_list;
}


//============================================================================================================
//=== ADDITIONAL MEMBER FUNCTIONS ============================================================================
//============================================================================================================


//------------------------------------------------------------------------------------------------------------
//! \brief  Returns the key used to search input lists for determining the time integrator type to construct.
//------------------------------------------------------------------------------------------------------------
std::string TimeIntegrator::GetInputKey( void ) {

    PRINT_STATUS( "Executing TimeIntegrator::%s.\n", __func__ )

    return "time_integrator";
}



//------------------------------------------------------------------------------------------------------------
//! \brief  Prints a summary of the integrator configuration to the logging interface.
//!
//! \param[in]  prefix  String prepended to each line of output, intended to be used for controlling
//!                     indentation of output.
//------------------------------------------------------------------------------------------------------------
void TimeIntegrator::Print (

    const std::string & prefix // = "  "

) const {

    PRINT_LOG( "%s%-*s  %.4e\n", prefix.c_str(), Global::col_width, "CFL:", this->CFL );
    PRINT_LOG( "%s%-*s  %.4e\n", prefix.c_str(), Global::col_width, "Step size:", this->max_step_size );

    if ( this->checkpoint_steps ) {

        PRINT_LOG( "%s%-*s  %" PRId64 " steps\n", prefix.c_str(),
                   Global::col_width, "Checkpoint:", this->checkpoint_steps );
    } else {

        PRINT_LOG( "%s%-*s  %s\n", prefix.c_str(), Global::col_width, "Checkpoint:", "no" );
    }
}


//------------------------------------------------------------------------------------------------------------
//! \brief  Outputs the solution at the current time.
//!
//! \param[in]  psi             OrdinateFlux containing solution (or its uncollided component).
//! \param[in]  psic            (optional) <br>
//!                             If not NULL, points an OrdinateFlux object containing the collided component
//!                             of the solution.
//! \param[in]  use_timestamp   Specifies whether the output filename should include the current time.
//------------------------------------------------------------------------------------------------------------
void TimeIntegrator::OutputSolution (

    const RKDG::OrdinateFlux & psi,
    const RKDG::OrdinateFlux * const psic,  // = nullptr
    const bool use_timestamp                // = false

) const {

    char time_buff[16] = "";

    if ( use_timestamp )
        std::sprintf( time_buff, "_%.4e", this->current_time );

    if ( this->output_sf ) {

        RKDG::ScalarFlux scalar_temp( psi, psi.DG_degree );
        TransportOperator::Pmv( 1.0, 0.0, psi, scalar_temp );

        if ( psic != nullptr )
            TransportOperator::Pmv( 1.0, 1.0, *psic, scalar_temp );

        scalar_temp.WriteToDisk( Global::file_prefix + std::string( time_buff )
                                                     + std::string( ".sd" )
                                                     + std::to_string( SPACE_DIMS ) );
    }

    if ( this->output_af ) {

        psi.WriteToDisk( Global::file_prefix + std::string( time_buff )
                                             + ( psic == nullptr ? std::string("") : ".u" )
                                             + std::string( ".af" )
                                             + std::to_string( SPACE_DIMS ) );
        if ( psic != nullptr ) {

            psic->WriteToDisk( Global::file_prefix + std::string( time_buff )
                                                   + std::string( ".c" )
                                                   + std::string( ".af" )
                                                   + std::to_string( SPACE_DIMS ) );
        }
    }
}


//============================================================================================================
//=== PROTECTED HELPER FUNCTIONS =============================================================================
//============================================================================================================


# if defined (ENABLE_MPI) || defined (DOXYCOMPILE)

//------------------------------------------------------------------------------------------------------------
//! \brief  Write values, one from each rank, contiguously to the given file pointer.
//!
//! Data is written to the file pointer \pp{fp} starting at its initial value. Upon return, \pp{fp} has been
//! incremented to the first position past the data written.
//!
//! \todo   Add check of return codes.
//!
//! \param[in,out]  fp          File pointer to write data to.
//! \param[in]      buff        Pointer to data to write.
//! \param[in]      count       Number of elements from buff written by each rank.
//! \param[in]      datatype    Type of data at buff to write.
//------------------------------------------------------------------------------------------------------------
static void EachWrite (

    MPI_File & fp,
    const void * const buff,
    int count,
    MPI_Datatype datatype
) {
    int size = 0;

    // Get size of type.
    MPI_Type_size( datatype, &size );

    // Get initial offset of file.
    MPI_Offset start_pos;
    MPI_File_get_position( fp, &start_pos );

    // Seek to write position for each rank.
    MPI_File_seek( fp,
                   start_pos + Global::MPI_rank * size * count,
                   MPI_SEEK_SET );

    // Perform write.
    MPI_File_write( fp, buff, count, datatype, MPI_STATUS_IGNORE );

    // Seek to end position.
    MPI_File_seek( fp,
                   start_pos + Global::MPI_num_ranks * size * count,
                   MPI_SEEK_SET );
}


//------------------------------------------------------------------------------------------------------------
//! \brief  Reads values, one to each rank, contiguously from the given file pointer.
//!
//! Data is read from the file pointer \pp{fp} starting at its initial value. Upon return, \pp{fp} has been
//! incremented to the first position past the data read.
//!
//! \todo   Add check of return codes.
//!
//! \param[in,out]  fp          File pointer to read data from.
//! \param[in]      buff        Pointer to store data into.
//! \param[in]      count       Number of elements from buff read by each rank.
//! \param[in]      datatype    Type of data at buff to read.
//------------------------------------------------------------------------------------------------------------
static int EachRead (

    MPI_File & fp,
    void * const buff,
    int count,
    MPI_Datatype datatype
) {
    int size = 0;

    // Get size of type.
    MPI_Type_size( datatype, &size );

    // Get initial offset of file.
    MPI_Offset start_pos;
    MPI_File_get_position( fp, &start_pos );

    // Seek to read position for each rank.
    MPI_File_seek( fp,
                   start_pos + Global::MPI_rank * size * count,
                   MPI_SEEK_SET );

    // Perform read.
    MPI_File_read( fp, buff, count, datatype, MPI_STATUS_IGNORE );

    // Seek to end position.
    MPI_File_seek( fp,
                   start_pos + Global::MPI_num_ranks * size * count,
                   MPI_SEEK_SET );

    return MPI_SUCCESS;
}

# endif // if defined (ENABLE_MPI)


//------------------------------------------------------------------------------------------------------------
//! \brief  Writes a checkpoint file for the current time.
//!
//! \param[in,out]  uncollided  The current solution (if relabeling is used) or the current uncollided
//!                             solution (if no relabeling).
//! \param[in,out]  collided    (optional) <br>
//!                             The current collided solution (if no relabeling). This value is not
//!                             dereferenced if it is null.
//------------------------------------------------------------------------------------------------------------
void TimeIntegrator::WriteCheckpoint (

    const RKDG::OrdinateFlux & uncollided,
    const RKDG::OrdinateFlux * const collided   // = nullptr

) const {

    char time_buff[16] = "";

    std::sprintf( time_buff, "%04" PRIX64, (uint64_t) this->current_step );

    const std::string filename =   Global::file_prefix
                                 + "."
                                 + std::string(time_buff)
                                 + ".chk";

    Global::TMR_io.Start();

# if defined (ENABLE_MPI)

    MPI_File fp;

    MPI_File_open( Global::MPI_cart_comm, filename.c_str(),
                   MPI_MODE_CREATE | MPI_MODE_RDWR,
                   MPI_INFO_NULL,
                   &fp );

    // Write current state of time integrator.
    if ( Global::MPI_rank == 0 ) {

        MPI_File_write( fp, &this->current_step, 1, MPI_INT64_T, MPI_STATUS_IGNORE );
        MPI_File_write( fp, &this->current_time, 1, MPI_DOUBLE,  MPI_STATUS_IGNORE );

    } else {

        MPI_Offset offset = sizeof(int64_t) + sizeof(double);
        MPI_File_seek( fp, offset, MPI_SEEK_CUR );
    }

    EachWrite( fp, &this->solve_time.GetElapsedTime(), 1, MPI_DOUBLE );

    // Write uncollided flux.
    uncollided.WriteToDisk( fp );

    // Write collided flux.
    if ( collided != nullptr )
        collided->WriteToDisk( fp );

    // Write timer data.
    EachWrite( fp, &Global::TMR_AF_Lsv.GetElapsedTime(),            1, MPI_DOUBLE );
    EachWrite( fp, &Global::TMR_AF_Lmv.GetElapsedTime(),            1, MPI_DOUBLE );
    EachWrite( fp, &Global::TMR_AF_Rmv.GetElapsedTime(),            1, MPI_DOUBLE );
    EachWrite( fp, &Global::TMR_Pmv.GetElapsedTime(),               1, MPI_DOUBLE );
    EachWrite( fp, &Global::TMR_AF_Smv.GetElapsedTime(),            1, MPI_DOUBLE );
    EachWrite( fp, &Global::TMR_AF_axpy.GetElapsedTime(),           1, MPI_DOUBLE );
    EachWrite( fp, &Global::TMR_AF_scal.GetElapsedTime(),           1, MPI_DOUBLE );
    EachWrite( fp, &Global::TMR_AF_copy.GetElapsedTime(),           1, MPI_DOUBLE );
    EachWrite( fp, &Global::TMR_AF_zero.GetElapsedTime(),           1, MPI_DOUBLE );
    EachWrite( fp, &Global::TMR_AF_boundary.GetElapsedTime(),       1, MPI_DOUBLE );
    EachWrite( fp, &Global::TMR_SD_zero.GetElapsedTime(),           1, MPI_DOUBLE );
    EachWrite( fp, &Global::TMR_SD_axpy.GetElapsedTime(),           1, MPI_DOUBLE );
    EachWrite( fp, &Global::TMR_SD_pack.GetElapsedTime(),           1, MPI_DOUBLE );
    EachWrite( fp, &Global::TMR_PETSc.GetElapsedTime(),             1, MPI_DOUBLE );
    EachWrite( fp, &Global::TMR_SynchronizeHalos.GetElapsedTime(),  1, MPI_DOUBLE );
    EachWrite( fp, &Global::TMR_DSA_Solve.GetElapsedTime(),         1, MPI_DOUBLE );
    EachWrite( fp, &Global::TMR_collided.GetElapsedTime(),          1, MPI_DOUBLE );
    EachWrite( fp, &Global::TMR_uncollided.GetElapsedTime(),        1, MPI_DOUBLE );
    EachWrite( fp, &Global::TMR_relabel.GetElapsedTime(),           1, MPI_DOUBLE );
    EachWrite( fp, &Global::TMR_io.GetElapsedTime(),                1, MPI_DOUBLE );

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

    // Write current state of time integrator.
    std::fwrite( &this->current_step, sizeof(int64_t), 1, fp );
    std::fwrite( &this->current_time, sizeof(double),  1, fp );

    std::fwrite( &this->solve_time.GetElapsedTime(), sizeof(double), 1, fp );

    // Write uncollided flux.
    uncollided.WriteToDisk( fp );

    // Write collided flux.
    if ( collided != nullptr )
        collided->WriteToDisk( fp );

    // Write timer data.
    std::fwrite( &Global::TMR_AF_Lsv.GetElapsedTime(),            sizeof(double), 1, fp );
    std::fwrite( &Global::TMR_AF_Lmv.GetElapsedTime(),            sizeof(double), 1, fp );
    std::fwrite( &Global::TMR_AF_Rmv.GetElapsedTime(),            sizeof(double), 1, fp );
    std::fwrite( &Global::TMR_Pmv.GetElapsedTime(),               sizeof(double), 1, fp );
    std::fwrite( &Global::TMR_AF_Smv.GetElapsedTime(),            sizeof(double), 1, fp );
    std::fwrite( &Global::TMR_AF_axpy.GetElapsedTime(),           sizeof(double), 1, fp );
    std::fwrite( &Global::TMR_AF_scal.GetElapsedTime(),           sizeof(double), 1, fp );
    std::fwrite( &Global::TMR_AF_copy.GetElapsedTime(),           sizeof(double), 1, fp );
    std::fwrite( &Global::TMR_AF_zero.GetElapsedTime(),           sizeof(double), 1, fp );
    std::fwrite( &Global::TMR_AF_boundary.GetElapsedTime(),       sizeof(double), 1, fp );
    std::fwrite( &Global::TMR_SD_zero.GetElapsedTime(),           sizeof(double), 1, fp );
    std::fwrite( &Global::TMR_SD_axpy.GetElapsedTime(),           sizeof(double), 1, fp );
    std::fwrite( &Global::TMR_SD_pack.GetElapsedTime(),           sizeof(double), 1, fp );
    std::fwrite( &Global::TMR_PETSc.GetElapsedTime(),             sizeof(double), 1, fp );
    std::fwrite( &Global::TMR_SynchronizeHalos.GetElapsedTime(),  sizeof(double), 1, fp );
    std::fwrite( &Global::TMR_DSA_Solve.GetElapsedTime(),         sizeof(double), 1, fp );
    std::fwrite( &Global::TMR_collided.GetElapsedTime(),          sizeof(double), 1, fp );
    std::fwrite( &Global::TMR_uncollided.GetElapsedTime(),        sizeof(double), 1, fp );
    std::fwrite( &Global::TMR_relabel.GetElapsedTime(),           sizeof(double), 1, fp );
    std::fwrite( &Global::TMR_io.GetElapsedTime(),                sizeof(double), 1, fp );

    std::fclose( fp );

# endif // if defined (ENABLE_MPI)

    Global::TMR_io.Stop();
}


//------------------------------------------------------------------------------------------------------------
//! \brief  Loads program state from a checkpoint file.
//!
//! \param[in,out]  uncollided  The current solution (if relabeling is used) or the current uncollided
//!                             solution (if no relabeling).
//! \param[in,out]  collided    (optional) <br>
//!                             The current collided solution (if no relabeling). This value is not
//!                             dereferenced if it is null.
//------------------------------------------------------------------------------------------------------------
void TimeIntegrator::ResumeCheckpoint (

    RKDG::OrdinateFlux & uncollided,
    RKDG::OrdinateFlux * const collided   // = nullptr
) {

    const std::string filename = Global::file_prefix + ".chk";

    DomainDecomposition spatial_params;
    std::string error_message = "";
    double temp;

    Global::TMR_io.Start();

# if defined (ENABLE_MPI)

//!
//! \brief  Macro to check return values from MPI commands.
//!
# define CHECK_READ \
    if ( ret != MPI_SUCCESS ) {  goto read_error;  }

    MPI_File fp;
    char err_str [MPI_MAX_ERROR_STRING];
    int err_len, ret = MPI_SUCCESS;


    // Try to open file to see if it exists.
    {
        std::FILE * std_fp = std::fopen( filename.c_str(), "rb" );

        if ( std_fp == nullptr ) {

                error_message =   "Failed to open file '"
                                + filename
                                + "' in '"
                                + std::string(__func__)
                                + "'.\n";
            goto read_error;

        } else {

            std::fclose( std_fp );
        }
    }


    ret = MPI_File_open( Global::MPI_cart_comm, filename.c_str(), MPI_MODE_RDONLY, MPI_INFO_NULL, &fp );
    CHECK_READ

    spatial_params = uncollided;

    // Read current state of time integrator.
    ret = MPI_File_read( fp, &this->current_step, 1, MPI_INT64_T, MPI_STATUS_IGNORE ); CHECK_READ
    ret = MPI_File_read( fp, &this->current_time, 1, MPI_DOUBLE,  MPI_STATUS_IGNORE ); CHECK_READ

    // Read time from this->solve time.
    EachRead( fp, &temp, 1, MPI_DOUBLE );   this->solve_time.SetElapsedTime( temp );

    // Read uncollided flux.
    uncollided.ReadFromDisk( fp );

    // Read collided flux.
    if ( collided != nullptr )
        collided->ReadFromDisk( fp );

    // Read timer data.
    EachRead( fp, &temp, 1, MPI_DOUBLE );   Global::TMR_AF_Lsv.             SetElapsedTime( temp );
    EachRead( fp, &temp, 1, MPI_DOUBLE );   Global::TMR_AF_Lmv.             SetElapsedTime( temp );
    EachRead( fp, &temp, 1, MPI_DOUBLE );   Global::TMR_AF_Rmv.             SetElapsedTime( temp );
    EachRead( fp, &temp, 1, MPI_DOUBLE );   Global::TMR_Pmv.                SetElapsedTime( temp );
    EachRead( fp, &temp, 1, MPI_DOUBLE );   Global::TMR_AF_Smv.             SetElapsedTime( temp );
    EachRead( fp, &temp, 1, MPI_DOUBLE );   Global::TMR_AF_axpy.            SetElapsedTime( temp );
    EachRead( fp, &temp, 1, MPI_DOUBLE );   Global::TMR_AF_scal.            SetElapsedTime( temp );
    EachRead( fp, &temp, 1, MPI_DOUBLE );   Global::TMR_AF_copy.            SetElapsedTime( temp );
    EachRead( fp, &temp, 1, MPI_DOUBLE );   Global::TMR_AF_zero.            SetElapsedTime( temp );
    EachRead( fp, &temp, 1, MPI_DOUBLE );   Global::TMR_AF_boundary.        SetElapsedTime( temp );
    EachRead( fp, &temp, 1, MPI_DOUBLE );   Global::TMR_SD_zero.            SetElapsedTime( temp );
    EachRead( fp, &temp, 1, MPI_DOUBLE );   Global::TMR_SD_axpy.            SetElapsedTime( temp );
    EachRead( fp, &temp, 1, MPI_DOUBLE );   Global::TMR_SD_pack.            SetElapsedTime( temp );
    EachRead( fp, &temp, 1, MPI_DOUBLE );   Global::TMR_PETSc.              SetElapsedTime( temp );
    EachRead( fp, &temp, 1, MPI_DOUBLE );   Global::TMR_SynchronizeHalos.   SetElapsedTime( temp );
    EachRead( fp, &temp, 1, MPI_DOUBLE );   Global::TMR_DSA_Solve.          SetElapsedTime( temp );
    EachRead( fp, &temp, 1, MPI_DOUBLE );   Global::TMR_collided.           SetElapsedTime( temp );
    EachRead( fp, &temp, 1, MPI_DOUBLE );   Global::TMR_uncollided.         SetElapsedTime( temp );
    EachRead( fp, &temp, 1, MPI_DOUBLE );   Global::TMR_relabel.            SetElapsedTime( temp );

    EachRead( fp, &temp, 1, MPI_DOUBLE );

    Global::TMR_io.Stop();
    Global::TMR_io.SetElapsedTime( Global::TMR_io.GetElapsedTime() + temp );

    /*  Handle offset for walltime based checkpoint.
     *
     *  When we resume, the value of this->solve_time is set from the checkpoint file. So this->solve_time
     *  stores cumulative solve time (over all restarts), while this->checkpoint_wtime as read from the input
     *  file should account for only the solve time of the current restart. Since checkpoint_wtime has already
     *  been set based on the input file, we update the value here with the current cumulative solve time at
     *  the time when we resume from a restart.
     */
    if ( this->checkpoint_wtime > 0.0 )
    {
        this->last_wtime = this->solve_time.GetElapsedTime();
        this->checkpoint_wtime += this->last_wtime;
    }

    return;


# undef CHECK_READ

read_error:

    MPI_Error_string( ret, err_str, &err_len );
    std::string error_string = std::string( err_str );
    std::replace( error_string.begin(), error_string.end(), '\n', ' ' );

    Global::TMR_io.Stop();

    error_message =   "Failed to resume from checkpoint file '"
                    + filename
                    + "' with error '"
                    + error_string
                    + "'.\n";

    throw std::runtime_error( error_message );

# else // if defined (ENABLE_MPI)

//!
//! \brief  Macro to check return values from \c std::fread for errors.
//!
# define CHECK_READ \
    success &= ((bool)ret); \
    if ( !success ) {  goto read_error;  }


    bool success = true;
    size_t ret;

    std::FILE * fp = std::fopen( filename.c_str(), "rb" );

    if ( fp == nullptr ) {

        error_message =   "Failed to open file '"
                        + filename
                        + "' in '"
                        + std::string(__func__)
                        + "'.\n";
        goto read_error;
    }

    spatial_params = uncollided;

    ret = std::fread( &this->current_step, sizeof(int64_t), 1, fp ); CHECK_READ
    ret = std::fread( &this->current_time, sizeof(double),  1, fp ); CHECK_READ

    ret = std::fread( &temp, sizeof(double), 1, fp ); CHECK_READ  this->solve_time.SetElapsedTime( temp );

    uncollided.ReadFromDisk( fp );

    if ( collided != nullptr )
        collided->ReadFromDisk( fp );

    ret = std::fread( &temp, sizeof(double), 1, fp ); CHECK_READ  Global::TMR_AF_Lsv.           SetElapsedTime( temp );
    ret = std::fread( &temp, sizeof(double), 1, fp ); CHECK_READ  Global::TMR_AF_Lmv.           SetElapsedTime( temp );
    ret = std::fread( &temp, sizeof(double), 1, fp ); CHECK_READ  Global::TMR_AF_Rmv.           SetElapsedTime( temp );
    ret = std::fread( &temp, sizeof(double), 1, fp ); CHECK_READ  Global::TMR_Pmv.              SetElapsedTime( temp );
    ret = std::fread( &temp, sizeof(double), 1, fp ); CHECK_READ  Global::TMR_AF_Smv.           SetElapsedTime( temp );
    ret = std::fread( &temp, sizeof(double), 1, fp ); CHECK_READ  Global::TMR_AF_axpy.          SetElapsedTime( temp );
    ret = std::fread( &temp, sizeof(double), 1, fp ); CHECK_READ  Global::TMR_AF_scal.          SetElapsedTime( temp );
    ret = std::fread( &temp, sizeof(double), 1, fp ); CHECK_READ  Global::TMR_AF_copy.          SetElapsedTime( temp );
    ret = std::fread( &temp, sizeof(double), 1, fp ); CHECK_READ  Global::TMR_AF_zero.          SetElapsedTime( temp );
    ret = std::fread( &temp, sizeof(double), 1, fp ); CHECK_READ  Global::TMR_AF_boundary.      SetElapsedTime( temp );
    ret = std::fread( &temp, sizeof(double), 1, fp ); CHECK_READ  Global::TMR_SD_zero.          SetElapsedTime( temp );
    ret = std::fread( &temp, sizeof(double), 1, fp ); CHECK_READ  Global::TMR_SD_axpy.          SetElapsedTime( temp );
    ret = std::fread( &temp, sizeof(double), 1, fp ); CHECK_READ  Global::TMR_SD_pack.          SetElapsedTime( temp );
    ret = std::fread( &temp, sizeof(double), 1, fp ); CHECK_READ  Global::TMR_PETSc.            SetElapsedTime( temp );
    ret = std::fread( &temp, sizeof(double), 1, fp ); CHECK_READ  Global::TMR_SynchronizeHalos. SetElapsedTime( temp );
    ret = std::fread( &temp, sizeof(double), 1, fp ); CHECK_READ  Global::TMR_DSA_Solve.        SetElapsedTime( temp );
    ret = std::fread( &temp, sizeof(double), 1, fp ); CHECK_READ  Global::TMR_collided.         SetElapsedTime( temp );
    ret = std::fread( &temp, sizeof(double), 1, fp ); CHECK_READ  Global::TMR_uncollided.       SetElapsedTime( temp );
    ret = std::fread( &temp, sizeof(double), 1, fp ); CHECK_READ  Global::TMR_relabel.          SetElapsedTime( temp );

    ret = std::fread( &temp, sizeof(double), 1, fp ); CHECK_READ

    Global::TMR_io.Stop();
    Global::TMR_io.SetElapsedTime( Global::TMR_io.GetElapsedTime() + temp );

    /*  Handle offset for walltime based checkpoint.
     *
     *  When we resume, the value of this->solve_time is set from the checkpoint file. So this->solve_time
     *  stores cumulative solve time (over all restarts), while this->checkpoint_wtime as read from the input
     *  file should account for only the solve time of the current restart. Since checkpoint_wtime has already
     *  been set based on the input file, we update the value here with the current cumulative solve time at
     *  the time when we resume from a restart.
     */
    if ( this->checkpoint_wtime > 0.0 )
    {
        this->last_wtime = this->solve_time.GetElapsedTime();
        this->checkpoint_wtime += this->last_wtime;
    }

    return;


# undef CHECK_READ

read_error:

    Global::TMR_io.Stop();

    if ( error_message == "" ) {

        error_message =   "Failed to resume from checkpoint file '"
                      + std::string( filename )
                      + "'.\n";
    }

    throw std::runtime_error( error_message );

# endif // if defined (ENABLE_MPI)
}


//------------------------------------------------------------------------------------------------------------
//! \brief  Sets a signal handler to checkpoint the program state using the given objects.
//!
//! \param[in]  uncollided  Object that will contain the current solution (if relabeling is used) or the
//!                         current uncollided solution (if no relabeling).
//! \param[in]  collided    (optional) <br>
//!                         The current collided solution (if no relabeling). This value is not dereferenced
//!                         if it is null.
//------------------------------------------------------------------------------------------------------------
void TimeIntegrator::SetSignalHandler (

    RKDG::OrdinateFlux & uncollided,
    RKDG::OrdinateFlux * const collided // = nullptr
) {

    SignalCheckpointData::integrator = this;
    SignalCheckpointData::uncollided = &uncollided;
    SignalCheckpointData::collided   = collided;

    std::signal( SIGINT,  CheckpointSignalHandler );
    std::signal( SIGTERM, CheckpointSignalHandler );
}


//------------------------------------------------------------------------------------------------------------
//! \brief  Checks walltime estimate to see if we need to checkpoint.
//!
//! Checkpoint timing is disabled if \pp{checkpoint_wtime} is zero.
//!
//! \return     Returns true if we need to checkpoint at the current time, and false otherwise.
//------------------------------------------------------------------------------------------------------------
bool TimeIntegrator::CheckWalltimeCheckpoint ( void ) {

    // Check if checkpoint timing is disabled.
    if ( this->checkpoint_wtime == 0.0 ) {  return false;  }

    // Initialize estimate after first step, otherwise update estimate with exponential weighting.
    if ( this->first_step )
    {
        this->first_step = false;
        this->est_wtime = this->solve_time.GetElapsedTime() - this->last_wtime;
    }
    else
    {
        const double step_wtime = this->solve_time.GetElapsedTime() - this->last_wtime;
        this->est_wtime = this->A_wtime * step_wtime + (1.0 - this->A_wtime) * this->est_wtime;
    }

    // Update the stored time.
    this->last_wtime = this->solve_time.GetElapsedTime();

    // Determine the result of the estimate.
    bool do_checkpoint = false;

    const double wtime_left = this->checkpoint_wtime - this->last_wtime;
    const int64_t steps_left = wtime_left / this->est_wtime;

    if (    ( wtime_left <= 0.0 )
         || ( steps_left <= 0 )
    ) {
        do_checkpoint = true;
    }

# if defined (ENABLE_MPI)
    if ( Global::MPI_rank == 0 ) {
# endif
        PRINT_NOTE( "Current walltime used: %.4e\n", this->last_wtime )
        PRINT_NOTE( "Checkpoint walltime:   %.4e\n", this->checkpoint_wtime )
        PRINT_NOTE( "Walltime remaining:    %.4e\n", wtime_left )
        PRINT_NOTE( "Current estimate:      %.4e\n", this->est_wtime )
        PRINT_NOTE( "Estimated steps left:  %" PRId64 "\n", steps_left )
# if defined (ENABLE_MPI)
    }
# endif

# if defined (ENABLE_MPI)
    MPI_Bcast( &do_checkpoint, 1, MPI_CXX_BOOL, 0, Global::MPI_cart_comm );
# endif

    if (    do_checkpoint
    # if defined (ENABLE_MPI)
         && Global::MPI_rank == 0
    # endif
    ) {
        PRINT_LOG( "Checkpointing current run based on walltime estimate:\n" )
        PRINT_LOG( "    Current walltime used: %.4e\n", this->last_wtime )
        PRINT_LOG( "    Checkpoint walltime:   %.4e\n", this->checkpoint_wtime )
        PRINT_LOG( "    Current estimate:      %.4e\n", this->est_wtime )
        PRINT_LOG( "\n" )
    }

    return do_checkpoint;
}
