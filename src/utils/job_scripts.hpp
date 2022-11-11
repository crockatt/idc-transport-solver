//------------------------------------------------------------------------------------------------------------
//! \file   utils/job_scripts.hpp
//! \brief  Defines strings which are used to output job scripts. See \ref main_deckmaker.cpp.
//!
//! \author Michael M. Crockatt
//! \date   February 2018
//------------------------------------------------------------------------------------------------------------

# include <string>


//! \brief  Placeholder for job name used in batch scripts.
const std::string jobname_placeholder = "<__JOBNAME__>";

//! \brief  Placeholder for account name used in batch scripts.
const std::string accountname_placeholder = "<__ACCOUNT__>";

//! \brief  Placeholder for number of spatial dimensions (for executable name).
const std::string spacedims_placeholder = "<__SPACE_DIMS__>";

//! \brief  Placeholder for setting OMP_NUM_THREADS environment variable.
const std::string ompnumthreads_placeholder = "<__OMP_NUM_THREADS__>";

//! \brief  Placeholder for setting OMP_PLACES environment variable.
const std::string ompplaces_placeholder = "<__OMP_PLACES__>";

//! \brief  Placeholder for setting npernode flag passed to mpiexec.
const std::string npernode_placeholder = "<__NPERNODE__>";


# include "utils/job_scripts/omp_parallel_restart.shh"

//------------------------------------------------------------------------------------------------------------
//! \brief  Contains script for running checkpointed runs with OpenMP thread parallelism on each node.
//!
//! Assumes that restart capability is provided by restarting after a fixed number of steps (to maintain
//! consistency of timings).
//!
//! Uses N OpenMP threads which are bound to the cores with ID > 0.
//------------------------------------------------------------------------------------------------------------
const std::string omp_parallel_restart(
        reinterpret_cast<const char *>( omp_parallel_restart_sh ),
        omp_parallel_restart_sh_len
    );


# include "utils/job_scripts/mpi_parallel_restart.shh"

//------------------------------------------------------------------------------------------------------------
//! \brief  Contains script for running checkpointed runs with MPI parallelism on each node.
//!
//! Assumes that restart capability is provided by restarting after a fixed number of steps (to maintain
//! consistency of timings).
//!
//! Uses N MPI ranks on each node bound to cores.
//------------------------------------------------------------------------------------------------------------
const std::string mpi_parallel_restart(
        reinterpret_cast<const char *>( mpi_parallel_restart_sh ),
        mpi_parallel_restart_sh_len
    );

