//------------------------------------------------------------------------------------------------------------
//! \file   defines.h
//! \brief  File containing descriptions of configuration macros.
//!
//! \author Michael Crockatt
//! \date   March 2017
//------------------------------------------------------------------------------------------------------------


//!
//! \brief  Enables ANSI color codes in output if defined during compilation.
//!
//! \see    \c __RED
//! \see    \c __YEL
//! \see    \c __GRN
//! \see    \c __CYN
//! \see    \c __BLU
//! \see    \c __MAG
//! \see    \c __RESET
//!
# define COLOR_TERM


//!
//! \brief  Sets the level of logging to use for status messages.
//!
//! Increasing this value increases logging verbosity.
//! Setting \c #LOGLEVEL=1 is usually a good default choice.
//!
//! \param     0    Prints:
//!                     - Problem configuration.
//!                     - Current timestep.
//!                     - Short summary of results.
//! \param     1    Prints all of the above plus:
//!                     - Number of iterations and final error after each iterative solve.
//!                     - Turns on warning messages.
//! \param     2    Prints all of the above plus:
//!                     - A list of all key-value pairs read from the input file.
//!                     - Turns on performance notes (warnings for things which may be detrimental to performance).
//!                     - Error at each iteration of iterative solves.
//!                     - Prepends each line of output with the file and line number which generated the message.
//! \param     3    Prints all of the above plus turns on verbose status messages.
//!                 Generally only useful for debugging purposes.
//!
# define LOGLEVEL


//!
//! \brief  Enables strict error checks if defined during compilation.
//!
//! Examples of the types of checks include checking for null pointers passed to functions, error messages
//! from library calls, and consistency checks between discretization parameters of angular flux and scalar
//! density objects, etc.
//!
# define STRICT_CHECK


//!
//! \brief  Enables PETSc library routines if defined during compilation.
//!
//! Requires linking a compatible PETSc library.
//! Required to enable the solver type Abstract::ImplicitSolver::SolveType::GMRES.
//!
//! \see    Abstract::ImplicitSolver::SolveType::GMRES
//!
# define ENABLE_PETSC


//!
//! \brief  Enables MPI domain decomposition if defined during compilation.
//!
//! Requires linking a compatible MPI library.
//!
# define ENABLE_MPI


//!
//! \brief  Enables the use of the hwloc library.
//!
//! Requires linking to the hwloc library.
//!
# define ENABLE_HWLOC


//!
//! \brief  Turns on use of hwloc's memory allocation routines to bind memory to NUMA nodes in desirable ways.
//!
//! Requires that hwloc is enabled by setting the \c #ENABLE_HWLOC macro during compilation.
//!
//! \see    \c #ENABLE_HWLOC
//!
# define USE_HWLOC_ALLOC


//!
//! \brief  Enables the use of \c aligned_alloc (added to standard library in C11) in certain locations.
//!
//! Requires linking a more recent version of libc, which may not be available on all systems.
//!
//! The alignment value used by \c aligned_alloc can be adjusted by setting \c #ALIGNED_ALLOC_ALIGNMENT.
//!
//! \see    \c #ALIGNED_ALLOC_ALIGNMENT
//!
# define USE_ALIGNED_ALLOC
