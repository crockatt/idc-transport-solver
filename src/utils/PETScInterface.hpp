//------------------------------------------------------------------------------------------------------------
//! \file   utils/PETScInterface.hpp
//! \brief  Header for common definitions used by solver implementations using PETSc packages.
//!
//! \author Michael Crockatt
//! \date   January 2018
//------------------------------------------------------------------------------------------------------------

# ifndef __PETSC_INTERFACE_HPP__
# define __PETSC_INTERFACE_HPP__

# if defined (ENABLE_PETSC)


# include <petscksp.h>


//!
//! \brief  Typedef of PETSc's Vec (vector) type.
//!
typedef Vec PETScVec;

//!
//! \brief  Typedef of PETSc's Mat (matrix) type.
//!
typedef Mat PETScMat;

//!
//! \brief  Typedef of PETSc's KSP (Krylov subspace solver) type.
//!
typedef KSP PETScKSP;

//!
//! \brief  Typedef of PETSc's PC (preconditioner) type.
//!
typedef PC PETScPC;


# endif // if defined (ENABLE_PETSC)
# endif // ifndef __PETSC_INTERFACE_HPP__
