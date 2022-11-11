//------------------------------------------------------------------------------------------------------------
//! \file   linear_solvers/RKDG/ImplicitSolver/PETScPeierlsSolver.hpp
//! \brief  Header file containing declarations of PETScPeierlsSolver member function specializations for
//!         RKDG::OrdinateFlux objects.
//!
//! \author Michael M. Crockatt
//! \date   March 2018
//------------------------------------------------------------------------------------------------------------

# ifndef __RKDG__PETSC_PEIERLS_SOLVER_HPP__
# define __RKDG__PETSC_PEIERLS_SOLVER_HPP__


# include "objects/RKDG/OrdinateFlux.hpp"


# if SPACE_DIMS == 1

template<>
void PETScPeierlsSolver<RKDG::OrdinateFlux>::DSA_ComputeMatrix(
    const CrossSection &,
    const CrossSection &,
    const double
);

template<>
void PETScPeierlsSolver<RKDG::OrdinateFlux>::DSA_ComputeRHS(
    const ScalarFlux &,
    ScalarFlux &
);

template<>
double PETScPeierlsSolver<RKDG::OrdinateFlux>::Kappa(
    const int64_t,
    const int64_t,
    const CrossSection &,
    const double
) const;

# endif // if SPACE_DIMS == 1


# endif // ifndef __RKDG__PETSC_PEIERLS_SOLVER_HPP__
