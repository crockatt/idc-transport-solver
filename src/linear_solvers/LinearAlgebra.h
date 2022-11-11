//------------------------------------------------------------------------------------------------------------
//! \file   linear_solvers/LinearAlgebra.h
//! \brief  Header file for basic linear algebra routines.
//!
//! \author Michael Crockatt
//! \date   May 2017
//------------------------------------------------------------------------------------------------------------

# ifndef __LINEAR_ALGEBRA_H__
# define __LINEAR_ALGEBRA_H__


# if defined (__cplusplus)
    # define restrict
    extern "C" {
# endif


//============================================================================================================
//=== DECLARATION OF LAPACK ROUTINES =========================================================================
//============================================================================================================


//------------------------------------------------------------------------------------------------------------
//! \brief See <a href="http://www.netlib.org/lapack/explore-html/d7/d15/group__double__blas__level2_gadd421a107a488d524859b4a64c1901a9.html#gadd421a107a488d524859b4a64c1901a9">netlib</a>.
//------------------------------------------------------------------------------------------------------------
int dgemv_( char *, int *, int *, double *, double *, int *, double *, int *, double *, double *, int * );

//------------------------------------------------------------------------------------------------------------
//! \brief See <a href="http://www.netlib.org/lapack/explore-html/d1/d54/group__double__blas__level3_gaeda3cbd99c8fb834a60a6412878226e1.html#gaeda3cbd99c8fb834a60a6412878226e1">netlib</a>.
//------------------------------------------------------------------------------------------------------------
int dgemm_( char *, char *, int *, int *, int *, double *, double *, int *, double *, int *, double *,
            double *, int *);

//------------------------------------------------------------------------------------------------------------
//! \brief See <a href="http://www.netlib.org/lapack/explore-html/dd/d9a/group__double_g_ecomputational_ga2611cc9dfdc84e2a08ec57a5dd6cdd2e.html#ga2611cc9dfdc84e2a08ec57a5dd6cdd2e">netlib</a>.
//------------------------------------------------------------------------------------------------------------
int dgehrd_( int *, int *, int *, double *, int *, double *, double *, int *, int * );

//------------------------------------------------------------------------------------------------------------
//! \brief See <a href="http://www.netlib.org/lapack/explore-html/d1/d7e/group__double_g_esing_gad8e0f1c83a78d3d4858eaaa88a1c5ab1.html#gad8e0f1c83a78d3d4858eaaa88a1c5ab1">netlib</a>.
//------------------------------------------------------------------------------------------------------------
int dgesdd_( char *, int *, int *, double *, int *, double *, double *, int *, double *, int *, double *,
             int *, int *, int * );

//------------------------------------------------------------------------------------------------------------
//! \brief See <a href="http://www.netlib.org/lapack/explore-html/d7/d3b/group__double_g_esolve_ga5ee879032a8365897c3ba91e3dc8d512.html#ga5ee879032a8365897c3ba91e3dc8d512">netlib</a>.
//------------------------------------------------------------------------------------------------------------
int dgesv_( int *, int *, double *, int *, int *, double *, int *, int * );

//------------------------------------------------------------------------------------------------------------
//! \brief See <a href="http://www.netlib.org/lapack/explore-html/d7/d43/group__aux_o_t_h_e_rauxiliary_ga7eb8731ffab2734378157c40964bf788.html#ga7eb8731ffab2734378157c40964bf788">netlib</a>.
//------------------------------------------------------------------------------------------------------------
int dlacpy_( char *, int *, int *, double *, int *, double *, int * );

//------------------------------------------------------------------------------------------------------------
//! \brief See <a href="http://www.netlib.org/lapack/explore-html/de/d39/group__double_g_eauxiliary_gaefa80dbd8cd1732740478618b8b622a1.html">netlib</a>.
//------------------------------------------------------------------------------------------------------------
double dlange_( char *, int *, int *, double *, int *, double * );

//------------------------------------------------------------------------------------------------------------
//! \brief See <a href="http://www.netlib.org/lapack/explore-html/dd/d8a/dorghr_8f_adacfe7750b7fbd625d8101c118174dec.html#adacfe7750b7fbd625d8101c118174dec">netlib</a>.
//------------------------------------------------------------------------------------------------------------
int dorghr_( int *, int *, int *, double *, int *, double *, double *, int *, int * );

//------------------------------------------------------------------------------------------------------------
//! \brief See <a href="http://www.netlib.org/lapack/explore-html/d3/d07/dormhr_8f_a61c5cd2bdf8252f282b30765f8ea8c31.html#a61c5cd2bdf8252f282b30765f8ea8c31">netlib</a>.
//------------------------------------------------------------------------------------------------------------
int dormhr_( char * , char *, int *, int *, int *, int *, double *,
             int *, double *, double *, int *, double *, int *, int * );


//============================================================================================================
//=== LOCAL IMPLEMENTATIONS ==================================================================================
//============================================================================================================


//------------------------------------------------------------------------------------------------------------
//! \brief  Solves the given linear system \f$ Ax = b \f$ using Gaussian elimination without pivoting.
//------------------------------------------------------------------------------------------------------------
void GaussianElimination_CM(
    const int N,
    double * const restrict A,
    double * const restrict b
);

//------------------------------------------------------------------------------------------------------------
//! \brief  Solves a linear system of the form \f$ Hx = b \f$ where \f$ H \f$ is upper Hessenberg.
//------------------------------------------------------------------------------------------------------------
void UpperHessenbergSolve_CM(
    const int N,
    double * const restrict H,
    double * const restrict b,
    const double * const restrict tau,
    int * const restrict piv
);

# if defined (__cplusplus)
    } // extern "C"
#endif


# endif // ifndef __LINEAR_ALGEBRA_H__
